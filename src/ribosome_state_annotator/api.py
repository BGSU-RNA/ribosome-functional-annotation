"""High-level library API (spec §14.1, §14.2).

This module is the outermost orchestrator. It glues together every other
module:

::

    rcsb_client.fetch_entry_payload
        → rcsb_client.parse_assemblies
            → classify.classify_assembly
                → bgsu_client.fetch_correspondence (per site)
                    → correspondence.build_correspondence_result
                        → coordinates.get_assembly_structure
                            → infer.assign_functional_chains
                                → infer.compute_trna_states
                                    → RibosomeAnnotation

The cache wrapping layer lives in this module: each external call goes
through ``_get_or_*`` helpers that check the :class:`Cache` before
falling through to the network/disk fetch. The per-module clients stay
pure HTTP — they don't know about the cache.

Per spec §14.1, three public entry points:

- :func:`annotate_pdb` — annotate every biological assembly in one PDB.
- :func:`annotate_assembly` — annotate one specific assembly.
- :func:`annotate_many` — batch wrapper over a list of PDB IDs.
"""

from __future__ import annotations

import logging
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any

import gemmi
import httpx

from ribosome_state_annotator import constants as C
from ribosome_state_annotator.bgsu_client import (
    DEFAULT_BGSU_DEPTH,
    DEFAULT_BGSU_RESOLUTION,
    DEFAULT_BGSU_SCOPE,
    fetch_correspondence,
    parse_alignment_response,
)
from ribosome_state_annotator.cache import Cache
from ribosome_state_annotator.classify import (
    ClassificationResult,
    classify_assembly,
    partition_rna_chains_by_role,
)
from ribosome_state_annotator.config import CompletenessThresholds
from ribosome_state_annotator.coordinates import (
    CoordinateSource,
    get_assembly_structure,
)
from ribosome_state_annotator.correspondence import build_correspondence_result
from ribosome_state_annotator.exceptions import (
    ApiRequestError,
    CoordinateDownloadError,
    CoordinateParsingError,
    CorrespondenceMappingError,
    RibosomeAnnotatorError,
)
from ribosome_state_annotator.infer import (
    ChainAssignments,
    TRNAStates,
    assign_functional_chains,
    compute_trna_states,
)
from ribosome_state_annotator.models import (
    AssemblyContext,
    ChainRef,
    CorrespondenceResult,
    LargeScaleMovements,
    LigandRef,
    RibosomeAnnotation,
)
from ribosome_state_annotator.multiribo import (
    detect_multi_ribosome,
    pair_ssu_lsu_by_centroid,
    partition_chains_by_ribosome,
)
from ribosome_state_annotator.raddb import (
    RADdbDataset,
    ensure_raddb_available,
    get_motion_metrics,
    load_raddb_dataset,
)
from ribosome_state_annotator.rcsb_client import (
    fetch_entry_payload,
    parse_assemblies,
)
from ribosome_state_annotator.taxonomy import aggregate_assembly_lineage
from ribosome_state_annotator.rfam_pdb_region import (
    RfamPdbRegionDataset,
    ensure_rfam_pdb_region_available,
    get_rfam_mapping_for_pdb,
    load_rfam_pdb_region_dataset,
)
from ribosome_state_annotator.trna_mrna import extract_trna_mrna_interactions

logger = logging.getLogger(__name__)

WARNING_MULTIPLE_CHAINS_FOR_LEGACY_CSV = "multiple_chains_for_legacy_csv_column"
"""Spec §15.3: emitted when >1 SSU or >1 LSU main rRNA chain is present
so a downstream legacy-CSV column is left blank."""


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def annotate_pdb(
    pdb_id: str,
    *,
    assembly_id: str | None = None,
    contact_cutoff_angstrom: float = C.DEFAULT_CONTACT_CUTOFF_ANGSTROM,
    cache: Cache | None = None,
    cache_dir: Path | None = None,
    no_cache: bool = False,
    strict_complete_check: bool = False,
    completeness_thresholds: CompletenessThresholds | None = None,
    coordinate_source: CoordinateSource = "auto",
    local_coordinate_path: Path | None = None,
    client: httpx.Client | None = None,
    raddb_dataset: RADdbDataset | None = None,
    refresh_raddb: bool = False,
    no_raddb: bool = False,
    rfam_dataset: RfamPdbRegionDataset | None = None,
    refresh_rfam: bool = False,
    no_rfam: bool = False,
    no_fr3d: bool = False,
) -> list[RibosomeAnnotation]:
    """Annotate every biological assembly in one PDB entry.

    Args:
        pdb_id: Four-character PDB accession (case-insensitive).
        assembly_id: Restrict to one biological-assembly ID. When ``None``
            every assembly in the entry is processed.
        contact_cutoff_angstrom: Gemmi neighbour-search cutoff (spec §10.2).
        cache: Pre-built :class:`Cache`. Takes precedence over
            ``cache_dir`` / ``no_cache``.
        cache_dir: Directory to back the cache, used when ``cache`` is None.
        no_cache: Disable caching entirely (no read or write).
        strict_complete_check: Spec §7.4 strict mode — skip assemblies
            below the ribosomal-protein-count threshold instead of warning.
        completeness_thresholds: Override the default §7.4 thresholds.
        coordinate_source: ``"auto"`` (download/cache) or ``"local"``.
        local_coordinate_path: Required when ``coordinate_source == "local"``.
        client: Optional shared :class:`httpx.Client`. When omitted, each
            sub-call opens and closes its own.
        raddb_dataset: Pre-loaded RADdb dataset. When omitted, the
            function loads it lazily on first use (and reuses it across
            assemblies in this call). Pass an explicit dataset to share
            it across many ``annotate_pdb`` calls in a batch.
        refresh_raddb: Force an online check for a newer RADdb release
            even if the local copy is fresh. No-op when ``raddb_dataset``
            is supplied explicitly.
        no_raddb: Skip RADdb integration entirely. The output JSON still
            contains a ``large_scale_movements`` block with ``rad_date=None``
            and null metrics so consumers see a stable schema.
        rfam_dataset: Pre-loaded Rfam pdb_full_region dataset. When
            omitted, loaded lazily on first use and reused across
            assemblies in this call.
        refresh_rfam: Force an online check for a newer Rfam
            pdb_full_region file even if the local copy is fresh.
        no_rfam: Skip the Rfam file augmentation entirely. RCSB-supplied
            Rfam tags (when present) are still used.
        no_fr3d: Skip the tRNA-mRNA codon/anticodon extraction (no FR3D
            call). The output JSON still contains
            ``trna_mrna_interactions`` as an empty list.

    Returns:
        A list of :class:`RibosomeAnnotation` — one per processed assembly
        (or per skipped/failed result). Entry-level skips (NMR) collapse
        to a single annotation with ``assembly_id=None`` per spec §9.1.
    """
    resolved_cache = _resolve_cache(cache, cache_dir, no_cache)
    pdb_id_upper = pdb_id.upper()
    logger.info("annotating %s%s", pdb_id_upper, f" (assembly {assembly_id})" if assembly_id else "")

    try:
        entry_payload = _get_or_fetch_entry(pdb_id_upper, resolved_cache, client)
    except ApiRequestError as exc:
        logger.error("RCSB fetch failed for %s: %s", pdb_id_upper, exc)
        return [
            RibosomeAnnotation(
                pdb_id=pdb_id_upper,
                assembly_id=None,
                status="failed",
                skip_reason=f"rcsb_fetch_failed: {exc}",
            )
        ]

    try:
        assemblies = parse_assemblies(entry_payload)
    except ApiRequestError as exc:
        logger.error("RCSB parse failed for %s: %s", pdb_id_upper, exc)
        return [
            RibosomeAnnotation(
                pdb_id=pdb_id_upper,
                assembly_id=None,
                status="failed",
                skip_reason=f"rcsb_parse_failed: {exc}",
            )
        ]

    # Entry-level NMR skip (§7.3) — one annotation, not one per assembly.
    methods = assemblies[0].experimental_methods if assemblies else []
    if any("NMR" in m.upper() for m in methods):
        return [
            RibosomeAnnotation(
                pdb_id=pdb_id_upper,
                assembly_id=None,
                status="skipped",
                skip_reason=C.SKIP_NMR,
            )
        ]

    if assembly_id is not None:
        assemblies = [a for a in assemblies if a.assembly_id == assembly_id]
        if not assemblies:
            logger.warning("no assembly with id %s in %s", assembly_id, pdb_id_upper)
            return []

    # Rfam augmentation: RCSB no longer supplies Rfam annotations for
    # rRNA chains. We pull single-best-score Rfam tags from the EBI
    # ``pdb_full_region.txt.gz`` flat file (locally cached, weekly
    # refresh), which selects one Rfam per chain by highest bit-score
    # and eliminates the multi-family over-annotation noise PDBe's REST
    # endpoint surfaces (the source of the historical MIXED rrna_core
    # edge case on entries like 9B0S).
    if no_rfam:
        resolved_rfam: RfamPdbRegionDataset | None = None
    elif rfam_dataset is not None:
        resolved_rfam = rfam_dataset
    else:
        resolved_rfam = _load_rfam_safely(client=client, force_refresh=refresh_rfam)

    rfam_by_chain = get_rfam_mapping_for_pdb(resolved_rfam, pdb_id_upper)
    if rfam_by_chain:
        for assembly in assemblies:
            _apply_rfam_pdb_region(assembly.rna_chains, rfam_by_chain)

    if no_raddb:
        resolved_raddb: RADdbDataset | None = None
    elif raddb_dataset is not None:
        resolved_raddb = raddb_dataset
    else:
        resolved_raddb = _load_raddb_safely(client=client, force_refresh=refresh_raddb)

    results: list[RibosomeAnnotation] = []
    for assembly in assemblies:
        results.extend(
            _annotate_one_assembly(
                assembly,
                cache=resolved_cache,
                cutoff=contact_cutoff_angstrom,
                strict_complete_check=strict_complete_check,
                completeness_thresholds=completeness_thresholds,
                coordinate_source=coordinate_source,
                local_coordinate_path=local_coordinate_path,
                client=client,
                raddb_dataset=resolved_raddb,
                no_fr3d=no_fr3d,
            )
        )
    return results


def annotate_assembly(
    pdb_id: str,
    assembly_id: str,
    **kwargs: Any,
) -> RibosomeAnnotation:
    """Annotate one biological assembly. Convenience wrapper.

    Returns the :class:`RibosomeAnnotation` for ``(pdb_id, assembly_id)``.

    Multi-ribosome bundles (e.g. ``8R3V`` assembly ``1``) are split into
    sub-annotations with suffixed assembly IDs (``"1-1"``, ``"1-2"``).
    ``assembly_id="1-1"`` returns that specific sub-ribosome;
    ``assembly_id="1"`` returns the first sub-ribosome (with a warning
    on the annotation). Callers that need every sub-ribosome should use
    :func:`annotate_pdb` instead.
    """
    base_assembly_id = assembly_id.split("-", 1)[0]
    results = annotate_pdb(pdb_id, assembly_id=base_assembly_id, **kwargs)
    if not results:
        return RibosomeAnnotation(
            pdb_id=pdb_id.upper(),
            assembly_id=assembly_id,
            status="failed",
            skip_reason=f"assembly_not_found: {assembly_id}",
        )
    for ann in results:
        if ann.assembly_id == assembly_id:
            return ann
    return results[0]


def annotate_many(
    pdb_ids: Iterable[str],
    *,
    continue_on_error: bool = True,
    refresh_raddb: bool = False,
    raddb_dataset: RADdbDataset | None = None,
    no_raddb: bool = False,
    refresh_rfam: bool = False,
    rfam_dataset: RfamPdbRegionDataset | None = None,
    no_rfam: bool = False,
    no_fr3d: bool = False,
    client: httpx.Client | None = None,
    **kwargs: Any,
) -> list[RibosomeAnnotation]:
    """Batch wrapper. Iterates ``pdb_ids`` and calls :func:`annotate_pdb` on each.

    With ``continue_on_error=True`` (the default), per-entry exceptions
    are caught and recorded as a ``status="failed"`` annotation rather
    than propagating. Pass ``continue_on_error=False`` to abort the
    batch on the first error instead.

    ``refresh_raddb`` and ``refresh_rfam`` are honoured once at the
    start of the batch; the resulting datasets (or the explicit
    ``raddb_dataset`` / ``rfam_dataset`` if supplied) are threaded
    through every per-entry call so each is loaded at most once per
    batch.
    """
    aggregated: list[RibosomeAnnotation] = []
    pdb_ids_list = list(pdb_ids)
    total = len(pdb_ids_list)
    if no_raddb:
        resolved_raddb: RADdbDataset | None = None
    elif raddb_dataset is not None:
        resolved_raddb = raddb_dataset
    else:
        resolved_raddb = _load_raddb_safely(client=client, force_refresh=refresh_raddb)

    if no_rfam:
        resolved_rfam: RfamPdbRegionDataset | None = None
    elif rfam_dataset is not None:
        resolved_rfam = rfam_dataset
    else:
        resolved_rfam = _load_rfam_safely(client=client, force_refresh=refresh_rfam)

    for index, pdb_id in enumerate(pdb_ids_list, start=1):
        logger.info("[batch %d/%d] %s", index, total, pdb_id.upper())
        try:
            aggregated.extend(
                annotate_pdb(
                    pdb_id,
                    client=client,
                    raddb_dataset=resolved_raddb,
                    no_raddb=no_raddb,
                    rfam_dataset=resolved_rfam,
                    no_rfam=no_rfam,
                    no_fr3d=no_fr3d,
                    **kwargs,
                )
            )
        except RibosomeAnnotatorError as exc:
            if not continue_on_error:
                raise
            logger.error("annotation failed for %s: %s", pdb_id, exc)
            aggregated.append(
                RibosomeAnnotation(
                    pdb_id=pdb_id.upper(),
                    assembly_id=None,
                    status="failed",
                    skip_reason=f"unhandled_error: {exc}",
                )
            )
    return aggregated


# ---------------------------------------------------------------------------
# Per-assembly orchestration
# ---------------------------------------------------------------------------


WARNING_MULTI_RIBOSOME_SPLIT = "multi_ribosome_bundle_split"
"""Emitted on each sub-ribosome annotation when an assembly was split
into multiple :class:`RibosomeAnnotation` results because it packed
multiple complete SSU+LSU pairs (e.g. 8R3V, 9O3L)."""


def _annotate_one_assembly(
    assembly: AssemblyContext,
    *,
    cache: Cache | None,
    cutoff: float,
    strict_complete_check: bool,
    completeness_thresholds: CompletenessThresholds | None,
    coordinate_source: CoordinateSource,
    local_coordinate_path: Path | None,
    client: httpx.Client | None,
    raddb_dataset: RADdbDataset | None = None,
    no_fr3d: bool = False,
) -> list[RibosomeAnnotation]:
    """Annotate one biological assembly, returning one or more annotations.

    Most assemblies return a single-element list. Multi-ribosome bundles
    (e.g. Mycobacterium 70S dimers like 8R3V) are split into one
    annotation per ribosome, each carrying a suffixed ``assembly_id``
    such as ``"1-1"`` / ``"1-2"`` so downstream CSV consumers see one
    row per ribosome.
    """
    pdb_id = assembly.pdb_id
    aid = assembly.assembly_id

    logger.info("classifying %s assembly %s", pdb_id, aid)
    classification_result = classify_assembly(
        assembly,
        thresholds=completeness_thresholds,
        strict_complete_check=strict_complete_check,
    )
    if classification_result.is_skip:
        return [_build_skip_annotation(assembly, classification_result)]

    by_role = partition_rna_chains_by_role(assembly.rna_chains)

    warnings: list[str] = list(classification_result.warnings)
    n_ribosomes = detect_multi_ribosome(by_role)
    n_ssu = len(by_role.get("ssu_main_rrna", []))
    n_lsu = len(by_role.get("lsu_main_rrna", []))
    if n_ribosomes == 1 and (n_ssu >= 2 or n_lsu >= 2) and n_ssu != n_lsu:
        # Fragmented rRNA: asymmetric SSU/LSU chain counts (e.g. one SSU
        # chain with the 28S split into LSUa/LSUb/SR2 fragments). The
        # canonical BGSU anchors don't transfer onto fragmented chains
        # because the anchor residue numbers are relative to a single
        # reference rRNA molecule, so contact-transfer can't resolve A/P/E
        # sites. Skip with a clear reason; batch callers continue to the
        # next entry.
        logger.info(
            "fragmented ribosome detected for %s assembly %s (SSU chains=%d, LSU chains=%d); skipping",
            pdb_id,
            aid,
            n_ssu,
            n_lsu,
        )
        return [
            RibosomeAnnotation(
                pdb_id=pdb_id,
                assembly_id=aid,
                status="skipped",
                skip_reason=C.SKIP_FRAGMENTED_RIBOSOME,
                ribosome_classification=classification_result.classification,
                topology=classification_result.topology,
                ssu_main_rrna_chains=by_role.get("ssu_main_rrna", []),
                lsu_main_rrna_chains=by_role.get("lsu_main_rrna", []),
                lsu_associated_rrna_chains=by_role.get("lsu_associated_rrna", []),
                classification_evidence=classification_result.evidence,
                warnings=warnings,
            )
        ]

    # Reference set selection per §6.4
    assert classification_result.classification is not None
    reference_units = _select_reference_units(classification_result.classification)

    assembly_chains_set = {
        chain.auth_asym_id for chain in (*assembly.rna_chains, *assembly.protein_chains)
    }

    # Same-organism fallback: when BGSU returns only NR equivalence-class
    # representatives (and not the target PDB itself), we substitute the
    # reference's SSU/LSU chain segments with the target's SSU/LSU chain
    # IDs. This works because same-organism deposits share residue
    # numbering convention (the dominant case for E. coli → another E. coli).
    chain_substitution = _build_chain_substitution(reference_units, by_role)

    # Per-subunit batching: concatenate all SSU site anchors into one BGSU
    # call and all LSU site anchors into another. The biological premise
    # (these anchors are textbook conserved tRNA-interacting nucleotides:
    # decoding-centre A1492/A1493, PTC G2252/G2253, etc.) is why the
    # intersection-semantic risk of larger batches is acceptable here.
    correspondence_by_site: dict[str, CorrespondenceResult] = {}
    for subunit in ("ssu", "lsu"):
        subunit_groups = {
            key: list(units)
            for key, units in reference_units.items()
            if key.startswith(f"{subunit}_") and units
        }
        if not subunit_groups:
            continue
        try:
            subunit_results = _get_or_fetch_subunit_correspondence(
                subunit_groups,
                target_pdb_id=pdb_id,
                assembly_chains=assembly_chains_set,
                chain_substitution=chain_substitution,
                cache=cache,
                client=client,
            )
        except (ApiRequestError, CorrespondenceMappingError) as exc:
            logger.warning(
                "correspondence fetch failed for %s subunit %s: %s",
                pdb_id,
                subunit,
                exc,
            )
            for site_key in subunit_groups:
                warnings.append(f"correspondence_fetch_failed_for_{site_key}")
            continue
        for site_key, result in subunit_results.items():
            correspondence_by_site[site_key] = result
            warnings.extend(result.warnings)

    try:
        logger.info("loading coordinates for %s assembly %s", pdb_id, aid)
        structure = get_assembly_structure(
            pdb_id,
            aid,
            source=coordinate_source,
            local_path=local_coordinate_path,
            cache=cache,
            client=client,
        )
    except (CoordinateDownloadError, CoordinateParsingError) as exc:
        logger.error("coordinate failure for %s assembly %s: %s", pdb_id, aid, exc)
        return [
            RibosomeAnnotation(
                pdb_id=pdb_id,
                assembly_id=aid,
                status="failed",
                skip_reason=f"coordinate_failure: {exc}",
                ribosome_classification=classification_result.classification,
                topology=classification_result.topology,
                ssu_main_rrna_chains=by_role.get("ssu_main_rrna", []),
                lsu_main_rrna_chains=by_role.get("lsu_main_rrna", []),
                lsu_associated_rrna_chains=by_role.get("lsu_associated_rrna", []),
                classification_evidence=classification_result.evidence,
                warnings=warnings,
            )
        ]

    # Multi-ribosome bundle: split into per-ribosome sub-contexts and
    # emit one annotation per ribosome.
    if n_ribosomes >= 2:
        return _annotate_multi_ribosome_bundle(
            assembly=assembly,
            structure=structure,
            by_role=by_role,
            reference_units=reference_units,
            correspondence_by_site=correspondence_by_site,
            classification_result=classification_result,
            warnings=warnings,
            cache=cache,
            client=client,
            cutoff=cutoff,
            raddb_dataset=raddb_dataset,
            no_fr3d=no_fr3d,
        )

    annotation = _run_assignment_for_assembly(
        assembly=assembly,
        structure=structure,
        by_role=by_role,
        correspondence_by_site=correspondence_by_site,
        classification_result=classification_result,
        warnings=warnings,
        cutoff=cutoff,
        cache=cache,
        client=client,
        raddb_dataset=raddb_dataset,
        no_fr3d=no_fr3d,
    )
    return [annotation]


def _run_assignment_for_assembly(
    *,
    assembly: AssemblyContext,
    structure: gemmi.Structure,
    by_role: dict[str, list[ChainRef]],
    correspondence_by_site: dict[str, CorrespondenceResult],
    classification_result: ClassificationResult,
    warnings: list[str],
    cutoff: float,
    cache: Cache | None,
    client: httpx.Client | None,
    raddb_dataset: RADdbDataset | None,
    no_fr3d: bool,
) -> RibosomeAnnotation:
    """Run the assignment + state + FR3D pipeline for one (sub-)assembly."""
    pdb_id = assembly.pdb_id
    aid = assembly.assembly_id
    site_warnings = list(warnings)

    logger.info("assigning functional chains for %s assembly %s", pdb_id, aid)
    assignments = assign_functional_chains(
        structure,
        assembly,
        by_role,
        correspondence_by_site,
        cutoff=cutoff,
        topology=classification_result.topology,
    )
    site_warnings.extend(assignments.warnings)

    logger.info("computing tRNA states for %s assembly %s", pdb_id, aid)
    states = compute_trna_states(
        structure,
        assembly,
        assignments,
        correspondence_by_site,
        cutoff=cutoff,
        topology=classification_result.topology,
    )
    site_warnings.extend(states.warnings)

    # Safeguard: a state of ``**/**`` means the chain is a too-short
    # fragment with no anchor contact on either subunit. The assignment
    # passes already require an anchor contact within cutoff, so this
    # should not appear from the normal pipeline — but defensive against
    # pathological inputs we demote the chain to unmapped_rna_chains.
    assignments, states = _demote_no_contact_fragments(assignments, states)

    annotation = _build_annotated_annotation(
        assembly=assembly,
        by_role=by_role,
        classification_result=classification_result,
        assignments=assignments,
        states=states,
        warnings=site_warnings,
        raddb_dataset=raddb_dataset,
    )

    # tRNA-mRNA codon/anticodon evidence. Best-effort — every failure
    # path returns an empty list, never raises, so the annotation
    # pipeline cannot be broken by an FR3D outage.
    if not no_fr3d:
        try:
            annotation.trna_mrna_interactions = extract_trna_mrna_interactions(
                annotation, structure, cache=cache, client=client
            )
        except Exception as exc:  # defensive — module already swallows known errors
            logger.warning(
                "tRNA-mRNA extraction failed for %s assembly %s: %s",
                pdb_id,
                aid,
                exc,
            )
            annotation.trna_mrna_interactions = []

    return annotation


def _annotate_multi_ribosome_bundle(
    *,
    assembly: AssemblyContext,
    structure: gemmi.Structure,
    by_role: dict[str, list[ChainRef]],
    reference_units: Mapping[str, tuple[str, ...]],
    correspondence_by_site: dict[str, CorrespondenceResult],
    classification_result: ClassificationResult,
    warnings: list[str],
    cache: Cache | None,
    client: httpx.Client | None,
    cutoff: float,
    raddb_dataset: RADdbDataset | None,
    no_fr3d: bool,
) -> list[RibosomeAnnotation]:
    """Split a multi-ribosome assembly and emit one annotation per ribosome.

    Pairing rule: greedy nearest-centroid SSU↔LSU. Each non-rRNA RNA
    chain is partitioned to the ribosome whose combined SSU+LSU centroid
    it is closest to. Non-ribosomal proteins are partitioned the same
    way (each ribosome carries its own bound factors).

    Per-ribosome correspondence resolution: BGSU's NR set typically maps
    each canonical anchor to **one** chain per subunit (e.g. for 8R3V
    only chains ``A1`` and ``72`` are mapped, even though the deposit
    has ``A1+71`` and ``A2+72`` as its two ribosomes). To rescue the
    other ribosome's tRNAs the bundle handler rebuilds a per-ribosome
    :class:`CorrespondenceResult` using the same-organism chain
    substitution fallback (§5.2.2), substituting the BGSU
    representative chain with this ribosome's SSU/LSU author IDs.
    """
    pairs = pair_ssu_lsu_by_centroid(
        structure,
        by_role.get("ssu_main_rrna", []),
        by_role.get("lsu_main_rrna", []),
    )
    if len(pairs) < 2:
        # Pairing failed (centroid lookup misses) — fall back to single-
        # ribosome path so the assembly still gets some annotation.
        logger.warning(
            "multi-ribosome split for %s assembly %s degenerated to %d pair(s); "
            "falling back to single-ribosome path",
            assembly.pdb_id,
            assembly.assembly_id,
            len(pairs),
        )
        single = _run_assignment_for_assembly(
            assembly=assembly,
            structure=structure,
            by_role=by_role,
            correspondence_by_site=correspondence_by_site,
            classification_result=classification_result,
            warnings=warnings,
            cutoff=cutoff,
            cache=cache,
            client=client,
            raddb_dataset=raddb_dataset,
            no_fr3d=no_fr3d,
        )
        return [single]

    # Partition non-main-rRNA RNA chains and protein chains by geometric
    # proximity to each ribosome.
    non_main_rna = [
        c for c in assembly.rna_chains
        if c not in by_role.get("ssu_main_rrna", []) and c not in by_role.get("lsu_main_rrna", [])
    ]
    rna_groups = partition_chains_by_ribosome(structure, pairs, non_main_rna)
    protein_groups = partition_chains_by_ribosome(structure, pairs, assembly.protein_chains)

    # Identify the reference SSU/LSU chain segments once — these are the
    # BGSU "representative" chains we'll substitute per ribosome.
    ssu_ref_chain, lsu_ref_chain = _reference_subunit_chains(reference_units)

    annotations: list[RibosomeAnnotation] = []
    for index, ((ssu, lsu), rna_chain_ids, protein_chain_ids) in enumerate(
        zip(pairs, rna_groups, protein_groups, strict=True), start=1
    ):
        sub_warnings = list(warnings)
        sub_warnings.append(WARNING_MULTI_RIBOSOME_SPLIT)

        # Carve out the per-ribosome chain set: this ribosome's SSU + LSU +
        # geometrically-proximate other RNA chains + per-ribosome proteins.
        sub_rna_chains = [ssu, lsu] + [c for c in non_main_rna if c.auth_asym_id in rna_chain_ids]
        sub_protein_chains = [c for c in assembly.protein_chains if c.auth_asym_id in protein_chain_ids]

        sub_by_role = partition_rna_chains_by_role(sub_rna_chains)
        sub_assembly = AssemblyContext(
            pdb_id=assembly.pdb_id,
            assembly_id=f"{assembly.assembly_id}-{index}",
            experimental_methods=list(assembly.experimental_methods),
            rna_chains=sub_rna_chains,
            protein_chains=sub_protein_chains,
            ligands=list(assembly.ligands),
            coordinate_path=assembly.coordinate_path,
        )

        # Per-ribosome correspondence: re-run the BGSU fetch (cache-resident)
        # with this ribosome's SSU/LSU as substitution targets so the
        # §5.2.2 fallback covers the anchors BGSU left pointing at the
        # *other* ribosome's chains.
        sub_correspondence = _rebuild_correspondence_for_ribosome(
            reference_units=reference_units,
            target_pdb_id=assembly.pdb_id,
            assembly_chains={c.auth_asym_id for c in (*sub_rna_chains, *sub_protein_chains)},
            ssu_ref_chain=ssu_ref_chain,
            lsu_ref_chain=lsu_ref_chain,
            ssu_target_chain=ssu.auth_asym_id,
            lsu_target_chain=lsu.auth_asym_id,
            cache=cache,
            client=client,
        )

        annotation = _run_assignment_for_assembly(
            assembly=sub_assembly,
            structure=structure,
            by_role=sub_by_role,
            correspondence_by_site=sub_correspondence,
            classification_result=classification_result,
            warnings=sub_warnings,
            cutoff=cutoff,
            cache=cache,
            client=client,
            raddb_dataset=raddb_dataset,
            no_fr3d=no_fr3d,
        )
        annotations.append(annotation)

    return annotations


_NO_CONTACT_FRAGMENT_STATE = "**/**"


def _demote_no_contact_fragments(
    assignments: ChainAssignments,
    states: TRNAStates,
) -> tuple[ChainAssignments, TRNAStates]:
    """Clear A or P chain assignments whose state is ``**/**``.

    A state of ``**/**`` means the chain is shorter than
    :data:`constants.ASL_FRAGMENT_MAX_LENGTH` and makes no anchor contact
    on either subunit. Such a chain shouldn't be claiming a tRNA role —
    let it fall through to ``unmapped_rna_chains`` instead. E-tRNA is
    not demoted here: its assignment already requires a ≤cutoff contact
    at the LSU E-site (or the SSU E-site for isolated_ssu), so a state
    like ``**/E`` represents a real partial-tRNA bound at the E-site
    (e.g. 7Q0R chain 10, 14 modeled residues out of a 76-nt entity).
    """
    chain_updates: dict[str, Any] = {}
    state_updates: dict[str, Any] = {}
    if states.aminoacyl_trna_state == _NO_CONTACT_FRAGMENT_STATE:
        chain_updates["aminoacyl_trna_chain"] = None
        state_updates["aminoacyl_trna_state"] = None
    if states.peptidyl_trna_state == _NO_CONTACT_FRAGMENT_STATE:
        chain_updates["peptidyl_trna_chain"] = None
        state_updates["peptidyl_trna_state"] = None
    if chain_updates:
        assignments = assignments.model_copy(update=chain_updates)
    if state_updates:
        states = states.model_copy(update=state_updates)
    return assignments, states


def _reference_subunit_chains(
    reference_units: Mapping[str, tuple[str, ...]],
) -> tuple[str | None, str | None]:
    """Return ``(ssu_ref_chain, lsu_ref_chain)`` parsed from the first SSU and
    LSU reference units. Both may be ``None`` if the corresponding set is
    empty (shouldn't happen for known classifications)."""
    from ribosome_state_annotator.correspondence import parse_unit_id

    ssu_ref: str | None = None
    lsu_ref: str | None = None
    for site_key, units in reference_units.items():
        if not units:
            continue
        try:
            parsed = parse_unit_id(units[0])
        except ValueError:
            continue
        if site_key.startswith("ssu") and ssu_ref is None:
            ssu_ref = parsed.chain
        elif site_key.startswith("lsu") and lsu_ref is None:
            lsu_ref = parsed.chain
    return ssu_ref, lsu_ref


def _rebuild_correspondence_for_ribosome(
    *,
    reference_units: Mapping[str, tuple[str, ...]],
    target_pdb_id: str,
    assembly_chains: set[str],
    ssu_ref_chain: str | None,
    lsu_ref_chain: str | None,
    ssu_target_chain: str,
    lsu_target_chain: str,
    cache: Cache | None,
    client: httpx.Client | None,
) -> dict[str, CorrespondenceResult]:
    """Rebuild per-site :class:`CorrespondenceResult` for one ribosome.

    Re-uses :func:`_get_or_fetch_subunit_correspondence` with this
    ribosome's SSU/LSU as the substitution targets. The BGSU fetch is
    cache-resident so this adds no network cost — the rebuild
    re-applies the §5.2.2 filter + chain-substitution fallback against
    the per-ribosome ``assembly_chains`` and ``chain_substitution``,
    which is what produces the correct anchor residues even when
    BGSU's NR response only mapped one of the two SSU (or LSU) chains.
    """
    chain_substitution: dict[str, str] = {}
    if ssu_ref_chain is not None:
        chain_substitution[ssu_ref_chain] = ssu_target_chain
    if lsu_ref_chain is not None:
        chain_substitution[lsu_ref_chain] = lsu_target_chain

    out: dict[str, CorrespondenceResult] = {}
    for subunit in ("ssu", "lsu"):
        subunit_groups = {
            key: list(units)
            for key, units in reference_units.items()
            if key.startswith(f"{subunit}_") and units
        }
        if not subunit_groups:
            continue
        try:
            subunit_results = _get_or_fetch_subunit_correspondence(
                subunit_groups,
                target_pdb_id=target_pdb_id,
                assembly_chains=assembly_chains,
                chain_substitution=chain_substitution,
                cache=cache,
                client=client,
            )
        except (ApiRequestError, CorrespondenceMappingError) as exc:
            # Shouldn't happen on the warm cache path that fed the
            # bundle-wide call, but degrade gracefully if it does.
            logger.warning(
                "per-ribosome correspondence rebuild failed for %s subunit %s: %s",
                target_pdb_id,
                subunit,
                exc,
            )
            for site_key in subunit_groups:
                out[site_key] = CorrespondenceResult(reference_key=site_key)
            continue
        out.update(subunit_results)
    return out


# ---------------------------------------------------------------------------
# Annotation builders
# ---------------------------------------------------------------------------


def _build_skip_annotation(
    assembly: AssemblyContext, classification_result: ClassificationResult
) -> RibosomeAnnotation:
    return RibosomeAnnotation(
        pdb_id=assembly.pdb_id,
        assembly_id=assembly.assembly_id,
        status="skipped",
        skip_reason=classification_result.skip_reason,
        topology=classification_result.topology,
        classification_evidence=classification_result.evidence,
        warnings=classification_result.warnings,
    )


def _build_annotated_annotation(
    *,
    assembly: AssemblyContext,
    by_role: Mapping[str, list[ChainRef]],
    classification_result: ClassificationResult,
    assignments: ChainAssignments,
    states: TRNAStates,
    warnings: list[str],
    raddb_dataset: RADdbDataset | None = None,
) -> RibosomeAnnotation:
    # Compute unmapped_rna_chains: RNA chains that weren't placed in any
    # OUTPUT bucket (rRNA roles + the four functional-chain slots).
    #
    # The "trna" partition bucket is intentionally NOT in this exclusion
    # list — it's an intermediate internal classification used by the
    # mRNA-pool exclusion in assign_functional_chains. A chain in the
    # trna bucket that doesn't get assigned to A/P/E by contact transfer
    # (e.g. mt-tRNA-Val in 3J9M, which sits at the central protuberance
    # as a structural 5S surrogate and so doesn't contact any canonical
    # tRNA-binding-site anchor) must still surface somewhere in the
    # output, otherwise it disappears from view entirely.
    placed_ifes: set[str] = set()
    for role in ("ssu_main_rrna", "lsu_main_rrna", "lsu_associated_rrna"):
        placed_ifes |= {chain.ife for chain in by_role.get(role, [])}
    for chain in (
        assignments.mrna_chain,
        assignments.aminoacyl_trna_chain,
        assignments.peptidyl_trna_chain,
        assignments.exit_trna_chain,
    ):
        if chain is not None:
            placed_ifes.add(chain.ife)
    unmapped_rna_chains = [chain for chain in assembly.rna_chains if chain.ife not in placed_ifes]

    non_ribosomal_proteins = [
        chain for chain in assembly.protein_chains if not chain.is_ribosomal_protein
    ]

    bound_ligands = [
        ligand for ligand in assembly.ligands if ligand.comp_id not in C.DEFAULT_LIGAND_EXCLUSIONS
    ]
    # Dedup by comp_id (the parser already dedups but we re-assert here so
    # this function is robust when used outside the standard pipeline).
    bound_ligands = _dedupe_by_comp_id(bound_ligands)

    large_scale_movements = _build_large_scale_movements(
        pdb_id=assembly.pdb_id,
        ssu_chains=list(by_role.get("ssu_main_rrna", [])),
        lsu_chains=list(by_role.get("lsu_main_rrna", [])),
        raddb_dataset=raddb_dataset,
        warnings=warnings,
    )

    # Per-assembly NCBI taxonomy (§34). Vote with rRNA chains only; if
    # none of them carries a lineage tag emit the dedicated warning.
    voting_rrna = (
        list(by_role.get("ssu_main_rrna", []))
        + list(by_role.get("lsu_main_rrna", []))
        + list(by_role.get("lsu_associated_rrna", []))
    )
    all_chains = list(assembly.rna_chains) + list(assembly.protein_chains)
    assembly_taxonomy = aggregate_assembly_lineage(voting_rrna, all_chains=all_chains)
    if assembly_taxonomy is None and voting_rrna:
        warnings.append("no_source_organism_taxonomy")

    return RibosomeAnnotation(
        pdb_id=assembly.pdb_id,
        assembly_id=assembly.assembly_id,
        status="annotated",
        ribosome_classification=classification_result.classification,
        topology=classification_result.topology,
        ssu_main_rrna_chains=list(by_role.get("ssu_main_rrna", [])),
        lsu_main_rrna_chains=list(by_role.get("lsu_main_rrna", [])),
        lsu_associated_rrna_chains=list(by_role.get("lsu_associated_rrna", [])),
        unmapped_rna_chains=unmapped_rna_chains,
        mrna_chain=assignments.mrna_chain,
        aminoacyl_trna_chain=assignments.aminoacyl_trna_chain,
        peptidyl_trna_chain=assignments.peptidyl_trna_chain,
        exit_trna_chain=assignments.exit_trna_chain,
        aminoacyl_trna_state=states.aminoacyl_trna_state,
        peptidyl_trna_state=states.peptidyl_trna_state,
        exit_trna_state=states.exit_trna_state,
        non_ribosomal_proteins=non_ribosomal_proteins,
        bound_ligands=bound_ligands,
        large_scale_movements=large_scale_movements,
        assembly_taxonomy=assembly_taxonomy,
        classification_evidence={
            **classification_result.evidence,
            **states.trna_state_evidence,
        },
        warnings=warnings,
    )


def _load_raddb_safely(
    *, client: httpx.Client | None, force_refresh: bool
) -> RADdbDataset | None:
    """Best-effort RADdb load. Returns ``None`` on any failure.

    The annotation pipeline must never crash because RADdb is missing or
    unreachable — when this returns ``None``, downstream code emits
    ``large_scale_movements`` with null metrics and ``rad_date=None``.
    """
    try:
        metadata = ensure_raddb_available(client=client, force_refresh=force_refresh)
    except Exception as exc:
        logger.warning("RADdb refresh check failed: %s", exc)
        return None
    if metadata is None:
        return None
    try:
        return load_raddb_dataset(metadata=metadata)
    except Exception as exc:
        logger.warning("RADdb dataset load failed: %s", exc)
        return None


def _load_rfam_safely(
    *, client: httpx.Client | None, force_refresh: bool
) -> RfamPdbRegionDataset | None:
    """Best-effort Rfam pdb_full_region load. Returns ``None`` on any failure.

    When this returns ``None``, the pipeline falls back to whatever Rfam
    tags RCSB supplies directly on the entry (which is typically nothing
    for rRNA chains, but the path still runs without error).
    """
    try:
        metadata = ensure_rfam_pdb_region_available(
            client=client, force_refresh=force_refresh
        )
    except Exception as exc:
        logger.warning("Rfam pdb_full_region refresh check failed: %s", exc)
        return None
    if metadata is None:
        return None
    try:
        return load_rfam_pdb_region_dataset(metadata=metadata)
    except Exception as exc:
        logger.warning("Rfam pdb_full_region dataset load failed: %s", exc)
        return None


def _apply_rfam_pdb_region(
    chains: list[ChainRef], rfam_by_chain: dict[str, list[str]]
) -> None:
    """Replace each chain's ``rfam_accessions`` with the single best-score
    Rfam from the EBI ``pdb_full_region`` file (when present).

    Chains with no entry in the file keep their existing
    ``rfam_accessions`` (which may be empty or RCSB-supplied). The file
    is authoritative when it has a hit — single-best-score selection
    eliminates the multi-family over-annotation pattern (RF00177 +
    RF01959 + RF01960 on the same chain) that PDBe's REST endpoint
    surfaces.
    """
    for chain in chains:
        file_rfam = rfam_by_chain.get(chain.auth_asym_id)
        if file_rfam:
            chain.rfam_accessions = list(file_rfam)


def _build_large_scale_movements(
    *,
    pdb_id: str,
    ssu_chains: list[ChainRef],
    lsu_chains: list[ChainRef],
    raddb_dataset: RADdbDataset | None,
    warnings: list[str],
) -> LargeScaleMovements:
    """Resolve RADdb metrics for one assembly and return the JSON-shape block.

    Multi-chain assemblies (>1 SSU or >1 LSU main rRNA chain): the
    cartesian product of LSU x SSU is tried and metrics are returned only
    when exactly one pair matches a RADdb row. Otherwise a warning is
    recorded and the metrics fields are null.
    """
    rad_date = raddb_dataset.metadata.rad_date if raddb_dataset is not None else None
    if raddb_dataset is None:
        return LargeScaleMovements(rad_date=None)

    matches: list[dict[str, Any]] = []
    for lsu in lsu_chains:
        for ssu in ssu_chains:
            metrics = get_motion_metrics(
                raddb_dataset, pdb_id, lsu.auth_asym_id, ssu.auth_asym_id
            )
            if metrics is not None:
                matches.append(metrics)

    if len(matches) == 1:
        return LargeScaleMovements(
            rad_date=rad_date,
            intersubunit_rotation=matches[0]["intersubunit_rotation"],
            ssu_head_rotation=matches[0]["ssu_head_rotation"],
        )
    if len(matches) > 1:
        warnings.append("raddb_ambiguous_chain_pair")
    return LargeScaleMovements(rad_date=rad_date)


def _dedupe_by_comp_id(ligands: list[LigandRef]) -> list[LigandRef]:
    seen: set[str] = set()
    out: list[LigandRef] = []
    for ligand in ligands:
        if ligand.comp_id in seen:
            continue
        seen.add(ligand.comp_id)
        out.append(ligand)
    return out


# ---------------------------------------------------------------------------
# Cache wrappers
# ---------------------------------------------------------------------------


def _resolve_cache(cache: Cache | None, cache_dir: Path | None, no_cache: bool) -> Cache | None:
    """Resolve the three cache-related parameters into a single ``Cache | None``."""
    if no_cache:
        return None
    if cache is not None:
        return cache
    if cache_dir is not None:
        return Cache(cache_dir)
    # Default: use the standard cache dir (spec §17).
    return Cache()


def _build_chain_substitution(
    reference_units: Mapping[str, tuple[str, ...]],
    by_role: Mapping[str, list[ChainRef]],
) -> dict[str, str]:
    """Build a ``{reference_chain: target_chain}`` map for the §5.2.2 fallback.

    See :func:`correspondence.build_correspondence_result` for the rationale.
    The map carries SSU and LSU rRNA chain substitutions inferred from the
    reference unit IDs (which tell us which chain in the reference holds
    SSU vs LSU anchors) and the target's Rfam-role partition.
    """
    from ribosome_state_annotator.correspondence import parse_unit_id

    ssu_ref_chain: str | None = None
    lsu_ref_chain: str | None = None
    for site_key, units in reference_units.items():
        if not units:
            continue
        try:
            parsed = parse_unit_id(units[0])
        except ValueError:
            continue
        if site_key.startswith("ssu") and ssu_ref_chain is None:
            ssu_ref_chain = parsed.chain
        elif site_key.startswith("lsu") and lsu_ref_chain is None:
            lsu_ref_chain = parsed.chain

    substitution: dict[str, str] = {}
    ssu_main = by_role.get("ssu_main_rrna") or []
    lsu_main = by_role.get("lsu_main_rrna") or []
    if ssu_ref_chain is not None and len(ssu_main) == 1:
        substitution[ssu_ref_chain] = ssu_main[0].auth_asym_id
    if lsu_ref_chain is not None and len(lsu_main) == 1:
        substitution[lsu_ref_chain] = lsu_main[0].auth_asym_id
    return substitution


def _get_or_fetch_entry(
    pdb_id: str, cache: Cache | None, client: httpx.Client | None
) -> dict[str, Any]:
    """RCSB entry fetch with cache wrap."""
    if cache is not None:
        cached = cache.get_rcsb_payload(pdb_id)
        if cached is not None:
            logger.debug("rcsb cache hit for %s", pdb_id)
            return cached
    logger.info("fetching RCSB entry payload for %s", pdb_id)
    payload = fetch_entry_payload(pdb_id, client=client)
    if cache is not None:
        cache.put_rcsb_payload(pdb_id, payload)
    return payload


def _bgsu_cache_key(units: list[str], scope: str, resolution: str, depth: int) -> str:
    """Deterministic key for BGSU responses. Independent of URL encoding.

    ``depth`` is part of the key because changing it materially changes
    the response (the number of equivalence-class members BGSU returns).
    Cache entries from a smaller depth would silently mislead a
    deeper-depth request.
    """
    return f"units={','.join(units)}|scope={scope}|resolution={resolution}|depth={depth}"


def _get_or_fetch_subunit_correspondence(
    site_groups: Mapping[str, list[str]],
    *,
    target_pdb_id: str,
    assembly_chains: set[str],
    cache: Cache | None,
    client: httpx.Client | None,
    chain_substitution: dict[str, str] | None = None,
    scope: str = DEFAULT_BGSU_SCOPE,
    resolution: str = DEFAULT_BGSU_RESOLUTION,
    depth: int = DEFAULT_BGSU_DEPTH,
) -> dict[str, CorrespondenceResult]:
    """Batched BGSU fetch for all sites in one subunit.

    Concatenates every site's anchor units into a single BGSU
    ``map_across_chains`` query, retrieves the response (cached on the
    full union of unit IDs), and slices the parsed alignment back into
    per-site :class:`CorrespondenceResult` objects so downstream
    consumers see exactly the same per-site shape as the
    unconsolidated path.

    Two BGSU calls per cold-cache organism (SSU + LSU) instead of
    seven. The biological premise — the anchor units are textbook
    conserved tRNA-interacting nucleotides — is what makes the
    larger-batch intersection-semantic risk acceptable. Anchors that
    were empirically weakly-conserved have already been curated out
    (e.g. yeast 25S position 2454; see REFERENCES.md and the spec
    v3.1 addendum).

    Returns a ``{site_key: CorrespondenceResult}`` dict with one entry
    per non-empty key in ``site_groups``. An empty ``site_groups``
    returns an empty dict.
    """
    site_unit_ranges: dict[str, tuple[int, int]] = {}
    all_units: list[str] = []
    for site_key, units in site_groups.items():
        if not units:
            continue
        start = len(all_units)
        all_units.extend(units)
        site_unit_ranges[site_key] = (start, len(all_units))

    if not all_units:
        return {
            site_key: CorrespondenceResult(reference_key=site_key)
            for site_key in site_groups
        }

    cache_key = _bgsu_cache_key(all_units, scope, resolution, depth)
    raw: dict[str, list[str]] | None = None
    if cache is not None:
        cached_payload = cache.get_bgsu_payload(cache_key)
        if cached_payload is not None:
            logger.debug("bgsu cache hit for subunit batch (%d units)", len(all_units))
            raw = _ensure_alignment_dict(cached_payload)

    if raw is None:
        logger.info(
            "fetching BGSU correspondence for subunit batch (%d units, %d sites)",
            len(all_units),
            len(site_unit_ranges),
        )
        raw = fetch_correspondence(
            all_units, scope=scope, resolution=resolution, depth=depth, client=client
        )
        if cache is not None:
            alignment_payload = {
                "alignment": [
                    {"reference_unit": ref, "mapped_units": list(mapped)}
                    for ref, mapped in raw.items()
                ]
            }
            cache.put_bgsu_payload(cache_key, alignment_payload)

    # Slice the union response back into per-site CorrespondenceResults.
    # The per-site §5.2.2 filter, chain-substitution fallback, and
    # missing-anchor warnings are preserved unchanged by routing each
    # site's slice through build_correspondence_result.
    results: dict[str, CorrespondenceResult] = {}
    for site_key, units in site_groups.items():
        if not units:
            results[site_key] = CorrespondenceResult(reference_key=site_key)
            continue
        site_units_list = list(units)
        site_raw = {unit: list(raw.get(unit, [])) for unit in site_units_list}
        results[site_key] = build_correspondence_result(
            site_key,
            site_units_list,
            site_raw,
            target_pdb_id=target_pdb_id,
            assembly_chains=assembly_chains,
            chain_substitution=chain_substitution,
        )
    return results


def _get_or_fetch_correspondence(
    reference_key: str,
    units: list[str],
    *,
    target_pdb_id: str,
    assembly_chains: set[str],
    cache: Cache | None,
    client: httpx.Client | None,
    chain_substitution: dict[str, str] | None = None,
    scope: str = DEFAULT_BGSU_SCOPE,
    resolution: str = DEFAULT_BGSU_RESOLUTION,
    depth: int = DEFAULT_BGSU_DEPTH,
) -> CorrespondenceResult:
    """BGSU correspondence fetch + filter with cache wrap.

    The cache stores the RAW alignment dict (post-parse, pre-filter)
    because the §5.2.2 PDB-prefix + assembly-chain filter depends on
    per-assembly state and would invalidate cache hits across assemblies
    if applied before storage.
    """
    if not units:
        return CorrespondenceResult(reference_key=reference_key)
    cache_key = _bgsu_cache_key(units, scope, resolution, depth)

    raw: dict[str, list[str]] | None = None
    if cache is not None:
        cached_payload = cache.get_bgsu_payload(cache_key)
        if cached_payload is not None:
            logger.debug("bgsu cache hit for %s", reference_key)
            # The cached payload is the raw {ref: mapped_units} dict. Parse-
            # alignment-response is a no-op here, but we re-validate to
            # tolerate any future cache-format additions.
            raw = _ensure_alignment_dict(cached_payload)

    if raw is None:
        logger.info(
            "fetching BGSU correspondence for site=%s (%d units)",
            reference_key,
            len(units),
        )
        raw = fetch_correspondence(
            units, scope=scope, resolution=resolution, depth=depth, client=client
        )
        if cache is not None:
            # Store as alignment-shaped JSON so future cache reads remain valid
            # under :func:`parse_alignment_response`.
            alignment_payload = {
                "alignment": [
                    {"reference_unit": ref, "mapped_units": list(mapped)}
                    for ref, mapped in raw.items()
                ]
            }
            cache.put_bgsu_payload(cache_key, alignment_payload)

    return build_correspondence_result(
        reference_key,
        units,
        raw,
        target_pdb_id=target_pdb_id,
        assembly_chains=assembly_chains,
        chain_substitution=chain_substitution,
    )


def _ensure_alignment_dict(payload: dict[str, Any]) -> dict[str, list[str]]:
    """Adapt a cached BGSU response back to the parsed-alignment shape.

    The cache stores BGSU responses in the canonical alignment shape
    (``{"alignment": [{"reference_unit": ..., "mapped_units": [...]}]}``),
    so we can call :func:`parse_alignment_response` to get the
    ``{ref_unit: [mapped, ...]}`` form back out. We tolerate the legacy
    flat-dict form too, in case a future cache version drops the
    wrapping.
    """
    if "alignment" in payload:
        return parse_alignment_response(payload)
    # Legacy / flat: assume {ref_unit: [mapped, ...]}.
    return {
        key: [val for val in values if isinstance(val, str)]
        for key, values in payload.items()
        if isinstance(values, list)
    }


# ---------------------------------------------------------------------------
# Reference selection (§6.4)
# ---------------------------------------------------------------------------


def _select_reference_units(
    classification: str,
) -> Mapping[str, tuple[str, ...]]:
    if classification == "bacterial_ribosome":
        return C.ECOLI_REFERENCE_UNITS
    if classification == "eukaryotic_organellar_ribosome":
        # Organellar rRNA cores are bacterial-like, so we reuse the E.
        # coli reference — but with the four mt-incompatible anchors
        # removed (see :data:`constants.ECOLI_REFERENCE_UNITS_ORGANELLAR`
        # for the rationale). BGSU's structural-alignment correspondence
        # cross-walks the remaining 16 anchors onto mt-12S / mt-16S
        # residues directly, regardless of the target deposit's
        # numbering scheme.
        return C.ECOLI_REFERENCE_UNITS_ORGANELLAR
    if classification == "eukaryotic_ribosome":
        return C.YEAST_REFERENCE_UNITS
    # Shouldn't reach here: the only callers feed in non-None
    # RibosomeClassification values. Defensive empty mapping prevents an
    # exception in the unlikely case.
    logger.warning("no reference units for classification %r", classification)
    return {}


__all__ = [
    "annotate_assembly",
    "annotate_many",
    "annotate_pdb",
]
