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
    LigandRef,
    RibosomeAnnotation,
)
from ribosome_state_annotator.pdbe_client import (
    fetch_rfam_mappings,
    parse_rfam_mappings,
)
from ribosome_state_annotator.rcsb_client import (
    fetch_entry_payload,
    parse_assemblies,
)

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

    # PDBe Rfam augmentation: RCSB no longer supplies Rfam annotations
    # for rRNA chains, so we cross-reference PDBe's
    # `/nucleic_mappings/rfam` endpoint and merge into each ChainRef
    # before classification runs.
    rfam_by_chain = _get_or_fetch_rfam_mappings(pdb_id_upper, resolved_cache, client)
    if rfam_by_chain:
        for assembly in assemblies:
            _augment_chain_rfam(assembly.rna_chains, rfam_by_chain)

    results: list[RibosomeAnnotation] = []
    for assembly in assemblies:
        results.append(
            _annotate_one_assembly(
                assembly,
                cache=resolved_cache,
                cutoff=contact_cutoff_angstrom,
                strict_complete_check=strict_complete_check,
                completeness_thresholds=completeness_thresholds,
                coordinate_source=coordinate_source,
                local_coordinate_path=local_coordinate_path,
                client=client,
            )
        )
    return results


def annotate_assembly(
    pdb_id: str,
    assembly_id: str,
    **kwargs: Any,
) -> RibosomeAnnotation:
    """Annotate one biological assembly. Convenience wrapper.

    Returns the single :class:`RibosomeAnnotation` for ``(pdb_id, assembly_id)``.
    """
    results = annotate_pdb(pdb_id, assembly_id=assembly_id, **kwargs)
    if not results:
        return RibosomeAnnotation(
            pdb_id=pdb_id.upper(),
            assembly_id=assembly_id,
            status="failed",
            skip_reason=f"assembly_not_found: {assembly_id}",
        )
    return results[0]


def annotate_many(
    pdb_ids: Iterable[str],
    *,
    continue_on_error: bool = True,
    **kwargs: Any,
) -> list[RibosomeAnnotation]:
    """Batch wrapper. Iterates ``pdb_ids`` and calls :func:`annotate_pdb` on each.

    With ``continue_on_error=True`` (the default), per-entry exceptions
    are caught and recorded as a ``status="failed"`` annotation rather
    than propagating. Pass ``continue_on_error=False`` to abort the
    batch on the first error instead.
    """
    aggregated: list[RibosomeAnnotation] = []
    pdb_ids_list = list(pdb_ids)
    total = len(pdb_ids_list)
    for index, pdb_id in enumerate(pdb_ids_list, start=1):
        logger.info("[batch %d/%d] %s", index, total, pdb_id.upper())
        try:
            aggregated.extend(annotate_pdb(pdb_id, **kwargs))
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
) -> RibosomeAnnotation:
    pdb_id = assembly.pdb_id
    aid = assembly.assembly_id

    logger.info("classifying %s assembly %s", pdb_id, aid)
    classification_result = classify_assembly(
        assembly,
        thresholds=completeness_thresholds,
        strict_complete_check=strict_complete_check,
    )
    if classification_result.is_skip:
        return _build_skip_annotation(assembly, classification_result)

    by_role = partition_rna_chains_by_role(assembly.rna_chains)

    warnings: list[str] = list(classification_result.warnings)
    if len(by_role.get("ssu_main_rrna", [])) > 1 or len(by_role.get("lsu_main_rrna", [])) > 1:
        warnings.append(WARNING_MULTIPLE_CHAINS_FOR_LEGACY_CSV)

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

    correspondence_by_site: dict[str, CorrespondenceResult] = {}
    for site_key, units in reference_units.items():
        try:
            result = _get_or_fetch_correspondence(
                site_key,
                list(units),
                target_pdb_id=pdb_id,
                assembly_chains=assembly_chains_set,
                chain_substitution=chain_substitution,
                cache=cache,
                client=client,
            )
        except (ApiRequestError, CorrespondenceMappingError) as exc:
            logger.warning(
                "correspondence fetch failed for %s site %s: %s",
                pdb_id,
                site_key,
                exc,
            )
            warnings.append(f"correspondence_fetch_failed_for_{site_key}")
            continue
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
        return RibosomeAnnotation(
            pdb_id=pdb_id,
            assembly_id=aid,
            status="failed",
            skip_reason=f"coordinate_failure: {exc}",
            ribosome_classification=classification_result.classification,
            ssu_main_rrna_chains=by_role.get("ssu_main_rrna", []),
            lsu_main_rrna_chains=by_role.get("lsu_main_rrna", []),
            lsu_associated_rrna_chains=by_role.get("lsu_associated_rrna", []),
            classification_evidence=classification_result.evidence,
            warnings=warnings,
        )

    logger.info("assigning functional chains for %s assembly %s", pdb_id, aid)
    assignments = assign_functional_chains(
        structure,
        assembly,
        by_role,
        correspondence_by_site,
        cutoff=cutoff,
    )
    warnings.extend(assignments.warnings)

    logger.info("computing tRNA states for %s assembly %s", pdb_id, aid)
    states = compute_trna_states(
        structure,
        assembly,
        assignments,
        correspondence_by_site,
        cutoff=cutoff,
    )
    warnings.extend(states.warnings)

    return _build_annotated_annotation(
        assembly=assembly,
        by_role=by_role,
        classification_result=classification_result,
        assignments=assignments,
        states=states,
        warnings=warnings,
    )


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
) -> RibosomeAnnotation:
    # Compute other_rna_chains: RNA chains that weren't placed in any role
    # AND weren't subsequently assigned as mRNA / A-tRNA / P-tRNA / E-tRNA.
    placed_ifes: set[str] = set()
    for role in ("ssu_main_rrna", "lsu_main_rrna", "lsu_associated_rrna", "trna"):
        placed_ifes |= {chain.ife for chain in by_role.get(role, [])}
    for chain in (
        assignments.mrna_chain,
        assignments.aminoacyl_trna_chain,
        assignments.peptidyl_trna_chain,
        assignments.exit_trna_chain,
    ):
        if chain is not None:
            placed_ifes.add(chain.ife)
    other_rna_chains = [chain for chain in assembly.rna_chains if chain.ife not in placed_ifes]

    non_ribosomal_proteins = [
        chain for chain in assembly.protein_chains if not chain.is_ribosomal_protein
    ]

    bound_ligands = [
        ligand for ligand in assembly.ligands if ligand.comp_id not in C.DEFAULT_LIGAND_EXCLUSIONS
    ]
    # Dedup by comp_id (the parser already dedups but we re-assert here so
    # this function is robust when used outside the standard pipeline).
    bound_ligands = _dedupe_by_comp_id(bound_ligands)

    return RibosomeAnnotation(
        pdb_id=assembly.pdb_id,
        assembly_id=assembly.assembly_id,
        status="annotated",
        ribosome_classification=classification_result.classification,
        ssu_main_rrna_chains=list(by_role.get("ssu_main_rrna", [])),
        lsu_main_rrna_chains=list(by_role.get("lsu_main_rrna", [])),
        lsu_associated_rrna_chains=list(by_role.get("lsu_associated_rrna", [])),
        other_rna_chains=other_rna_chains,
        mrna_chain=assignments.mrna_chain,
        aminoacyl_trna_chain=assignments.aminoacyl_trna_chain,
        peptidyl_trna_chain=assignments.peptidyl_trna_chain,
        exit_trna_chain=assignments.exit_trna_chain,
        aminoacyl_trna_state=states.aminoacyl_trna_state,
        peptidyl_trna_state=states.peptidyl_trna_state,
        exit_trna_state=states.exit_trna_state,
        non_ribosomal_proteins=non_ribosomal_proteins,
        bound_ligands=bound_ligands,
        classification_evidence={
            **classification_result.evidence,
            **states.trna_state_evidence,
        },
        warnings=warnings,
    )


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


def _augment_chain_rfam(chains: list[ChainRef], rfam_by_chain: dict[str, list[str]]) -> None:
    """Merge PDBe-supplied Rfam accessions into each chain's
    :attr:`ChainRef.rfam_accessions` list, preserving any RCSB-supplied
    accessions already present."""
    for chain in chains:
        for accession in rfam_by_chain.get(chain.auth_asym_id, []):
            if accession not in chain.rfam_accessions:
                chain.rfam_accessions.append(accession)


def _get_or_fetch_rfam_mappings(
    pdb_id: str, cache: Cache | None, client: httpx.Client | None
) -> dict[str, list[str]]:
    """PDBe Rfam fetch with cache wrap.

    Returns an empty dict on any failure (logs a warning) rather than
    propagating — Rfam augmentation is best-effort; classification will
    simply skip assemblies that end up with no detectable rRNA.
    """
    if cache is not None:
        cached = cache.get_pdbe_payload(pdb_id)
        if cached is not None:
            logger.debug("pdbe cache hit for %s", pdb_id)
            return parse_rfam_mappings(cached, pdb_id=pdb_id)
    logger.info("fetching PDBe Rfam mappings for %s", pdb_id)
    try:
        result = fetch_rfam_mappings(pdb_id, client=client)
    except ApiRequestError as exc:
        logger.warning("PDBe Rfam fetch failed for %s: %s", pdb_id, exc)
        return {}
    if cache is not None:
        # Persist the parsed-but-still-canonical form so a cache hit on
        # the next call returns the same result without re-parsing the
        # raw PDBe response.
        cache.put_pdbe_payload(
            pdb_id,
            {pdb_id.lower(): {"Rfam": _rebuild_rfam_block(result)}},
        )
    return result


def _rebuild_rfam_block(mapping: dict[str, list[str]]) -> dict[str, Any]:
    """Build a minimal PDBe-shaped Rfam block from a chain → accessions map.

    Symmetric with :func:`pdbe_client.parse_rfam_mappings` so a
    round-trip through the cache is lossless for the only field we
    actually use (chain → Rfam-accession).
    """
    out: dict[str, Any] = {}
    for chain_id, accessions in mapping.items():
        for accession in accessions:
            block = out.setdefault(accession, {"mappings": []})
            block["mappings"].append({"chain_id": chain_id})
    return out


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
        # Spec §6.4: organellar rRNA cores are bacterial-like, so we use the
        # E. coli reference set for v1.
        return C.ECOLI_REFERENCE_UNITS
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
