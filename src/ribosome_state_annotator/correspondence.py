"""BGSU unit-ID parsing and §5.2.2 correspondence filtering.

This module owns:

- :class:`UnitId` and :func:`parse_unit_id` / :func:`try_parse_unit_id`,
  which parse BGSU unit IDs per spec §10.3 (handling modified-nucleotide
  residue names like ``4OC`` and ``OMG``);
- the §5.2.2 two-step filter (:func:`filter_mapped_units`) and the
  :class:`CorrespondenceResult` builder (:func:`build_correspondence_result`);
- the orchestrator :func:`fetch_and_filter_correspondence` that bridges
  :mod:`.bgsu_client` with the filter.

Step 8's neighbour-search layer will reuse :class:`UnitId` to look up
the corresponding residue in a parsed Gemmi structure.
"""

from __future__ import annotations

import logging
from collections.abc import Iterable

import httpx
from pydantic import BaseModel, ConfigDict, computed_field

from ribosome_state_annotator.bgsu_client import fetch_correspondence
from ribosome_state_annotator.models import CorrespondenceResult

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# UnitId parser (spec §10.3)
# ---------------------------------------------------------------------------


class UnitId(BaseModel):
    """Parsed BGSU unit ID (spec §10.3).

    Format: ``PDB|model|chain|residue_name|residue_number`` plus optional
    extended segments (``|atom|alt|insertion|symop``) that v1 ignores. The
    first 5 segments are required and non-empty; ``residue_number`` must
    be an integer (may be negative).

    The model is frozen because a unit ID is an identity — once parsed,
    none of its fields should be mutated. Equality is structural across
    all 5 fields.
    """

    model_config = ConfigDict(frozen=True)

    pdb_id: str
    model: str
    chain: str
    residue_name: str
    residue_number: int

    @computed_field  # type: ignore[prop-decorator]
    @property
    def base_unit_id(self) -> str:
        """Canonical 5-segment string form, e.g. ``"5J7L|1|AA|G|926"``."""
        return f"{self.pdb_id}|{self.model}|{self.chain}|{self.residue_name}|{self.residue_number}"


def parse_unit_id(unit_id: str) -> UnitId:
    """Parse a BGSU unit ID string into a :class:`UnitId`.

    Accepts the canonical 5-segment form as well as the extended 6-to-9
    segment form (BGSU sometimes appends ``|atom|alt|insertion|symop``);
    only the first 5 segments are interpreted in v1.

    Modified-nucleotide residue names (``4OC``, ``OMG``, etc.) are
    accepted as-is — no normalisation, no validation against PDB's
    chemical component dictionary. The §10.4 Gemmi lookup tolerates
    name differences and matches on residue number instead.

    Raises:
        ValueError: if ``unit_id`` is empty, has fewer than 5 segments,
            has an empty required field, or has a non-integer
            residue_number.
    """
    if not unit_id:
        raise ValueError("unit ID is empty")
    parts = unit_id.split("|")
    if len(parts) < 5:
        raise ValueError(
            f"unit ID {unit_id!r} has {len(parts)} segments; need at least 5 "
            "(PDB|model|chain|residue_name|residue_number)"
        )
    pdb_id, model, chain, residue_name, residue_number_str = parts[:5]
    if not pdb_id:
        raise ValueError(f"unit ID {unit_id!r} has empty PDB ID")
    if not model:
        raise ValueError(f"unit ID {unit_id!r} has empty model segment")
    if not chain:
        raise ValueError(f"unit ID {unit_id!r} has empty chain segment")
    if not residue_name:
        raise ValueError(f"unit ID {unit_id!r} has empty residue_name segment")
    if not residue_number_str:
        raise ValueError(f"unit ID {unit_id!r} has empty residue_number segment")
    try:
        residue_number = int(residue_number_str)
    except ValueError as exc:
        raise ValueError(
            f"unit ID {unit_id!r} has non-integer residue_number {residue_number_str!r}"
        ) from exc
    return UnitId(
        pdb_id=pdb_id,
        model=model,
        chain=chain,
        residue_name=residue_name,
        residue_number=residue_number,
    )


def try_parse_unit_id(unit_id: str) -> UnitId | None:
    """Forgiving wrapper around :func:`parse_unit_id`.

    Returns ``None`` (and logs at DEBUG) on any parse failure. Used by
    the §5.2.2 filter so that one malformed unit ID in a BGSU response
    doesn't take down the whole pipeline.
    """
    try:
        return parse_unit_id(unit_id)
    except ValueError as exc:
        logger.debug("skipping malformed unit ID %r: %s", unit_id, exc)
        return None


# ---------------------------------------------------------------------------
# §5.2.2 filter
# ---------------------------------------------------------------------------


def filter_mapped_units(
    mapped_units: Iterable[str],
    *,
    target_pdb_id: str,
    assembly_chains: set[str],
) -> list[str]:
    """Apply the spec §5.2.2 two-step filter to a list of mapped unit IDs.

    Step 1: keep only unit IDs whose PDB segment matches ``target_pdb_id``
    (case-insensitive). Discards rows from any other PDB entry in the
    Rfam alignment.

    Step 2: keep only unit IDs whose chain segment appears in
    ``assembly_chains``. This is essential because a PDB entry can carry
    multiple biological assemblies and the BGSU correspondence response
    crosses all of them.

    Returns the surviving unit IDs in input order, deduplicated.
    Malformed unit IDs (fewer than three pipe-separated segments) are
    silently skipped.
    """
    target_lower = target_pdb_id.lower()
    seen: set[str] = set()
    out: list[str] = []
    for unit in mapped_units:
        parsed = try_parse_unit_id(unit)
        if parsed is None:
            continue
        if parsed.pdb_id.lower() != target_lower:
            continue
        if parsed.chain not in assembly_chains:
            continue
        if unit in seen:
            continue
        seen.add(unit)
        out.append(unit)
    return out


def group_by_chain(mapped_units: Iterable[str]) -> dict[str, list[str]]:
    """Group mapped unit IDs by their chain segment, preserving input order."""
    by_chain: dict[str, list[str]] = {}
    for unit in mapped_units:
        parsed = try_parse_unit_id(unit)
        if parsed is None:
            continue
        by_chain.setdefault(parsed.chain, []).append(unit)
    return by_chain


def _try_chain_substitution(
    reference_unit: str,
    *,
    target_pdb_id: str,
    assembly_chains: set[str],
    chain_substitution: dict[str, str],
) -> str | None:
    """Apply ``chain_substitution`` to ``reference_unit`` and return the
    target-PDB unit ID, or ``None`` if the substitution doesn't apply or
    the resulting chain isn't in the assembly.

    This is the same-organism fallback for when BGSU's correspondence
    response doesn't include the target PDB (BGSU returns only NR
    equivalence-class representatives).
    """
    parsed = try_parse_unit_id(reference_unit)
    if parsed is None:
        return None
    target_chain = chain_substitution.get(parsed.chain)
    if target_chain is None or target_chain not in assembly_chains:
        return None
    return (
        f"{target_pdb_id}|{parsed.model}|{target_chain}|"
        f"{parsed.residue_name}|{parsed.residue_number}"
    )


# ---------------------------------------------------------------------------
# CorrespondenceResult builder
# ---------------------------------------------------------------------------


def build_correspondence_result(
    reference_key: str,
    reference_units: list[str],
    raw_alignment: dict[str, list[str]],
    *,
    target_pdb_id: str,
    assembly_chains: set[str],
    chain_substitution: dict[str, str] | None = None,
) -> CorrespondenceResult:
    """Assemble a :class:`CorrespondenceResult` for one functional-site key.

    For each reference unit in ``reference_units``, filter the alignment's
    candidates through :func:`filter_mapped_units` and accumulate the
    survivors. Reference units that map to **zero** in-assembly units
    contribute a per-anchor warning using the exact §5.2.3 format
    ``"correspondence_missing_for_{reference_key}_{unit_id}"`` and do not
    skip the whole site (per §5.2.3, "a site … with at least one mapped
    nucleotide in the assembly is sufficient").

    The resulting ``mapped_units`` list is deduplicated across all anchors
    in the site (the same residue is sometimes anchored by adjacent
    reference units). ``mapped_units_by_chain`` groups by author chain ID
    so the Gemmi layer (step 8) can iterate per chain.

    ``chain_substitution`` is an optional ``{reference_chain: target_chain}``
    fallback used when BGSU returns no mappings for the target PDB —
    common in the **same-organism** case because BGSU only returns
    non-redundant equivalence-class representatives, never the query PDB
    itself. When the reference and target share residue numbering
    convention (true for E. coli reference → another E. coli structure),
    substituting the chain segment of the reference unit produces a
    valid in-target unit ID without needing BGSU's help.
    """
    all_mapped: list[str] = []
    warnings: list[str] = []
    for ref_unit in reference_units:
        candidates = raw_alignment.get(ref_unit, [])
        filtered = filter_mapped_units(
            candidates,
            target_pdb_id=target_pdb_id,
            assembly_chains=assembly_chains,
        )
        if not filtered and chain_substitution:
            substituted = _try_chain_substitution(
                ref_unit,
                target_pdb_id=target_pdb_id,
                assembly_chains=assembly_chains,
                chain_substitution=chain_substitution,
            )
            if substituted is not None:
                filtered = [substituted]
        if not filtered:
            warnings.append(f"correspondence_missing_for_{reference_key}_{ref_unit}")
            continue
        all_mapped.extend(filtered)

    # Deduplicate across reference units within the same site while
    # preserving insertion order.
    dedup_mapped = list(dict.fromkeys(all_mapped))
    by_chain = group_by_chain(dedup_mapped)

    return CorrespondenceResult(
        reference_key=reference_key,
        reference_units=list(reference_units),
        mapped_units=dedup_mapped,
        mapped_units_by_chain=by_chain,
        warnings=warnings,
    )


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------


def fetch_and_filter_correspondence(
    reference_key: str,
    reference_units: list[str],
    *,
    target_pdb_id: str,
    assembly_chains: set[str],
    client: httpx.Client | None = None,
) -> CorrespondenceResult:
    """High-level: fetch BGSU mappings for one site and apply §5.2.2 filter.

    Wraps :func:`.bgsu_client.fetch_correspondence` and
    :func:`build_correspondence_result`. Always returns a
    :class:`CorrespondenceResult`; per-anchor failures become warnings on
    the result, not exceptions. Whole-fetch failures
    (:class:`~.exceptions.ApiRequestError` / :class:`~.exceptions.CorrespondenceMappingError`)
    propagate to the caller.
    """
    if not reference_units:
        return CorrespondenceResult(reference_key=reference_key)
    raw = fetch_correspondence(reference_units, client=client)
    return build_correspondence_result(
        reference_key,
        reference_units,
        raw,
        target_pdb_id=target_pdb_id,
        assembly_chains=assembly_chains,
    )
