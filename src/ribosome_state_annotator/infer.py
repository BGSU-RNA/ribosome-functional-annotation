"""Functional-chain assignment and tRNA state inference (spec §11, §12).

This module is the package's biological core. It carries the most
prototype-derived logic; **rewrite from this spec, do not copy from the
prototype's** ``process_annotation.py``. The relationship to the prototype:

- :func:`assign_functional_chains` corresponds to the prototype's
  ``infer_interacting_chain`` series (``infer_atrna_ife_new``,
  ``infer_etrna_ife_new``, plus the pool-removal ``remove_ife`` logic).
- :func:`compute_trna_states` corresponds to the prototype's
  ``infer_tRNA_state``. The same conceptual contact-pattern truth tables
  (A/A vs A/P vs ap/AP vs **/* etc.) are implemented here, but every
  primitive — chain length, residue lookup, neighbour search — goes
  through :mod:`.gemmi_contacts` rather than a SQLAlchemy session.

Inputs are framework objects (:class:`gemmi.Structure`,
:class:`AssemblyContext`, :class:`CorrespondenceResult` per site). No
network, no filesystem.
"""

from __future__ import annotations

import logging
from collections.abc import Mapping
from typing import Any, cast

import gemmi
from pydantic import BaseModel, Field

from ribosome_state_annotator.constants import (
    ASL_FRAGMENT_MAX_LENGTH,
    DEFAULT_CONTACT_CUTOFF_ANGSTROM,
)
from ribosome_state_annotator.correspondence import parse_unit_id
from ribosome_state_annotator.gemmi_contacts import (
    find_neighbouring_chains,
    find_residues,
)
from ribosome_state_annotator.models import (
    AssemblyContext,
    ChainRef,
    CorrespondenceResult,
)

logger = logging.getLogger(__name__)

# Site keys for the per-site correspondence dict. Mirrors REFERENCE_SITE_KEYS
# in :mod:`constants` but kept here as separate names so the inference layer
# documents its dependencies explicitly.
SSU_MRNA = "ssu_mrna"
SSU_ATRNA = "ssu_atrna"
SSU_PTRNA = "ssu_ptrna"
SSU_ETRNA = "ssu_etrna"
LSU_ATRNA = "lsu_atrna"
LSU_PTRNA = "lsu_ptrna"
LSU_ETRNA = "lsu_etrna"


# ---------------------------------------------------------------------------
# Result models
# ---------------------------------------------------------------------------


class ChainAssignments(BaseModel):
    """Output of :func:`assign_functional_chains` (spec §11)."""

    mrna_chain: ChainRef | None = None
    aminoacyl_trna_chain: ChainRef | None = None
    peptidyl_trna_chain: ChainRef | None = None
    exit_trna_chain: ChainRef | None = None
    warnings: list[str] = Field(default_factory=list)


class TRNAStates(BaseModel):
    """Output of :func:`compute_trna_states` (spec §12).

    State strings have the form ``"<SSU>/<LSU>"`` for A-tRNA / P-tRNA
    and are exactly ``"E/E"`` for E-tRNA when assigned. Per the
    Zhou *et al.* 2014 convention (see REFERENCES.md):

    - **Lowercase letters denote SSU contacts**; uppercase letters
      denote LSU contacts. A doubled letter on either side means the
      tRNA simultaneously contacts two adjacent sites on that subunit
      (chimeric on that subunit).
    - The SSU half is ``"A"``, ``"P"``, ``"E"``, ``"ap"``, ``"pe"``,
      or ``"*"``.
    - The LSU half is ``"A"``, ``"P"``, ``"E"``, ``"AP"``, ``"PE"``,
      ``"*"``, ``"**"`` (contact found but unlabelable), or the name
      of a protein factor near the tRNA's CCA end (§12.4) — preferring
      the chain's ``uniprot_name`` and falling back to
      ``pdbx_description`` when no UniProt name is available.
    """

    aminoacyl_trna_state: str | None = None
    peptidyl_trna_state: str | None = None
    exit_trna_state: str | None = None
    trna_state_evidence: dict[str, Any] = Field(default_factory=dict)
    warnings: list[str] = Field(default_factory=list)


_RELAXED_ETRNA_CUTOFF_ANGSTROM = 8.0
"""Relaxed cutoff for the leftover-tRNA E-site heuristic in isolated SSU
structures (§33). The canonical 5 Å cutoff catches tight base-pair-style
contacts at SSU anchors; the E-site of an isolated 30S is often more
loosely engaged, so we widen the search to 8 Å but only when exactly
one tRNA-like chain is still in the candidate pool."""


# ---------------------------------------------------------------------------
# Chain assignment (§11)
# ---------------------------------------------------------------------------


def assign_functional_chains(
    structure: gemmi.Structure,
    assembly: AssemblyContext,
    by_role: Mapping[str, list[ChainRef]],
    correspondence_by_site: Mapping[str, CorrespondenceResult],
    *,
    cutoff: float = DEFAULT_CONTACT_CUTOFF_ANGSTROM,
    topology: str = "complete",
) -> ChainAssignments:
    """Assign the mRNA, A-tRNA, P-tRNA, and E-tRNA chains.

    Implements spec §11.1-§11.5 in order. Each successful assignment
    removes the chain from the candidate pool so the same chain can't be
    assigned to two roles (mirrors the prototype's ``remove_ife``).

    Args:
        structure: Parsed biological-assembly :class:`gemmi.Structure`.
        assembly: :class:`AssemblyContext` carrying all chains and ligands.
        by_role: Output of :func:`classify.partition_rna_chains_by_role`.
            The rRNA buckets identify which chains are NOT candidates.
        correspondence_by_site: Per-site :class:`CorrespondenceResult`
            from :func:`correspondence.fetch_and_filter_correspondence`.
            Sites without a result simply skip the corresponding
            assignment.
        cutoff: Contact cutoff in angstroms (§10.2 default 5.0).
        topology: ``"complete"`` (default), ``"isolated_ssu"``, or
            ``"isolated_lsu"`` — see §33. For ``isolated_ssu`` the LSU
            anchor passes are skipped (E-tRNA falls back to
            ``ssu_etrna`` and a leftover-tRNA heuristic). For
            ``isolated_lsu`` the SSU anchor passes are skipped (A and P
            are assigned via ``lsu_atrna`` / ``lsu_ptrna`` directly; E
            uses ``lsu_etrna``).
    """
    warnings: list[str] = []

    # Candidate pool: every RNA chain that is NOT in any rRNA role.
    rrna_chain_ifes: set[str] = set()
    for role in ("ssu_main_rrna", "lsu_main_rrna", "lsu_associated_rrna"):
        rrna_chain_ifes |= {c.ife for c in by_role.get(role, [])}
    candidate_pool: list[ChainRef] = [
        chain for chain in assembly.rna_chains if chain.ife not in rrna_chain_ifes
    ]

    has_ssu = topology in ("complete", "isolated_ssu")
    has_lsu = topology in ("complete", "isolated_lsu")

    # Pre-compute target residues for every site we'll use. Misses (residues
    # the BGSU correspondence pointed at but the mmCIF doesn't have) surface
    # as warnings here.
    site_targets: dict[str, dict[str, list[gemmi.Residue]]] = {}
    sites_to_load: tuple[str, ...]
    if topology == "isolated_ssu":
        sites_to_load = (SSU_MRNA, SSU_PTRNA, SSU_ATRNA, SSU_ETRNA)
    elif topology == "isolated_lsu":
        sites_to_load = (LSU_ATRNA, LSU_PTRNA, LSU_ETRNA)
    else:  # complete
        sites_to_load = (SSU_MRNA, SSU_PTRNA, SSU_ATRNA, LSU_ATRNA, LSU_PTRNA, LSU_ETRNA)
    for site_key in sites_to_load:
        residues_by_chain, missing = _lookup_target_residues(
            structure, correspondence_by_site.get(site_key)
        )
        site_targets[site_key] = residues_by_chain
        for unit_id in missing:
            warnings.append(f"missing_residue_in_coords_{site_key}_{unit_id}")

    trna_rfam_ifes = {c.ife for c in by_role.get("trna", [])}
    mrna_chain: ChainRef | None = None
    peptidyl_chain: ChainRef | None = None
    aminoacyl_chain: ChainRef | None = None
    exit_chain: ChainRef | None = None

    if has_ssu:
        # 1. mRNA (§11.2)
        #
        # tRNAs (chains annotated as Rfam ``RF00005``) cannot also be mRNA —
        # the decoding-loop SSU anchors (A1492 / A1493 / etc.) contact tRNAs
        # too. Exclude known-tRNA chains from the mRNA candidate pool only;
        # they're still eligible for A-tRNA / P-tRNA / E-tRNA assignment.
        # Description-based fallback catches mitoribosome A-tRNAs and
        # synthetic tRNA analogs that lack the Rfam RF00005 tag.
        mrna_candidates = [
            c for c in candidate_pool
            if c.ife not in trna_rfam_ifes and not _looks_like_trna(c)
        ]
        mrna_chain = _pick_chain_by_min_distance(
            structure,
            mrna_candidates,
            site_targets[SSU_MRNA],
            cutoff=cutoff,
        )
        if mrna_chain is not None:
            candidate_pool = _without(candidate_pool, mrna_chain)

        # 2. P-tRNA (§11.3) — uses ssu_ptrna anchors.
        peptidyl_chain = _pick_chain_by_min_distance(
            structure,
            candidate_pool,
            site_targets[SSU_PTRNA],
            cutoff=cutoff,
        )
        if peptidyl_chain is not None:
            candidate_pool = _without(candidate_pool, peptidyl_chain)

        # 3. A-tRNA (§11.4) — uses ssu_atrna anchors.
        aminoacyl_chain = _pick_chain_by_min_distance(
            structure,
            candidate_pool,
            site_targets[SSU_ATRNA],
            cutoff=cutoff,
        )
        if aminoacyl_chain is not None:
            candidate_pool = _without(candidate_pool, aminoacyl_chain)

    # 4. E-tRNA — anchor source depends on topology.
    #
    # For complete / isolated_lsu: use lsu_etrna (canonical).
    # For isolated_ssu: use ssu_etrna. If no chain contacts the SSU E-site
    # anchor within cutoff, fall back to the "leftover tRNA = likely E"
    # heuristic — if exactly one tRNA-like chain is still in the
    # candidate pool and contacts ssu_etrna within a relaxed 8 Å cutoff,
    # assign it. This catches typical 30S structures with 3 bound tRNAs
    # (A, P, E) where the third tRNA sits at the SSU exit site but
    # doesn't tightly contact the canonical SSU E anchors (G693 / A694
    # in E. coli).
    if has_lsu:
        exit_chain = _pick_chain_by_min_distance(
            structure,
            candidate_pool,
            site_targets[LSU_ETRNA],
            cutoff=cutoff,
        )
    else:  # isolated_ssu
        exit_chain = _pick_chain_by_min_distance(
            structure,
            candidate_pool,
            site_targets[SSU_ETRNA],
            cutoff=cutoff,
        )
        if exit_chain is None:
            leftover_trnas = [
                c
                for c in candidate_pool
                if c.ife in trna_rfam_ifes or _looks_like_trna(c)
            ]
            if len(leftover_trnas) == 1:
                relaxed = _RELAXED_ETRNA_CUTOFF_ANGSTROM
                dists = _compute_chain_distances(
                    structure,
                    leftover_trnas,
                    site_targets[SSU_ETRNA],
                    cutoff=relaxed,
                )
                if dists:
                    exit_chain = leftover_trnas[0]
    if exit_chain is not None:
        candidate_pool = _without(candidate_pool, exit_chain)

    # 5. LSU-based fallback for unfilled A and P slots.
    #
    # Synthetic tRNA analogs (acceptor-end-only fragments in 7RQA, 8T8C),
    # real pre-accommodation tRNAs (3JAG, 7O7Z, 7OSM, 7UG7), and
    # mitoribosome P-tRNAs (where the SSU decoding centre anchors don't
    # map cleanly) contact the PTC tightly but are far from the SSU
    # anchors. After the canonical pass fills what it can, fall back to
    # assigning each remaining candidate to whichever LSU site
    # (lsu_atrna or lsu_ptrna) it is **closer** to — not just to A
    # first.
    #
    # The closer-site choice is essential for mitoribosomal P-tRNAs
    # which incidentally contact both PTC anchor clusters: with strict
    # A-before-P ordering 8OIR's P-Met-tRNA wins the A slot before
    # P-fallback runs, even though it's closer to the P anchors. The
    # closer-site choice picks the correct site without changing the
    # canonical-SSU happy path (where A-tRNA is picked via ssu_atrna
    # before the fallback fires at all).
    if has_lsu and (aminoacyl_chain is None or peptidyl_chain is None):
        lsu_a_distances = _compute_chain_distances(
            structure, candidate_pool, site_targets[LSU_ATRNA], cutoff=cutoff
        )
        lsu_p_distances = _compute_chain_distances(
            structure, candidate_pool, site_targets[LSU_PTRNA], cutoff=cutoff
        )
        # Build (chain, d_a, d_p) tuples for candidates contacting either site.
        triples: list[tuple[ChainRef, float | None, float | None]] = []
        for chain in candidate_pool:
            da = lsu_a_distances.get(chain.auth_asym_id)
            dp = lsu_p_distances.get(chain.auth_asym_id)
            if da is None and dp is None:
                continue
            triples.append((chain, da, dp))

        # Best-first ordering: candidates with the smallest LSU-side
        # contact get assigned first. Ties go to A.
        def _best_distance(t: tuple[ChainRef, float | None, float | None]) -> float:
            distances = [d for d in (t[1], t[2]) if d is not None]
            return min(distances) if distances else float("inf")

        triples.sort(key=_best_distance)

        for chain, da, dp in triples:
            if aminoacyl_chain is not None and peptidyl_chain is not None:
                break
            # Preferred site = the LSU anchor cluster this chain sits
            # closer to. Ties (or single-site contact) collapse to that
            # single site. If the preferred site is already filled, fall
            # back to the other site (only if this chain also makes a
            # ≤cutoff contact there).
            if da is not None and (dp is None or da <= dp):
                preferred: tuple[str, float | None] = ("A", da)
                alternative: tuple[str, float | None] = ("P", dp)
            else:
                preferred = ("P", dp)
                alternative = ("A", da)
            for site, _d in (preferred, alternative):
                if site == "A" and aminoacyl_chain is None:
                    aminoacyl_chain = chain
                    candidate_pool = _without(candidate_pool, chain)
                    break
                if site == "P" and peptidyl_chain is None:
                    peptidyl_chain = chain
                    candidate_pool = _without(candidate_pool, chain)
                    break

    return ChainAssignments(
        mrna_chain=mrna_chain,
        aminoacyl_trna_chain=aminoacyl_chain,
        peptidyl_trna_chain=peptidyl_chain,
        exit_trna_chain=exit_chain,
        warnings=warnings,
    )


def _pick_chain_by_min_distance(
    structure: gemmi.Structure,
    candidates: list[ChainRef],
    target_residues_by_chain: dict[str, list[gemmi.Residue]],
    *,
    cutoff: float,
) -> ChainRef | None:
    """Return the candidate ChainRef with the smallest min-distance contact
    to the anchor residues, or ``None`` if no candidate makes contact.

    Ties (same min distance to multiple candidates) are broken by the
    candidate's position in the input list — deterministic but otherwise
    arbitrary; in practice ties are vanishingly rare.
    """
    if not candidates or not target_residues_by_chain:
        return None
    by_auth = {c.auth_asym_id: c for c in candidates}
    distances = find_neighbouring_chains(
        structure,
        target_residues_by_chain,
        cutoff=cutoff,
        candidate_chains=by_auth.keys(),
    )
    if not distances:
        return None
    best_chain_name = min(distances, key=lambda name: distances[name])
    return by_auth.get(best_chain_name)


def _without(pool: list[ChainRef], to_remove: ChainRef) -> list[ChainRef]:
    """Return ``pool`` minus ``to_remove`` (compared by IFE — the canonical
    chain identity used throughout the package)."""
    return [c for c in pool if c.ife != to_remove.ife]


def _compute_chain_distances(
    structure: gemmi.Structure,
    candidates: list[ChainRef],
    target_residues_by_chain: dict[str, list[gemmi.Residue]],
    *,
    cutoff: float,
) -> dict[str, float]:
    """Return ``{auth_asym_id: min_distance}`` for each candidate that has at
    least one atom within ``cutoff`` of any anchor residue. Same neighbour
    search as :func:`_pick_chain_by_min_distance`, just returning the whole
    distance map rather than the winner."""
    if not candidates or not target_residues_by_chain:
        return {}
    by_auth = {c.auth_asym_id: c for c in candidates}
    return find_neighbouring_chains(
        structure,
        target_residues_by_chain,
        cutoff=cutoff,
        candidate_chains=by_auth.keys(),
    )


def _looks_like_trna(chain: ChainRef) -> bool:
    """Description-based tRNA detection for the mRNA-pool exclusion.

    Catches tRNAs that lack the Rfam ``RF00005`` tag but whose deposit
    description still says "tRNA" (e.g. mitoribosome A-tRNAs in 7QI5
    chain Aw, or peptidyl-tRNA analogs in 7RQA / 8T8C). Case-insensitive
    substring match on ``description``.
    """
    desc = chain.description or ""
    return "trna" in desc.lower()


def _lookup_target_residues(
    structure: gemmi.Structure,
    correspondence: CorrespondenceResult | None,
) -> tuple[dict[str, list[gemmi.Residue]], list[str]]:
    """Resolve a :class:`CorrespondenceResult`'s mapped unit IDs to gemmi
    residues, grouped by chain. Returns ``(residues_by_chain, missing_ids)``.
    """
    if correspondence is None or not correspondence.mapped_units:
        return {}, []
    units = [parse_unit_id(u) for u in correspondence.mapped_units]
    return find_residues(structure, units)


# ---------------------------------------------------------------------------
# tRNA state computation (§12)
# ---------------------------------------------------------------------------


def compute_trna_states(
    structure: gemmi.Structure,
    assembly: AssemblyContext,
    assignments: ChainAssignments,
    correspondence_by_site: Mapping[str, CorrespondenceResult],
    *,
    cutoff: float = DEFAULT_CONTACT_CUTOFF_ANGSTROM,
    topology: str = "complete",
) -> TRNAStates:
    """Compute the SSU/LSU functional state strings for the assigned tRNAs.

    Spec §12.1 / §12.2 / §12.3 (+ §33 topology extension):

    - A-tRNA: ``<SSU state>/<LSU state>`` with SSU ∈ {``A``, ``ap``, ``*``,
      ``**``, ``-``} and LSU ∈ {``A``, ``P``, ``AP``, ``**``, ``*``,
      ``-``, <factor description>}.
    - P-tRNA: SSU ∈ {``P``, ``pe``, ``*``, ``**``, ``-``} and LSU ∈
      {``P``, ``E``, ``PE``, ``**``, ``*``, ``-``, <factor>}.
    - E-tRNA: SSU ∈ {``E``, ``*``, ``**``, ``-``} and LSU ∈ {``E``, ``-``}.

    Half-state vocabulary:

      - ``A``, ``P``, ``E`` — canonical contact at the corresponding site.
      - ``ap`` / ``pe`` — SSU hybrid.
      - ``AP`` / ``PE`` — LSU hybrid.
      - ``*`` — full-length polymer (≥30 nt) with no canonical contact at
        this subunit; positionally displaced.
      - ``**`` — chain shorter than :data:`ASL_FRAGMENT_MAX_LENGTH`; can't
        physically reach this subunit.
      - ``-`` — subunit not present in the assembly (topology=isolated_*).
      - ``<factor>`` — LSU-side label when the tRNA's CCA end is engaged
        by a non-ribosomal protein factor (§12.4).

    Evidence keys populated on ``trna_state_evidence`` (when applicable):
    ``aminoacyl_trna_factor_*`` and ``peptidyl_trna_factor_*`` per §12.4.
    """
    warnings: list[str] = []
    evidence: dict[str, Any] = {}

    has_ssu = topology in ("complete", "isolated_ssu")
    has_lsu = topology in ("complete", "isolated_lsu")

    # Pre-compute every site's target residues.
    state_sites = (SSU_ATRNA, SSU_PTRNA, SSU_ETRNA, LSU_ATRNA, LSU_PTRNA, LSU_ETRNA)
    site_targets: dict[str, dict[str, list[gemmi.Residue]]] = {}
    for site_key in state_sites:
        residues_by_chain, missing = _lookup_target_residues(
            structure, correspondence_by_site.get(site_key)
        )
        site_targets[site_key] = residues_by_chain
        for unit_id in missing:
            warnings.append(f"missing_residue_in_coords_{site_key}_{unit_id}")

    aminoacyl_state = _compute_atrna_state(
        structure,
        assignments.aminoacyl_trna_chain,
        site_targets,
        assembly.protein_chains,
        cutoff=cutoff,
        evidence=evidence,
        has_ssu=has_ssu,
        has_lsu=has_lsu,
    )
    peptidyl_state = _compute_ptrna_state(
        structure,
        assignments.peptidyl_trna_chain,
        site_targets,
        assembly.protein_chains,
        cutoff=cutoff,
        evidence=evidence,
        has_ssu=has_ssu,
        has_lsu=has_lsu,
    )
    exit_state = _compute_etrna_state(
        structure,
        assignments.exit_trna_chain,
        site_targets,
        cutoff=cutoff,
        has_ssu=has_ssu,
        has_lsu=has_lsu,
    )

    return TRNAStates(
        aminoacyl_trna_state=aminoacyl_state,
        peptidyl_trna_state=peptidyl_state,
        exit_trna_state=exit_state,
        trna_state_evidence=evidence,
        warnings=warnings,
    )


def _compute_atrna_state(
    structure: gemmi.Structure,
    atrna_chain: ChainRef | None,
    site_targets: Mapping[str, dict[str, list[gemmi.Residue]]],
    protein_chains: list[ChainRef],
    *,
    cutoff: float,
    evidence: dict[str, Any],
    has_ssu: bool = True,
    has_lsu: bool = True,
) -> str | None:
    """Implement the §12.1 A-tRNA state rules."""
    if atrna_chain is None:
        return None
    chain_name = atrna_chain.auth_asym_id

    if has_ssu:
        contacts_ssu_atrna = _chain_contacts_site(
            structure, chain_name, site_targets[SSU_ATRNA], cutoff=cutoff
        )
        contacts_ssu_ptrna = _chain_contacts_site(
            structure, chain_name, site_targets[SSU_PTRNA], cutoff=cutoff
        )
        # SSU state. When the chain doesn't engage the SSU decoding centre
        # at all, distinguish a structurally-too-short fragment (CCA-end
        # analogs, < 30 nt) from a full-length tRNA that's positionally
        # displaced. Convention: ``**`` = fragment, ``*`` = displaced full
        # polymer.
        if not contacts_ssu_atrna:
            chain_length = _get_chain_length(structure, chain_name)
            ssu_state = "**" if chain_length < ASL_FRAGMENT_MAX_LENGTH else "*"
        elif contacts_ssu_ptrna:
            ssu_state = "ap"
        else:
            ssu_state = "A"
    else:
        ssu_state = "-"

    if has_lsu:
        contacts_lsu_atrna = _chain_contacts_site(
            structure, chain_name, site_targets[LSU_ATRNA], cutoff=cutoff
        )
        contacts_lsu_ptrna = _chain_contacts_site(
            structure, chain_name, site_targets[LSU_PTRNA], cutoff=cutoff
        )
        lsu_state = _resolve_lsu_state(
            structure,
            atrna_chain,
            protein_chains,
            contacts_primary=contacts_lsu_atrna,
            contacts_secondary=contacts_lsu_ptrna,
            primary_label="A",
            secondary_label="P",
            combined_label="AP",
            cutoff=cutoff,
            evidence_prefix="aminoacyl_trna",
            evidence=evidence,
        )
    else:
        lsu_state = "-"
    return f"{ssu_state}/{lsu_state}"


def _compute_ptrna_state(
    structure: gemmi.Structure,
    ptrna_chain: ChainRef | None,
    site_targets: Mapping[str, dict[str, list[gemmi.Residue]]],
    protein_chains: list[ChainRef],
    *,
    cutoff: float,
    evidence: dict[str, Any],
    has_ssu: bool = True,
    has_lsu: bool = True,
) -> str | None:
    """Implement the §12.2 P-tRNA state rules."""
    if ptrna_chain is None:
        return None
    chain_name = ptrna_chain.auth_asym_id

    if has_ssu:
        contacts_ssu_ptrna = _chain_contacts_site(
            structure, chain_name, site_targets[SSU_PTRNA], cutoff=cutoff
        )
        contacts_ssu_etrna = _chain_contacts_site(
            structure, chain_name, site_targets[SSU_ETRNA], cutoff=cutoff
        )
        # Same fragment-vs-displaced convention as the A-tRNA path.
        if not contacts_ssu_ptrna:
            chain_length = _get_chain_length(structure, chain_name)
            ssu_state = "**" if chain_length < ASL_FRAGMENT_MAX_LENGTH else "*"
        elif contacts_ssu_etrna:
            ssu_state = "pe"
        else:
            ssu_state = "P"
    else:
        ssu_state = "-"

    if has_lsu:
        contacts_lsu_ptrna = _chain_contacts_site(
            structure, chain_name, site_targets[LSU_PTRNA], cutoff=cutoff
        )
        contacts_lsu_etrna = _chain_contacts_site(
            structure, chain_name, site_targets[LSU_ETRNA], cutoff=cutoff
        )
        lsu_state = _resolve_lsu_state(
            structure,
            ptrna_chain,
            protein_chains,
            contacts_primary=contacts_lsu_ptrna,
            contacts_secondary=contacts_lsu_etrna,
            primary_label="P",
            secondary_label="E",
            combined_label="PE",
            cutoff=cutoff,
            evidence_prefix="peptidyl_trna",
            evidence=evidence,
        )
    else:
        lsu_state = "-"
    return f"{ssu_state}/{lsu_state}"


def _compute_etrna_state(
    structure: gemmi.Structure,
    etrna_chain: ChainRef | None,
    site_targets: Mapping[str, dict[str, list[gemmi.Residue]]],
    *,
    cutoff: float,
    has_ssu: bool = True,
    has_lsu: bool = True,
) -> str | None:
    """Implement the §12.3 E-tRNA state rules.

    Canonical full-length E-tRNA produces ``E/E``. A chain that doesn't
    contact the SSU exit-site anchor (``ssu_etrna``) gets ``*/E`` if
    full-length, ``**/E`` if shorter than
    :data:`ASL_FRAGMENT_MAX_LENGTH`. For isolated topologies the absent
    subunit's half-state is ``-``.
    """
    if etrna_chain is None:
        return None
    chain_name = etrna_chain.auth_asym_id

    if has_ssu:
        contacts_ssu_etrna = _chain_contacts_site(
            structure, chain_name, site_targets[SSU_ETRNA], cutoff=cutoff
        )
        if contacts_ssu_etrna:
            ssu_state = "E"
        else:
            chain_length = _get_chain_length(structure, chain_name)
            ssu_state = "**" if chain_length < ASL_FRAGMENT_MAX_LENGTH else "*"
    else:
        ssu_state = "-"
    lsu_state = "E" if has_lsu else "-"
    return f"{ssu_state}/{lsu_state}"


def _resolve_lsu_state(
    structure: gemmi.Structure,
    trna_chain: ChainRef,
    protein_chains: list[ChainRef],
    *,
    contacts_primary: bool,
    contacts_secondary: bool,
    primary_label: str,
    secondary_label: str,
    combined_label: str,
    cutoff: float,
    evidence_prefix: str,
    evidence: dict[str, Any],
) -> str:
    """Shared §12.1/§12.2 LSU-side truth table.

    - contacts both → combined (``AP`` / ``PE``)
    - contacts primary only → primary (``A`` / ``P``)
    - contacts secondary only → secondary (``P`` / ``E``)
    - contacts neither → ASL fragment (``**``) or §12.4 protein-factor /
      ``*``
    """
    if contacts_primary and contacts_secondary:
        return combined_label
    if contacts_primary:
        return primary_label
    if contacts_secondary:
        return secondary_label

    chain_length = _get_chain_length(structure, trna_chain.auth_asym_id)
    if chain_length < ASL_FRAGMENT_MAX_LENGTH:
        return "**"
    factor = _compute_protein_factor_label(structure, trna_chain, protein_chains, cutoff=cutoff)
    if factor is None:
        return "*"
    factor_chain, factor_distance = factor
    evidence[f"{evidence_prefix}_factor_chain"] = factor_chain.ife
    evidence[f"{evidence_prefix}_factor_description"] = factor_chain.description
    evidence[f"{evidence_prefix}_factor_uniprot_name"] = factor_chain.uniprot_name
    evidence[f"{evidence_prefix}_factor_distance"] = factor_distance
    # Prefer the UniProt protein name in the state label — it is the
    # canonical cross-deposit identifier and is more stable than the
    # pdbx_description (which often carries copy suffixes like "Tu 2").
    # Fall back to pdbx_description, then to "*" when neither is set.
    return factor_chain.uniprot_name or factor_chain.description or "*"


# ---------------------------------------------------------------------------
# §12.4 protein-factor LSU label
# ---------------------------------------------------------------------------


def _compute_protein_factor_label(
    structure: gemmi.Structure,
    trna_chain: ChainRef,
    protein_chains: list[ChainRef],
    *,
    cutoff: float,
) -> tuple[ChainRef, float] | None:
    """Spec §12.4: identify the nearest non-ribosomal protein at the
    tRNA's CCA end.

    Algorithm (verbatim from §12.4):

    1. Take the **last three residues** of the tRNA chain by highest
       author seqid (approximate CCA end).
    2. Run a Gemmi neighbour search around all atoms of those residues
       at ``cutoff``.
    3. Among protein chains in the assembly with ≥1 atom inside the
       cutoff, **exclude chains classified as ribosomal protein**
       (the §13.1 narrow rule, surfaced as ``ChainRef.is_ribosomal_protein``).
    4. Pick the chain with the minimum atom-atom distance to any CCA-end
       atom.

    Returns ``(chain, min_distance)`` or ``None`` if no qualifying
    protein chain was found.
    """
    cca_residues = _get_cca_end_residues(structure, trna_chain.auth_asym_id, n=3)
    if not cca_residues:
        return None
    non_ribosomal = [p for p in protein_chains if not p.is_ribosomal_protein]
    if not non_ribosomal:
        return None
    by_auth = {p.auth_asym_id: p for p in non_ribosomal}
    distances = find_neighbouring_chains(
        structure,
        {trna_chain.auth_asym_id: cca_residues},
        cutoff=cutoff,
        candidate_chains=by_auth.keys(),
    )
    if not distances:
        return None
    best_chain_name = min(distances, key=lambda name: distances[name])
    best_chain = by_auth.get(best_chain_name)
    if best_chain is None:
        return None
    return (best_chain, distances[best_chain_name])


# ---------------------------------------------------------------------------
# Gemmi helpers
# ---------------------------------------------------------------------------


def _chain_contacts_site(
    structure: gemmi.Structure,
    chain_name: str,
    target_residues_by_chain: Mapping[str, list[gemmi.Residue]],
    *,
    cutoff: float,
) -> bool:
    """Return True iff any atom on ``chain_name`` is within ``cutoff`` of
    any atom on the target residues."""
    if not target_residues_by_chain:
        return False
    if chain_name in target_residues_by_chain:
        # The site's anchor residues are on the same chain as the query
        # chain — a target rRNA can't "contact" itself meaningfully for
        # our purposes.
        return False
    distances = find_neighbouring_chains(
        structure,
        target_residues_by_chain,
        cutoff=cutoff,
        candidate_chains={chain_name},
    )
    return chain_name in distances


def _get_chain_length(structure: gemmi.Structure, chain_name: str) -> int:
    """Return the number of resolved residues for ``chain_name`` in model 0."""
    if len(structure) == 0:
        return 0
    # gemmi 0.7.5 stubs declare find_chain → Chain but it returns None at
    # runtime when the chain is missing. Cast to Optional so mypy lets us
    # short-circuit.
    chain = cast(gemmi.Chain | None, structure[0].find_chain(chain_name))
    if chain is None:
        return 0
    return sum(1 for _ in chain)


def _get_cca_end_residues(
    structure: gemmi.Structure, chain_name: str, *, n: int = 3
) -> list[gemmi.Residue]:
    """Return the ``n`` polymer residues with the highest author seqid on
    ``chain_name``.

    Per spec §12.4 this approximates the universally-conserved 3'-CCA end
    of a mature tRNA. Returning fewer than ``n`` is acceptable when the
    chain has fewer resolved polymer residues.

    Filters out non-polymer residues (waters, ions, ligands) that mmCIF
    files often park under the same chain ID as the macromolecule. A
    tRNA chain typically has its real CCA end at residue ~76 followed by
    several HETATM cofactors at much higher seqids; if we naively took
    the last three by seqid we'd pick up `MG` / `HOH` and the §12.4
    factor search would silently miss the actual elongation factor (the
    distance from a bound water to the factor protein is meaningless).
    Modified nucleotides (`4SU`, `PSU`, `MIA`, …) remain `EntityType.Polymer`
    so they're correctly retained as candidate CCA-end residues.
    """
    if len(structure) == 0:
        return []
    chain = cast(gemmi.Chain | None, structure[0].find_chain(chain_name))
    if chain is None:
        return []
    polymer_residues = [r for r in chain if r.entity_type == gemmi.EntityType.Polymer]
    # gemmi typing declares Residue.seqid.num as int | None; in practice every
    # residue we iterate has a resolved seqid. Treat None as -1 so the sort
    # is total-ordered even in the pathological case.
    polymer_residues.sort(key=lambda r: r.seqid.num if r.seqid.num is not None else -1)
    return polymer_residues[-n:]


__all__ = [
    "ChainAssignments",
    "TRNAStates",
    "assign_functional_chains",
    "compute_trna_states",
]
