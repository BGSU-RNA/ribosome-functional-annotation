"""Multi-ribosome assembly splitting (spec §11 addendum).

Some PDB deposits pack **two complete ribosomes** into a single
biological assembly (e.g. Mycobacterium 70S dimer entries 8R3V, 8RCL,
9GFT, 9O3L). The default contact-transfer in :mod:`infer` fills only
one A/P/E slot per assembly, leaving the second ribosome's mRNA + tRNAs
in ``unmapped_rna_chains``.

This module detects multi-ribosome bundles and partitions the assembly
into per-ribosome sub-contexts:

  1. :func:`detect_multi_ribosome` — gate test on the role-partition
     output (``len(ssu_main) == len(lsu_main) >= 2``).
  2. :func:`pair_ssu_lsu_by_centroid` — greedy nearest-centroid matching
     of SSU↔LSU chains using atom-coordinate centroids.
  3. :func:`partition_chains_by_ribosome` — assigns each non-rRNA RNA
     chain to the ribosome (pair) whose combined centroid it is closest
     to.
  4. :func:`filter_correspondence_for_chains` — restricts an existing
     :class:`CorrespondenceResult` to a per-ribosome chain subset so
     contact-transfer only sees that ribosome's anchor residues.

The orchestrator (:mod:`api`) uses these helpers to emit one
:class:`RibosomeAnnotation` per detected ribosome with an
``assembly_id`` suffix like ``"1-1"`` / ``"1-2"`` while keeping
single-ribosome assemblies unchanged.
"""

from __future__ import annotations

import logging
from collections.abc import Iterable

import gemmi

from ribosome_state_annotator.correspondence import group_by_chain, try_parse_unit_id
from ribosome_state_annotator.models import ChainRef, CorrespondenceResult

logger = logging.getLogger(__name__)


def detect_multi_ribosome(by_role: dict[str, list[ChainRef]]) -> int:
    """Return the number of ribosomes if this is a multi-ribosome bundle, else 1.

    A bundle qualifies when there are at least two SSU main rRNA chains
    AND at least two LSU main rRNA chains AND the two counts are equal —
    i.e. the deposit packs N complete SSU + LSU pairs. Equal counts is a
    deliberate guard: unequal counts (e.g. 2 SSU + 1 LSU) imply some
    fragmentation rather than a clean dimer and is handled by the
    fragmented-LSU path instead.
    """
    n_ssu = len(by_role.get("ssu_main_rrna", []))
    n_lsu = len(by_role.get("lsu_main_rrna", []))
    if n_ssu >= 2 and n_lsu >= 2 and n_ssu == n_lsu:
        return n_ssu
    return 1


def _chain_centroid(structure: gemmi.Structure, chain_id: str) -> gemmi.Position | None:
    """Return the centroid of all atom positions for ``chain_id`` in model 0.

    Returns ``None`` if the chain is absent or has no atoms.
    """
    model = structure[0]
    n = 0
    sx = sy = sz = 0.0
    for chain in model:
        if chain.name != chain_id:
            continue
        for res in chain:
            for atom in res:
                sx += atom.pos.x
                sy += atom.pos.y
                sz += atom.pos.z
                n += 1
    if n == 0:
        return None
    return gemmi.Position(sx / n, sy / n, sz / n)


def pair_ssu_lsu_by_centroid(
    structure: gemmi.Structure,
    ssu_chains: list[ChainRef],
    lsu_chains: list[ChainRef],
) -> list[tuple[ChainRef, ChainRef]]:
    """Greedy nearest-centroid pairing of SSU and LSU chains.

    Each iteration picks the (SSU, LSU) pair with the smallest centroid
    distance from the remaining pools and removes both from contention.
    The result is deterministic given the input order. Chains whose
    centroid can't be computed (absent or empty) are skipped.

    Returns a list of ``(ssu_chain, lsu_chain)`` tuples of length
    ``min(len(ssu_chains), len(lsu_chains))`` minus any centroid-missing
    chains.
    """
    ssu_centroids: dict[str, gemmi.Position] = {}
    lsu_centroids: dict[str, gemmi.Position] = {}
    for c in ssu_chains:
        pt = _chain_centroid(structure, c.auth_asym_id)
        if pt is not None:
            ssu_centroids[c.auth_asym_id] = pt
    for c in lsu_chains:
        pt = _chain_centroid(structure, c.auth_asym_id)
        if pt is not None:
            lsu_centroids[c.auth_asym_id] = pt

    remaining_ssu = [c for c in ssu_chains if c.auth_asym_id in ssu_centroids]
    remaining_lsu = [c for c in lsu_chains if c.auth_asym_id in lsu_centroids]

    pairs: list[tuple[ChainRef, ChainRef]] = []
    while remaining_ssu and remaining_lsu:
        best: tuple[ChainRef, ChainRef] | None = None
        best_d = float("inf")
        for s in remaining_ssu:
            sc = ssu_centroids[s.auth_asym_id]
            for lc in remaining_lsu:
                lcent = lsu_centroids[lc.auth_asym_id]
                d = sc.dist(lcent)
                if d < best_d:
                    best_d = d
                    best = (s, lc)
        if best is None:
            break
        pairs.append(best)
        remaining_ssu = [c for c in remaining_ssu if c.auth_asym_id != best[0].auth_asym_id]
        remaining_lsu = [c for c in remaining_lsu if c.auth_asym_id != best[1].auth_asym_id]
    return pairs


def _ribosome_centroids(
    structure: gemmi.Structure,
    pairs: list[tuple[ChainRef, ChainRef]],
) -> list[gemmi.Position]:
    """Compute one combined centroid per ribosome (mean of SSU + LSU centroids)."""
    out: list[gemmi.Position] = []
    for ssu, lsu in pairs:
        sc = _chain_centroid(structure, ssu.auth_asym_id)
        lc = _chain_centroid(structure, lsu.auth_asym_id)
        if sc is None and lc is None:
            # Should not happen because pair_ssu_lsu_by_centroid filtered
            # these out, but guard against caller misuse.
            continue
        if sc is None:
            assert lc is not None
            out.append(lc)
            continue
        if lc is None:
            out.append(sc)
            continue
        out.append(gemmi.Position((sc.x + lc.x) / 2, (sc.y + lc.y) / 2, (sc.z + lc.z) / 2))
    return out


def partition_chains_by_ribosome(
    structure: gemmi.Structure,
    pairs: list[tuple[ChainRef, ChainRef]],
    chains: Iterable[ChainRef],
) -> list[set[str]]:
    """Assign each chain to one ribosome by centroid proximity.

    Returns a list of ``set[auth_asym_id]`` of the same length as
    ``pairs``. A chain whose centroid can't be computed is dropped from
    the assignment (it ends up in no ribosome).
    """
    rib_centroids = _ribosome_centroids(structure, pairs)
    if not rib_centroids:
        return [set() for _ in pairs]

    groups: list[set[str]] = [set() for _ in pairs]
    for chain in chains:
        pt = _chain_centroid(structure, chain.auth_asym_id)
        if pt is None:
            continue
        best_i = 0
        best_d = float("inf")
        for i, rc in enumerate(rib_centroids):
            d = pt.dist(rc)
            if d < best_d:
                best_d = d
                best_i = i
        groups[best_i].add(chain.auth_asym_id)
    return groups


def filter_correspondence_for_chains(
    correspondence_by_site: dict[str, CorrespondenceResult],
    allowed_chains: set[str],
) -> dict[str, CorrespondenceResult]:
    """Restrict each site's mapped units to those whose chain is in ``allowed_chains``.

    The returned :class:`CorrespondenceResult` instances preserve the
    original ``reference_units`` and ``warnings`` but their
    ``mapped_units`` / ``mapped_units_by_chain`` lists are filtered.
    Sites whose mapped_units become empty after filtering remain in the
    dict — downstream contact-transfer treats them as "no anchors" the
    same way it would for a top-level empty result.
    """
    out: dict[str, CorrespondenceResult] = {}
    for site_key, corr in correspondence_by_site.items():
        filtered: list[str] = []
        for unit in corr.mapped_units:
            parsed = try_parse_unit_id(unit)
            if parsed is None:
                continue
            if parsed.chain in allowed_chains:
                filtered.append(unit)
        out[site_key] = CorrespondenceResult(
            reference_key=corr.reference_key,
            reference_units=list(corr.reference_units),
            mapped_units=filtered,
            mapped_units_by_chain=group_by_chain(filtered),
            warnings=list(corr.warnings),
        )
    return out
