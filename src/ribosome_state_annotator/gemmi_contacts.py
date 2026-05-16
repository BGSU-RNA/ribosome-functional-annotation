"""Gemmi-based residue lookup and neighbour search (spec §10).

This module owns three concerns:

- :func:`find_residue` / :func:`find_residues`: look up the
  :class:`gemmi.Residue` for one or more :class:`UnitId` instances by
  walking model 0 of the structure. Per spec §10.4, only the first model
  is read; multi-model structures log a warning and use model 0 anyway.
- :func:`find_neighbouring_chains`: the core contact-detection function.
  Given a set of target residues (grouped by their owning chain),
  return ``{chain_name: min_distance}`` for every other chain in the
  structure that has at least one atom within ``cutoff`` of any target
  atom (spec §10.2).

The package's higher layers (step 9 :mod:`infer`) wrap this with
biological semantics — "which RNA chain neighbours the SSU mRNA anchor?"
— but everything biological is the caller's concern. This module is
pure geometry.
"""

from __future__ import annotations

import logging
from collections.abc import Iterable, Mapping
from typing import cast

import gemmi

from ribosome_state_annotator.constants import DEFAULT_CONTACT_CUTOFF_ANGSTROM
from ribosome_state_annotator.correspondence import UnitId

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Residue lookup
# ---------------------------------------------------------------------------


def find_residue(
    structure: gemmi.Structure,
    unit_id: UnitId,
) -> tuple[str, gemmi.Residue] | None:
    """Locate the :class:`gemmi.Residue` for one :class:`UnitId`.

    Returns ``(chain_name, residue)`` on success, ``None`` if the chain
    or residue is absent.

    Per spec §10.4:

    - Only model 0 is read (biological-assembly mmCIFs use chain
      duplication, not model duplication, since the May 2022 wwPDB
      format change).
    - Chains are looked up by author asym ID (the chain segment of the
      BGSU unit ID).
    - Residues are matched on author sequence ID. Residue-name
      mismatches are tolerated (modified-nucleotide annotations sometimes
      differ between BGSU and the mmCIF chemical component dictionary);
      we log at DEBUG and return the residue anyway.
    - If the structure has more than one model, a WARNING is logged and
      we use model 0 only.
    """
    if len(structure) == 0:
        return None
    if len(structure) > 1:
        logger.warning(
            "structure %r has %d models; using only model 0 per spec §10.4",
            structure.name,
            len(structure),
        )

    model = structure[0]
    # gemmi's stub declares find_chain → Chain but at runtime it returns
    # None when the chain is absent. Cast to Optional so mypy lets us
    # check for None.
    chain = cast(gemmi.Chain | None, model.find_chain(unit_id.chain))
    if chain is None:
        return None

    for residue in chain:
        if residue.seqid.num == unit_id.residue_number:
            if unit_id.residue_name and residue.name != unit_id.residue_name:
                logger.debug(
                    "residue %s|%d name mismatch: BGSU=%s mmCIF=%s",
                    unit_id.chain,
                    unit_id.residue_number,
                    unit_id.residue_name,
                    residue.name,
                )
            return (chain.name, residue)
    return None


def find_residues(
    structure: gemmi.Structure,
    unit_ids: Iterable[UnitId],
) -> tuple[dict[str, list[gemmi.Residue]], list[str]]:
    """Bulk lookup. Returns ``(residues_by_chain, missing_base_unit_ids)``.

    ``residues_by_chain`` is keyed by author chain name and maps to the
    list of residues found in that chain, preserving input order.
    ``missing_base_unit_ids`` is the list of canonical 5-segment unit ID
    strings that could not be located — per spec §10.4, missing residues
    are not errors; the caller surfaces them as warnings.
    """
    residues_by_chain: dict[str, list[gemmi.Residue]] = {}
    missing: list[str] = []
    for unit_id in unit_ids:
        result = find_residue(structure, unit_id)
        if result is None:
            missing.append(unit_id.base_unit_id)
            continue
        chain_name, residue = result
        residues_by_chain.setdefault(chain_name, []).append(residue)
    return residues_by_chain, missing


# ---------------------------------------------------------------------------
# Neighbour search
# ---------------------------------------------------------------------------


def find_neighbouring_chains(
    structure: gemmi.Structure,
    target_residues_by_chain: Mapping[str, Iterable[gemmi.Residue]],
    *,
    cutoff: float = DEFAULT_CONTACT_CUTOFF_ANGSTROM,
    candidate_chains: Iterable[str] | None = None,
) -> dict[str, float]:
    """Find every other chain with an atom within ``cutoff`` Å of any target atom.

    Returns ``{chain_name: min_distance_angstrom}``. Chains in
    ``target_residues_by_chain`` are automatically **excluded** from the
    result — the intent is "which OTHER chains contact our targets?".

    Args:
        structure: A parsed :class:`gemmi.Structure`. Must have a unit cell
            (RCSB mmCIFs always do; the test fixtures set a P1 cell).
        target_residues_by_chain: Map of owning chain name to the residues
            on that chain whose neighbourhood we want.
        cutoff: Distance cutoff in angstroms (spec §10.2 default 5.0).
        candidate_chains: Optional whitelist of chains to consider as
            neighbours. When given, every chain in the result will be in
            this set (per spec §10.5: restrict to the current assembly
            even if the loaded mmCIF has extra chains). When omitted, all
            non-target chains in model 0 are considered.

    Returns:
        Empty dict if there are no targets, no atoms within ``cutoff``,
        or all neighbours fall into the excluded target chains.
    """
    if not target_residues_by_chain:
        return {}
    if len(structure) == 0:
        return {}
    model = structure[0]

    exclude_chains = set(target_residues_by_chain.keys())
    allowed_chains: set[str] | None = (
        set(candidate_chains) if candidate_chains is not None else None
    )

    neighbour_search = gemmi.NeighborSearch(model, structure.cell, cutoff).populate()

    min_distances: dict[str, float] = {}
    for residues in target_residues_by_chain.values():
        for residue in residues:
            for atom in residue:
                for mark in neighbour_search.find_atoms(atom.pos, "\0", radius=cutoff):
                    chain = model[mark.chain_idx]
                    chain_name = chain.name
                    if chain_name in exclude_chains:
                        continue
                    if allowed_chains is not None and chain_name not in allowed_chains:
                        continue
                    neighbour_atom = chain[mark.residue_idx][mark.atom_idx]
                    distance = atom.pos.dist(neighbour_atom.pos)
                    # PBC / search-radius rounding can return hits whose
                    # straight-line distance is marginally above the cutoff.
                    if distance > cutoff:
                        continue
                    prev = min_distances.get(chain_name)
                    if prev is None or distance < prev:
                        min_distances[chain_name] = distance
    return min_distances


__all__ = [
    "find_neighbouring_chains",
    "find_residue",
    "find_residues",
]
