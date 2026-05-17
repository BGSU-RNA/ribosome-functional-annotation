"""tRNA-mRNA codon/anticodon base-pair extraction (codon-anticodon spec).

For each A/P/E-site tRNA in an annotated assembly, this module pulls the
FR3D-curated base-pair interactions from
``rna.bgsu.edu/rna3dhub/pdb/<pdb>/interactions/fr3d/basepairs/csv`` and
extracts the codon ↔ anticodon evidence:

- The three anticodon residues (tRNA biological positions 34, 35, 36,
  identified by polymer-sequence index in the mmCIF — so this works
  even for organellar mt-tRNAs whose author residue numbers don't run
  1..76).
- The mRNA codon residues that pair with them under FR3D, falling back
  to mmCIF-polymer-order reconstruction for codon positions that FR3D
  doesn't observe, then to mRNA-frame inference for whole sites whose
  codon couldn't be assigned directly.
- The raw FR3D interaction label (``cWW``, ``tHS``, etc.) for each
  observed pair.

The module is **evidence-only**: no cognate / near-cognate / non-cognate
classification, no Watson-Crick interpretation. Downstream consumers
read the raw pairs and decide.
"""

from __future__ import annotations

import csv
import io
import itertools
import logging
from dataclasses import dataclass
from typing import Literal

import gemmi
import httpx

from ribosome_state_annotator.cache import Cache
from ribosome_state_annotator.ccd_client import fetch_chem_comp
from ribosome_state_annotator.models import (
    Anticodon,
    AnticodonResidue,
    BasePair,
    ChainRef,
    Codon,
    CodonResidue,
    RibosomeAnnotation,
    TRNAmRNAInteraction,
)

logger = logging.getLogger(__name__)

FR3D_BASEPAIRS_URL_TEMPLATE = (
    "http://rna.bgsu.edu/rna3dhub/pdb/{pdb_id}/interactions/fr3d/basepairs/csv"
)
"""FR3D basepair CSV endpoint, by lowercase PDB ID."""

DEFAULT_HTTP_TIMEOUT = 30.0

# Biological convention for codon ↔ anticodon antiparallel pairing:
#   anticodon 34 ↔ codon position 3 (wobble)
#   anticodon 35 ↔ codon position 2
#   anticodon 36 ↔ codon position 1
_TRNA_POSITION_TO_CODON_POSITION = {34: 3, 35: 2, 36: 1}
_CODON_POSITION_TO_TRNA_POSITION = {3: 34, 2: 35, 1: 36}

SITES: tuple[Literal["A", "P", "E"], ...] = ("A", "P", "E")

# Minimum number of consecutive mRNA polymer residues needed before we
# attempt to infer P-site / E-site codons by stepping the A-site frame
# upstream. 3 codons * 3 nt = 9 residues end-to-end.
_MIN_MRNA_RUN_FOR_FRAME_INFERENCE = 9


# ---------------------------------------------------------------------------
# Public dataclasses for internal pipeline use
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class _ParsedFr3dRow:
    unit_id_1: str
    interaction: str
    unit_id_2: str


@dataclass(frozen=True)
class _AnticodonHit:
    """One picked anticodon residue: biological position + mmCIF residue."""

    trna_position: int  # 34 / 35 / 36
    residue: gemmi.Residue


# ---------------------------------------------------------------------------
# FR3D fetch + parse
# ---------------------------------------------------------------------------


def fetch_fr3d_basepairs(
    pdb_id: str,
    *,
    cache: Cache | None = None,
    client: httpx.Client | None = None,
) -> list[_ParsedFr3dRow] | None:
    """Return the FR3D basepair rows for ``pdb_id``, or ``None`` on failure.

    Cached on disk under the ``fr3d/`` namespace as the raw CSV bytes.
    Failures (network, non-200, parse) log a warning and return ``None``
    — the annotation pipeline degrades to an empty
    ``trna_mrna_interactions`` list.
    """
    if cache is not None:
        cached = cache.get_fr3d_csv(pdb_id)
        if cached is not None:
            logger.debug("fr3d cache hit for %s", pdb_id)
            return _parse_fr3d_csv(cached)
    url = FR3D_BASEPAIRS_URL_TEMPLATE.format(pdb_id=pdb_id.lower())
    http = client if client is not None else httpx.Client(timeout=DEFAULT_HTTP_TIMEOUT)
    try:
        logger.info("fetching FR3D basepairs for %s", pdb_id)
        try:
            response = http.get(url, follow_redirects=True)
        except httpx.HTTPError as exc:
            logger.warning("FR3D fetch failed for %s: %s", pdb_id, exc)
            return None
        if response.status_code != 200:
            logger.warning("FR3D returned HTTP %d for %s", response.status_code, pdb_id)
            return None
        body = response.content
    finally:
        if client is None:
            http.close()
    if cache is not None:
        cache.put_fr3d_csv(pdb_id, body)
    return _parse_fr3d_csv(body)


def _parse_fr3d_csv(data: bytes) -> list[_ParsedFr3dRow]:
    """Parse the 3-column FR3D CSV (``unit1, interaction, unit2``)."""
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as exc:
        logger.warning("FR3D CSV is not valid UTF-8 (%s); treating as empty", exc)
        return []
    reader = csv.reader(io.StringIO(text))
    rows: list[_ParsedFr3dRow] = []
    for raw in reader:
        if len(raw) < 3:
            continue
        u1, interaction, u2 = (cell.strip() for cell in raw[:3])
        if not u1 or not u2 or not interaction:
            continue
        rows.append(_ParsedFr3dRow(unit_id_1=u1, interaction=interaction, unit_id_2=u2))
    return rows


# ---------------------------------------------------------------------------
# Anticodon residue identification (polymer-sequence-index 34/35/36)
# ---------------------------------------------------------------------------


def _polymer_residues_in_order(structure: gemmi.Structure, chain_id: str) -> list[gemmi.Residue]:
    """Return the polymer residues of ``chain_id`` in mmCIF polymer order.

    HETATM ligands (Mg, HOH, modified-nt-as-ligand) are filtered out via
    ``EntityType.Polymer`` so they don't displace real polymer residues
    when counting positions 34/35/36. Sorted by ``seqid.num`` —
    insertion-coded residues (very rare in tRNAs) sort by their numeric
    base, which matches the way mmCIF polymer order is conventionally
    interpreted.
    """
    model = structure[0]
    try:
        chain = model[chain_id]
    except (KeyError, IndexError):
        return []
    polymer = [r for r in chain if r.entity_type == gemmi.EntityType.Polymer]
    polymer.sort(key=lambda r: (r.seqid.num if r.seqid.num is not None else -1, r.seqid.icode or ""))
    return polymer


def _pick_anticodon_residues(
    structure: gemmi.Structure,
    chain_id: str,
) -> list[_AnticodonHit] | None:
    """Return the polymer residues at auth_seq_id 34, 35, 36 in ``chain_id``.

    The pick is anchored on the canonical first residue (auth_seq_id 1)
    rather than on the 34th element of the polymer list. This matters
    when a deposited chain has a pre-residue numbered ``0`` (or negative)
    at the 5' end — e.g. 5UYM chain W, whose polymer runs 0..75, where
    counting "the 34th polymer residue" lands on auth_seq_id 33 (one
    too early) and misses the true Sprinzl-34 wobble residue.

    Returns ``None`` when the chain does not have a polymer residue at
    every one of auth_seq_id 34, 35, 36 (anticodon-stem-loop fragments,
    organellar mt-tRNAs that renumber, severely truncated chains).
    """
    residues = _polymer_residues_in_order(structure, chain_id)
    by_seqid = {r.seqid.num: r for r in residues if r.seqid.num is not None}
    hits: list[_AnticodonHit] = []
    for biological_position in (34, 35, 36):
        residue = by_seqid.get(biological_position)
        if residue is None:
            return None
        hits.append(_AnticodonHit(trna_position=biological_position, residue=residue))
    return hits


def _build_unit_id(pdb_id: str, chain_id: str, residue: gemmi.Residue) -> str:
    """Construct a BGSU/FR3D-style unit ID for a residue.

    Uses the residue's **author sequence number** (``seqid.num``) — FR3D
    unit IDs embed deposited residue numbers, not polymer indices. The
    chem_comp_id is the residue's observed name (so modified-nucleotide
    codes like ``PSU`` / ``5MU`` are preserved verbatim).
    """
    return f"{pdb_id.upper()}|1|{chain_id}|{residue.name}|{residue.seqid.num}"


def _parent_base_info(
    comp_id: str,
    *,
    cache: Cache | None = None,
    client: httpx.Client | None = None,
) -> tuple[str, bool]:
    """Return ``(parent_base_uppercase, is_modified)`` for a CCD code.

    Resolution order:

    1. **Gemmi tabulated dictionary** (in-process; covers
       canonical A/C/G/U and most common modifications: ``PSU`` → ``"u"``,
       ``5MU`` → ``"u"``, ``MIA`` → ``"a"``, etc.).
    2. **Per-component CCD fetch** when Gemmi returns blank / missing.
       The authoritative ``mon_nstd_parent_comp_id`` (or
       ``one_letter_code``) from the PDB CCD entry — covers unusual
       modifications that ship without a Gemmi one-letter code (e.g.
       ``U8U`` = 5-methylaminomethyl-2-thiouridine-5'-monophosphate,
       a *T. thermophilus* tRNA wobble modification). Cached under the
       ``ccd/`` namespace; one network call per unrecognized comp_id
       per cache lifetime.
    3. **First-character heuristic** as last resort: take
       ``comp_id[0]`` and flag ``is_modified=True`` if the comp_id is
       longer than one character. Works because PDB CCD names for
       modified nts almost universally start with the parent base
       letter (``PSU`` → U, ``7MG`` → G).

    Passing ``cache=None`` / ``client=None`` skips step 2 — useful
    for tests that don't want to hit the network.
    """
    info = gemmi.find_tabulated_residue(comp_id)
    letter = info.one_letter_code if info is not None else None
    if letter and letter.strip():
        return letter.upper(), letter.islower()

    if cache is not None or client is not None:
        ccd_info = fetch_chem_comp(comp_id, cache=cache, client=client)
        if ccd_info is not None:
            parent = ccd_info.parent_comp_id or ccd_info.one_letter_code
            if parent and parent.strip():
                clean = parent.strip()
                # If the parent string is a single canonical base, that base
                # is the parent and we're modified relative to it. If the
                # CCD's parent reference is itself a multi-char code (rare),
                # fall back to its first character.
                parent_letter = clean[:1].upper()
                is_modified = (
                    ccd_info.parent_comp_id is not None
                    and ccd_info.parent_comp_id.strip().upper() != comp_id.upper()
                )
                # Also treat a lowercase one_letter_code as modified, mirroring
                # Gemmi's convention.
                if ccd_info.parent_comp_id is None and clean.islower():
                    is_modified = True
                return parent_letter, is_modified

    if not comp_id:
        return "", False
    return comp_id[:1].upper(), len(comp_id) > 1


# ---------------------------------------------------------------------------
# Pair filtering + codon assembly
# ---------------------------------------------------------------------------


def _filter_cross_chain_pairs(
    fr3d_rows: list[_ParsedFr3dRow],
    *,
    mrna_chain_id: str,
    trna_chain_id: str,
    anticodon_auth_seq_ids: set[int],
) -> list[tuple[str, str, str]]:
    """Return deduplicated mRNA ↔ tRNA-anticodon pairs from FR3D.

    Each FR3D pair appears twice (once per direction). We canonicalise by
    sorting the two unit IDs so the same physical pair only appears once
    in the returned list. The returned tuples are always
    ``(mrna_unit_id, interaction, trna_unit_id)`` regardless of which
    direction the row had originally.
    """
    seen: set[tuple[str, str, str]] = set()
    out: list[tuple[str, str, str]] = []
    for row in fr3d_rows:
        sides = _classify_pair_sides(row, mrna_chain_id, trna_chain_id, anticodon_auth_seq_ids)
        if sides is None:
            continue
        mrna_uid, trna_uid = sides
        key = (mrna_uid, row.interaction, trna_uid)
        if key in seen:
            continue
        seen.add(key)
        out.append(key)
    return out


def _classify_pair_sides(
    row: _ParsedFr3dRow,
    mrna_chain_id: str,
    trna_chain_id: str,
    anticodon_auth_seq_ids: set[int],
) -> tuple[str, str] | None:
    """Return ``(mrna_unit_id, trna_unit_id)`` if the row crosses
    mRNA ↔ this-tRNA at one of the anticodon residues; else ``None``."""
    from ribosome_state_annotator.correspondence import parse_unit_id

    try:
        u1 = parse_unit_id(row.unit_id_1)
        u2 = parse_unit_id(row.unit_id_2)
    except ValueError:
        return None
    if u1.chain == mrna_chain_id and u2.chain == trna_chain_id and u2.residue_number in anticodon_auth_seq_ids:
        return row.unit_id_1, row.unit_id_2
    if u2.chain == mrna_chain_id and u1.chain == trna_chain_id and u1.residue_number in anticodon_auth_seq_ids:
        return row.unit_id_2, row.unit_id_1
    return None


def _assign_codon_from_pairs(
    cross_pairs: list[tuple[str, str, str]],
    anticodon_by_seq_id: dict[int, _AnticodonHit],
) -> tuple[dict[int, str], dict[int, str], dict[int, str], set[int]]:
    """For each codon position, pick the mRNA unit ID and FR3D label.

    Returns:
        codon_unit_id_by_position:  ``{codon_position: mrna_unit_id}``
        interaction_by_position:    ``{codon_position: fr3d_label}``
        codon_base_by_position:     ``{codon_position: mrna_base_letter}``
        ambiguous_positions:        codon positions where more than one
                                    mRNA residue paired and we couldn't
                                    pick a single cWW representative.

    Resolution rule (spec §5):
      - If multiple mRNA residues pair with the same anticodon position,
        prefer the cWW pairing.
      - If still ambiguous (no cWW, or multiple cWW), mark the codon
        position as ambiguous; assignment is left to the per-residue
        candidate the caller picks (we keep all candidates surfaced via
        the BasePair list, so consumers can still see them).
    """
    from ribosome_state_annotator.correspondence import parse_unit_id

    # Group pairs by anticodon position.
    pairs_by_trna_position: dict[int, list[tuple[str, str, str]]] = {}
    for mrna_uid, interaction, trna_uid in cross_pairs:
        try:
            trna_parsed = parse_unit_id(trna_uid)
        except ValueError:
            continue
        hit = anticodon_by_seq_id.get(trna_parsed.residue_number)
        if hit is None:
            continue
        pairs_by_trna_position.setdefault(hit.trna_position, []).append(
            (mrna_uid, interaction, trna_uid)
        )

    codon_unit_id_by_position: dict[int, str] = {}
    interaction_by_position: dict[int, str] = {}
    codon_base_by_position: dict[int, str] = {}
    ambiguous_positions: set[int] = set()

    for trna_position, candidates in pairs_by_trna_position.items():
        codon_position = _TRNA_POSITION_TO_CODON_POSITION[trna_position]
        if len(candidates) == 1:
            chosen = candidates[0]
        else:
            cww = [c for c in candidates if c[1] == "cWW"]
            if len(cww) == 1:
                chosen = cww[0]
            else:
                ambiguous_positions.add(codon_position)
                chosen = candidates[0]  # surface something so the pair list isn't empty
        mrna_uid, interaction, _ = chosen
        codon_unit_id_by_position[codon_position] = mrna_uid
        interaction_by_position[codon_position] = interaction
        try:
            codon_base_by_position[codon_position] = parse_unit_id(mrna_uid).residue_name
        except ValueError:
            continue

    return codon_unit_id_by_position, interaction_by_position, codon_base_by_position, ambiguous_positions


# ---------------------------------------------------------------------------
# Codon reconstruction (mmCIF polymer order) + frame inference
# ---------------------------------------------------------------------------


def _mrna_polymer_index(structure: gemmi.Structure, mrna_chain_id: str) -> tuple[
    list[gemmi.Residue], dict[int, int]
]:
    """Return ``(polymer_residues, auth_seq_id_to_index)`` for the mRNA chain."""
    residues = _polymer_residues_in_order(structure, mrna_chain_id)
    index = {r.seqid.num: i for i, r in enumerate(residues) if r.seqid.num is not None}
    return residues, index


def _reconstruct_codon_locally(
    pdb_id: str,
    mrna_chain_id: str,
    fr3d_codon_unit_ids: dict[int, str],
    mrna_residues: list[gemmi.Residue],
    mrna_index_by_seq_id: dict[int, int],
) -> dict[int, CodonResidue]:
    """Fill missing codon positions from mmCIF polymer order.

    If any codon position is FR3D-observed, the other two positions of
    the same codon triplet are inferred from the mmCIF polymer order:
    position 1 is two upstream of position 3, position 2 is one upstream
    of position 3, etc.

    Returns a ``{codon_position: CodonResidue}`` dict that contains:
      - The FR3D-observed residues (source="fr3d_observed")
      - Any reconstructed residues (source="mmcif_reconstructed")
    Positions that couldn't be filled are absent from the dict.
    """
    from ribosome_state_annotator.correspondence import parse_unit_id

    if not fr3d_codon_unit_ids:
        return {}

    by_position: dict[int, CodonResidue] = {}
    # First seed observed positions.
    for codon_position, unit_id in fr3d_codon_unit_ids.items():
        try:
            parsed = parse_unit_id(unit_id)
        except ValueError:
            continue
        by_position[codon_position] = CodonResidue(
            codon_position=codon_position,
            unit_id=unit_id,
            base=parsed.residue_name,
            source="fr3d_observed",
        )

    if not by_position:
        return {}

    # Pick the highest-codon-position seed (rightmost in the triplet) and
    # walk outward in polymer order. The mRNA is antiparallel to the
    # anticodon, so codon position 1 sits 5' of position 3. In mmCIF
    # polymer order, that means position 1 is 2 residues *before*
    # position 3 (smaller seq index).
    seed_codon_position = max(by_position.keys())
    seed_uid = by_position[seed_codon_position].unit_id
    seed_resnum = parse_unit_id(seed_uid).residue_number
    seed_polymer_idx = mrna_index_by_seq_id.get(seed_resnum)
    if seed_polymer_idx is None:
        return by_position

    # The polymer-index of codon position 1 relative to seed:
    #   codon_pos 1 ↔ seed_polymer_idx - (seed_codon_position - 1)
    base_polymer_idx = seed_polymer_idx - (seed_codon_position - 1)
    for codon_position in (1, 2, 3):
        if codon_position in by_position:
            continue
        polymer_idx = base_polymer_idx + (codon_position - 1)
        if not (0 <= polymer_idx < len(mrna_residues)):
            continue
        residue = mrna_residues[polymer_idx]
        unit_id = _build_unit_id(pdb_id, mrna_chain_id, residue)
        by_position[codon_position] = CodonResidue(
            codon_position=codon_position,
            unit_id=unit_id,
            base=residue.name,
            source="mmcif_reconstructed",
        )
    return by_position


# Offset in codons between sites along the mRNA, read 5' → 3'.
# E is most upstream (smallest auth_seq_id), then P, then A (most downstream).
# So inferring P from A is one codon upstream (offset -3 residues);
# inferring A from E is two codons downstream (offset +6 residues); etc.
_SITE_FRAME_INDEX: dict[Literal["A", "P", "E"], int] = {"E": 0, "P": 1, "A": 2}


def _frame_inferred_codon(
    pdb_id: str,
    mrna_chain_id: str,
    anchor_site: Literal["A", "P", "E"],
    anchor_codon_residues: dict[int, CodonResidue],
    target_site: Literal["A", "P", "E"],
    mrna_residues: list[gemmi.Residue],
    mrna_index_by_seq_id: dict[int, int],
) -> dict[int, CodonResidue] | None:
    """Infer the ``target_site`` codon by stepping the mRNA frame from
    a fully-resolved ``anchor_site`` codon.

    Returns ``None`` if the anchor codon is incomplete, the target is
    the same as the anchor, or the mRNA polymer doesn't have enough
    upstream / downstream residues to span both codons as a single
    gap-free run.

    Inference geometry: codon position 1 of any site sits at the same
    mRNA polymer position as codon position 1 of any other site, offset
    by 3 residues per codon. With E < P < A along the mRNA (5' → 3'),
    the polymer-index offset from anchor codon-position-1 to target
    codon-position-1 is ``3 * (target_index - anchor_index)``.
    """
    from ribosome_state_annotator.correspondence import parse_unit_id

    if target_site == anchor_site:
        return None
    if len(anchor_codon_residues) < 3:
        return None
    anchor_pos1 = anchor_codon_residues.get(1)
    if anchor_pos1 is None:
        return None
    try:
        anchor_pos1_seqnum = parse_unit_id(anchor_pos1.unit_id).residue_number
    except ValueError:
        return None
    anchor_pos1_polymer_idx = mrna_index_by_seq_id.get(anchor_pos1_seqnum)
    if anchor_pos1_polymer_idx is None:
        return None

    codon_offset = _SITE_FRAME_INDEX[target_site] - _SITE_FRAME_INDEX[anchor_site]
    base_polymer_idx = anchor_pos1_polymer_idx + 3 * codon_offset
    if base_polymer_idx < 0:
        return None
    if base_polymer_idx + 2 >= len(mrna_residues):
        return None

    # Verify the span from the leftmost involved residue to the rightmost
    # is one consecutive polymer run by auth_seq_id. If there's a gap,
    # frame inference is unreliable.
    span_first = min(base_polymer_idx, anchor_pos1_polymer_idx)
    span_last = max(base_polymer_idx + 2, anchor_pos1_polymer_idx + 2)
    span_residues = mrna_residues[span_first : span_last + 1]
    if not _is_consecutive_run(span_residues):
        return None

    out: dict[int, CodonResidue] = {}
    for codon_position in (1, 2, 3):
        polymer_idx = base_polymer_idx + (codon_position - 1)
        residue = mrna_residues[polymer_idx]
        unit_id = _build_unit_id(pdb_id, mrna_chain_id, residue)
        out[codon_position] = CodonResidue(
            codon_position=codon_position,
            unit_id=unit_id,
            base=residue.name,
            source="mrna_frame_inference",
        )
    return out


def _is_consecutive_run(residues: list[gemmi.Residue]) -> bool:
    """True iff residues form a gap-free auth_seq_id+1 run."""
    if len(residues) <= 1:
        return True
    for prev, curr in itertools.pairwise(residues):
        if curr.seqid.num is None or prev.seqid.num is None:
            return False
        if curr.seqid.num - prev.seqid.num != 1:
            return False
    return True


# ---------------------------------------------------------------------------
# Public extraction entry point
# ---------------------------------------------------------------------------


_RFAM_TRNA = "RF00005"

_SITE_ORDER: tuple[Literal["A", "P", "E"], ...] = ("A", "P", "E")


def fr3d_codon_pairing_fallback(
    annotation: RibosomeAnnotation,
    structure: gemmi.Structure,
    fr3d_rows: list[_ParsedFr3dRow],
) -> None:
    """Fill empty A / P / E tRNA slots from FR3D codon-pairing evidence.

    Contact-transfer (:func:`infer.assign_functional_chains`) compares
    candidate tRNA chains against the canonical SSU decoding-centre
    monitor bases (E. coli A1492/A1493 → yeast 18S A1824/A1825 → etc.).
    Those monitor bases only flip out and contact the codon-anticodon
    duplex once the aa-tRNA has been **accommodated** — in
    pre-accommodation states (e.g. eEF1A-GTP-aa-tRNA decoding
    intermediates such as 3JAG) the tRNA's anticodon is base-paired to
    the mRNA codon but the rRNA fingerprint hasn't formed yet, so
    contact-transfer correctly refuses to claim the chain.

    This fallback closes the gap. For every unassigned tRNA-Rfam chain
    that makes at least one cWW pair to the mRNA at anticodon positions
    34/35/36, we determine **which mRNA codon** it engages (by reading
    the mRNA residues' ``auth_seq_id``). The mRNA is translated 5' → 3':

    - A-site codon = most 3' (largest auth_seq_id)
    - P-site codon = middle (one codon = 3 residues 5' of A)
    - E-site codon = most 5' (two codons 5' of A)

    So candidates are sorted by mRNA codon position descending, and
    assigned in that order to the remaining empty A / P / E slots. When
    other sites are already filled by contact-transfer their codon
    positions (if any FR3D pairs exist) are read first as frame anchors,
    so the fallback respects the existing frame.

    The annotation's ``other_rna_chains`` is updated and warnings record
    each FR3D-derived assignment so consumers can distinguish it from
    canonical contact-transfer assignments.

    Mutates ``annotation`` in place; no return value.
    """
    if annotation.mrna_chain is None:
        return

    empty_slots: list[Literal["A", "P", "E"]] = []
    if annotation.aminoacyl_trna_chain is None:
        empty_slots.append("A")
    if annotation.peptidyl_trna_chain is None:
        empty_slots.append("P")
    if annotation.exit_trna_chain is None:
        empty_slots.append("E")
    if not empty_slots:
        return

    mrna_chain_id = annotation.mrna_chain.auth_asym_id
    assigned_ifes: set[str] = {
        chain.ife
        for chain in (
            annotation.aminoacyl_trna_chain,
            annotation.peptidyl_trna_chain,
            annotation.exit_trna_chain,
        )
        if chain is not None
    }
    # Candidates: tRNA-Rfam chains currently surfaced as `other_rna_chains`
    # (i.e. not already placed in any role).
    trna_candidates: list[ChainRef] = [
        chain
        for chain in annotation.other_rna_chains
        if _RFAM_TRNA in chain.rfam_accessions and chain.ife not in assigned_ifes
    ]
    if not trna_candidates:
        return

    @dataclass
    class _Candidate:
        chain: ChainRef
        codon_position: float  # mean auth_seq_id across mRNA codon residues
        cww_count: int

    scored: list[_Candidate] = []
    for cand in trna_candidates:
        codon_position, cww_count = _score_candidate_for_fallback(
            cand.auth_asym_id, mrna_chain_id, structure, fr3d_rows
        )
        if codon_position is None:
            continue  # no cWW pair to mRNA — skip
        scored.append(_Candidate(chain=cand, codon_position=codon_position, cww_count=cww_count))

    if not scored:
        return

    # Sort A-side first: most-downstream mRNA codon position → most likely A-site.
    scored.sort(key=lambda c: -c.codon_position)
    # Empty slots in canonical A → P → E order so the most-downstream
    # candidate goes to the most-A-side empty slot.
    empty_slots.sort(key=_SITE_ORDER.index)

    newly_assigned: list[tuple[Literal["A", "P", "E"], ChainRef, int]] = []
    for cand_entry, site in zip(scored, empty_slots, strict=False):
        if site == "A":
            annotation.aminoacyl_trna_chain = cand_entry.chain
        elif site == "P":
            annotation.peptidyl_trna_chain = cand_entry.chain
        else:
            annotation.exit_trna_chain = cand_entry.chain
        newly_assigned.append((site, cand_entry.chain, cand_entry.cww_count))
        logger.info(
            "FR3D fallback assigned %s-tRNA = chain %s "
            "(mRNA codon centroid auth_seq_id %.1f, %d cWW pair(s))",
            site,
            cand_entry.chain.auth_asym_id,
            cand_entry.codon_position,
            cand_entry.cww_count,
        )

    if not newly_assigned:
        return

    promoted_ifes = {chain.ife for _, chain, _ in newly_assigned}
    annotation.other_rna_chains = [
        chain for chain in annotation.other_rna_chains if chain.ife not in promoted_ifes
    ]
    for site, chain, _ in newly_assigned:
        annotation.warnings.append(
            f"{site.lower()}trna_assigned_from_fr3d_codon_pairing_{chain.auth_asym_id}"
        )
    annotation.classification_evidence.setdefault("fr3d_codon_pairing_fallback", {})
    for site, chain, _ in newly_assigned:
        annotation.classification_evidence["fr3d_codon_pairing_fallback"][site] = chain.auth_asym_id


def _score_candidate_for_fallback(
    trna_chain_id: str,
    mrna_chain_id: str,
    structure: gemmi.Structure,
    fr3d_rows: list[_ParsedFr3dRow],
) -> tuple[float | None, int]:
    """Return ``(mRNA codon position, cWW pair count)`` for a candidate.

    Returns ``(None, 0)`` if the chain has fewer than 36 polymer
    residues or makes no cWW pair to the mRNA.

    "Codon position" is the mean ``auth_seq_id`` of the mRNA residues
    that pair (under any FR3D interaction, not just cWW) with the
    candidate's anticodon residues. Used to disambiguate A vs P vs E
    when multiple candidates need to be placed: 5' → 3' on the mRNA
    corresponds to E → P → A in the ribosome.
    """
    from ribosome_state_annotator.correspondence import parse_unit_id

    anticodon_hits = _pick_anticodon_residues(structure, trna_chain_id)
    if anticodon_hits is None:
        return None, 0
    anticodon_auth_seq_ids = {
        hit.residue.seqid.num for hit in anticodon_hits if hit.residue.seqid.num is not None
    }
    cross_pairs = _filter_cross_chain_pairs(
        fr3d_rows,
        mrna_chain_id=mrna_chain_id,
        trna_chain_id=trna_chain_id,
        anticodon_auth_seq_ids=anticodon_auth_seq_ids,
    )
    cww_count = sum(1 for _, interaction, _ in cross_pairs if interaction == "cWW")
    if cww_count == 0:
        return None, 0

    mrna_residue_numbers: list[int] = []
    for mrna_uid, _, _ in cross_pairs:
        try:
            mrna_residue_numbers.append(parse_unit_id(mrna_uid).residue_number)
        except ValueError:
            continue
    if not mrna_residue_numbers:
        return None, 0
    return sum(mrna_residue_numbers) / len(mrna_residue_numbers), cww_count


def extract_trna_mrna_interactions(
    annotation: RibosomeAnnotation,
    structure: gemmi.Structure,
    *,
    cache: Cache | None = None,
    client: httpx.Client | None = None,
) -> list[TRNAmRNAInteraction]:
    """Extract codon/anticodon evidence for every populated A/P/E site.

    Before the per-site extraction, this also runs the FR3D
    codon-pairing fallback: if a site (A or E) was left unassigned by
    contact-transfer but FR3D shows an unclaimed tRNA-Rfam chain making
    cWW codon-anticodon pairs to the mRNA, that chain is assigned to
    the empty site (see :func:`fr3d_codon_pairing_fallback`). This
    rescues pre-accommodation A-site states (e.g. 3JAG, where the
    decoding-centre monitor bases haven't flipped out and so the
    canonical SSU anchors don't contact the aa-tRNA).

    Returns an empty list if the run condition fails (no mRNA chain, or
    no A/P/E tRNA assigned and no fallback hit), if FR3D is unreachable,
    or if every site fails to produce evidence. Per-site errors degrade
    to a per-site warning rather than raising.
    """
    if annotation.mrna_chain is None:
        logger.info("skipping codon/anticodon extraction: no mRNA chain")
        return []

    pdb_id = annotation.pdb_id

    # Skip FR3D fetch entirely if there's nothing to extract or to rescue:
    # no A/P/E tRNA already assigned and no unassigned tRNA-Rfam candidate
    # for the codon-pairing fallback.
    has_assigned_trna = any(
        chain is not None
        for chain in (
            annotation.aminoacyl_trna_chain,
            annotation.peptidyl_trna_chain,
            annotation.exit_trna_chain,
        )
    )
    has_trna_rfam_candidate = any(
        _RFAM_TRNA in chain.rfam_accessions for chain in annotation.other_rna_chains
    )
    if not (has_assigned_trna or has_trna_rfam_candidate):
        logger.info(
            "skipping codon/anticodon extraction: no A/P/E-site tRNA assigned "
            "and no tRNA-Rfam candidate for fallback"
        )
        return []

    fr3d_rows = fetch_fr3d_basepairs(pdb_id, cache=cache, client=client)
    if fr3d_rows is None:
        logger.warning("FR3D unavailable for %s; codon/anticodon evidence will be empty", pdb_id)
        return []

    # Pre-extraction fallback: fill in any A/P/E site that contact-transfer
    # missed, using FR3D codon-anticodon evidence.
    fr3d_codon_pairing_fallback(annotation, structure, fr3d_rows)

    site_trna: dict[Literal["A", "P", "E"], ChainRef | None] = {
        "A": annotation.aminoacyl_trna_chain,
        "P": annotation.peptidyl_trna_chain,
        "E": annotation.exit_trna_chain,
    }
    if not any(v is not None for v in site_trna.values()):
        logger.info("skipping codon/anticodon extraction: no A/P/E-site tRNA assigned")
        return []

    mrna_chain_id = annotation.mrna_chain.auth_asym_id
    present_sites = [site for site, ref in site_trna.items() if ref is not None]
    logger.info(
        "codon/anticodon extraction: mRNA=%s, tRNA sites=%s",
        mrna_chain_id,
        ",".join(present_sites),
    )

    mrna_residues, mrna_index_by_seq_id = _mrna_polymer_index(structure, mrna_chain_id)
    mrna_consecutive_runs = _consecutive_run_lengths(mrna_residues)
    has_long_run = any(length >= _MIN_MRNA_RUN_FOR_FRAME_INFERENCE for length in mrna_consecutive_runs)

    # First pass: extract every site directly from FR3D + local reconstruction.
    per_site_interactions: dict[Literal["A", "P", "E"], TRNAmRNAInteraction] = {}
    for site in SITES:
        trna_ref = site_trna[site]
        if trna_ref is None:
            continue
        interaction = _extract_one_site(
            site=site,
            pdb_id=pdb_id,
            mrna_chain_id=mrna_chain_id,
            trna_chain_id=trna_ref.auth_asym_id,
            structure=structure,
            fr3d_rows=fr3d_rows,
            mrna_residues=mrna_residues,
            mrna_index_by_seq_id=mrna_index_by_seq_id,
            cache=cache,
            client=client,
        )
        if interaction is not None:
            per_site_interactions[site] = interaction

    # Second pass: frame-inference for any site whose codon is still
    # incomplete, anchored on whichever site IS complete. Priority for
    # the anchor: A → P → E. This generalises the A-anchored inference
    # so cases like 5UYN (only P-site has FR3D evidence; A is held by
    # EF-Tu pre-accommodation, E is post-translocation) still resolve
    # the A and E codons by stepping the mRNA frame from the resolved
    # P-codon.
    anchor_site: Literal["A", "P", "E"] | None = None
    for candidate_site in SITES:
        interaction = per_site_interactions.get(candidate_site)
        if interaction is not None and interaction.codon.assignment_status == "complete":
            anchor_site = candidate_site
            break

    if anchor_site is not None and has_long_run:
        anchor_codon_residues = {
            r.codon_position: r
            for r in per_site_interactions[anchor_site].codon.residues
        }
        for site in SITES:
            if site == anchor_site:
                continue
            existing = per_site_interactions.get(site)
            if existing is None or existing.codon.assignment_status == "complete":
                continue
            inferred = _frame_inferred_codon(
                pdb_id,
                mrna_chain_id,
                anchor_site,
                anchor_codon_residues,
                site,
                mrna_residues,
                mrna_index_by_seq_id,
            )
            if inferred is None:
                continue
            existing.codon.residues = _sorted_codon_residues(inferred)
            existing.codon.assignment_status = "complete"
            existing.codon.sequence = _codon_sequence(existing.codon.residues)
    elif not has_long_run:
        # If frame inference would have helped but the mRNA is too short,
        # surface that as a warning on any site we couldn't complete.
        for site in SITES:
            existing = per_site_interactions.get(site)
            if existing is not None and existing.codon.assignment_status != "complete":
                existing.warnings.append("insufficient_mrna_for_frame_inference")

    # Preserve A → P → E order.
    return [per_site_interactions[site] for site in SITES if site in per_site_interactions]


def _extract_one_site(
    *,
    site: Literal["A", "P", "E"],
    pdb_id: str,
    mrna_chain_id: str,
    trna_chain_id: str,
    structure: gemmi.Structure,
    fr3d_rows: list[_ParsedFr3dRow],
    mrna_residues: list[gemmi.Residue],
    mrna_index_by_seq_id: dict[int, int],
    cache: Cache | None = None,
    client: httpx.Client | None = None,
) -> TRNAmRNAInteraction | None:
    """Build one site's :class:`TRNAmRNAInteraction`. Returns ``None`` if
    the anticodon residues cannot be picked (tRNA chain too short)."""
    anticodon_hits = _pick_anticodon_residues(structure, trna_chain_id)
    if anticodon_hits is None:
        logger.warning(
            "tRNA chain %s for site %s has fewer than 36 polymer residues; skipping",
            trna_chain_id,
            site,
        )
        return None

    anticodon_residues: list[AnticodonResidue] = []
    anticodon_by_seq_id: dict[int, _AnticodonHit] = {}
    anticodon_auth_seq_ids: set[int] = set()
    parent_bases_by_position: dict[int, str] = {}
    is_modified_by_position: dict[int, bool] = {}
    chem_comp_id_by_position: dict[int, str] = {}
    for hit in anticodon_hits:
        residue = hit.residue
        parent, is_modified = _parent_base_info(residue.name, cache=cache, client=client)
        parent_bases_by_position[hit.trna_position] = parent
        is_modified_by_position[hit.trna_position] = is_modified
        chem_comp_id_by_position[hit.trna_position] = residue.name
        unit_id = _build_unit_id(pdb_id, trna_chain_id, residue)
        anticodon_residues.append(
            AnticodonResidue(
                trna_position=hit.trna_position,
                unit_id=unit_id,
                parent_base=parent,
                trna_chem_comp_id=residue.name,
                is_modified=is_modified,
            )
        )
        if residue.seqid.num is not None:
            anticodon_by_seq_id[residue.seqid.num] = hit
            anticodon_auth_seq_ids.add(residue.seqid.num)

    sequence_parent = "".join(parent_bases_by_position.get(p, "?") for p in (34, 35, 36))
    anticodon = Anticodon(
        sequence_parent=sequence_parent if "?" not in sequence_parent else None,
        residues=anticodon_residues,
    )

    cross_pairs = _filter_cross_chain_pairs(
        fr3d_rows,
        mrna_chain_id=mrna_chain_id,
        trna_chain_id=trna_chain_id,
        anticodon_auth_seq_ids=anticodon_auth_seq_ids,
    )
    logger.info(
        "site %s (tRNA %s): %d mRNA-anticodon FR3D pair(s) retained",
        site,
        trna_chain_id,
        len(cross_pairs),
    )

    # interaction_by_position / codon_base_by_position aren't consumed here —
    # the BasePair list below carries the same information per-pair. We keep
    # the four-tuple return shape because tests assert on it.
    (
        codon_unit_id_by_position,
        _,
        _,
        ambiguous_positions,
    ) = _assign_codon_from_pairs(cross_pairs, anticodon_by_seq_id)

    codon_residues_by_position = _reconstruct_codon_locally(
        pdb_id,
        mrna_chain_id,
        codon_unit_id_by_position,
        mrna_residues,
        mrna_index_by_seq_id,
    )

    # BasePair objects from FR3D-observed pairs.
    pairs = _build_basepairs(
        cross_pairs,
        anticodon_by_seq_id,
        parent_bases_by_position,
        is_modified_by_position,
        chem_comp_id_by_position,
        ambiguous_positions,
    )

    codon_residues_sorted = _sorted_codon_residues(codon_residues_by_position)
    assignment_status = _codon_assignment_status(codon_residues_sorted)
    sequence = _codon_sequence(codon_residues_sorted) if assignment_status == "complete" else None
    codon = Codon(
        sequence=sequence,
        assignment_status=assignment_status,
        residues=codon_residues_sorted,
    )

    warnings: list[str] = []
    if ambiguous_positions:
        for pos in sorted(ambiguous_positions):
            warnings.append(f"ambiguous_codon_assignment_position_{pos}")
            logger.warning(
                "site %s codon position %d is ambiguous: multiple mRNA FR3D pairs", site, pos
            )
    if assignment_status == "missing":
        warnings.append("no_fr3d_codon_interactions")
        logger.warning("site %s: no mRNA-anticodon FR3D interactions found", site)

    return TRNAmRNAInteraction(
        site=site,
        mrna_chain_id=mrna_chain_id,
        trna_chain_id=trna_chain_id,
        anticodon_position_source="polymer_sequence_index",
        codon=codon,
        anticodon=anticodon,
        pairs=pairs,
        warnings=warnings,
    )


def _build_basepairs(
    cross_pairs: list[tuple[str, str, str]],
    anticodon_by_seq_id: dict[int, _AnticodonHit],
    parent_bases_by_position: dict[int, str],
    is_modified_by_position: dict[int, bool],
    chem_comp_id_by_position: dict[int, str],
    ambiguous_positions: set[int],
) -> list[BasePair]:
    from ribosome_state_annotator.correspondence import parse_unit_id

    out: list[BasePair] = []
    for mrna_uid, interaction, trna_uid in cross_pairs:
        try:
            trna_parsed = parse_unit_id(trna_uid)
            mrna_parsed = parse_unit_id(mrna_uid)
        except ValueError:
            continue
        hit = anticodon_by_seq_id.get(trna_parsed.residue_number)
        if hit is None:
            continue
        trna_position = hit.trna_position
        codon_position = _TRNA_POSITION_TO_CODON_POSITION[trna_position]
        codon_base = mrna_parsed.residue_name
        parent_base = parent_bases_by_position.get(trna_position, "")
        assignment_status: Literal["assigned", "ambiguous"] = (
            "ambiguous" if codon_position in ambiguous_positions else "assigned"
        )
        out.append(
            BasePair(
                codon_position=codon_position,
                trna_position=trna_position,
                codon_unit_id=mrna_uid,
                trna_unit_id=trna_uid,
                codon_base=codon_base,
                trna_parent_base=parent_base,
                trna_chem_comp_id=chem_comp_id_by_position.get(trna_position, ""),
                trna_is_modified=is_modified_by_position.get(trna_position, False),
                fr3d_interaction=interaction,
                basepair=f"{codon_base}-{parent_base}",
                is_wobble_position=(trna_position == 34),
                assignment_status=assignment_status,
            )
        )
    # Sort by (codon_position, trna_position, interaction) for stable JSON.
    out.sort(key=lambda p: (p.codon_position, p.trna_position, p.fr3d_interaction))
    return out


def _sorted_codon_residues(by_position: dict[int, CodonResidue]) -> list[CodonResidue]:
    return [by_position[p] for p in sorted(by_position.keys())]


def _codon_assignment_status(residues: list[CodonResidue]) -> Literal["complete", "partial", "missing"]:
    positions = {r.codon_position for r in residues}
    if positions == {1, 2, 3}:
        return "complete"
    if not positions:
        return "missing"
    return "partial"


def _codon_sequence(residues: list[CodonResidue]) -> str:
    """3-letter codon string, parent bases. Caller must ensure complete."""
    by_position = {r.codon_position: r for r in residues}
    out_chars: list[str] = []
    for pos in (1, 2, 3):
        res = by_position.get(pos)
        if res is None:
            return ""
        parent, _ = _parent_base_info(res.base)
        out_chars.append(parent or res.base[:1].upper())
    return "".join(out_chars)


def _consecutive_run_lengths(residues: list[gemmi.Residue]) -> list[int]:
    """Length of each gap-free auth_seq_id run in the polymer."""
    if not residues:
        return []
    runs: list[int] = []
    current = 1
    for prev, curr in itertools.pairwise(residues):
        if (
            prev.seqid.num is not None
            and curr.seqid.num is not None
            and curr.seqid.num - prev.seqid.num == 1
        ):
            current += 1
        else:
            runs.append(current)
            current = 1
    runs.append(current)
    return runs


__all__ = [
    "FR3D_BASEPAIRS_URL_TEMPLATE",
    "extract_trna_mrna_interactions",
    "fetch_fr3d_basepairs",
    "fr3d_codon_pairing_fallback",
]
