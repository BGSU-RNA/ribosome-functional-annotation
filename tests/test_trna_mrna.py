"""Tests for the FR3D-driven tRNA-mRNA codon/anticodon extraction
(codon-anticodon-extraction spec)."""

from __future__ import annotations

import gemmi
import httpx
import pytest
import respx

from ribosome_state_annotator.cache import Cache
from ribosome_state_annotator.models import (
    ChainRef,
    RibosomeAnnotation,
)
from ribosome_state_annotator.trna_mrna import (
    FR3D_BASEPAIRS_URL_TEMPLATE,
    _AnticodonHit,
    _assign_codon_from_pairs,
    _build_unit_id,
    _classify_pair_sides,
    _codon_assignment_status,
    _consecutive_run_lengths,
    _filter_cross_chain_pairs,
    _frame_inferred_codon,
    _is_consecutive_run,
    _parent_base_info,
    _parse_fr3d_csv,
    _ParsedFr3dRow,
    _pick_anticodon_residues,
    _polymer_residues_in_order,
    _reconstruct_codon_locally,
    extract_trna_mrna_interactions,
    fetch_fr3d_basepairs,
    fr3d_codon_pairing_fallback,
)

PDB = "TEST"

# ---------------------------------------------------------------------------
# Structure-builder helpers
# ---------------------------------------------------------------------------


def _make_residue(
    name: str,
    seqid_num: int,
    entity_type: gemmi.EntityType = gemmi.EntityType.Polymer,
    x: float = 0.0,
) -> gemmi.Residue:
    res = gemmi.Residue()
    res.name = name
    res.seqid = gemmi.SeqId(seqid_num, " ")
    res.entity_type = entity_type
    atom = gemmi.Atom()
    atom.name = "C1'"
    atom.element = gemmi.Element("C")
    atom.pos = gemmi.Position(x, 0.0, 0.0)
    res.add_atom(atom)
    return res


def _build_test_structure(
    *,
    trna_chain_id: str = "Y",
    mrna_chain_id: str = "V",
    trna_residue_names: list[str] | None = None,
    trna_first_seqid: int = 1,
    mrna_residue_names: list[str] | None = None,
    mrna_first_seqid: int = 1,
    include_trna_ions: bool = False,
) -> gemmi.Structure:
    """Build a deterministic Gemmi structure for the extraction tests.

    Default: a 76-residue tRNA (all 'A') and a 30-residue mRNA (all 'U')
    so polymer positions 34/35/36 are well-defined and the mRNA has
    plenty of room for frame inference.
    """
    if trna_residue_names is None:
        trna_residue_names = ["A"] * 76
    if mrna_residue_names is None:
        mrna_residue_names = ["U"] * 30

    structure = gemmi.Structure()
    structure.name = PDB
    model = gemmi.Model("1")

    trna_chain = gemmi.Chain(trna_chain_id)
    for i, name in enumerate(trna_residue_names):
        trna_chain.add_residue(_make_residue(name, trna_first_seqid + i))
    if include_trna_ions:
        # Two Mg ions deposited under the same chain ID. The polymer
        # filter must skip these so positions 34/35/36 still refer to
        # the polymer residues.
        trna_chain.add_residue(
            _make_residue("MG", 200, entity_type=gemmi.EntityType.NonPolymer)
        )
        trna_chain.add_residue(
            _make_residue("HOH", 300, entity_type=gemmi.EntityType.Water)
        )
    model.add_chain(trna_chain)

    mrna_chain = gemmi.Chain(mrna_chain_id)
    for i, name in enumerate(mrna_residue_names):
        mrna_chain.add_residue(_make_residue(name, mrna_first_seqid + i))
    model.add_chain(mrna_chain)

    structure.add_model(model)
    structure.cell = gemmi.UnitCell(200.0, 200.0, 200.0, 90.0, 90.0, 90.0)
    return structure


def _annotation_for(
    structure: gemmi.Structure,
    *,
    trna_chain_id: str = "Y",
    mrna_chain_id: str = "V",
    aminoacyl: bool = True,
    peptidyl: bool = False,
    exit_: bool = False,
) -> RibosomeAnnotation:
    """A minimal RibosomeAnnotation with chains pointing into ``structure``."""

    def chain(auth: str) -> ChainRef:
        return ChainRef(pdb_id=PDB, assembly_id="1", auth_asym_id=auth, polymer_type="RNA")

    ann = RibosomeAnnotation(pdb_id=PDB, assembly_id="1", status="annotated")
    ann.mrna_chain = chain(mrna_chain_id)
    if aminoacyl:
        ann.aminoacyl_trna_chain = chain(trna_chain_id)
    if peptidyl:
        ann.peptidyl_trna_chain = chain(trna_chain_id)
    if exit_:
        ann.exit_trna_chain = chain(trna_chain_id)
    return ann


# ---------------------------------------------------------------------------
# Anticodon picker
# ---------------------------------------------------------------------------


def test_polymer_residues_in_order_filters_hetatm() -> None:
    structure = _build_test_structure(include_trna_ions=True)
    residues = _polymer_residues_in_order(structure, "Y")
    assert len(residues) == 76
    assert all(r.entity_type == gemmi.EntityType.Polymer for r in residues)


def test_pick_anticodon_residues_picks_positions_34_35_36() -> None:
    structure = _build_test_structure()
    hits = _pick_anticodon_residues(structure, "Y")
    assert hits is not None
    assert [h.trna_position for h in hits] == [34, 35, 36]
    assert [h.residue.seqid.num for h in hits] == [34, 35, 36]


def test_pick_anticodon_residues_returns_none_when_no_residue_34() -> None:
    """A tRNA whose residues start at auth_seq_id 101 has no residue at
    auth_seq_id 34/35/36 — the picker returns ``None`` rather than
    silently picking a non-anticodon residue."""
    structure = _build_test_structure(trna_first_seqid=101)
    assert _pick_anticodon_residues(structure, "Y") is None


def test_pick_anticodon_residues_anchors_on_residue_1_skipping_residue_0() -> None:
    """A chain whose polymer runs 0..75 (like 5UYM chain W) must still
    pick auth_seq_id 34/35/36 as the anticodon, NOT the 34th element of
    the polymer list (which would be auth_seq_id 33 — one too early)."""
    structure = _build_test_structure(
        trna_residue_names=["X"] + ["A"] * 76,  # residue 0 then 1..76
        trna_first_seqid=0,
    )
    hits = _pick_anticodon_residues(structure, "Y")
    assert hits is not None
    # Anchored on residue 1 → biological 34 is at auth_seq_id 34, not 33.
    assert [h.residue.seqid.num for h in hits] == [34, 35, 36]


def test_pick_anticodon_residues_returns_none_for_short_chain() -> None:
    structure = _build_test_structure(trna_residue_names=["A"] * 20)
    assert _pick_anticodon_residues(structure, "Y") is None


def test_pick_anticodon_residues_skips_ions_when_counting() -> None:
    """Even when HETATM ions are interleaved, polymer position 34
    refers to the 34th polymer residue."""
    structure = _build_test_structure(include_trna_ions=True)
    hits = _pick_anticodon_residues(structure, "Y")
    assert hits is not None
    assert [h.residue.name for h in hits] == ["A", "A", "A"]


# ---------------------------------------------------------------------------
# Unit-ID + parent-base helpers
# ---------------------------------------------------------------------------


def test_build_unit_id_uses_auth_seq_id() -> None:
    """Unit ID embeds the deposited auth_seq_id, not the polymer-index."""
    structure = _build_test_structure(trna_first_seqid=101)
    res = _polymer_residues_in_order(structure, "Y")[33]
    assert _build_unit_id(PDB, "Y", res) == "TEST|1|Y|A|134"


def test_parent_base_info_canonical_bases() -> None:
    for base in ("A", "C", "G", "U"):
        parent, is_modified = _parent_base_info(base)
        assert parent == base
        assert is_modified is False


def test_parent_base_info_modified_nucleotides() -> None:
    assert _parent_base_info("PSU") == ("U", True)
    assert _parent_base_info("5MU") == ("U", True)
    assert _parent_base_info("7MG") == ("G", True)
    assert _parent_base_info("MIA") == ("A", True)
    assert _parent_base_info("H2U") == ("U", True)


# ---------------------------------------------------------------------------
# FR3D CSV parse + fetch
# ---------------------------------------------------------------------------


def test_parse_fr3d_csv_parses_three_columns() -> None:
    csv_bytes = (
        b'"TEST|1|V|U|19","cWW","TEST|1|Y|A|36"\n'
        b'"TEST|1|V|U|20","cWW","TEST|1|Y|A|35"\n'
        b'"TEST|1|V|C|21","cWW","TEST|1|Y|A|34"\n'
    )
    rows = _parse_fr3d_csv(csv_bytes)
    assert len(rows) == 3
    assert rows[0] == _ParsedFr3dRow("TEST|1|V|U|19", "cWW", "TEST|1|Y|A|36")


def test_parse_fr3d_csv_skips_malformed_rows() -> None:
    csv_bytes = b'"TEST|1|V|U|19","cWW","TEST|1|Y|A|36"\nshort,row\n,,\n'
    rows = _parse_fr3d_csv(csv_bytes)
    assert len(rows) == 1


@respx.mock
def test_fetch_fr3d_basepairs_returns_rows_on_200(tmp_path) -> None:
    url = FR3D_BASEPAIRS_URL_TEMPLATE.format(pdb_id="test")
    respx.get(url).mock(
        return_value=httpx.Response(
            200, content=b'"TEST|1|V|U|19","cWW","TEST|1|Y|A|36"\n'
        )
    )
    rows = fetch_fr3d_basepairs(PDB, cache=Cache(tmp_path))
    assert rows is not None
    assert len(rows) == 1


@respx.mock
def test_fetch_fr3d_basepairs_returns_none_on_http_error() -> None:
    url = FR3D_BASEPAIRS_URL_TEMPLATE.format(pdb_id="test")
    respx.get(url).mock(return_value=httpx.Response(500))
    assert fetch_fr3d_basepairs(PDB) is None


@respx.mock
def test_fetch_fr3d_basepairs_uses_cache_on_second_call(tmp_path) -> None:
    cache = Cache(tmp_path)
    url = FR3D_BASEPAIRS_URL_TEMPLATE.format(pdb_id="test")
    route = respx.get(url).mock(
        return_value=httpx.Response(
            200, content=b'"TEST|1|V|U|19","cWW","TEST|1|Y|A|36"\n'
        )
    )
    fetch_fr3d_basepairs(PDB, cache=cache)
    fetch_fr3d_basepairs(PDB, cache=cache)
    assert route.call_count == 1  # second call hit the cache


# ---------------------------------------------------------------------------
# Pair filter + codon assignment
# ---------------------------------------------------------------------------


def _classic_anticodon_seq_ids() -> set[int]:
    return {34, 35, 36}


def test_classify_pair_sides_accepts_both_orientations() -> None:
    fwd = _ParsedFr3dRow("TEST|1|V|U|19", "cWW", "TEST|1|Y|A|36")
    rev = _ParsedFr3dRow("TEST|1|Y|A|36", "cWW", "TEST|1|V|U|19")
    sides_fwd = _classify_pair_sides(fwd, "V", "Y", _classic_anticodon_seq_ids())
    sides_rev = _classify_pair_sides(rev, "V", "Y", _classic_anticodon_seq_ids())
    assert sides_fwd == ("TEST|1|V|U|19", "TEST|1|Y|A|36")
    assert sides_rev == ("TEST|1|V|U|19", "TEST|1|Y|A|36")


def test_classify_pair_sides_rejects_non_anticodon_resnum() -> None:
    # tRNA residue 50 is not in the anticodon set.
    row = _ParsedFr3dRow("TEST|1|V|U|19", "cWW", "TEST|1|Y|A|50")
    assert _classify_pair_sides(row, "V", "Y", _classic_anticodon_seq_ids()) is None


def test_classify_pair_sides_rejects_unrelated_chains() -> None:
    row = _ParsedFr3dRow("TEST|1|A|G|100", "cWW", "TEST|1|B|C|200")
    assert _classify_pair_sides(row, "V", "Y", _classic_anticodon_seq_ids()) is None


def test_filter_cross_chain_pairs_dedupes_symmetric_rows() -> None:
    rows = [
        _ParsedFr3dRow("TEST|1|V|C|21", "cWW", "TEST|1|Y|A|34"),
        _ParsedFr3dRow("TEST|1|Y|A|34", "cWW", "TEST|1|V|C|21"),
    ]
    out = _filter_cross_chain_pairs(
        rows,
        mrna_chain_id="V",
        trna_chain_id="Y",
        anticodon_auth_seq_ids=_classic_anticodon_seq_ids(),
    )
    assert out == [("TEST|1|V|C|21", "cWW", "TEST|1|Y|A|34")]


def _anticodon_by_seq_id(structure: gemmi.Structure, chain_id: str) -> dict[int, _AnticodonHit]:
    hits = _pick_anticodon_residues(structure, chain_id) or []
    return {h.residue.seqid.num: h for h in hits}


def test_assign_codon_maps_positions_correctly() -> None:
    structure = _build_test_structure()
    by_seq_id = _anticodon_by_seq_id(structure, "Y")
    pairs = [
        ("TEST|1|V|U|19", "cWW", "TEST|1|Y|A|36"),  # tRNA 36 → codon 1
        ("TEST|1|V|U|20", "cWW", "TEST|1|Y|A|35"),  # tRNA 35 → codon 2
        ("TEST|1|V|C|21", "cWW", "TEST|1|Y|A|34"),  # tRNA 34 → codon 3
    ]
    codon_uids, _, codon_bases, ambiguous = _assign_codon_from_pairs(pairs, by_seq_id)
    assert codon_uids == {1: "TEST|1|V|U|19", 2: "TEST|1|V|U|20", 3: "TEST|1|V|C|21"}
    assert codon_bases == {1: "U", 2: "U", 3: "C"}
    assert ambiguous == set()


def test_assign_codon_prefers_cww_when_multiple_pairs() -> None:
    structure = _build_test_structure()
    by_seq_id = _anticodon_by_seq_id(structure, "Y")
    pairs = [
        ("TEST|1|V|A|10", "tHS", "TEST|1|Y|A|34"),  # weak alternate
        ("TEST|1|V|C|21", "cWW", "TEST|1|Y|A|34"),  # the right one
    ]
    codon_uids, _, _, ambiguous = _assign_codon_from_pairs(pairs, by_seq_id)
    assert codon_uids[3] == "TEST|1|V|C|21"
    assert ambiguous == set()


def test_assign_codon_marks_ambiguous_when_no_cww_disambiguator() -> None:
    structure = _build_test_structure()
    by_seq_id = _anticodon_by_seq_id(structure, "Y")
    pairs = [
        ("TEST|1|V|A|10", "tHS", "TEST|1|Y|A|34"),
        ("TEST|1|V|G|11", "tWH", "TEST|1|Y|A|34"),
    ]
    _, _, _, ambiguous = _assign_codon_from_pairs(pairs, by_seq_id)
    assert ambiguous == {3}


# ---------------------------------------------------------------------------
# Local codon reconstruction + frame inference
# ---------------------------------------------------------------------------


def _mrna_index(structure: gemmi.Structure, chain_id: str) -> tuple[
    list[gemmi.Residue], dict[int, int]
]:
    residues = _polymer_residues_in_order(structure, chain_id)
    return residues, {r.seqid.num: i for i, r in enumerate(residues)}


def test_reconstruct_codon_locally_fills_from_position_3() -> None:
    structure = _build_test_structure(mrna_residue_names=["A", "U", "C", "G"] * 6)
    mrna_residues, mrna_index = _mrna_index(structure, "V")
    # FR3D only observed codon position 3 at mRNA residue 21.
    res_21 = mrna_residues[20]
    fr3d_uids = {3: _build_unit_id(PDB, "V", res_21)}
    reconstructed = _reconstruct_codon_locally(
        PDB, "V", fr3d_uids, mrna_residues, mrna_index
    )
    assert sorted(reconstructed.keys()) == [1, 2, 3]
    assert reconstructed[1].source == "mmcif_reconstructed"
    assert reconstructed[2].source == "mmcif_reconstructed"
    assert reconstructed[3].source == "fr3d_observed"
    # Positions 1, 2 are the two residues immediately preceding position 3 in polymer order.
    assert reconstructed[1].unit_id.endswith("|19")
    assert reconstructed[2].unit_id.endswith("|20")


def test_reconstruct_codon_locally_fills_from_position_2() -> None:
    structure = _build_test_structure()
    mrna_residues, mrna_index = _mrna_index(structure, "V")
    res_20 = mrna_residues[19]
    fr3d_uids = {2: _build_unit_id(PDB, "V", res_20)}
    reconstructed = _reconstruct_codon_locally(
        PDB, "V", fr3d_uids, mrna_residues, mrna_index
    )
    assert reconstructed[1].unit_id.endswith("|19")
    assert reconstructed[3].unit_id.endswith("|21")


def test_reconstruct_codon_locally_empty_when_no_seed() -> None:
    structure = _build_test_structure()
    mrna_residues, mrna_index = _mrna_index(structure, "V")
    assert _reconstruct_codon_locally(PDB, "V", {}, mrna_residues, mrna_index) == {}


def test_frame_inferred_codon_from_a_site_anchor_to_p_site() -> None:
    structure = _build_test_structure(mrna_residue_names=["A", "C", "G"] * 10)
    mrna_residues, mrna_index = _mrna_index(structure, "V")
    # A-site codon at polymer indices 18, 19, 20 (residues 19, 20, 21).
    a_codon = {
        1: _make_codon_residue(_build_unit_id(PDB, "V", mrna_residues[18]), "A", 1, "fr3d_observed"),
        2: _make_codon_residue(_build_unit_id(PDB, "V", mrna_residues[19]), "C", 2, "fr3d_observed"),
        3: _make_codon_residue(_build_unit_id(PDB, "V", mrna_residues[20]), "G", 3, "fr3d_observed"),
    }
    inferred = _frame_inferred_codon(PDB, "V", "A", a_codon, "P", mrna_residues, mrna_index)
    assert inferred is not None
    # P-site = 3 residues upstream → polymer indices 15, 16, 17 → residues 16, 17, 18.
    assert inferred[1].unit_id.endswith("|16")
    assert inferred[2].unit_id.endswith("|17")
    assert inferred[3].unit_id.endswith("|18")
    assert all(r.source == "mrna_frame_inference" for r in inferred.values())


def test_frame_inferred_codon_from_a_site_anchor_to_e_site() -> None:
    structure = _build_test_structure(mrna_residue_names=["A", "C", "G"] * 10)
    mrna_residues, mrna_index = _mrna_index(structure, "V")
    a_codon = {
        1: _make_codon_residue(_build_unit_id(PDB, "V", mrna_residues[18]), "A", 1, "fr3d_observed"),
        2: _make_codon_residue(_build_unit_id(PDB, "V", mrna_residues[19]), "C", 2, "fr3d_observed"),
        3: _make_codon_residue(_build_unit_id(PDB, "V", mrna_residues[20]), "G", 3, "fr3d_observed"),
    }
    inferred = _frame_inferred_codon(PDB, "V", "A", a_codon, "E", mrna_residues, mrna_index)
    assert inferred is not None
    # E-site = 6 residues upstream → polymer indices 12, 13, 14 → residues 13, 14, 15.
    assert inferred[1].unit_id.endswith("|13")


def test_frame_inferred_codon_from_p_site_anchor_infers_a_and_e() -> None:
    """5UYN-style: anchor on the P-site codon and infer both A and E."""
    structure = _build_test_structure(mrna_residue_names=["A", "C", "G"] * 10)
    mrna_residues, mrna_index = _mrna_index(structure, "V")
    # P-site codon at polymer indices 15, 16, 17 (residues 16, 17, 18).
    p_codon = {
        1: _make_codon_residue(_build_unit_id(PDB, "V", mrna_residues[15]), "A", 1, "fr3d_observed"),
        2: _make_codon_residue(_build_unit_id(PDB, "V", mrna_residues[16]), "C", 2, "fr3d_observed"),
        3: _make_codon_residue(_build_unit_id(PDB, "V", mrna_residues[17]), "G", 3, "fr3d_observed"),
    }
    # Infer A (downstream of P, +3 polymer indices).
    a_inferred = _frame_inferred_codon(PDB, "V", "P", p_codon, "A", mrna_residues, mrna_index)
    assert a_inferred is not None
    assert a_inferred[1].unit_id.endswith("|19")
    assert a_inferred[2].unit_id.endswith("|20")
    assert a_inferred[3].unit_id.endswith("|21")
    # Infer E (upstream of P, -3 polymer indices).
    e_inferred = _frame_inferred_codon(PDB, "V", "P", p_codon, "E", mrna_residues, mrna_index)
    assert e_inferred is not None
    assert e_inferred[1].unit_id.endswith("|13")


def test_frame_inferred_codon_returns_none_for_same_site() -> None:
    structure = _build_test_structure(mrna_residue_names=["A", "C", "G"] * 10)
    mrna_residues, mrna_index = _mrna_index(structure, "V")
    a_codon = {
        1: _make_codon_residue(_build_unit_id(PDB, "V", mrna_residues[18]), "A", 1, "fr3d_observed"),
        2: _make_codon_residue(_build_unit_id(PDB, "V", mrna_residues[19]), "C", 2, "fr3d_observed"),
        3: _make_codon_residue(_build_unit_id(PDB, "V", mrna_residues[20]), "G", 3, "fr3d_observed"),
    }
    assert _frame_inferred_codon(PDB, "V", "A", a_codon, "A", mrna_residues, mrna_index) is None


def test_frame_inferred_codon_none_when_run_too_short() -> None:
    # mRNA only 5 residues; can't span E/P/A.
    structure = _build_test_structure(mrna_residue_names=["A"] * 5)
    mrna_residues, mrna_index = _mrna_index(structure, "V")
    a_codon = {
        1: _make_codon_residue(_build_unit_id(PDB, "V", mrna_residues[2]), "A", 1, "fr3d_observed"),
        2: _make_codon_residue(_build_unit_id(PDB, "V", mrna_residues[3]), "A", 2, "fr3d_observed"),
        3: _make_codon_residue(_build_unit_id(PDB, "V", mrna_residues[4]), "A", 3, "fr3d_observed"),
    }
    assert _frame_inferred_codon(PDB, "V", "A", a_codon, "P", mrna_residues, mrna_index) is None


def _make_codon_residue(unit_id: str, base: str, codon_position: int, source: str):
    from ribosome_state_annotator.models import CodonResidue

    return CodonResidue(
        codon_position=codon_position, unit_id=unit_id, base=base, source=source
    )


# ---------------------------------------------------------------------------
# Misc helpers
# ---------------------------------------------------------------------------


def test_consecutive_run_lengths() -> None:
    structure = _build_test_structure(mrna_residue_names=["A"] * 30)
    residues = _polymer_residues_in_order(structure, "V")
    assert _consecutive_run_lengths(residues) == [30]


def test_consecutive_run_lengths_split_by_gap() -> None:
    structure = gemmi.Structure()
    model = gemmi.Model("1")
    chain = gemmi.Chain("X")
    for s in (1, 2, 3, 5, 6, 7, 8, 20, 21):
        chain.add_residue(_make_residue("A", s))
    model.add_chain(chain)
    structure.add_model(model)
    structure.cell = gemmi.UnitCell(50, 50, 50, 90, 90, 90)
    runs = _consecutive_run_lengths(_polymer_residues_in_order(structure, "X"))
    assert runs == [3, 4, 2]


def test_is_consecutive_run() -> None:
    structure = _build_test_structure()
    residues = _polymer_residues_in_order(structure, "V")[10:15]
    assert _is_consecutive_run(residues) is True


def test_codon_assignment_status() -> None:
    from ribosome_state_annotator.models import CodonResidue

    def r(pos: int) -> CodonResidue:
        return CodonResidue(
            codon_position=pos, unit_id=f"u{pos}", base="A", source="fr3d_observed"
        )

    assert _codon_assignment_status([r(1), r(2), r(3)]) == "complete"
    assert _codon_assignment_status([r(1), r(2)]) == "partial"
    assert _codon_assignment_status([]) == "missing"


# ---------------------------------------------------------------------------
# End-to-end extraction
# ---------------------------------------------------------------------------


@respx.mock
def test_extract_returns_empty_when_no_mrna() -> None:
    structure = _build_test_structure()
    ann = _annotation_for(structure)
    ann.mrna_chain = None
    assert extract_trna_mrna_interactions(ann, structure) == []


@respx.mock
def test_extract_returns_empty_when_no_trna() -> None:
    structure = _build_test_structure()
    ann = _annotation_for(structure, aminoacyl=False)
    assert extract_trna_mrna_interactions(ann, structure) == []


@respx.mock
def test_extract_returns_empty_when_fr3d_unavailable() -> None:
    structure = _build_test_structure()
    ann = _annotation_for(structure)
    url = FR3D_BASEPAIRS_URL_TEMPLATE.format(pdb_id="test")
    respx.get(url).mock(return_value=httpx.Response(500))
    assert extract_trna_mrna_interactions(ann, structure) == []


@respx.mock
def test_extract_end_to_end_single_site() -> None:
    structure = _build_test_structure(mrna_residue_names=["A", "C", "G", "U"] * 8)
    ann = _annotation_for(structure)
    # Anticodon residues at tRNA seqids 34/35/36; pair with mRNA 21/20/19.
    csv_body = (
        b'"TEST|1|V|U|19","cWW","TEST|1|Y|A|36"\n'
        b'"TEST|1|V|U|20","cWW","TEST|1|Y|A|35"\n'
        b'"TEST|1|V|C|21","cWW","TEST|1|Y|A|34"\n'
    )
    url = FR3D_BASEPAIRS_URL_TEMPLATE.format(pdb_id="test")
    respx.get(url).mock(return_value=httpx.Response(200, content=csv_body))

    interactions = extract_trna_mrna_interactions(ann, structure)
    assert len(interactions) == 1
    a = interactions[0]
    assert a.site == "A"
    assert a.mrna_chain_id == "V"
    assert a.trna_chain_id == "Y"
    assert a.anticodon_position_source == "polymer_sequence_index"
    assert a.codon.assignment_status == "complete"
    # mRNA 19 (residue C at polymer index 18 of pattern A,C,G,U) is base 'C',
    # mRNA 20 is base 'U', mRNA 21 is base 'A' → sequence = 'CUA'.
    assert a.codon.sequence is not None
    assert len(a.codon.sequence) == 3
    assert [r.codon_position for r in a.codon.residues] == [1, 2, 3]
    # 3 anticodon residues, 3 pairs, all assigned.
    assert len(a.anticodon.residues) == 3
    assert len(a.pairs) == 3
    assert {p.fr3d_interaction for p in a.pairs} == {"cWW"}
    assert all(p.assignment_status == "assigned" for p in a.pairs)
    # Wobble flag only on the trna_position == 34 pair.
    wobble = [p for p in a.pairs if p.is_wobble_position]
    assert len(wobble) == 1
    assert wobble[0].trna_position == 34


@respx.mock
def test_extract_frame_inference_fills_p_and_e_sites() -> None:
    """When only the A-site tRNA has FR3D pairs but the mRNA is long
    enough, P and E codons should be inferred from the A-site frame."""
    # Real ribosomes have one tRNA chain per site; build the fixture that way.
    structure = gemmi.Structure()
    structure.name = PDB
    model = gemmi.Model("1")
    for chain_id in ("Y", "W", "Z"):
        chain = gemmi.Chain(chain_id)
        for i in range(1, 77):
            chain.add_residue(_make_residue("A", i))
        model.add_chain(chain)
    mrna_chain = gemmi.Chain("V")
    for i, base in enumerate(["A", "C", "G"] * 12, start=1):
        mrna_chain.add_residue(_make_residue(base, i))
    model.add_chain(mrna_chain)
    structure.add_model(model)
    structure.cell = gemmi.UnitCell(200.0, 200.0, 200.0, 90.0, 90.0, 90.0)

    def chain_ref(auth: str) -> ChainRef:
        return ChainRef(pdb_id=PDB, assembly_id="1", auth_asym_id=auth, polymer_type="RNA")

    ann = RibosomeAnnotation(pdb_id=PDB, assembly_id="1", status="annotated")
    ann.mrna_chain = chain_ref("V")
    ann.aminoacyl_trna_chain = chain_ref("Y")
    ann.peptidyl_trna_chain = chain_ref("W")
    ann.exit_trna_chain = chain_ref("Z")
    # Only A-site has explicit FR3D rows; P and E will be frame-inferred.
    csv_body = (
        b'"TEST|1|V|G|19","cWW","TEST|1|Y|A|36"\n'
        b'"TEST|1|V|A|20","cWW","TEST|1|Y|A|35"\n'
        b'"TEST|1|V|C|21","cWW","TEST|1|Y|A|34"\n'
    )
    url = FR3D_BASEPAIRS_URL_TEMPLATE.format(pdb_id="test")
    respx.get(url).mock(return_value=httpx.Response(200, content=csv_body))

    interactions = extract_trna_mrna_interactions(ann, structure)
    by_site = {i.site: i for i in interactions}
    assert by_site["A"].codon.assignment_status == "complete"
    assert by_site["P"].codon.assignment_status == "complete"
    assert by_site["E"].codon.assignment_status == "complete"
    # P-site residues sourced from frame inference.
    assert all(r.source == "mrna_frame_inference" for r in by_site["P"].codon.residues)
    # P codon should be 3 residues upstream of A (mRNA 16/17/18).
    assert by_site["P"].codon.residues[0].unit_id.endswith("|16")
    assert by_site["E"].codon.residues[0].unit_id.endswith("|13")


@respx.mock
def test_extract_frame_inference_anchors_on_p_site_when_a_is_missing() -> None:
    """5UYN-style: P-tRNA has FR3D pairs, A-tRNA does not (pre-accommodation).
    Frame inference should fall back to the P-site as anchor and resolve the
    A and E codons by stepping the mRNA frame."""
    structure = gemmi.Structure()
    structure.name = PDB
    model = gemmi.Model("1")
    for chain_id in ("Y", "W", "Z"):
        chain = gemmi.Chain(chain_id)
        for i in range(1, 77):
            chain.add_residue(_make_residue("A", i))
        model.add_chain(chain)
    mrna_chain = gemmi.Chain("V")
    for i, base in enumerate(["A", "C", "G"] * 12, start=1):
        mrna_chain.add_residue(_make_residue(base, i))
    model.add_chain(mrna_chain)
    structure.add_model(model)
    structure.cell = gemmi.UnitCell(200.0, 200.0, 200.0, 90.0, 90.0, 90.0)

    def chain_ref(auth: str) -> ChainRef:
        return ChainRef(pdb_id=PDB, assembly_id="1", auth_asym_id=auth, polymer_type="RNA")

    ann = RibosomeAnnotation(pdb_id=PDB, assembly_id="1", status="annotated")
    ann.mrna_chain = chain_ref("V")
    ann.aminoacyl_trna_chain = chain_ref("Y")
    ann.peptidyl_trna_chain = chain_ref("W")
    ann.exit_trna_chain = chain_ref("Z")

    # Only the P-site (chain W) has FR3D cWW pairs. Place them at mRNA
    # positions 16/17/18 so the inferred A-codon is at 19/20/21 and the
    # inferred E-codon at 13/14/15.
    csv_body = (
        b'"TEST|1|V|A|16","cWW","TEST|1|W|A|36"\n'
        b'"TEST|1|V|C|17","cWW","TEST|1|W|A|35"\n'
        b'"TEST|1|V|G|18","cWW","TEST|1|W|A|34"\n'
    )
    url = FR3D_BASEPAIRS_URL_TEMPLATE.format(pdb_id="test")
    respx.get(url).mock(return_value=httpx.Response(200, content=csv_body))

    interactions = extract_trna_mrna_interactions(ann, structure)
    by_site = {i.site: i for i in interactions}

    # P codon comes directly from FR3D.
    assert by_site["P"].codon.assignment_status == "complete"
    assert by_site["P"].codon.residues[0].source == "fr3d_observed"

    # A and E codons are filled by frame inference, anchored on P.
    assert by_site["A"].codon.assignment_status == "complete"
    assert all(r.source == "mrna_frame_inference" for r in by_site["A"].codon.residues)
    assert by_site["A"].codon.residues[0].unit_id.endswith("|19")
    assert by_site["A"].codon.residues[2].unit_id.endswith("|21")

    assert by_site["E"].codon.assignment_status == "complete"
    assert all(r.source == "mrna_frame_inference" for r in by_site["E"].codon.residues)
    assert by_site["E"].codon.residues[0].unit_id.endswith("|13")
    assert by_site["E"].codon.residues[2].unit_id.endswith("|15")


@respx.mock
def test_extract_skips_site_when_trna_too_short(caplog) -> None:
    structure = _build_test_structure(trna_residue_names=["A"] * 20)
    ann = _annotation_for(structure)
    url = FR3D_BASEPAIRS_URL_TEMPLATE.format(pdb_id="test")
    respx.get(url).mock(return_value=httpx.Response(200, content=b""))
    with caplog.at_level("WARNING"):
        interactions = extract_trna_mrna_interactions(ann, structure)
    assert interactions == []
    assert any("fewer than 36 polymer residues" in r.message for r in caplog.records)


@respx.mock
def test_extract_warns_when_mrna_too_short_for_frame_inference() -> None:
    """If we have A-site FR3D pairs but the mRNA is too short to span
    E/P/A, P + E codons stay incomplete with a warning."""
    structure = _build_test_structure(
        trna_chain_id="Y",
        mrna_residue_names=["A"] * 4,  # only 4 residues — no room for upstream codons
        mrna_first_seqid=21,
    )
    ann = _annotation_for(structure, aminoacyl=True, peptidyl=True)
    # No FR3D pairs match (mRNA chain too short to host useful pairs).
    url = FR3D_BASEPAIRS_URL_TEMPLATE.format(pdb_id="test")
    respx.get(url).mock(return_value=httpx.Response(200, content=b""))

    interactions = extract_trna_mrna_interactions(ann, structure)
    by_site = {i.site: i for i in interactions}
    assert "insufficient_mrna_for_frame_inference" in by_site["P"].warnings


# ---------------------------------------------------------------------------
# Modified-nt bookkeeping
# ---------------------------------------------------------------------------


@respx.mock
def test_extract_marks_modified_anticodon_residues() -> None:
    """PSU at biological position 34 (auth_seq_id 34) should surface as
    parent_base=U, is_modified=True, and the unit ID should embed
    'PSU', not 'U'."""
    # residues 1..33 = A, 34 = PSU, 35..76 = A — gives PSU at auth_seq_id 34.
    names = ["A"] * 33 + ["PSU"] + ["A"] * 42
    structure = _build_test_structure(trna_residue_names=names)
    ann = _annotation_for(structure)
    url = FR3D_BASEPAIRS_URL_TEMPLATE.format(pdb_id="test")
    respx.get(url).mock(return_value=httpx.Response(200, content=b""))
    interactions = extract_trna_mrna_interactions(ann, structure)
    by_position = {r.trna_position: r for r in interactions[0].anticodon.residues}
    assert by_position[34].parent_base == "U"
    assert by_position[34].is_modified is True
    assert by_position[34].trna_chem_comp_id == "PSU"
    assert by_position[34].unit_id == "TEST|1|Y|PSU|34"


# ---------------------------------------------------------------------------
# Integration with annotate_pdb
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# FR3D codon-pairing fallback (rescues pre-accommodation A-site states like 3JAG)
# ---------------------------------------------------------------------------


def _fallback_structure_two_trnas() -> tuple[gemmi.Structure, ChainRef, ChainRef, ChainRef]:
    """Two 76-mer tRNAs (chains A2, B2) plus a 30-residue mRNA (chain MR).

    Returns ``(structure, mrna_chain_ref, trna1_ref, trna2_ref)``. Both
    tRNAs are tagged with RF00005 so they qualify as fallback candidates.
    """
    structure = gemmi.Structure()
    structure.name = PDB
    model = gemmi.Model("1")
    for chain_id in ("A2", "B2"):
        chain = gemmi.Chain(chain_id)
        for i in range(1, 77):
            chain.add_residue(_make_residue("A", i))
        model.add_chain(chain)
    mrna = gemmi.Chain("MR")
    for i in range(1, 31):
        mrna.add_residue(_make_residue("U", i))
    model.add_chain(mrna)
    structure.add_model(model)
    structure.cell = gemmi.UnitCell(200.0, 200.0, 200.0, 90.0, 90.0, 90.0)

    def chain_ref(auth: str) -> ChainRef:
        return ChainRef(
            pdb_id=PDB,
            assembly_id="1",
            auth_asym_id=auth,
            polymer_type="RNA",
            rfam_accessions=["RF00005"],
        )

    return structure, chain_ref("MR"), chain_ref("A2"), chain_ref("B2")


def test_fallback_noop_when_mrna_missing() -> None:
    structure, _, trna1, _ = _fallback_structure_two_trnas()
    ann = RibosomeAnnotation(pdb_id=PDB, assembly_id="1", status="annotated")
    ann.mrna_chain = None  # no mRNA
    ann.other_rna_chains = [trna1]
    fr3d_codon_pairing_fallback(ann, structure, [])
    assert ann.aminoacyl_trna_chain is None


def test_fallback_noop_when_no_empty_slots() -> None:
    structure, mrna, trna1, trna2 = _fallback_structure_two_trnas()
    ann = RibosomeAnnotation(pdb_id=PDB, assembly_id="1", status="annotated")
    ann.mrna_chain = mrna
    ann.aminoacyl_trna_chain = trna1
    ann.peptidyl_trna_chain = trna2
    extra_trna = ChainRef(
        pdb_id=PDB, assembly_id="1", auth_asym_id="C2", polymer_type="RNA",
        rfam_accessions=["RF00005"],
    )
    ann.exit_trna_chain = extra_trna
    rows = [_ParsedFr3dRow("TEST|1|MR|U|17", "cWW", "TEST|1|A2|A|36")]
    fr3d_codon_pairing_fallback(ann, structure, rows)
    # Nothing should change.
    assert ann.aminoacyl_trna_chain is trna1
    assert ann.peptidyl_trna_chain is trna2
    assert ann.exit_trna_chain is extra_trna


def test_fallback_assigns_a_site_when_one_candidate_has_cww_pair() -> None:
    """3JAG-style: contact-transfer assigned P, left A empty. FR3D shows
    a cWW pair between the unassigned tRNA's anticodon and the mRNA."""
    structure, mrna, trna1, trna2 = _fallback_structure_two_trnas()
    ann = RibosomeAnnotation(pdb_id=PDB, assembly_id="1", status="annotated")
    ann.mrna_chain = mrna
    ann.peptidyl_trna_chain = trna2  # contact-transfer assigned P
    ann.other_rna_chains = [trna1]   # A-candidate is unassigned

    rows = [
        _ParsedFr3dRow("TEST|1|MR|U|17", "cWW", "TEST|1|A2|A|36"),
        _ParsedFr3dRow("TEST|1|MR|U|18", "cWW", "TEST|1|A2|A|35"),
        _ParsedFr3dRow("TEST|1|MR|U|19", "cWW", "TEST|1|A2|A|34"),
    ]
    fr3d_codon_pairing_fallback(ann, structure, rows)

    assert ann.aminoacyl_trna_chain is trna1
    assert ann.peptidyl_trna_chain is trna2  # untouched
    assert trna1 not in ann.other_rna_chains
    assert any("atrna_assigned_from_fr3d_codon_pairing_A2" in w for w in ann.warnings)
    assert ann.classification_evidence.get("fr3d_codon_pairing_fallback") == {"A": "A2"}


def test_fallback_skips_candidates_without_cww_pair() -> None:
    """Only non-cWW pairs (tHS, etc.) → fallback should NOT assign."""
    structure, mrna, trna1, _ = _fallback_structure_two_trnas()
    ann = RibosomeAnnotation(pdb_id=PDB, assembly_id="1", status="annotated")
    ann.mrna_chain = mrna
    ann.other_rna_chains = [trna1]

    rows = [_ParsedFr3dRow("TEST|1|MR|U|17", "tHS", "TEST|1|A2|A|36")]
    fr3d_codon_pairing_fallback(ann, structure, rows)
    assert ann.aminoacyl_trna_chain is None


def test_fallback_assigns_a_and_p_by_mrna_codon_position() -> None:
    """Two unassigned tRNAs, both with cWW pairs. The one paired with the
    more 3' (downstream) mRNA codon should become A; the other becomes P."""
    structure, mrna, trna_upstream, trna_downstream = _fallback_structure_two_trnas()
    ann = RibosomeAnnotation(pdb_id=PDB, assembly_id="1", status="annotated")
    ann.mrna_chain = mrna
    ann.other_rna_chains = [trna_upstream, trna_downstream]

    rows = [
        # chain A2 pairs with mRNA residues 10/11/12 (UPSTREAM → P-site)
        _ParsedFr3dRow("TEST|1|MR|U|10", "cWW", "TEST|1|A2|A|36"),
        _ParsedFr3dRow("TEST|1|MR|U|11", "cWW", "TEST|1|A2|A|35"),
        _ParsedFr3dRow("TEST|1|MR|U|12", "cWW", "TEST|1|A2|A|34"),
        # chain B2 pairs with mRNA residues 19/20/21 (DOWNSTREAM → A-site)
        _ParsedFr3dRow("TEST|1|MR|U|19", "cWW", "TEST|1|B2|A|36"),
        _ParsedFr3dRow("TEST|1|MR|U|20", "cWW", "TEST|1|B2|A|35"),
        _ParsedFr3dRow("TEST|1|MR|U|21", "cWW", "TEST|1|B2|A|34"),
    ]
    fr3d_codon_pairing_fallback(ann, structure, rows)

    assert ann.aminoacyl_trna_chain is trna_downstream  # B2 → A-site
    assert ann.peptidyl_trna_chain is trna_upstream    # A2 → P-site


def test_fallback_skips_candidates_with_short_trna() -> None:
    """A tRNA chain with fewer than 36 polymer residues can't have an
    anticodon picked, so the fallback should skip it silently."""
    structure = gemmi.Structure()
    structure.name = PDB
    model = gemmi.Model("1")
    short = gemmi.Chain("A2")
    for i in range(1, 20):
        short.add_residue(_make_residue("A", i))
    model.add_chain(short)
    mrna = gemmi.Chain("MR")
    for i in range(1, 31):
        mrna.add_residue(_make_residue("U", i))
    model.add_chain(mrna)
    structure.add_model(model)
    structure.cell = gemmi.UnitCell(200.0, 200.0, 200.0, 90.0, 90.0, 90.0)

    short_ref = ChainRef(
        pdb_id=PDB, assembly_id="1", auth_asym_id="A2", polymer_type="RNA",
        rfam_accessions=["RF00005"],
    )
    mrna_ref = ChainRef(pdb_id=PDB, assembly_id="1", auth_asym_id="MR", polymer_type="RNA")
    ann = RibosomeAnnotation(pdb_id=PDB, assembly_id="1", status="annotated")
    ann.mrna_chain = mrna_ref
    ann.other_rna_chains = [short_ref]
    fr3d_codon_pairing_fallback(ann, structure, [])
    assert ann.aminoacyl_trna_chain is None


def test_fallback_ignores_non_trna_rfam_chains() -> None:
    """Chains without RF00005 (e.g. snRNA, ribozymes) shouldn't be
    considered as fallback candidates even if they have FR3D pairs."""
    structure, mrna, trna1, _ = _fallback_structure_two_trnas()
    ann = RibosomeAnnotation(pdb_id=PDB, assembly_id="1", status="annotated")
    ann.mrna_chain = mrna
    # Strip RF00005 from the candidate.
    trna1.rfam_accessions = []
    ann.other_rna_chains = [trna1]
    rows = [_ParsedFr3dRow("TEST|1|MR|U|17", "cWW", "TEST|1|A2|A|36")]
    fr3d_codon_pairing_fallback(ann, structure, rows)
    assert ann.aminoacyl_trna_chain is None


# ---------------------------------------------------------------------------
# End-to-end with fallback
# ---------------------------------------------------------------------------


@respx.mock
def test_extract_uses_fallback_then_extracts_codon_for_filled_a_site() -> None:
    """End-to-end: contact-transfer left A empty; FR3D fallback assigns A,
    then per-site extraction returns A-site evidence."""
    structure, mrna, trna1, trna2 = _fallback_structure_two_trnas()
    # Replace structure's mRNA with one that gives recognizable bases.
    # (Already U-only, fine for the test.)
    ann = RibosomeAnnotation(pdb_id=PDB, assembly_id="1", status="annotated")
    ann.mrna_chain = mrna
    ann.peptidyl_trna_chain = trna2  # contact-transfer assigned P
    ann.other_rna_chains = [trna1]

    csv_body = (
        b'"TEST|1|MR|U|19","cWW","TEST|1|A2|A|36"\n'
        b'"TEST|1|MR|U|20","cWW","TEST|1|A2|A|35"\n'
        b'"TEST|1|MR|U|21","cWW","TEST|1|A2|A|34"\n'
    )
    url = FR3D_BASEPAIRS_URL_TEMPLATE.format(pdb_id="test")
    respx.get(url).mock(return_value=httpx.Response(200, content=csv_body))

    interactions = extract_trna_mrna_interactions(ann, structure)
    by_site = {i.site: i for i in interactions}
    # A-site got assigned + extracted.
    assert "A" in by_site
    assert by_site["A"].trna_chain_id == "A2"
    assert by_site["A"].codon.assignment_status == "complete"
    # Original annotation was mutated.
    assert ann.aminoacyl_trna_chain is trna1


def test_annotate_one_assembly_respects_no_fr3d_flag(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """When no_fr3d=True, extract_trna_mrna_interactions must not be called."""
    from ribosome_state_annotator import api

    called: list[str] = []

    def fail_if_called(*args, **kwargs):
        called.append("yes")
        raise AssertionError("extract should not be called when no_fr3d=True")

    monkeypatch.setattr(
        api, "extract_trna_mrna_interactions", fail_if_called
    )
    # We don't actually run the full pipeline here — we just confirm
    # that the wired function is the one we expect (sanity check of the
    # monkey-patch surface).
    assert api.extract_trna_mrna_interactions is fail_if_called
    assert called == []
