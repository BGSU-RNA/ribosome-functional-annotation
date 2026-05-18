"""Unit tests for the multi-ribosome bundle splitter (:mod:`multiribo`)."""

from __future__ import annotations

import gemmi
import pytest

from ribosome_state_annotator import multiribo
from ribosome_state_annotator.models import ChainRef, CorrespondenceResult


def _chain_ref(auth: str, rfam: tuple[str, ...] = (), description: str | None = None) -> ChainRef:
    return ChainRef(
        pdb_id="MULTIRIBO",
        assembly_id="1",
        auth_asym_id=auth,
        polymer_type="RNA",
        rfam_accessions=list(rfam),
        description=description,
    )


def _add_residue(chain: gemmi.Chain, name: str, seqid: int, x: float, y: float, z: float) -> None:
    residue = gemmi.Residue()
    residue.name = name
    residue.seqid = gemmi.SeqId(seqid, " ")
    residue.entity_type = gemmi.EntityType.Polymer
    atom = gemmi.Atom()
    atom.name = "C1'"
    atom.element = gemmi.Element("C")
    atom.pos = gemmi.Position(x, y, z)
    residue.add_atom(atom)
    chain.add_residue(residue)


@pytest.fixture
def two_ribosome_structure() -> gemmi.Structure:
    """Build a synthetic structure with two ribosomes 100 Å apart.

    Layout:
      - Ribosome 1 (x ≈ 0..10):
          SSU chain ``S1`` at (0, 0, 0) — 3 residues
          LSU chain ``L1`` at (5, 0, 0) — 3 residues
          tRNA chain ``T1`` at (3, 0, 0) — 1 residue (between SSU & LSU)
      - Ribosome 2 (x ≈ 100..110):
          SSU chain ``S2`` at (100, 0, 0) — 3 residues
          LSU chain ``L2`` at (105, 0, 0) — 3 residues
          tRNA chain ``T2`` at (103, 0, 0) — 1 residue
    """
    structure = gemmi.Structure()
    structure.name = "MULTIRIBO"
    structure.cell = gemmi.UnitCell(500.0, 500.0, 500.0, 90.0, 90.0, 90.0)

    model = gemmi.Model("1")
    for chain_name, base_x in [("S1", 0.0), ("L1", 5.0)]:
        c = gemmi.Chain(chain_name)
        for i in range(3):
            _add_residue(c, "G", seqid=i + 1, x=base_x + i, y=0.0, z=0.0)
        model.add_chain(c)
    t1 = gemmi.Chain("T1")
    _add_residue(t1, "G", seqid=1, x=3.0, y=0.0, z=0.0)
    model.add_chain(t1)

    for chain_name, base_x in [("S2", 100.0), ("L2", 105.0)]:
        c = gemmi.Chain(chain_name)
        for i in range(3):
            _add_residue(c, "G", seqid=i + 1, x=base_x + i, y=0.0, z=0.0)
        model.add_chain(c)
    t2 = gemmi.Chain("T2")
    _add_residue(t2, "G", seqid=1, x=103.0, y=0.0, z=0.0)
    model.add_chain(t2)

    structure.add_model(model)
    return structure


def test_detect_multi_ribosome_returns_count_for_paired_bundle() -> None:
    by_role = {
        "ssu_main_rrna": [_chain_ref("S1"), _chain_ref("S2")],
        "lsu_main_rrna": [_chain_ref("L1"), _chain_ref("L2")],
    }
    assert multiribo.detect_multi_ribosome(by_role) == 2


def test_detect_multi_ribosome_returns_one_for_single_ribosome() -> None:
    by_role = {
        "ssu_main_rrna": [_chain_ref("S1")],
        "lsu_main_rrna": [_chain_ref("L1")],
    }
    assert multiribo.detect_multi_ribosome(by_role) == 1


def test_detect_multi_ribosome_returns_one_for_unequal_counts() -> None:
    """Asymmetric counts (e.g. fragmented LSU = 1 SSU + 3 LSU fragments)
    are not multi-ribosome bundles. These get routed to the fragmented-
    ribosome skip path in the api orchestrator (see SKIP_FRAGMENTED_RIBOSOME)."""
    by_role = {
        "ssu_main_rrna": [_chain_ref("S1")],
        "lsu_main_rrna": [_chain_ref("L1a"), _chain_ref("L1b"), _chain_ref("L1c")],
    }
    assert multiribo.detect_multi_ribosome(by_role) == 1


def test_detect_multi_ribosome_returns_one_for_missing_rrna() -> None:
    by_role: dict[str, list[ChainRef]] = {"ssu_main_rrna": [], "lsu_main_rrna": []}
    assert multiribo.detect_multi_ribosome(by_role) == 1


def test_pair_ssu_lsu_by_centroid_groups_nearby_chains(
    two_ribosome_structure: gemmi.Structure,
) -> None:
    """S1↔L1 (centroids 5 Å apart) should pair before S1↔L2 (~105 Å)."""
    ssu = [_chain_ref("S1"), _chain_ref("S2")]
    lsu = [_chain_ref("L1"), _chain_ref("L2")]
    pairs = multiribo.pair_ssu_lsu_by_centroid(two_ribosome_structure, ssu, lsu)
    assert len(pairs) == 2
    paired_ids = {(s.auth_asym_id, large.auth_asym_id) for s, large in pairs}
    assert paired_ids == {("S1", "L1"), ("S2", "L2")}


def test_pair_ssu_lsu_skips_chains_with_no_atoms(
    two_ribosome_structure: gemmi.Structure,
) -> None:
    """A ChainRef whose auth_asym_id isn't in the structure has no centroid
    and should be silently dropped from the pairing pool."""
    ssu = [_chain_ref("S1"), _chain_ref("PHANTOM")]
    lsu = [_chain_ref("L1"), _chain_ref("L2")]
    pairs = multiribo.pair_ssu_lsu_by_centroid(two_ribosome_structure, ssu, lsu)
    assert len(pairs) == 1
    assert pairs[0][0].auth_asym_id == "S1"
    assert pairs[0][1].auth_asym_id == "L1"


def test_partition_chains_by_ribosome_assigns_by_centroid(
    two_ribosome_structure: gemmi.Structure,
) -> None:
    """T1 (x=3) should go to ribosome 1, T2 (x=103) to ribosome 2."""
    pairs = [
        (_chain_ref("S1"), _chain_ref("L1")),
        (_chain_ref("S2"), _chain_ref("L2")),
    ]
    chains = [_chain_ref("T1"), _chain_ref("T2")]
    groups = multiribo.partition_chains_by_ribosome(two_ribosome_structure, pairs, chains)
    assert len(groups) == 2
    assert groups[0] == {"T1"}
    assert groups[1] == {"T2"}


def test_partition_chains_skips_phantom_chains(
    two_ribosome_structure: gemmi.Structure,
) -> None:
    pairs = [
        (_chain_ref("S1"), _chain_ref("L1")),
        (_chain_ref("S2"), _chain_ref("L2")),
    ]
    chains = [_chain_ref("T1"), _chain_ref("PHANTOM")]
    groups = multiribo.partition_chains_by_ribosome(two_ribosome_structure, pairs, chains)
    # T1 goes to ribosome 1; PHANTOM has no centroid and lands nowhere.
    assert groups[0] == {"T1"}
    assert groups[1] == set()


def test_filter_correspondence_keeps_only_allowed_chains() -> None:
    corr = CorrespondenceResult(
        reference_key="ssu_atrna",
        reference_units=["REF|1|A|G|530"],
        mapped_units=[
            "MULTIRIBO|1|S1|G|530",
            "MULTIRIBO|1|S2|G|530",
        ],
    )
    filtered = multiribo.filter_correspondence_for_chains(
        {"ssu_atrna": corr}, allowed_chains={"S1"}
    )
    assert filtered["ssu_atrna"].mapped_units == ["MULTIRIBO|1|S1|G|530"]
    assert filtered["ssu_atrna"].mapped_units_by_chain == {"S1": ["MULTIRIBO|1|S1|G|530"]}


def test_filter_correspondence_drops_malformed_units() -> None:
    corr = CorrespondenceResult(
        reference_key="ssu_atrna",
        reference_units=["REF|1|A|G|530"],
        mapped_units=["MULTIRIBO|1|S1|G|530", "malformed"],
    )
    filtered = multiribo.filter_correspondence_for_chains(
        {"ssu_atrna": corr}, allowed_chains={"S1"}
    )
    assert filtered["ssu_atrna"].mapped_units == ["MULTIRIBO|1|S1|G|530"]


def test_filter_correspondence_empty_result_when_no_match() -> None:
    corr = CorrespondenceResult(
        reference_key="lsu_etrna",
        reference_units=["REF|1|L|G|2793"],
        mapped_units=["MULTIRIBO|1|L1|G|2793"],
        warnings=["upstream_warning"],
    )
    filtered = multiribo.filter_correspondence_for_chains(
        {"lsu_etrna": corr}, allowed_chains={"L2"}
    )
    # The site key is preserved with empty mapped_units so downstream code
    # sees the same shape as a never-mapped anchor.
    assert filtered["lsu_etrna"].mapped_units == []
    assert filtered["lsu_etrna"].mapped_units_by_chain == {}
    # Warnings carry over so the per-ribosome annotation still surfaces them.
    assert filtered["lsu_etrna"].warnings == ["upstream_warning"]
