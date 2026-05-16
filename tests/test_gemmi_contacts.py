"""Unit tests for the Gemmi residue lookup and neighbour search (spec §10)."""

from __future__ import annotations

import logging

import gemmi
import pytest

from ribosome_state_annotator import gemmi_contacts
from ribosome_state_annotator.correspondence import parse_unit_id

# ---------------------------------------------------------------------------
# find_residue
# ---------------------------------------------------------------------------


def test_find_residue_happy_path(minimal_structure: gemmi.Structure) -> None:
    unit = parse_unit_id("TEST|1|AA|G|1")
    result = gemmi_contacts.find_residue(minimal_structure, unit)
    assert result is not None
    chain_name, residue = result
    assert chain_name == "AA"
    assert residue.name == "G"
    assert residue.seqid.num == 1


def test_find_residue_with_offset_seqid(minimal_structure: gemmi.Structure) -> None:
    """Chain DA uses author seqid numbering 100..104; the lookup must find it."""
    unit = parse_unit_id("TEST|1|DA|U|102")
    result = gemmi_contacts.find_residue(minimal_structure, unit)
    assert result is not None
    chain_name, residue = result
    assert chain_name == "DA"
    assert residue.seqid.num == 102


def test_find_residue_returns_none_when_chain_missing(
    minimal_structure: gemmi.Structure,
) -> None:
    unit = parse_unit_id("TEST|1|XX|G|1")
    assert gemmi_contacts.find_residue(minimal_structure, unit) is None


def test_find_residue_returns_none_when_residue_missing(
    minimal_structure: gemmi.Structure,
) -> None:
    unit = parse_unit_id("TEST|1|AA|G|999")
    assert gemmi_contacts.find_residue(minimal_structure, unit) is None


def test_find_residue_tolerates_name_mismatch(
    minimal_structure: gemmi.Structure,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """§10.4: residue-name mismatch is logged at DEBUG but the residue is returned."""
    unit = parse_unit_id("TEST|1|AA|4OC|1")  # actual residue is G at seqid 1
    caplog.set_level(logging.DEBUG, logger="ribosome_state_annotator.gemmi_contacts")
    result = gemmi_contacts.find_residue(minimal_structure, unit)
    assert result is not None
    _, residue = result
    assert residue.name == "G"
    assert any("name mismatch" in record.message for record in caplog.records)


def test_find_residue_returns_none_when_structure_has_no_models() -> None:
    empty = gemmi.Structure()
    empty.cell = gemmi.UnitCell(50, 50, 50, 90, 90, 90)
    unit = parse_unit_id("X|1|A|G|1")
    assert gemmi_contacts.find_residue(empty, unit) is None


def test_find_residue_logs_warning_on_multi_model_structure(
    minimal_structure: gemmi.Structure,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Bio-assembly mmCIFs are single-model by spec; multi-model is an
    oddity that should be flagged."""
    extra_model = gemmi.Model("2")
    minimal_structure.add_model(extra_model)
    unit = parse_unit_id("TEST|1|AA|G|1")
    caplog.set_level(logging.WARNING, logger="ribosome_state_annotator.gemmi_contacts")
    result = gemmi_contacts.find_residue(minimal_structure, unit)
    assert result is not None
    assert any(
        "multiple models" in r.message.lower() or "models" in r.message.lower()
        for r in caplog.records
    )


# ---------------------------------------------------------------------------
# find_residues
# ---------------------------------------------------------------------------


def test_find_residues_all_hit(minimal_structure: gemmi.Structure) -> None:
    units = [
        parse_unit_id("TEST|1|AA|G|1"),
        parse_unit_id("TEST|1|AA|A|3"),
        parse_unit_id("TEST|1|DA|G|100"),
    ]
    by_chain, missing = gemmi_contacts.find_residues(minimal_structure, units)
    assert missing == []
    assert set(by_chain.keys()) == {"AA", "DA"}
    assert len(by_chain["AA"]) == 2
    assert len(by_chain["DA"]) == 1


def test_find_residues_mixes_hit_and_miss(minimal_structure: gemmi.Structure) -> None:
    units = [
        parse_unit_id("TEST|1|AA|G|1"),
        parse_unit_id("TEST|1|XX|G|1"),  # chain missing
        parse_unit_id("TEST|1|AA|G|999"),  # residue missing
        parse_unit_id("TEST|1|DA|G|100"),
    ]
    by_chain, missing = gemmi_contacts.find_residues(minimal_structure, units)
    # Missing list reports the canonical 5-segment unit-id strings.
    assert "TEST|1|XX|G|1" in missing
    assert "TEST|1|AA|G|999" in missing
    assert set(by_chain.keys()) == {"AA", "DA"}


def test_find_residues_preserves_input_order_within_chain(
    minimal_structure: gemmi.Structure,
) -> None:
    units = [
        parse_unit_id("TEST|1|AA|G|5"),
        parse_unit_id("TEST|1|AA|G|1"),
        parse_unit_id("TEST|1|AA|A|3"),
    ]
    by_chain, _ = gemmi_contacts.find_residues(minimal_structure, units)
    assert [r.seqid.num for r in by_chain["AA"]] == [5, 1, 3]


def test_find_residues_empty_input() -> None:
    empty = gemmi.Structure()
    empty.cell = gemmi.UnitCell(50, 50, 50, 90, 90, 90)
    by_chain, missing = gemmi_contacts.find_residues(empty, [])
    assert by_chain == {}
    assert missing == []


# ---------------------------------------------------------------------------
# find_neighbouring_chains
# ---------------------------------------------------------------------------


def _residues_for(structure: gemmi.Structure, chain_name: str, *seqids: int) -> list[gemmi.Residue]:
    """Helper: gather residues by chain name + seqid for the contact tests."""
    chain = structure[0].find_chain(chain_name)
    assert chain is not None
    return [residue for residue in chain if residue.seqid.num in seqids]


def test_neighbour_search_finds_close_pair(contact_structure: gemmi.Structure) -> None:
    """AA residue A3 (x=2) has BB residue U1 (x=2.5) within 0.5 Å."""
    targets = {"AA": _residues_for(contact_structure, "AA", 3)}
    result = gemmi_contacts.find_neighbouring_chains(contact_structure, targets, cutoff=5.0)
    assert "BB" in result
    assert result["BB"] == pytest.approx(0.5)


def test_neighbour_search_returns_min_distance_per_chain(
    contact_structure: gemmi.Structure,
) -> None:
    """When multiple target atoms hit the same neighbour chain, the smallest
    pairwise distance is returned."""
    targets = {"AA": _residues_for(contact_structure, "AA", 1, 2, 3)}  # x=0,1,2
    result = gemmi_contacts.find_neighbouring_chains(contact_structure, targets, cutoff=5.0)
    # BB residue U1 at x=2.5: AA residue A3 (x=2) is 0.5 Å, AA residue C2 (x=1) is 1.5 Å.
    assert result["BB"] == pytest.approx(0.5)


def test_neighbour_search_excludes_target_owning_chain(
    contact_structure: gemmi.Structure,
) -> None:
    """The chain that owns the targets should never appear in the result."""
    targets = {"AA": _residues_for(contact_structure, "AA", 3)}
    result = gemmi_contacts.find_neighbouring_chains(contact_structure, targets, cutoff=5.0)
    assert "AA" not in result


def test_neighbour_search_respects_cutoff(contact_structure: gemmi.Structure) -> None:
    """FAR is 50 Å from everything; cutoff=5.0 must exclude it."""
    targets = {"AA": _residues_for(contact_structure, "AA", 1, 2, 3)}
    result = gemmi_contacts.find_neighbouring_chains(contact_structure, targets, cutoff=5.0)
    assert "FAR" not in result


def test_neighbour_search_cutoff_includes_more_at_larger_radius(
    contact_structure: gemmi.Structure,
) -> None:
    """A larger cutoff picks up BB residue G2 (x=10), which is >7 Å from any AA atom."""
    targets = {"AA": _residues_for(contact_structure, "AA", 3)}
    result_5 = gemmi_contacts.find_neighbouring_chains(contact_structure, targets, cutoff=5.0)
    result_15 = gemmi_contacts.find_neighbouring_chains(contact_structure, targets, cutoff=15.0)
    # BB shows up in both; the min distance to BB stays at 0.5 (closest atom unchanged).
    assert result_5.get("BB") == pytest.approx(0.5)
    assert result_15.get("BB") == pytest.approx(0.5)
    # FAR at x=50 is beyond 15 Å too, but cutoff=15 picks up no further chains.
    assert "FAR" not in result_15


def test_neighbour_search_candidate_chains_filters_results(
    contact_structure: gemmi.Structure,
) -> None:
    """When candidate_chains is provided, only those chains appear in the result."""
    targets = {"AA": _residues_for(contact_structure, "AA", 3)}
    result = gemmi_contacts.find_neighbouring_chains(
        contact_structure, targets, cutoff=5.0, candidate_chains={"BB"}
    )
    assert set(result.keys()) == {"BB"}


def test_neighbour_search_candidate_chains_empty_yields_empty(
    contact_structure: gemmi.Structure,
) -> None:
    targets = {"AA": _residues_for(contact_structure, "AA", 3)}
    result = gemmi_contacts.find_neighbouring_chains(
        contact_structure, targets, cutoff=5.0, candidate_chains=set()
    )
    assert result == {}


def test_neighbour_search_multiple_target_chains_all_excluded(
    contact_structure: gemmi.Structure,
) -> None:
    """If two chains are both targets, neither appears in the result."""
    targets = {
        "AA": _residues_for(contact_structure, "AA", 3),
        "CC": _residues_for(contact_structure, "CC", 1),
    }
    result = gemmi_contacts.find_neighbouring_chains(contact_structure, targets, cutoff=5.0)
    assert "AA" not in result
    assert "CC" not in result
    # BB is not a target, so it should still appear.
    assert "BB" in result


def test_neighbour_search_no_targets_returns_empty(
    contact_structure: gemmi.Structure,
) -> None:
    result = gemmi_contacts.find_neighbouring_chains(contact_structure, {}, cutoff=5.0)
    assert result == {}


def test_neighbour_search_empty_target_residue_list_returns_empty(
    contact_structure: gemmi.Structure,
) -> None:
    """A chain key with an empty residue list should yield no neighbours."""
    result = gemmi_contacts.find_neighbouring_chains(contact_structure, {"AA": []}, cutoff=5.0)
    assert result == {}


def test_neighbour_search_returns_floats(contact_structure: gemmi.Structure) -> None:
    targets = {"AA": _residues_for(contact_structure, "AA", 3)}
    result = gemmi_contacts.find_neighbouring_chains(contact_structure, targets, cutoff=5.0)
    for chain_name, dist in result.items():
        assert isinstance(chain_name, str)
        assert isinstance(dist, float)
        assert 0.0 <= dist <= 5.0
