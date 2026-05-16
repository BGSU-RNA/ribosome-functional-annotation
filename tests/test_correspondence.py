"""Unit tests for the §5.2.2 two-step filter and CorrespondenceResult builder."""

from __future__ import annotations

import httpx
import respx

from ribosome_state_annotator import bgsu_client, correspondence
from ribosome_state_annotator.models import CorrespondenceResult

# Coverage for the UnitId parser used by the filter below lives in
# ``tests/test_unit_id_parser.py``; this file focuses on the filter and
# CorrespondenceResult behaviour.


# ---------------------------------------------------------------------------
# filter_mapped_units (§5.2.2)
# ---------------------------------------------------------------------------


def test_filter_keeps_only_target_pdb_case_insensitive() -> None:
    mapped = [
        "5J7L|1|AA|G|926",  # keep
        "7K00|1|a|G|926",  # drop (other PDB)
        "5j7l|1|AA|U|927",  # keep (lowercase form of target)
        "4V51|1|AA|G|926",  # drop
    ]
    result = correspondence.filter_mapped_units(
        mapped,
        target_pdb_id="5J7L",
        assembly_chains={"AA"},
    )
    assert result == ["5J7L|1|AA|G|926", "5j7l|1|AA|U|927"]


def test_filter_drops_chains_not_in_assembly() -> None:
    """A PDB with multiple assemblies — only chains for assembly N may pass."""
    mapped = [
        "5FDV|1|AA|G|530",  # assembly 1 chain — keep
        "5FDV|1|CA|G|530",  # assembly 2 chain — drop
        "5FDV|1|BB|G|530",  # not in assembly 1 — drop
    ]
    result = correspondence.filter_mapped_units(
        mapped,
        target_pdb_id="5FDV",
        assembly_chains={"AA"},
    )
    assert result == ["5FDV|1|AA|G|530"]


def test_filter_dedups_preserving_input_order() -> None:
    mapped = [
        "5J7L|1|AA|G|926",
        "5J7L|1|AA|G|926",  # dup
        "5J7L|1|AA|U|927",
        "5J7L|1|AA|G|926",  # dup
    ]
    result = correspondence.filter_mapped_units(
        mapped, target_pdb_id="5J7L", assembly_chains={"AA"}
    )
    assert result == ["5J7L|1|AA|G|926", "5J7L|1|AA|U|927"]


def test_filter_skips_malformed_unit_ids() -> None:
    mapped = ["5J7L|1|AA|G|926", "garbage", "5J7L|1", "5J7L|1|AA|U|927"]
    result = correspondence.filter_mapped_units(
        mapped, target_pdb_id="5J7L", assembly_chains={"AA"}
    )
    assert result == ["5J7L|1|AA|G|926", "5J7L|1|AA|U|927"]


def test_filter_empty_input_yields_empty_output() -> None:
    assert (
        correspondence.filter_mapped_units([], target_pdb_id="5J7L", assembly_chains={"AA"}) == []
    )


def test_filter_empty_assembly_chains_drops_everything() -> None:
    mapped = ["5J7L|1|AA|G|926", "5J7L|1|BB|G|927"]
    assert (
        correspondence.filter_mapped_units(mapped, target_pdb_id="5J7L", assembly_chains=set())
        == []
    )


# ---------------------------------------------------------------------------
# group_by_chain
# ---------------------------------------------------------------------------


def test_group_by_chain_basic() -> None:
    units = [
        "5J7L|1|AA|G|530",
        "5J7L|1|AA|A|1492",
        "5J7L|1|DA|G|2553",
    ]
    grouped = correspondence.group_by_chain(units)
    assert grouped == {
        "AA": ["5J7L|1|AA|G|530", "5J7L|1|AA|A|1492"],
        "DA": ["5J7L|1|DA|G|2553"],
    }


def test_group_by_chain_preserves_input_order_within_chain() -> None:
    units = [
        "5J7L|1|AA|U|1",
        "5J7L|1|AA|G|2",
        "5J7L|1|AA|C|3",
    ]
    grouped = correspondence.group_by_chain(units)
    assert grouped["AA"] == units


def test_group_by_chain_skips_malformed() -> None:
    units = ["5J7L|1|AA|G|530", "garbage", "5J7L|1|DA|G|2553"]
    grouped = correspondence.group_by_chain(units)
    assert set(grouped.keys()) == {"AA", "DA"}


# ---------------------------------------------------------------------------
# build_correspondence_result
# ---------------------------------------------------------------------------


def test_build_result_happy_path() -> None:
    raw = {
        "5J7L|1|AA|G|530": ["5J7L|1|AA|G|530", "7K00|1|a|G|530"],
        "5J7L|1|AA|A|1492": ["5J7L|1|AA|A|1492"],
    }
    result = correspondence.build_correspondence_result(
        "ssu_atrna",
        ["5J7L|1|AA|G|530", "5J7L|1|AA|A|1492"],
        raw,
        target_pdb_id="5J7L",
        assembly_chains={"AA"},
    )
    assert isinstance(result, CorrespondenceResult)
    assert result.reference_key == "ssu_atrna"
    assert result.reference_units == ["5J7L|1|AA|G|530", "5J7L|1|AA|A|1492"]
    # 7K00 mapping dropped by PDB filter
    assert result.mapped_units == ["5J7L|1|AA|G|530", "5J7L|1|AA|A|1492"]
    assert result.mapped_units_by_chain == {"AA": ["5J7L|1|AA|G|530", "5J7L|1|AA|A|1492"]}
    assert result.warnings == []


def test_build_result_missing_anchor_yields_warning() -> None:
    """An anchor with no in-assembly mapping warns but does not skip the site."""
    raw = {
        "5J7L|1|AA|G|530": ["7K00|1|a|G|530"],  # different PDB — filtered out
        "5J7L|1|AA|A|1492": ["5J7L|1|AA|A|1492"],
    }
    result = correspondence.build_correspondence_result(
        "ssu_atrna",
        ["5J7L|1|AA|G|530", "5J7L|1|AA|A|1492"],
        raw,
        target_pdb_id="5J7L",
        assembly_chains={"AA"},
    )
    assert result.mapped_units == ["5J7L|1|AA|A|1492"]
    assert result.warnings == ["correspondence_missing_for_ssu_atrna_5J7L|1|AA|G|530"]


def test_build_result_warning_format_exact_spec_5_2_3() -> None:
    """Spec §5.2.3 specifies the exact warning format."""
    raw: dict[str, list[str]] = {"5J7L|1|AA|G|926": []}
    result = correspondence.build_correspondence_result(
        "ssu_mrna",
        ["5J7L|1|AA|G|926"],
        raw,
        target_pdb_id="5J7L",
        assembly_chains={"AA"},
    )
    assert result.warnings == ["correspondence_missing_for_ssu_mrna_5J7L|1|AA|G|926"]


def test_build_result_all_anchors_missing_yields_all_warnings() -> None:
    raw: dict[str, list[str]] = {}
    result = correspondence.build_correspondence_result(
        "ssu_etrna",
        ["5J7L|1|AA|G|693", "5J7L|1|AA|A|694"],
        raw,
        target_pdb_id="5J7L",
        assembly_chains={"AA"},
    )
    assert result.mapped_units == []
    assert result.mapped_units_by_chain == {}
    assert result.warnings == [
        "correspondence_missing_for_ssu_etrna_5J7L|1|AA|G|693",
        "correspondence_missing_for_ssu_etrna_5J7L|1|AA|A|694",
    ]


def test_build_result_dedups_across_anchors() -> None:
    """Adjacent reference units sometimes pull the same mapped nt — the
    result.mapped_units list dedupes across anchors."""
    raw = {
        "5J7L|1|AA|G|1338": ["5J7L|1|AA|G|1338", "5J7L|1|AA|A|1339"],
        "5J7L|1|AA|A|1339": ["5J7L|1|AA|A|1339", "5J7L|1|AA|C|1400"],
    }
    result = correspondence.build_correspondence_result(
        "ssu_ptrna",
        ["5J7L|1|AA|G|1338", "5J7L|1|AA|A|1339"],
        raw,
        target_pdb_id="5J7L",
        assembly_chains={"AA"},
    )
    assert result.mapped_units == [
        "5J7L|1|AA|G|1338",
        "5J7L|1|AA|A|1339",
        "5J7L|1|AA|C|1400",
    ]


def test_build_result_groups_by_chain_across_ssu_and_lsu() -> None:
    """Reference units from one site may map to multiple chains in the query
    assembly (rare for SSU/LSU but the grouping must still be correct)."""
    raw = {
        "5J7L|1|AA|G|926": ["5J7L|1|AA|G|926", "5J7L|1|DA|G|2553"],
    }
    result = correspondence.build_correspondence_result(
        "test_site",
        ["5J7L|1|AA|G|926"],
        raw,
        target_pdb_id="5J7L",
        assembly_chains={"AA", "DA"},
    )
    assert set(result.mapped_units_by_chain.keys()) == {"AA", "DA"}


# ---------------------------------------------------------------------------
# fetch_and_filter_correspondence — orchestrator (mocked HTTP)
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Chain-substitution fallback (multi-assembly / non-NR-representative case)
# ---------------------------------------------------------------------------


def test_try_chain_substitution_translates_chain_segment() -> None:
    sub = correspondence._try_chain_substitution(
        "5J7L|1|AA|G|926",
        target_pdb_id="5J7L",
        assembly_chains={"BA", "CA"},
        chain_substitution={"AA": "BA"},
    )
    assert sub == "5J7L|1|BA|G|926"


def test_try_chain_substitution_returns_none_when_no_mapping() -> None:
    sub = correspondence._try_chain_substitution(
        "5J7L|1|AA|G|926",
        target_pdb_id="5J7L",
        assembly_chains={"BA"},
        chain_substitution={"DA": "CA"},  # no entry for AA
    )
    assert sub is None


def test_try_chain_substitution_returns_none_when_target_chain_not_in_assembly() -> None:
    """Defensive: never produce a substituted unit ID whose chain isn't
    actually in the assembly being processed."""
    sub = correspondence._try_chain_substitution(
        "5J7L|1|AA|G|926",
        target_pdb_id="5J7L",
        assembly_chains={"CA"},  # BA is mapped, but not in the assembly
        chain_substitution={"AA": "BA"},
    )
    assert sub is None


def test_build_correspondence_result_substitution_kicks_in_on_empty_filter() -> None:
    """When BGSU returns only wrong-chain rows (e.g., BGSU's NR
    representative is on a different chain than the one we're asking about),
    chain_substitution rescues the anchor."""
    raw = {
        "5J7L|1|AA|G|926": ["5J7L|1|AA|G|926"],  # BGSU returned AA, but we want BA
    }
    result = correspondence.build_correspondence_result(
        "ssu_mrna",
        ["5J7L|1|AA|G|926"],
        raw,
        target_pdb_id="5J7L",
        assembly_chains={"BA"},  # assembly 2 chains, AA not present
        chain_substitution={"AA": "BA"},
    )
    assert result.mapped_units == ["5J7L|1|BA|G|926"]
    assert result.warnings == []  # no missing-anchor warning


def test_build_correspondence_result_substitution_does_not_override_valid_hits() -> None:
    """If BGSU's response already has an in-assembly hit, the substitution
    fallback should NOT replace or augment it."""
    raw = {
        "5J7L|1|AA|G|926": ["5J7L|1|AA|G|926"],
    }
    result = correspondence.build_correspondence_result(
        "ssu_mrna",
        ["5J7L|1|AA|G|926"],
        raw,
        target_pdb_id="5J7L",
        assembly_chains={"AA"},
        chain_substitution={"AA": "BA"},  # would map AA→BA, but AA already hit
    )
    assert result.mapped_units == ["5J7L|1|AA|G|926"]


@respx.mock
def test_fetch_and_filter_end_to_end() -> None:
    """Hits the §5.1.2 example response shape, filters down to the target
    assembly, returns a fully-populated CorrespondenceResult."""
    respx.get(bgsu_client.BGSU_CORRESPONDENCE_URL).mock(
        return_value=httpx.Response(
            200,
            json={
                "alignment": [
                    {
                        "reference_unit": "5J7L|1|AA|G|530",
                        "mapped_units": [
                            "5J7L|1|AA|G|530",
                            "7K00|1|a|G|530",  # other PDB → dropped
                            "5J7L|1|BB|G|530",  # other chain → dropped
                        ],
                    }
                ]
            },
        )
    )
    result = correspondence.fetch_and_filter_correspondence(
        "ssu_atrna",
        ["5J7L|1|AA|G|530"],
        target_pdb_id="5J7L",
        assembly_chains={"AA"},
    )
    assert result.reference_key == "ssu_atrna"
    assert result.mapped_units == ["5J7L|1|AA|G|530"]
    assert result.warnings == []


def test_fetch_and_filter_empty_reference_units_returns_empty_result() -> None:
    """Defensive: an empty site (zero anchors) returns a placeholder
    CorrespondenceResult rather than calling BGSU with no input."""
    result = correspondence.fetch_and_filter_correspondence(
        "ssu_atrna",
        [],
        target_pdb_id="5J7L",
        assembly_chains={"AA"},
    )
    assert result.mapped_units == []
    assert result.warnings == []
