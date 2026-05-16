"""Unit tests for the BGSU unit-ID parser (spec §10.3)."""

from __future__ import annotations

import pytest
from pydantic import ValidationError

from ribosome_state_annotator.correspondence import (
    UnitId,
    parse_unit_id,
    try_parse_unit_id,
)

# ---------------------------------------------------------------------------
# parse_unit_id — happy path
# ---------------------------------------------------------------------------


def test_parse_basic_unit_id() -> None:
    u = parse_unit_id("5J7L|1|AA|G|926")
    assert u.pdb_id == "5J7L"
    assert u.model == "1"
    assert u.chain == "AA"
    assert u.residue_name == "G"
    assert u.residue_number == 926


@pytest.mark.parametrize(
    ("unit_id", "expected_residue_name"),
    [
        ("5J7L|1|AA|4OC|1402", "4OC"),
        ("5J7L|1|DA|OMG|2251", "OMG"),
        ("5TBW|1|1|G|2922", "G"),
    ],
)
def test_parse_modified_nucleotide(unit_id: str, expected_residue_name: str) -> None:
    """Spec §10.3: the parser must handle modified nucleotides like 4OC/OMG."""
    u = parse_unit_id(unit_id)
    assert u.residue_name == expected_residue_name


def test_parse_numeric_chain_name() -> None:
    """Yeast LSU chains are named "1" (a numeric string), not 1 (an int)."""
    u = parse_unit_id("5TBW|1|1|G|2922")
    assert u.chain == "1"
    assert isinstance(u.chain, str)


def test_parse_multi_character_chain() -> None:
    """Human ribosome (6ZMI) uses multi-character chain names like S2 / L11."""
    assert parse_unit_id("6ZMI|1|S2|G|1207").chain == "S2"
    assert parse_unit_id("6ZMI|1|L11|U|123").chain == "L11"


def test_parse_lowercase_chain() -> None:
    """7K00 and other entries use lowercase chain identifiers; preserve case."""
    u = parse_unit_id("7K00|1|a|G|926")
    assert u.chain == "a"


def test_parse_negative_residue_number() -> None:
    """Some entries have negative author-sequence numbering (e.g. signal peptides)."""
    u = parse_unit_id("XXXX|1|A|G|-3")
    assert u.residue_number == -3


def test_parse_extended_form_ignores_trailing_segments() -> None:
    """BGSU's full unit-ID format has up to 9 segments; v1 only interprets the
    first 5 and silently drops the rest."""
    extended = "5J7L|1|AA|G|926|C1'|||1_555"
    u = parse_unit_id(extended)
    assert u.residue_name == "G"
    assert u.residue_number == 926
    # The base_unit_id strips the extended segments — useful as a canonical form.
    assert u.base_unit_id == "5J7L|1|AA|G|926"


def test_parse_six_segment_form() -> None:
    """6 segments (e.g. with atom suffix) parses fine; extra dropped."""
    u = parse_unit_id("5J7L|1|AA|G|926|C1'")
    assert u.residue_number == 926


# ---------------------------------------------------------------------------
# Spec §6.3 anchors — every reference unit must parse
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "unit_id",
    [
        # E. coli (5J7L) anchors from spec §6.3
        "5J7L|1|AA|G|926",
        "5J7L|1|AA|4OC|1402",
        "5J7L|1|AA|C|1403",
        "5J7L|1|AA|G|1338",
        "5J7L|1|AA|A|1492",
        "5J7L|1|AA|A|1493",
        "5J7L|1|DA|G|2553",
        "5J7L|1|DA|OMG|2251",
        # Yeast (5TBW) anchors with numeric LSU chain "1"
        "5TBW|1|A|A|1150",
        "5TBW|1|1|G|2922",
        "5TBW|1|1|G|2619",
    ],
)
def test_every_spec_6_3_anchor_parses(unit_id: str) -> None:
    """Regression: every functional-site anchor declared in spec §6.3 must
    parse with the v1 parser, otherwise the inference pipeline can't even
    start."""
    parsed = parse_unit_id(unit_id)
    assert parsed.base_unit_id == unit_id


# ---------------------------------------------------------------------------
# base_unit_id round-trips
# ---------------------------------------------------------------------------


def test_base_unit_id_roundtrips_basic() -> None:
    parsed = parse_unit_id("5J7L|1|AA|G|926")
    assert parsed.base_unit_id == "5J7L|1|AA|G|926"
    assert parse_unit_id(parsed.base_unit_id) == parsed


def test_base_unit_id_dropped_extended_segments() -> None:
    """parse(base_unit_id(parse(extended))) == parse(extended) modulo the
    extended segments — base_unit_id is the canonical 5-segment form."""
    extended = parse_unit_id("5J7L|1|AA|G|926|C1'|.||1_555")
    canonical = parse_unit_id(extended.base_unit_id)
    assert canonical == extended


# ---------------------------------------------------------------------------
# parse_unit_id — malformed inputs raise ValueError
# ---------------------------------------------------------------------------


def test_parse_empty_input_raises() -> None:
    with pytest.raises(ValueError, match="empty"):
        parse_unit_id("")


@pytest.mark.parametrize(
    "unit_id",
    [
        "5J7L",
        "5J7L|1",
        "5J7L|1|AA",
        "5J7L|1|AA|G",
    ],
)
def test_parse_too_few_segments_raises(unit_id: str) -> None:
    with pytest.raises(ValueError, match="need at least 5"):
        parse_unit_id(unit_id)


@pytest.mark.parametrize(
    ("unit_id", "match"),
    [
        ("|1|AA|G|926", "empty PDB"),
        ("5J7L||AA|G|926", "empty model"),
        ("5J7L|1||G|926", "empty chain"),
        ("5J7L|1|AA||926", "empty residue_name"),
        ("5J7L|1|AA|G|", "empty residue_number"),
    ],
)
def test_parse_empty_required_field_raises(unit_id: str, match: str) -> None:
    with pytest.raises(ValueError, match=match):
        parse_unit_id(unit_id)


def test_parse_non_integer_residue_number_raises() -> None:
    with pytest.raises(ValueError, match="non-integer residue_number"):
        parse_unit_id("5J7L|1|AA|G|926A")


def test_parse_residue_number_with_insertion_code_raises() -> None:
    """BGSU expects insertion codes in a separate segment, not concatenated."""
    with pytest.raises(ValueError, match="non-integer residue_number"):
        parse_unit_id("5J7L|1|AA|G|100A")


# ---------------------------------------------------------------------------
# try_parse_unit_id — forgiving wrapper
# ---------------------------------------------------------------------------


def test_try_parse_returns_unit_id_on_valid_input() -> None:
    result = try_parse_unit_id("5J7L|1|AA|G|926")
    assert result is not None
    assert result.pdb_id == "5J7L"


@pytest.mark.parametrize(
    "bad",
    [
        "",
        "garbage",
        "5J7L|1|AA",  # too few segments
        "|1|AA|G|926",  # empty PDB
        "5J7L|1|AA|G|abc",  # non-integer residue_number
    ],
)
def test_try_parse_returns_none_on_malformed(bad: str) -> None:
    assert try_parse_unit_id(bad) is None


# ---------------------------------------------------------------------------
# UnitId model behaviour
# ---------------------------------------------------------------------------


def test_unit_id_is_frozen() -> None:
    """UnitIds are identity objects — they should not be mutable after creation."""
    u = parse_unit_id("5J7L|1|AA|G|926")
    with pytest.raises(ValidationError):
        u.residue_number = 999  # type: ignore[misc]


def test_unit_id_equality_is_structural() -> None:
    a = parse_unit_id("5J7L|1|AA|G|926")
    b = parse_unit_id("5J7L|1|AA|G|926")
    c = parse_unit_id("5J7L|1|AA|G|927")
    assert a == b
    assert a != c


def test_unit_id_is_hashable() -> None:
    """Frozen Pydantic models are hashable — useful for set/dict membership."""
    a = parse_unit_id("5J7L|1|AA|G|926")
    b = parse_unit_id("5J7L|1|AA|G|926")
    s = {a, b}
    assert len(s) == 1


def test_unit_id_model_dump_includes_base_unit_id() -> None:
    """The base_unit_id computed_field shows up in model_dump output (so a
    serialised UnitId can be reconstructed without recomputing)."""
    u = parse_unit_id("5J7L|1|AA|G|926")
    dump = u.model_dump()
    assert dump["pdb_id"] == "5J7L"
    assert dump["base_unit_id"] == "5J7L|1|AA|G|926"


def test_direct_construction_rejects_non_integer_residue_number() -> None:
    """Pydantic field validation catches construction with a wrong-type residue_number."""
    with pytest.raises(ValidationError):
        UnitId(
            pdb_id="X",
            model="1",
            chain="A",
            residue_name="G",
            residue_number="not-an-int",  # type: ignore[arg-type]
        )
