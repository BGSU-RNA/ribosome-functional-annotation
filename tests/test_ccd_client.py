"""Tests for the per-component CCD client (parent-base lookup)."""

from __future__ import annotations

from pathlib import Path

import httpx
import pytest
import respx

from ribosome_state_annotator.cache import Cache
from ribosome_state_annotator.ccd_client import (
    CCD_URL_TEMPLATE,
    _parse_ccd_cif,
    fetch_chem_comp,
)

# ---------------------------------------------------------------------------
# Realistic CCD CIF samples
# ---------------------------------------------------------------------------


_U8U_CCD = b"""\
data_U8U
#
_chem_comp.id                                    U8U
_chem_comp.name                                  "5-METHYLAMINOMETHYL-2-THIOURIDINE-5'-MONOPHOSPHATE"
_chem_comp.type                                  "RNA LINKING"
_chem_comp.formula                               "C11 H18 N3 O8 P S"
_chem_comp.mon_nstd_parent_comp_id               U
_chem_comp.formula_weight                        383.315
_chem_comp.one_letter_code                       U
_chem_comp.three_letter_code                     U8U
#
"""

_OMG_CCD = b"""\
data_OMG
#
_chem_comp.id                                    OMG
_chem_comp.name                                  "O2'-METHYLGUANOSINE-5'-MONOPHOSPHATE"
_chem_comp.type                                  "RNA LINKING"
_chem_comp.mon_nstd_parent_comp_id               G
_chem_comp.one_letter_code                       g
_chem_comp.three_letter_code                     OMG
#
"""

# Canonical U: parent_comp_id is "?" sentinel (= no non-standard parent)
_U_CCD = b"""\
data_U
#
_chem_comp.id                                    U
_chem_comp.name                                  "URIDINE-5'-MONOPHOSPHATE"
_chem_comp.type                                  "RNA LINKING"
_chem_comp.mon_nstd_parent_comp_id               ?
_chem_comp.one_letter_code                       U
_chem_comp.three_letter_code                     U
#
"""


# ---------------------------------------------------------------------------
# Parser unit tests
# ---------------------------------------------------------------------------


def test_parse_ccd_cif_extracts_parent_and_one_letter() -> None:
    info = _parse_ccd_cif("U8U", _U8U_CCD)
    assert info is not None
    assert info.id == "U8U"
    assert info.parent_comp_id == "U"
    assert info.one_letter_code == "U"
    assert info.type == "RNA LINKING"


def test_parse_ccd_cif_lowercase_one_letter() -> None:
    """OMG ships parent G + lowercase one-letter code 'g' (Gemmi convention)."""
    info = _parse_ccd_cif("OMG", _OMG_CCD)
    assert info is not None
    assert info.parent_comp_id == "G"
    assert info.one_letter_code == "g"


def test_parse_ccd_cif_canonical_residue_has_no_parent_via_question_mark() -> None:
    """Canonical bases like U use the '?' CIF sentinel for parent_comp_id."""
    info = _parse_ccd_cif("U", _U_CCD)
    assert info is not None
    assert info.parent_comp_id is None  # '?' becomes None
    assert info.one_letter_code == "U"


def test_parse_ccd_cif_rejects_garbage() -> None:
    assert _parse_ccd_cif("XXX", b"this is not cif") is None


# ---------------------------------------------------------------------------
# fetch_chem_comp HTTP + cache
# ---------------------------------------------------------------------------


@respx.mock
def test_fetch_chem_comp_returns_info_on_200() -> None:
    respx.get(CCD_URL_TEMPLATE.format(comp_id="U8U")).mock(
        return_value=httpx.Response(200, content=_U8U_CCD)
    )
    info = fetch_chem_comp("U8U")
    assert info is not None
    assert info.parent_comp_id == "U"


@respx.mock
def test_fetch_chem_comp_uppercases_comp_id() -> None:
    """Lowercase input should still hit the uppercase-URL endpoint."""
    respx.get(CCD_URL_TEMPLATE.format(comp_id="U8U")).mock(
        return_value=httpx.Response(200, content=_U8U_CCD)
    )
    info = fetch_chem_comp("u8u")
    assert info is not None
    assert info.id == "U8U"


@respx.mock
def test_fetch_chem_comp_returns_none_on_http_error() -> None:
    respx.get(CCD_URL_TEMPLATE.format(comp_id="ZZZ")).mock(return_value=httpx.Response(500))
    assert fetch_chem_comp("ZZZ") is None


@respx.mock
def test_fetch_chem_comp_returns_none_on_404() -> None:
    respx.get(CCD_URL_TEMPLATE.format(comp_id="QQQ")).mock(return_value=httpx.Response(404))
    assert fetch_chem_comp("QQQ") is None


@respx.mock
def test_fetch_chem_comp_returns_none_on_network_failure() -> None:
    respx.get(CCD_URL_TEMPLATE.format(comp_id="ZZZ")).mock(
        side_effect=httpx.ConnectError("offline")
    )
    assert fetch_chem_comp("ZZZ") is None


def test_fetch_chem_comp_empty_comp_id_returns_none() -> None:
    assert fetch_chem_comp("") is None
    assert fetch_chem_comp("   ") is None


@respx.mock
def test_fetch_chem_comp_caches_on_disk(tmp_path: Path) -> None:
    """Second call should hit the cache, not the network."""
    cache = Cache(tmp_path)
    route = respx.get(CCD_URL_TEMPLATE.format(comp_id="U8U")).mock(
        return_value=httpx.Response(200, content=_U8U_CCD)
    )

    info1 = fetch_chem_comp("U8U", cache=cache)
    info2 = fetch_chem_comp("U8U", cache=cache)
    assert info1 == info2
    assert route.call_count == 1
    assert cache.ccd_path("U8U").is_file()


@respx.mock
def test_fetch_chem_comp_loads_corrupt_cache_recovers(tmp_path: Path) -> None:
    """A corrupt cached CIF should yield None at parse time, but the
    function shouldn't crash."""
    cache = Cache(tmp_path)
    cache.put_ccd_cif("U8U", b"not valid cif data")
    # No respx mock — if the function tried to fetch, the test would fail.
    info = fetch_chem_comp("U8U", cache=cache)
    assert info is None


# ---------------------------------------------------------------------------
# Integration with _parent_base_info
# ---------------------------------------------------------------------------


@respx.mock
def test_parent_base_info_falls_through_to_ccd_for_unknown_gemmi_residue(
    tmp_path: Path,
) -> None:
    """U8U: Gemmi has no one_letter_code; CCD says parent=U."""
    from ribosome_state_annotator.trna_mrna import _parent_base_info

    respx.get(CCD_URL_TEMPLATE.format(comp_id="U8U")).mock(
        return_value=httpx.Response(200, content=_U8U_CCD)
    )
    parent, is_modified = _parent_base_info("U8U", cache=Cache(tmp_path))
    assert parent == "U"
    assert is_modified is True


@respx.mock
def test_parent_base_info_skips_ccd_when_no_cache_or_client(
    tmp_path: Path,
) -> None:
    """When cache + client are both None, no CCD fetch happens — even for
    Gemmi-unrecognized residues — and the heuristic is used."""
    from ribosome_state_annotator.trna_mrna import _parent_base_info

    # Don't mock anything; if a fetch happened it would raise.
    parent, is_modified = _parent_base_info("U8U")  # no cache, no client
    assert parent == "U"  # heuristic: first char of comp_id
    assert is_modified is True


def test_parent_base_info_canonical_skips_ccd_fetch(tmp_path: Path) -> None:
    """For canonical bases, Gemmi already returns a clean letter — no CCD
    fetch should be triggered. Test by passing a cache and confirming no
    file is written."""
    from ribosome_state_annotator.trna_mrna import _parent_base_info

    cache = Cache(tmp_path)
    parent, is_modified = _parent_base_info("U", cache=cache)
    assert parent == "U"
    assert is_modified is False
    assert not cache.ccd_path("U").is_file()


@respx.mock
def test_parent_base_info_handles_ccd_fetch_failure_gracefully(
    tmp_path: Path,
) -> None:
    """When the CCD fetch fails, fall through to the first-char heuristic."""
    from ribosome_state_annotator.trna_mrna import _parent_base_info

    respx.get(CCD_URL_TEMPLATE.format(comp_id="U8U")).mock(
        return_value=httpx.Response(500)
    )
    parent, is_modified = _parent_base_info("U8U", cache=Cache(tmp_path))
    assert parent == "U"  # heuristic
    assert is_modified is True


@pytest.mark.parametrize(
    "comp_id, expected_parent, expected_modified",
    [
        ("A", "A", False),
        ("C", "C", False),
        ("G", "G", False),
        ("U", "U", False),
        # These all go through Gemmi's lowercase-letter convention:
        ("PSU", "U", True),
        ("5MU", "U", True),
        ("7MG", "G", True),
        ("H2U", "U", True),
        ("MIA", "A", True),
    ],
)
def test_parent_base_info_gemmi_known_residues(
    comp_id: str, expected_parent: str, expected_modified: bool, tmp_path: Path
) -> None:
    """Sanity check: residues Gemmi knows about don't trigger the CCD path."""
    from ribosome_state_annotator.trna_mrna import _parent_base_info

    cache = Cache(tmp_path)
    parent, is_modified = _parent_base_info(comp_id, cache=cache)
    assert parent == expected_parent
    assert is_modified is expected_modified
    # None of these should have triggered a CCD cache write.
    assert not cache.ccd_path(comp_id).is_file()
