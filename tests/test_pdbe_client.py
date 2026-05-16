"""Unit tests for the PDBe Rfam-mapping client."""

from __future__ import annotations

from typing import Any

import httpx
import pytest
import respx

from ribosome_state_annotator import pdbe_client
from ribosome_state_annotator.exceptions import ApiRequestError

# Synthesised after the live 5J7L response shape.
SAMPLE_5J7L_PAYLOAD: dict[str, Any] = {
    "5j7l": {
        "Rfam": {
            "RF00177": {
                "identifier": "Bacterial small subunit ribosomal RNA",
                "family": "SSU_rRNA_bacteria",
                "mappings": [
                    {"chain_id": "AA", "entity_id": 1},
                    {"chain_id": "BA", "entity_id": 1},
                ],
            },
            "RF02541": {
                "identifier": "Bacterial large subunit ribosomal RNA",
                "family": "LSU_rRNA_bacteria",
                "mappings": [
                    {"chain_id": "CA", "entity_id": 31},
                    {"chain_id": "DA", "entity_id": 54},
                ],
            },
            "RF00001": {
                "identifier": "5S ribosomal RNA",
                "family": "5S_rRNA",
                "mappings": [{"chain_id": "CB", "entity_id": 28}],
            },
        }
    }
}


# ---------------------------------------------------------------------------
# parse_rfam_mappings — pure
# ---------------------------------------------------------------------------


def test_parse_flattens_chain_to_accessions() -> None:
    mapping = pdbe_client.parse_rfam_mappings(SAMPLE_5J7L_PAYLOAD, pdb_id="5J7L")
    assert mapping == {
        "AA": ["RF00177"],
        "BA": ["RF00177"],
        "CA": ["RF02541"],
        "DA": ["RF02541"],
        "CB": ["RF00001"],
    }


def test_parse_uses_lowercased_pdb_id_for_lookup() -> None:
    """Top-level key in PDBe is lowercased; passing uppercase pdb_id still works."""
    mapping = pdbe_client.parse_rfam_mappings(SAMPLE_5J7L_PAYLOAD, pdb_id="5J7L")
    assert "AA" in mapping


def test_parse_missing_pdb_block_returns_empty() -> None:
    assert pdbe_client.parse_rfam_mappings({"someother": {}}, pdb_id="5J7L") == {}


def test_parse_missing_rfam_block_returns_empty() -> None:
    payload = {"5j7l": {"Pfam": {}}}
    assert pdbe_client.parse_rfam_mappings(payload, pdb_id="5J7L") == {}


def test_parse_skips_malformed_mapping_entries() -> None:
    payload = {
        "5j7l": {
            "Rfam": {
                "RF00177": {"mappings": ["not-a-dict", {"chain_id": "AA"}, {}]},
                "RF02541": "not-an-object",
            }
        }
    }
    mapping = pdbe_client.parse_rfam_mappings(payload, pdb_id="5J7L")
    assert mapping == {"AA": ["RF00177"]}


def test_parse_dedups_accessions_per_chain() -> None:
    """If a chain appears under the same Rfam family twice, only one
    accession lands in the output list."""
    payload = {
        "5j7l": {
            "Rfam": {
                "RF00177": {
                    "mappings": [
                        {"chain_id": "AA"},
                        {"chain_id": "AA"},  # dup
                    ]
                }
            }
        }
    }
    mapping = pdbe_client.parse_rfam_mappings(payload, pdb_id="5J7L")
    assert mapping == {"AA": ["RF00177"]}


def test_parse_chain_with_multiple_rfam_families() -> None:
    # Top-level key is lowercased to match the live PDBe response shape.
    payload = {
        "x": {
            "Rfam": {
                "RF00177": {"mappings": [{"chain_id": "AA"}]},
                "RF00001": {"mappings": [{"chain_id": "AA"}]},
            }
        }
    }
    mapping = pdbe_client.parse_rfam_mappings(payload, pdb_id="X")
    assert set(mapping["AA"]) == {"RF00177", "RF00001"}


# ---------------------------------------------------------------------------
# fetch_rfam_mappings — HTTP (mocked)
# ---------------------------------------------------------------------------


@respx.mock
def test_fetch_uses_lowercased_pdb_id_in_url() -> None:
    route = respx.get(pdbe_client.PDBE_RFAM_URL_TEMPLATE.format(pdb_id="5j7l")).mock(
        return_value=httpx.Response(200, json=SAMPLE_5J7L_PAYLOAD)
    )
    result = pdbe_client.fetch_rfam_mappings("5J7L")
    assert route.called
    assert result["AA"] == ["RF00177"]


@respx.mock
def test_fetch_404_is_not_an_error() -> None:
    """PDBe returns 404 for entries it doesn't recognise — treated as
    "no Rfam mappings" rather than a hard failure."""
    respx.get(pdbe_client.PDBE_RFAM_URL_TEMPLATE.format(pdb_id="xxxx")).mock(
        return_value=httpx.Response(404)
    )
    assert pdbe_client.fetch_rfam_mappings("XXXX") == {}


@respx.mock
def test_fetch_500_raises() -> None:
    respx.get(pdbe_client.PDBE_RFAM_URL_TEMPLATE.format(pdb_id="5j7l")).mock(
        return_value=httpx.Response(500)
    )
    with pytest.raises(ApiRequestError, match="HTTP 500"):
        pdbe_client.fetch_rfam_mappings("5J7L")


@respx.mock
def test_fetch_transport_error_raises() -> None:
    respx.get(pdbe_client.PDBE_RFAM_URL_TEMPLATE.format(pdb_id="5j7l")).mock(
        side_effect=httpx.ConnectError("down")
    )
    with pytest.raises(ApiRequestError, match="request failed"):
        pdbe_client.fetch_rfam_mappings("5J7L")


@respx.mock
def test_fetch_non_json_body_raises() -> None:
    respx.get(pdbe_client.PDBE_RFAM_URL_TEMPLATE.format(pdb_id="5j7l")).mock(
        return_value=httpx.Response(200, text="<html>oh no</html>")
    )
    with pytest.raises(ApiRequestError, match="non-JSON"):
        pdbe_client.fetch_rfam_mappings("5J7L")


@respx.mock
def test_fetch_non_object_json_raises() -> None:
    respx.get(pdbe_client.PDBE_RFAM_URL_TEMPLATE.format(pdb_id="5j7l")).mock(
        return_value=httpx.Response(200, json=[1, 2, 3])
    )
    with pytest.raises(ApiRequestError, match="non-object JSON"):
        pdbe_client.fetch_rfam_mappings("5J7L")


@respx.mock
def test_fetch_accepts_caller_supplied_client() -> None:
    respx.get(pdbe_client.PDBE_RFAM_URL_TEMPLATE.format(pdb_id="5j7l")).mock(
        return_value=httpx.Response(200, json=SAMPLE_5J7L_PAYLOAD)
    )
    with httpx.Client(timeout=5.0) as client:
        result = pdbe_client.fetch_rfam_mappings("5J7L", client=client)
    assert result["AA"] == ["RF00177"]
