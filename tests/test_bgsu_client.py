"""Unit tests for the BGSU correspondence HTTP client (spec §5.2)."""

from __future__ import annotations

from typing import Any

import httpx
import pytest
import respx

from ribosome_state_annotator import bgsu_client
from ribosome_state_annotator.exceptions import (
    ApiRequestError,
    CorrespondenceMappingError,
)

# The §5.2.1 example response, verbatim.
SPEC_5_2_1_RESPONSE: dict[str, Any] = {
    "query": ["5J7L|1|AA|G|926", "5J7L|1|AA|4OC|1402", "5J7L|1|AA|C|1403"],
    "alignment": [
        {
            "reference_unit": "5J7L|1|AA|G|926",
            "mapped_units": [
                "5J7L|1|AA|G|926",
                "7K00|1|a|G|926",
                "4V51|1|AA|G|926",
                "6ZMI|1|S2|G|1207",
                "5TBW|1|A|A|1150",
            ],
        },
        {
            "reference_unit": "5J7L|1|AA|4OC|1402",
            "mapped_units": [
                "5J7L|1|AA|4OC|1402",
                "7K00|1|a|C|1402",
                "6ZMI|1|S2|C|1703",
            ],
        },
    ],
}


# ---------------------------------------------------------------------------
# fetch_correspondence — HTTP layer
# ---------------------------------------------------------------------------


@respx.mock
def test_fetch_correspondence_sends_correct_url_and_params() -> None:
    route = respx.get(bgsu_client.BGSU_CORRESPONDENCE_URL).mock(
        return_value=httpx.Response(200, json=SPEC_5_2_1_RESPONSE)
    )
    result = bgsu_client.fetch_correspondence(
        ["5J7L|1|AA|G|926", "5J7L|1|AA|4OC|1402", "5J7L|1|AA|C|1403"]
    )
    assert route.called
    request = route.calls.last.request
    assert request.url.params["id"] == "5J7L|1|AA|G|926,5J7L|1|AA|4OC|1402,5J7L|1|AA|C|1403"
    assert request.url.params["scope"] == "Rfam"
    assert request.url.params["resolution"] == "20.0A"
    assert request.url.params["depth"] == str(bgsu_client.DEFAULT_BGSU_DEPTH)
    assert request.url.params["format"] == "json"
    assert request.headers["accept"] == "application/json"
    assert "5J7L|1|AA|G|926" in result
    assert "7K00|1|a|G|926" in result["5J7L|1|AA|G|926"]


@respx.mock
def test_fetch_correspondence_honours_custom_depth() -> None:
    respx.get(bgsu_client.BGSU_CORRESPONDENCE_URL).mock(
        return_value=httpx.Response(200, json={"alignment": []})
    )
    bgsu_client.fetch_correspondence(["5J7L|1|AA|G|926"], depth=1500)
    request = respx.routes[0].calls.last.request
    assert request.url.params["depth"] == "1500"


@respx.mock
def test_fetch_correspondence_honours_custom_scope_and_resolution() -> None:
    respx.get(bgsu_client.BGSU_CORRESPONDENCE_URL).mock(
        return_value=httpx.Response(200, json={"alignment": []})
    )
    bgsu_client.fetch_correspondence(
        ["5J7L|1|AA|G|926"],
        scope="EC",
        resolution="3.0A",
    )
    request = respx.routes[0].calls.last.request
    assert request.url.params["scope"] == "EC"
    assert request.url.params["resolution"] == "3.0A"


def test_fetch_correspondence_empty_units_raises_value_error() -> None:
    with pytest.raises(ValueError, match="at least one unit ID"):
        bgsu_client.fetch_correspondence([])


@respx.mock
def test_fetch_correspondence_raises_on_http_500() -> None:
    respx.get(bgsu_client.BGSU_CORRESPONDENCE_URL).mock(return_value=httpx.Response(500))
    with pytest.raises(ApiRequestError, match="HTTP 500"):
        bgsu_client.fetch_correspondence(["5J7L|1|AA|G|926"])


@respx.mock
def test_fetch_correspondence_raises_on_network_error() -> None:
    respx.get(bgsu_client.BGSU_CORRESPONDENCE_URL).mock(
        side_effect=httpx.ConnectError("connection refused")
    )
    with pytest.raises(ApiRequestError, match="request failed"):
        bgsu_client.fetch_correspondence(["5J7L|1|AA|G|926"])


@respx.mock
def test_fetch_correspondence_raises_on_non_json_body() -> None:
    respx.get(bgsu_client.BGSU_CORRESPONDENCE_URL).mock(
        return_value=httpx.Response(200, text="<html>oh no</html>")
    )
    with pytest.raises(ApiRequestError, match="non-JSON"):
        bgsu_client.fetch_correspondence(["5J7L|1|AA|G|926"])


@respx.mock
def test_fetch_correspondence_raises_on_non_object_json() -> None:
    respx.get(bgsu_client.BGSU_CORRESPONDENCE_URL).mock(
        return_value=httpx.Response(200, json=[1, 2, 3])
    )
    with pytest.raises(ApiRequestError, match="non-object JSON"):
        bgsu_client.fetch_correspondence(["5J7L|1|AA|G|926"])


@respx.mock
def test_fetch_correspondence_accepts_caller_supplied_client() -> None:
    respx.get(bgsu_client.BGSU_CORRESPONDENCE_URL).mock(
        return_value=httpx.Response(200, json=SPEC_5_2_1_RESPONSE)
    )
    with httpx.Client(timeout=5.0) as client:
        result = bgsu_client.fetch_correspondence(["5J7L|1|AA|G|926"], client=client)
    assert "5J7L|1|AA|G|926" in result


# ---------------------------------------------------------------------------
# parse_alignment_response — pure
# ---------------------------------------------------------------------------


def test_parse_spec_5_2_1_example() -> None:
    parsed = bgsu_client.parse_alignment_response(SPEC_5_2_1_RESPONSE)
    assert set(parsed.keys()) == {"5J7L|1|AA|G|926", "5J7L|1|AA|4OC|1402"}
    assert parsed["5J7L|1|AA|G|926"] == [
        "5J7L|1|AA|G|926",
        "7K00|1|a|G|926",
        "4V51|1|AA|G|926",
        "6ZMI|1|S2|G|1207",
        "5TBW|1|A|A|1150",
    ]


def test_parse_empty_alignment_is_not_an_error() -> None:
    parsed = bgsu_client.parse_alignment_response({"alignment": []})
    assert parsed == {}


def test_parse_missing_alignment_raises_correspondence_error() -> None:
    with pytest.raises(
        CorrespondenceMappingError,
        match="missing required 'mappings' or 'alignment'",
    ):
        bgsu_client.parse_alignment_response({"query": []})


def test_parse_non_dict_payload_raises_correspondence_error() -> None:
    with pytest.raises(CorrespondenceMappingError, match="not a JSON object"):
        bgsu_client.parse_alignment_response("string-payload")  # type: ignore[arg-type]


def test_parse_tolerates_extra_top_level_fields() -> None:
    """Future BGSU responses may carry extra top-level fields — those must not
    break the parser."""
    payload = {
        "alignment": [{"reference_unit": "5J7L|1|AA|G|926", "mapped_units": ["7K00|1|a|G|926"]}],
        "computed_at": "2024-01-01",
        "version": "v2",
    }
    parsed = bgsu_client.parse_alignment_response(payload)
    assert parsed == {"5J7L|1|AA|G|926": ["7K00|1|a|G|926"]}


def test_parse_tolerates_extra_per_row_fields() -> None:
    payload = {
        "alignment": [
            {
                "reference_unit": "5J7L|1|AA|G|926",
                "mapped_units": ["7K00|1|a|G|926"],
                "alignment_score": 0.98,
                "rfam": "RF00177",
            }
        ]
    }
    assert bgsu_client.parse_alignment_response(payload) == {"5J7L|1|AA|G|926": ["7K00|1|a|G|926"]}


def test_parse_accepts_aligned_units_alias() -> None:
    """If BGSU renames mapped_units → aligned_units, the parser still works."""
    payload = {
        "alignment": [
            {
                "reference_unit": "5J7L|1|AA|G|926",
                "aligned_units": ["7K00|1|a|G|926"],
            }
        ]
    }
    assert bgsu_client.parse_alignment_response(payload) == {"5J7L|1|AA|G|926": ["7K00|1|a|G|926"]}


def test_parse_merges_duplicate_reference_unit_rows() -> None:
    payload = {
        "alignment": [
            {"reference_unit": "5J7L|1|AA|G|926", "mapped_units": ["A|1|x|G|1"]},
            {"reference_unit": "5J7L|1|AA|G|926", "mapped_units": ["B|1|y|G|2"]},
        ]
    }
    parsed = bgsu_client.parse_alignment_response(payload)
    assert parsed["5J7L|1|AA|G|926"] == ["A|1|x|G|1", "B|1|y|G|2"]


def test_parse_skips_malformed_rows_without_failing() -> None:
    payload = {
        "alignment": [
            "not-a-dict",
            {"mapped_units": ["A|1|x|G|1"]},  # missing reference_unit
            {"reference_unit": "", "mapped_units": []},  # empty ref
            {"reference_unit": "5J7L|1|AA|G|926", "mapped_units": "not-a-list"},
            {"reference_unit": "5J7L|1|AA|G|926", "mapped_units": ["7K00|1|a|G|926"]},
        ]
    }
    parsed = bgsu_client.parse_alignment_response(payload)
    assert parsed == {"5J7L|1|AA|G|926": ["7K00|1|a|G|926"]}


def test_parse_filters_non_string_mapped_units() -> None:
    payload = {
        "alignment": [
            {
                "reference_unit": "5J7L|1|AA|G|926",
                "mapped_units": ["7K00|1|a|G|926", 42, None, ""],
            }
        ]
    }
    assert bgsu_client.parse_alignment_response(payload) == {"5J7L|1|AA|G|926": ["7K00|1|a|G|926"]}


# ---------------------------------------------------------------------------
# fetch_correspondence — timeout + retry policy
# ---------------------------------------------------------------------------


@respx.mock
def test_fetch_correspondence_retries_after_timeout_then_succeeds() -> None:
    """A cold BGSU query can exceed the timeout; the second attempt wins."""
    route = respx.get(bgsu_client.BGSU_CORRESPONDENCE_URL).mock(
        side_effect=[
            httpx.ReadTimeout("slow"),
            httpx.Response(200, json=SPEC_5_2_1_RESPONSE),
        ]
    )
    result = bgsu_client.fetch_correspondence(["5J7L|1|AA|G|926"], retries=2)
    assert route.call_count == 2
    assert "5J7L|1|AA|G|926" in result


@respx.mock
def test_fetch_correspondence_retries_on_http_5xx() -> None:
    route = respx.get(bgsu_client.BGSU_CORRESPONDENCE_URL).mock(
        side_effect=[
            httpx.Response(503),
            httpx.Response(200, json={"alignment": []}),
        ]
    )
    bgsu_client.fetch_correspondence(["5J7L|1|AA|G|926"], retries=1)
    assert route.call_count == 2


@respx.mock
def test_fetch_correspondence_exhausted_retries_raises_transient_error() -> None:
    route = respx.get(bgsu_client.BGSU_CORRESPONDENCE_URL).mock(
        side_effect=httpx.ReadTimeout("slow")
    )
    with pytest.raises(bgsu_client.BgsuTransientError, match="timed out after 7 s"):
        bgsu_client.fetch_correspondence(["5J7L|1|AA|G|926"], retries=2, timeout=7.0)
    # first attempt + 2 retries
    assert route.call_count == 3


@respx.mock
def test_fetch_correspondence_does_not_retry_non_transient_errors() -> None:
    route = respx.get(bgsu_client.BGSU_CORRESPONDENCE_URL).mock(return_value=httpx.Response(404))
    with pytest.raises(ApiRequestError, match="HTTP 404"):
        bgsu_client.fetch_correspondence(["5J7L|1|AA|G|926"], retries=3)
    assert route.call_count == 1


@respx.mock
def test_fetch_correspondence_retries_zero_means_single_attempt() -> None:
    route = respx.get(bgsu_client.BGSU_CORRESPONDENCE_URL).mock(
        side_effect=httpx.ConnectError("refused")
    )
    with pytest.raises(ApiRequestError, match="request failed"):
        bgsu_client.fetch_correspondence(["5J7L|1|AA|G|926"], retries=0)
    assert route.call_count == 1


def test_transient_error_is_an_api_request_error() -> None:
    """Existing ``except ApiRequestError`` handlers must keep catching it."""
    assert issubclass(bgsu_client.BgsuTransientError, ApiRequestError)


@respx.mock
def test_fetch_correspondence_applies_timeout_per_request_even_with_caller_client() -> None:
    """The timeout is applied on the request itself, so a caller-supplied
    client with a short default doesn't silently override it."""
    route = respx.get(bgsu_client.BGSU_CORRESPONDENCE_URL).mock(
        return_value=httpx.Response(200, json={"alignment": []})
    )
    with httpx.Client(timeout=1.0) as client:
        bgsu_client.fetch_correspondence(["5J7L|1|AA|G|926"], client=client, timeout=123.0)
    request_timeout = route.calls.last.request.extensions["timeout"]
    assert request_timeout["read"] == 123.0


def test_default_timeout_exceeds_measured_cold_query_latency() -> None:
    """Guard against regressing the default below the ~90 s cold-query
    latency measured against live BGSU (September 2026)."""
    assert bgsu_client.DEFAULT_BGSU_TIMEOUT >= 120.0
    assert bgsu_client.DEFAULT_BGSU_RETRIES >= 1
