"""Tests for the RADdb integration (raddb spec)."""

from __future__ import annotations

import json
from datetime import datetime, timedelta, timezone
from pathlib import Path

import httpx
import pytest
import respx

from ribosome_state_annotator.raddb import (
    RADDB_URL_TEMPLATE,
    RADdbDataset,
    RADdbMetadata,
    build_raddb_lookup,
    download_raddb_csv,
    ensure_raddb_available,
    find_latest_raddb_url,
    get_local_raddb_csv_path,
    get_local_raddb_metadata_path,
    get_motion_metrics,
    load_raddb_dataset,
    load_raddb_metadata,
    load_raddb_table,
    parse_rad_date_from_url,
    save_raddb_metadata,
)

# ---------------------------------------------------------------------------
# Test fixtures
# ---------------------------------------------------------------------------

# Minimal CSV with all required columns. The header line uses the exact
# RADdb column names — including the trailing periods on `body rot.` and
# `head rot.`, which the lookup builder reads verbatim.
_SAMPLE_CSV = (
    "RCSB,RCSB version,notes,LSU chain ID,SSU chain ID,body rot.,body tilt,"
    "body tilt dir.,body trans.,head rot.,head tilt,head tilt dir.,head trans.,"
    "LSU core res.,body core res.,head core res.,LSU RMSD,body RMSD,head RMSD\n"
    "5J7L,3.2,,DA,AA,5.8,1.3,106.7,1.349,9.4,1.8,43.9,0.489,2663,917,394,0.64,1.12,0.78\n"
    "5J7L,3.2,,CA,BA,-0.6,0.9,42.4,0.651,5.9,2.6,56.9,0.746,2768,1029,435,0.68,0.71,0.86\n"
    "7K00,1.0,,A,B,2.1,0.5,12.0,0.300,3.4,0.7,-10.0,0.200,2600,900,400,0.5,0.5,0.5\n"
)

_FROZEN_NOW = datetime(2026, 5, 17, 12, 0, 0, tzinfo=timezone.utc)


# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------


def test_get_local_raddb_csv_path_uses_default_cache_root(tmp_path: Path) -> None:
    # The autouse _isolate_raddb_cache fixture already redirects
    # DEFAULT_CACHE_ROOT to a per-test tmp dir.
    csv_path = get_local_raddb_csv_path()
    assert csv_path.name == "RADdb.LSUSSU.csv"
    assert csv_path.parent.name == "raddb"


def test_get_local_raddb_csv_path_respects_explicit_root(tmp_path: Path) -> None:
    csv_path = get_local_raddb_csv_path(cache_root=tmp_path)
    assert csv_path == tmp_path / "raddb" / "RADdb.LSUSSU.csv"


# ---------------------------------------------------------------------------
# Metadata round-trip
# ---------------------------------------------------------------------------


def test_save_and_load_raddb_metadata_roundtrips(tmp_path: Path) -> None:
    meta = RADdbMetadata(
        source_url="https://example.test/RADdb.20260508.LSUSSU.csv",
        downloaded_at=datetime(2026, 5, 12, 9, 30, 0, tzinfo=timezone.utc),
        rad_date="20260508",
    )
    save_raddb_metadata(meta, cache_root=tmp_path)
    loaded = load_raddb_metadata(cache_root=tmp_path)
    assert loaded == meta


def test_load_raddb_metadata_returns_none_when_missing(tmp_path: Path) -> None:
    assert load_raddb_metadata(cache_root=tmp_path) is None


def test_load_raddb_metadata_tolerates_bad_json(tmp_path: Path) -> None:
    meta_path = get_local_raddb_metadata_path(cache_root=tmp_path)
    meta_path.parent.mkdir(parents=True, exist_ok=True)
    meta_path.write_text("not-json{", encoding="utf-8")
    assert load_raddb_metadata(cache_root=tmp_path) is None


# ---------------------------------------------------------------------------
# CSV parsing + lookup
# ---------------------------------------------------------------------------


def _write_csv(tmp_path: Path, body: str = _SAMPLE_CSV) -> Path:
    csv_path = get_local_raddb_csv_path(cache_root=tmp_path)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    csv_path.write_text(body, encoding="utf-8")
    return csv_path


def test_load_raddb_table_reads_rows(tmp_path: Path) -> None:
    _write_csv(tmp_path)
    rows = load_raddb_table(cache_root=tmp_path)
    assert rows is not None
    assert len(rows) == 3
    assert rows[0]["RCSB"] == "5J7L"
    assert rows[0]["LSU chain ID"] == "DA"
    assert rows[0]["SSU chain ID"] == "AA"
    assert rows[0]["body rot."] == "5.8"
    assert rows[0]["head rot."] == "9.4"


def test_load_raddb_table_returns_none_when_missing(tmp_path: Path) -> None:
    assert load_raddb_table(cache_root=tmp_path) is None


def test_build_raddb_lookup_maps_columns_to_public_names() -> None:
    rows = [
        {
            "RCSB": "5J7L",
            "LSU chain ID": "DA",
            "SSU chain ID": "AA",
            "body rot.": "5.8",
            "head rot.": "9.4",
        }
    ]
    lookup, dupes = build_raddb_lookup(rows)
    assert dupes == frozenset()
    assert lookup[("5J7L", "DA", "AA")]["intersubunit_rotation"] == pytest.approx(5.8)
    assert lookup[("5J7L", "DA", "AA")]["ssu_head_rotation"] == pytest.approx(9.4)


def test_build_raddb_lookup_uppercases_pdb_id_but_keeps_chain_case() -> None:
    rows = [
        {
            "RCSB": "11bq",  # lowercase in source
            "LSU chain ID": "a",
            "SSU chain ID": "A",
            "body rot.": "-1.1",
            "head rot.": "1.2",
        }
    ]
    lookup, _ = build_raddb_lookup(rows)
    assert ("11BQ", "a", "A") in lookup
    assert ("11BQ", "A", "A") not in lookup  # chain case preserved


def test_build_raddb_lookup_flags_duplicate_keys() -> None:
    rows = [
        {
            "RCSB": "5J7L",
            "LSU chain ID": "DA",
            "SSU chain ID": "AA",
            "body rot.": "5.8",
            "head rot.": "9.4",
        },
        {
            "RCSB": "5J7L",
            "LSU chain ID": "DA",
            "SSU chain ID": "AA",
            "body rot.": "6.0",
            "head rot.": "9.5",
        },
    ]
    _, dupes = build_raddb_lookup(rows)
    assert ("5J7L", "DA", "AA") in dupes


def test_build_raddb_lookup_skips_rows_missing_pdb_or_chain() -> None:
    rows = [
        {"RCSB": "", "LSU chain ID": "A", "SSU chain ID": "B", "body rot.": "1.0", "head rot.": "2.0"},
        {"RCSB": "5J7L", "LSU chain ID": "", "SSU chain ID": "B", "body rot.": "1.0", "head rot.": "2.0"},
        {"RCSB": "5J7L", "LSU chain ID": "A", "SSU chain ID": "", "body rot.": "1.0", "head rot.": "2.0"},
    ]
    lookup, _ = build_raddb_lookup(rows)
    assert lookup == {}


def test_build_raddb_lookup_handles_blank_numeric_values() -> None:
    rows = [
        {
            "RCSB": "5J7L",
            "LSU chain ID": "DA",
            "SSU chain ID": "AA",
            "body rot.": "",  # blank
            "head rot.": "garbage",  # unparseable
        }
    ]
    lookup, _ = build_raddb_lookup(rows)
    entry = lookup[("5J7L", "DA", "AA")]
    assert entry["intersubunit_rotation"] is None
    assert entry["ssu_head_rotation"] is None


# ---------------------------------------------------------------------------
# load_raddb_dataset + missing-columns warning
# ---------------------------------------------------------------------------


def test_load_raddb_dataset_returns_dataset(tmp_path: Path) -> None:
    _write_csv(tmp_path)
    meta = RADdbMetadata(
        source_url="https://example.test/RADdb.20260508.LSUSSU.csv",
        downloaded_at=datetime(2026, 5, 8, tzinfo=timezone.utc),
        rad_date="20260508",
    )
    save_raddb_metadata(meta, cache_root=tmp_path)
    dataset = load_raddb_dataset(cache_root=tmp_path)
    assert dataset is not None
    assert dataset.metadata.rad_date == "20260508"
    assert ("5J7L", "DA", "AA") in dataset.lookup


def test_load_raddb_dataset_warns_on_missing_columns(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    # CSV missing the body rot. / head rot. columns entirely.
    bad_csv = "RCSB,LSU chain ID,SSU chain ID\n5J7L,DA,AA\n"
    _write_csv(tmp_path, body=bad_csv)
    save_raddb_metadata(
        RADdbMetadata(
            source_url="https://example.test/RADdb.20260508.LSUSSU.csv",
            downloaded_at=datetime(2026, 5, 8, tzinfo=timezone.utc),
            rad_date="20260508",
        ),
        cache_root=tmp_path,
    )
    with caplog.at_level("WARNING"):
        dataset = load_raddb_dataset(cache_root=tmp_path)
    assert dataset is None
    assert any("missing required columns" in r.message for r in caplog.records)


# ---------------------------------------------------------------------------
# get_motion_metrics
# ---------------------------------------------------------------------------


def _dataset_for(rows: list[dict[str, str]]) -> RADdbDataset:
    lookup, dupes = build_raddb_lookup(rows)
    meta = RADdbMetadata(
        source_url="x",
        downloaded_at=datetime(2026, 1, 1, tzinfo=timezone.utc),
        rad_date="20260101",
    )
    return RADdbDataset(metadata=meta, lookup=lookup, duplicate_keys=dupes)


def test_get_motion_metrics_returns_row_on_match() -> None:
    dataset = _dataset_for(
        [
            {
                "RCSB": "5J7L",
                "LSU chain ID": "DA",
                "SSU chain ID": "AA",
                "body rot.": "5.8",
                "head rot.": "9.4",
            }
        ]
    )
    metrics = get_motion_metrics(dataset, "5j7l", "DA", "AA")  # lowercase input PDB
    assert metrics == {"intersubunit_rotation": pytest.approx(5.8), "ssu_head_rotation": pytest.approx(9.4)}


def test_get_motion_metrics_returns_none_on_miss() -> None:
    dataset = _dataset_for(
        [
            {
                "RCSB": "5J7L",
                "LSU chain ID": "DA",
                "SSU chain ID": "AA",
                "body rot.": "5.8",
                "head rot.": "9.4",
            }
        ]
    )
    assert get_motion_metrics(dataset, "5J7L", "ZZ", "AA") is None


def test_get_motion_metrics_returns_none_when_dataset_none() -> None:
    assert get_motion_metrics(None, "5J7L", "A", "B") is None


def test_get_motion_metrics_returns_none_on_duplicate_keys(
    caplog: pytest.LogCaptureFixture,
) -> None:
    dataset = _dataset_for(
        [
            {
                "RCSB": "5J7L",
                "LSU chain ID": "DA",
                "SSU chain ID": "AA",
                "body rot.": "5.8",
                "head rot.": "9.4",
            },
            {
                "RCSB": "5J7L",
                "LSU chain ID": "DA",
                "SSU chain ID": "AA",
                "body rot.": "6.0",
                "head rot.": "9.5",
            },
        ]
    )
    with caplog.at_level("WARNING"):
        metrics = get_motion_metrics(dataset, "5J7L", "DA", "AA")
    assert metrics is None
    assert any("ambiguous" in r.message for r in caplog.records)


# ---------------------------------------------------------------------------
# find_latest_raddb_url (HEAD walk back)
# ---------------------------------------------------------------------------


@respx.mock
def test_find_latest_raddb_url_picks_first_200() -> None:
    today_url = RADDB_URL_TEMPLATE.format(rad_date="20260517")
    yesterday_url = RADDB_URL_TEMPLATE.format(rad_date="20260516")
    older_url = RADDB_URL_TEMPLATE.format(rad_date="20260508")

    respx.head(today_url).mock(return_value=httpx.Response(404))
    respx.head(yesterday_url).mock(return_value=httpx.Response(404))
    respx.head(older_url).mock(return_value=httpx.Response(200))

    # Catch-all for the days in between (all 404).
    respx.route(method="HEAD", host="radtool.rc.northeastern.edu").mock(
        return_value=httpx.Response(404)
    )

    url, rad_date = find_latest_raddb_url(today=_FROZEN_NOW) or (None, None)
    assert rad_date == "20260508"
    assert url == older_url


@respx.mock
def test_find_latest_raddb_url_returns_none_if_nothing_in_window() -> None:
    respx.route(method="HEAD", host="radtool.rc.northeastern.edu").mock(
        return_value=httpx.Response(404)
    )
    assert find_latest_raddb_url(today=_FROZEN_NOW, lookback_days=3) is None


# ---------------------------------------------------------------------------
# download_raddb_csv
# ---------------------------------------------------------------------------


@respx.mock
def test_download_raddb_csv_writes_file_and_metadata(tmp_path: Path) -> None:
    url = RADDB_URL_TEMPLATE.format(rad_date="20260508")
    respx.get(url).mock(return_value=httpx.Response(200, content=_SAMPLE_CSV.encode("utf-8")))

    csv_path, metadata = download_raddb_csv(url, "20260508", cache_root=tmp_path)
    assert csv_path == get_local_raddb_csv_path(cache_root=tmp_path)
    assert csv_path.read_text(encoding="utf-8") == _SAMPLE_CSV
    assert metadata.rad_date == "20260508"
    assert metadata.source_url == url

    # Metadata persisted on disk.
    meta_json = json.loads(get_local_raddb_metadata_path(cache_root=tmp_path).read_text())
    assert meta_json["rad_date"] == "20260508"


@respx.mock
def test_download_raddb_csv_raises_on_http_error(tmp_path: Path) -> None:
    url = RADDB_URL_TEMPLATE.format(rad_date="20260508")
    respx.get(url).mock(return_value=httpx.Response(500))
    with pytest.raises(httpx.HTTPError):
        download_raddb_csv(url, "20260508", cache_root=tmp_path)


# ---------------------------------------------------------------------------
# ensure_raddb_available (staleness + refresh orchestration)
# ---------------------------------------------------------------------------


@respx.mock
def test_ensure_raddb_available_downloads_when_no_local_file(tmp_path: Path) -> None:
    today_url = RADDB_URL_TEMPLATE.format(rad_date="20260517")
    respx.head(today_url).mock(return_value=httpx.Response(200))
    respx.get(today_url).mock(return_value=httpx.Response(200, content=_SAMPLE_CSV.encode("utf-8")))
    respx.route(method="HEAD", host="radtool.rc.northeastern.edu").mock(
        return_value=httpx.Response(404)
    )

    metadata = ensure_raddb_available(cache_root=tmp_path, now=_FROZEN_NOW)
    assert metadata is not None
    assert metadata.rad_date == "20260517"
    assert get_local_raddb_csv_path(cache_root=tmp_path).is_file()


@respx.mock
def test_ensure_raddb_available_uses_cache_when_fresh(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    _write_csv(tmp_path)
    meta = RADdbMetadata(
        source_url="https://example.test/RADdb.20260514.LSUSSU.csv",
        downloaded_at=_FROZEN_NOW - timedelta(days=2),
        rad_date="20260514",
    )
    save_raddb_metadata(meta, cache_root=tmp_path)

    with caplog.at_level("INFO"):
        result = ensure_raddb_available(cache_root=tmp_path, now=_FROZEN_NOW)
    assert result == meta
    # No HTTP calls when fresh.
    assert all(not call.url for call in respx.calls)
    assert any("using cached RADdb" in r.message for r in caplog.records)


@respx.mock
def test_ensure_raddb_available_refreshes_when_stale(tmp_path: Path) -> None:
    _write_csv(tmp_path)
    stale_meta = RADdbMetadata(
        source_url="https://example.test/RADdb.20260501.LSUSSU.csv",
        downloaded_at=_FROZEN_NOW - timedelta(days=14),
        rad_date="20260501",
    )
    save_raddb_metadata(stale_meta, cache_root=tmp_path)

    new_url = RADDB_URL_TEMPLATE.format(rad_date="20260517")
    new_csv = "RCSB,LSU chain ID,SSU chain ID,body rot.,head rot.\n5J7L,DA,AA,5.8,9.4\n"
    respx.head(new_url).mock(return_value=httpx.Response(200))
    respx.get(new_url).mock(return_value=httpx.Response(200, content=new_csv.encode("utf-8")))
    respx.route(method="HEAD", host="radtool.rc.northeastern.edu").mock(
        return_value=httpx.Response(404)
    )

    refreshed = ensure_raddb_available(cache_root=tmp_path, now=_FROZEN_NOW)
    assert refreshed is not None
    assert refreshed.rad_date == "20260517"


@respx.mock
def test_ensure_raddb_available_keeps_cache_when_refresh_fails(tmp_path: Path) -> None:
    _write_csv(tmp_path)
    stale_meta = RADdbMetadata(
        source_url="https://example.test/RADdb.20260501.LSUSSU.csv",
        downloaded_at=_FROZEN_NOW - timedelta(days=14),
        rad_date="20260501",
    )
    save_raddb_metadata(stale_meta, cache_root=tmp_path)
    # Online probe finds nothing.
    respx.route(method="HEAD", host="radtool.rc.northeastern.edu").mock(
        return_value=httpx.Response(404)
    )

    result = ensure_raddb_available(cache_root=tmp_path, now=_FROZEN_NOW)
    assert result == stale_meta  # cached metadata returned unchanged


@respx.mock
def test_ensure_raddb_available_returns_none_when_no_cache_and_no_network(
    tmp_path: Path,
) -> None:
    respx.route(method="HEAD", host="radtool.rc.northeastern.edu").mock(
        side_effect=httpx.ConnectError("offline")
    )
    assert ensure_raddb_available(cache_root=tmp_path, now=_FROZEN_NOW) is None


@respx.mock
def test_ensure_raddb_available_force_refresh_skips_staleness_check(tmp_path: Path) -> None:
    _write_csv(tmp_path)
    fresh_meta = RADdbMetadata(
        source_url="https://example.test/RADdb.20260514.LSUSSU.csv",
        downloaded_at=_FROZEN_NOW - timedelta(hours=1),
        rad_date="20260514",
    )
    save_raddb_metadata(fresh_meta, cache_root=tmp_path)

    new_url = RADDB_URL_TEMPLATE.format(rad_date="20260517")
    new_csv = "RCSB,LSU chain ID,SSU chain ID,body rot.,head rot.\n5J7L,DA,AA,5.8,9.4\n"
    respx.head(new_url).mock(return_value=httpx.Response(200))
    respx.get(new_url).mock(return_value=httpx.Response(200, content=new_csv.encode("utf-8")))
    respx.route(method="HEAD", host="radtool.rc.northeastern.edu").mock(
        return_value=httpx.Response(404)
    )

    refreshed = ensure_raddb_available(
        cache_root=tmp_path, now=_FROZEN_NOW, force_refresh=True
    )
    assert refreshed is not None
    assert refreshed.rad_date == "20260517"


# ---------------------------------------------------------------------------
# Misc helpers
# ---------------------------------------------------------------------------


def test_parse_rad_date_from_url() -> None:
    assert (
        parse_rad_date_from_url(
            "https://radtool.rc.northeastern.edu/public_database/RADdb.20260508.LSUSSU.csv"
        )
        == "20260508"
    )


def test_parse_rad_date_from_url_returns_none_on_no_match() -> None:
    assert parse_rad_date_from_url("https://example.test/nothing") is None
