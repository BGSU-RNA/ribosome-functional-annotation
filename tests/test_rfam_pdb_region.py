"""Unit tests for :mod:`rfam_pdb_region`."""

from __future__ import annotations

import gzip
from datetime import datetime, timedelta, timezone
from pathlib import Path

import httpx
import respx

from ribosome_state_annotator import rfam_pdb_region as rfam

# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------


def test_get_rfam_cache_dir_defaults_to_package_cache(tmp_path: Path) -> None:
    assert rfam.get_rfam_cache_dir(tmp_path) == tmp_path / "rfam"


def test_get_local_rfam_file_path(tmp_path: Path) -> None:
    assert rfam.get_local_rfam_file_path(tmp_path) == tmp_path / "rfam" / rfam.RFAM_PDB_REGION_FILENAME


# ---------------------------------------------------------------------------
# Metadata I/O
# ---------------------------------------------------------------------------


def test_metadata_round_trip(tmp_path: Path) -> None:
    metadata = rfam.RfamPdbRegionMetadata(
        source_url=rfam.RFAM_PDB_REGION_URL,
        downloaded_at=datetime(2026, 5, 19, 12, 0, tzinfo=timezone.utc),
        last_modified="Mon, 18 May 2026 09:37:40 GMT",
    )
    rfam.save_rfam_metadata(metadata, cache_root=tmp_path)
    loaded = rfam.load_rfam_metadata(cache_root=tmp_path)
    assert loaded == metadata


def test_load_rfam_metadata_missing_returns_none(tmp_path: Path) -> None:
    assert rfam.load_rfam_metadata(cache_root=tmp_path) is None


def test_load_rfam_metadata_malformed_returns_none(tmp_path: Path) -> None:
    path = rfam.get_local_rfam_metadata_path(cache_root=tmp_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("not json", encoding="utf-8")
    assert rfam.load_rfam_metadata(cache_root=tmp_path) is None


# ---------------------------------------------------------------------------
# Parsing + best-score selection
# ---------------------------------------------------------------------------


def _make_local_rfam_file(tmp_path: Path, body: bytes) -> None:
    """Write a fake gzipped Rfam file + metadata into the cache root."""
    cache_dir = rfam.get_rfam_cache_dir(tmp_path)
    cache_dir.mkdir(parents=True, exist_ok=True)
    (cache_dir / rfam.RFAM_PDB_REGION_FILENAME).write_bytes(body)
    metadata = rfam.RfamPdbRegionMetadata(
        source_url=rfam.RFAM_PDB_REGION_URL,
        downloaded_at=datetime.now(timezone.utc),
        last_modified="fixture",
    )
    rfam.save_rfam_metadata(metadata, cache_root=tmp_path)


def _gzip_lines(rows: list[str]) -> bytes:
    return gzip.compress(("\n".join(rows) + "\n").encode("utf-8"))


def test_load_rfam_pdb_region_table_parses_rows(tmp_path: Path) -> None:
    """Tab-separated rows parsed into (rfam_acc, pdb_id, chain, bit_score)."""
    rows = [
        "RF00177\t5j7l\tAA\t1\t1534\t1552.7\t0\t1\t1533\t0064f4\t1",
        "RF01959\t5j7l\tAA\t7\t1534\t1002.4\t5.1e-297\t1\t1478\t0064f4\t0",
        # Malformed row (missing fields) — silently skipped:
        "RF00005\t8any",
        # Malformed bit_score — silently skipped:
        "RF00005\t8any\tAx\t1\t67\tnot-a-float\t4.8e-16\t1\t71\t0064f4\t1",
    ]
    _make_local_rfam_file(tmp_path, _gzip_lines(rows))
    parsed = rfam.load_rfam_pdb_region_table(cache_root=tmp_path)
    assert parsed == [
        ("RF00177", "5j7l", "AA", 1552.7),
        ("RF01959", "5j7l", "AA", 1002.4),
    ]


def test_build_rfam_pdb_region_lookup_picks_highest_bit_score() -> None:
    """When a chain has multiple Rfam hits, the highest bit_score wins."""
    # 9B0S chain s2: PDBe-style over-annotation. RF01960 has the highest
    # bit_score (eukaryotic 18S), so it should win.
    rows = [
        ("RF00177", "9b0s", "s2", 456.3),
        ("RF01959", "9b0s", "s2", 581.4),
        ("RF01960", "9b0s", "s2", 1669.0),
        ("RF02542", "9b0s", "s2", 822.0),
    ]
    lookup = rfam.build_rfam_pdb_region_lookup(rows)
    assert lookup[("9b0s", "s2")] == "RF01960"


def test_build_rfam_pdb_region_lookup_picks_one_per_chain() -> None:
    """Single best Rfam per (pdb_id, chain). Other chains independent."""
    rows = [
        ("RF00177", "5j7l", "AA", 1552.7),
        ("RF01959", "5j7l", "AA", 1002.4),
        ("RF02540", "5j7l", "DA", 1860.2),
        ("RF02541", "5j7l", "DA", 2876.2),
    ]
    lookup = rfam.build_rfam_pdb_region_lookup(rows)
    assert lookup == {
        ("5j7l", "AA"): "RF00177",  # 1552.7 > 1002.4
        ("5j7l", "DA"): "RF02541",  # 2876.2 > 1860.2
    }


def test_get_rfam_for_chain_lookups_case_insensitive_pdb(tmp_path: Path) -> None:
    rows = ["RF00177\t5j7l\tAA\t1\t1534\t1552.7\t0\t1\t1533\t0064f4\t1"]
    _make_local_rfam_file(tmp_path, _gzip_lines(rows))
    dataset = rfam.load_rfam_pdb_region_dataset(cache_root=tmp_path)
    assert dataset is not None
    # PDB ID uppercased on query — should still find it.
    assert rfam.get_rfam_for_chain(dataset, "5J7L", "AA") == "RF00177"
    # Chain is case-sensitive — different casing should miss.
    assert rfam.get_rfam_for_chain(dataset, "5J7L", "aa") is None


def test_get_rfam_for_chain_returns_none_for_unknown_pdb(tmp_path: Path) -> None:
    rows = ["RF00177\t5j7l\tAA\t1\t1534\t1552.7\t0\t1\t1533\t0064f4\t1"]
    _make_local_rfam_file(tmp_path, _gzip_lines(rows))
    dataset = rfam.load_rfam_pdb_region_dataset(cache_root=tmp_path)
    assert rfam.get_rfam_for_chain(dataset, "0000", "AA") is None


def test_get_rfam_mapping_for_pdb_returns_pdbe_compatible_shape(
    tmp_path: Path,
) -> None:
    """The mapping output shape matches what pdbe_client.parse_rfam_mappings
    used to produce — drop-in compatible with api._apply_rfam_pdb_region."""
    rows = [
        "RF00177\t5j7l\tAA\t1\t1534\t1552.7\t0\t1\t1533\t0064f4\t1",
        "RF01959\t5j7l\tAA\t7\t1534\t1002.4\t5.1e-297\t1\t1478\t0064f4\t0",
        "RF02541\t5j7l\tDA\t1\t2904\t2876.2\t0\t1\t2925\t0064f4\t1",
    ]
    _make_local_rfam_file(tmp_path, _gzip_lines(rows))
    dataset = rfam.load_rfam_pdb_region_dataset(cache_root=tmp_path)
    mapping = rfam.get_rfam_mapping_for_pdb(dataset, "5J7L")
    assert mapping == {"AA": ["RF00177"], "DA": ["RF02541"]}


# ---------------------------------------------------------------------------
# Refresh orchestration
# ---------------------------------------------------------------------------


@respx.mock
def test_ensure_available_downloads_when_missing(tmp_path: Path) -> None:
    """No cached file → download from EBI."""
    body = _gzip_lines(["RF00177\t5j7l\tAA\t1\t1534\t1552.7\t0\t1\t1533\t0064f4\t1"])
    respx.head(rfam.RFAM_PDB_REGION_URL).mock(
        return_value=httpx.Response(200, headers={"Last-Modified": "Mon, 18 May 2026 09:37:40 GMT"})
    )
    respx.get(rfam.RFAM_PDB_REGION_URL).mock(
        return_value=httpx.Response(
            200,
            content=body,
            headers={"Last-Modified": "Mon, 18 May 2026 09:37:40 GMT"},
        )
    )
    metadata = rfam.ensure_rfam_pdb_region_available(cache_root=tmp_path)
    assert metadata is not None
    assert metadata.last_modified == "Mon, 18 May 2026 09:37:40 GMT"
    assert rfam.get_local_rfam_file_path(tmp_path).is_file()


@respx.mock
def test_ensure_available_keeps_cache_when_unchanged(tmp_path: Path) -> None:
    """Local file present, fresh, and upstream Last-Modified matches → no re-download."""
    last_modified = "Mon, 18 May 2026 09:37:40 GMT"
    _make_local_rfam_file(tmp_path, _gzip_lines(["RF00177\t5j7l\tAA\t1\t1534\t1552.7\t0\t1\t1533\t0064f4\t1"]))
    # Override metadata to set fresh download timestamp + matching Last-Modified.
    rfam.save_rfam_metadata(
        rfam.RfamPdbRegionMetadata(
            source_url=rfam.RFAM_PDB_REGION_URL,
            downloaded_at=datetime.now(timezone.utc),
            last_modified=last_modified,
        ),
        cache_root=tmp_path,
    )
    # Force-refresh so the HEAD probe fires even though the cache is fresh.
    respx.head(rfam.RFAM_PDB_REGION_URL).mock(
        return_value=httpx.Response(200, headers={"Last-Modified": last_modified})
    )
    get_route = respx.get(rfam.RFAM_PDB_REGION_URL).mock(
        return_value=httpx.Response(200, content=b"should not be downloaded")
    )
    metadata = rfam.ensure_rfam_pdb_region_available(
        cache_root=tmp_path, force_refresh=True
    )
    assert metadata is not None
    assert metadata.last_modified == last_modified
    assert get_route.call_count == 0  # cached version reused


@respx.mock
def test_ensure_available_refreshes_when_stale(tmp_path: Path) -> None:
    """Local file older than STALENESS_DAYS → check upstream; if changed, re-download."""
    old_last_modified = "Sun, 10 May 2026 00:00:00 GMT"
    new_last_modified = "Mon, 18 May 2026 09:37:40 GMT"
    _make_local_rfam_file(tmp_path, _gzip_lines(["RF00001\t1abc\tA\t1\t100\t100.0\t0\t1\t120\t0064f4\t1"]))
    # Set a fake old downloaded_at timestamp so the cache appears stale.
    old_dt = datetime.now(timezone.utc) - timedelta(days=rfam.STALENESS_DAYS + 1)
    rfam.save_rfam_metadata(
        rfam.RfamPdbRegionMetadata(
            source_url=rfam.RFAM_PDB_REGION_URL,
            downloaded_at=old_dt,
            last_modified=old_last_modified,
        ),
        cache_root=tmp_path,
    )
    respx.head(rfam.RFAM_PDB_REGION_URL).mock(
        return_value=httpx.Response(200, headers={"Last-Modified": new_last_modified})
    )
    new_body = _gzip_lines(["RF00177\t5j7l\tAA\t1\t1534\t1552.7\t0\t1\t1533\t0064f4\t1"])
    respx.get(rfam.RFAM_PDB_REGION_URL).mock(
        return_value=httpx.Response(200, content=new_body, headers={"Last-Modified": new_last_modified})
    )
    metadata = rfam.ensure_rfam_pdb_region_available(cache_root=tmp_path)
    assert metadata is not None
    assert metadata.last_modified == new_last_modified


@respx.mock
def test_ensure_available_returns_none_when_no_cache_and_offline(tmp_path: Path) -> None:
    respx.head(rfam.RFAM_PDB_REGION_URL).mock(
        side_effect=httpx.ConnectError("offline")
    )
    respx.get(rfam.RFAM_PDB_REGION_URL).mock(
        side_effect=httpx.ConnectError("offline")
    )
    metadata = rfam.ensure_rfam_pdb_region_available(cache_root=tmp_path)
    assert metadata is None
