"""Unit tests for :mod:`bgsu_nr_list`."""

from __future__ import annotations

import gzip
from datetime import datetime, timedelta, timezone
from pathlib import Path

import httpx
import respx

from ribosome_state_annotator import bgsu_nr_list as nr

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

_HEADER = (
    '"ec_id","ife_id","assembly_id","pdb_resolution","na_type","ec_rank","rfam",'
    '"standardized_name","pdb_release_date"'
)

_ROWS = [
    # 28UI: released after EBI's last Rfam scan; every chain must resolve.
    '"NR_all_56726.173","28UI|1|A","1","2.1","RNA","18","RF00177","SSU rRNA","2026-05-20"',
    '"NR_all_83717.188","28UI|1|a","1","2.1","RNA","22","RF02541","LSU rRNA","2026-05-20"',
    '"NR_all_10157.186","28UI|1|b","1","2.1","RNA","15","RF00001","5S rRNA","2026-05-20"',
    '"NR_all_35542.170","28UI|1|Z","1","2.1","RNA","14","RF00005","tRNA","2026-05-20"',
    # Composite IFE: eukaryotic 25S + 5.8S joined; families zip positionally.
    '"NR_all_11111.1","4V88|1|A3+4V88|1|A4","1","3.0","RNA","1","RF02543+RF00002","LSU rRNA","2014-07-09"',
    # NA family: chain present but unmapped — must not appear in the lookup.
    '"NR_all_22222.1","9Q3Q|1|v","1","2.01","RNA","1","NA","mRNA","2026-08-26"',
    # Composite with NA members: only the RF00005 chain is kept.
    '"NR_all_33333.1","8XYZ|1|Q+8XYZ|1|R+8XYZ|1|S","1","3.2","RNA","1","RF00005+NA+NA","tRNA","2025-01-01"',
    # Starred family (BGSU marks proposed mappings with *): asterisk stripped.
    '"NR_all_44444.1","7ABC|1|X","1","2.8","RNA","1","RF02348*","other","2021-05-05"',
    # Mismatched chain/family counts: skipped.
    '"NR_all_55555.1","6BAD|1|A+6BAD|1|B","1","3.5","RNA","1","RF00177","broken","2020-01-01"',
]


def _csv_text(rows: list[str] | None = None) -> str:
    return "\n".join([_HEADER, *(rows if rows is not None else _ROWS)]) + "\n"


def _make_local_file(tmp_path: Path, text: str, *, release: str = "4.57") -> nr.BgsuNrListMetadata:
    cache_dir = nr.get_bgsu_nr_cache_dir(tmp_path)
    cache_dir.mkdir(parents=True, exist_ok=True)
    (cache_dir / nr.BGSU_NR_FILENAME).write_bytes(gzip.compress(text.encode("utf-8")))
    metadata = nr.BgsuNrListMetadata(
        source_url=nr.bgsu_nr_download_url(release),
        downloaded_at=datetime.now(timezone.utc).replace(microsecond=0),
        release=release,
        release_date="2026-09-16",
    )
    nr.save_bgsu_nr_metadata(metadata, cache_root=tmp_path)
    return metadata


_RELEASE_PAGE = """
<html><head><title>Representative set 4.57</title></head>
<body><h1>Release 4.57, 2026-09-16</h1>
<a href="https://rna.bgsu.edu/rna3dhub/nrlist/download/rna/4.57/all/csv">csv</a>
</body></html>
"""


# ---------------------------------------------------------------------------
# Path helpers + metadata I/O
# ---------------------------------------------------------------------------


def test_cache_paths(tmp_path: Path) -> None:
    assert nr.get_bgsu_nr_cache_dir(tmp_path) == tmp_path / "bgsu_nr"
    assert nr.get_local_bgsu_nr_file_path(tmp_path) == tmp_path / "bgsu_nr" / nr.BGSU_NR_FILENAME
    assert nr.bgsu_nr_download_url("4.57") == (
        "https://rna.bgsu.edu/rna3dhub/nrlist/download/rna/4.57/all/csv/full"
    )


def test_metadata_round_trip(tmp_path: Path) -> None:
    metadata = nr.BgsuNrListMetadata(
        source_url=nr.bgsu_nr_download_url("4.57"),
        downloaded_at=datetime(2026, 9, 17, 9, 0, tzinfo=timezone.utc),
        release="4.57",
        release_date="2026-09-16",
    )
    nr.save_bgsu_nr_metadata(metadata, cache_root=tmp_path)
    assert nr.load_bgsu_nr_metadata(cache_root=tmp_path) == metadata


def test_load_metadata_missing_or_malformed(tmp_path: Path) -> None:
    assert nr.load_bgsu_nr_metadata(cache_root=tmp_path) is None
    path = nr.get_local_bgsu_nr_metadata_path(cache_root=tmp_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("not json", encoding="utf-8")
    assert nr.load_bgsu_nr_metadata(cache_root=tmp_path) is None
    path.write_text('{"source_url": "x"}', encoding="utf-8")  # missing release
    assert nr.load_bgsu_nr_metadata(cache_root=tmp_path) is None


# ---------------------------------------------------------------------------
# Release-page parsing
# ---------------------------------------------------------------------------


def test_parse_current_release_from_title_and_date() -> None:
    assert nr.parse_current_release(_RELEASE_PAGE) == ("4.57", "2026-09-16")


def test_parse_current_release_falls_back_to_download_links() -> None:
    html = '<a href="/nrlist/release/rna/4.58/all">all</a>'
    assert nr.parse_current_release(html) == ("4.58", None)


def test_parse_current_release_unparseable() -> None:
    assert nr.parse_current_release("<html>nothing here</html>") is None


# ---------------------------------------------------------------------------
# CSV parsing
# ---------------------------------------------------------------------------


def test_parse_rows_simple_and_composite_and_na() -> None:
    lookup = nr.parse_bgsu_nr_rows(_csv_text())
    assert lookup[("28ui", "A")] == "RF00177"
    assert lookup[("28ui", "a")] == "RF02541"
    assert lookup[("28ui", "b")] == "RF00001"
    assert lookup[("28ui", "Z")] == "RF00005"
    # Composite IFE zipped positionally.
    assert lookup[("4v88", "A3")] == "RF02543"
    assert lookup[("4v88", "A4")] == "RF00002"
    # NA dropped, including inside composites.
    assert ("9q3q", "v") not in lookup
    assert lookup[("8xyz", "Q")] == "RF00005"
    assert ("8xyz", "R") not in lookup
    assert ("8xyz", "S") not in lookup
    # Asterisk stripped.
    assert lookup[("7abc", "X")] == "RF02348"
    # Mismatched counts skipped entirely.
    assert ("6bad", "A") not in lookup
    assert ("6bad", "B") not in lookup


def test_parse_rows_missing_columns_returns_empty() -> None:
    assert nr.parse_bgsu_nr_rows('"ec_id","ife_id"\n"x","28UI|1|A"\n') == {}


def test_load_dataset_and_lookups(tmp_path: Path) -> None:
    metadata = _make_local_file(tmp_path, _csv_text())
    dataset = nr.load_bgsu_nr_list_dataset(cache_root=tmp_path)
    assert dataset is not None
    assert dataset.metadata == metadata
    assert nr.get_bgsu_rfam_for_chain(dataset, "28ui", "A") == "RF00177"
    assert nr.get_bgsu_rfam_for_chain(dataset, "28UI", "X") is None
    assert nr.get_bgsu_rfam_mapping_for_pdb(dataset, "28UI") == {
        "A": ["RF00177"],
        "a": ["RF02541"],
        "b": ["RF00001"],
        "Z": ["RF00005"],
    }
    assert nr.get_bgsu_rfam_mapping_for_pdb(dataset, "9Q3Q") == {}
    assert nr.get_bgsu_rfam_mapping_for_pdb(None, "28UI") == {}
    assert nr.get_bgsu_rfam_for_chain(None, "28UI", "A") is None


def test_load_dataset_without_cache_returns_none(tmp_path: Path) -> None:
    assert nr.load_bgsu_nr_list_dataset(cache_root=tmp_path) is None


def test_load_dataset_file_missing_returns_none(tmp_path: Path) -> None:
    metadata = _make_local_file(tmp_path, _csv_text())
    nr.get_local_bgsu_nr_file_path(tmp_path).unlink()
    assert nr.load_bgsu_nr_list_dataset(cache_root=tmp_path, metadata=metadata) is None


# ---------------------------------------------------------------------------
# Download + refresh orchestration (network mocked)
# ---------------------------------------------------------------------------


@respx.mock
def test_ensure_downloads_when_missing(tmp_path: Path) -> None:
    respx.get(nr.BGSU_NR_CURRENT_RELEASE_URL).mock(
        return_value=httpx.Response(200, text=_RELEASE_PAGE)
    )
    download = respx.get(nr.bgsu_nr_download_url("4.57")).mock(
        return_value=httpx.Response(200, text=_csv_text())
    )
    metadata = nr.ensure_bgsu_nr_list_available(cache_root=tmp_path)
    assert metadata is not None
    assert metadata.release == "4.57"
    assert metadata.release_date == "2026-09-16"
    assert download.called
    dataset = nr.load_bgsu_nr_list_dataset(cache_root=tmp_path)
    assert dataset is not None
    assert dataset.lookup[("28ui", "A")] == "RF00177"
    assert nr.list_bgsu_nr_files(tmp_path) > 0


@respx.mock
def test_ensure_fresh_cache_makes_no_network_call(tmp_path: Path) -> None:
    _make_local_file(tmp_path, _csv_text())
    probe = respx.get(nr.BGSU_NR_CURRENT_RELEASE_URL).mock(
        return_value=httpx.Response(200, text=_RELEASE_PAGE)
    )
    metadata = nr.ensure_bgsu_nr_list_available(cache_root=tmp_path)
    assert metadata is not None and metadata.release == "4.57"
    assert not probe.called


@respx.mock
def test_ensure_stale_same_release_touches_metadata_only(tmp_path: Path) -> None:
    _make_local_file(tmp_path, _csv_text())
    probe = respx.get(nr.BGSU_NR_CURRENT_RELEASE_URL).mock(
        return_value=httpx.Response(200, text=_RELEASE_PAGE)
    )
    download = respx.get(nr.bgsu_nr_download_url("4.57")).mock(
        return_value=httpx.Response(200, text=_csv_text())
    )
    later = datetime.now(timezone.utc) + timedelta(days=nr.STALENESS_DAYS + 1)
    metadata = nr.ensure_bgsu_nr_list_available(cache_root=tmp_path, now=later)
    assert metadata is not None and metadata.release == "4.57"
    assert probe.called
    assert not download.called
    # downloaded_at was bumped so the clock resets.
    reloaded = nr.load_bgsu_nr_metadata(cache_root=tmp_path)
    assert reloaded is not None
    assert reloaded.downloaded_at >= metadata.downloaded_at - timedelta(seconds=1)


@respx.mock
def test_ensure_stale_new_release_downloads(tmp_path: Path) -> None:
    _make_local_file(tmp_path, _csv_text(), release="4.56")
    respx.get(nr.BGSU_NR_CURRENT_RELEASE_URL).mock(
        return_value=httpx.Response(200, text=_RELEASE_PAGE)
    )
    download = respx.get(nr.bgsu_nr_download_url("4.57")).mock(
        return_value=httpx.Response(200, text=_csv_text())
    )
    metadata = nr.ensure_bgsu_nr_list_available(cache_root=tmp_path, force_refresh=True)
    assert metadata is not None and metadata.release == "4.57"
    assert download.called


@respx.mock
def test_ensure_probe_failure_keeps_cached_file(tmp_path: Path) -> None:
    cached = _make_local_file(tmp_path, _csv_text())
    respx.get(nr.BGSU_NR_CURRENT_RELEASE_URL).mock(return_value=httpx.Response(503))
    metadata = nr.ensure_bgsu_nr_list_available(cache_root=tmp_path, force_refresh=True)
    assert metadata == cached


@respx.mock
def test_ensure_probe_failure_without_cache_returns_none(tmp_path: Path) -> None:
    respx.get(nr.BGSU_NR_CURRENT_RELEASE_URL).mock(side_effect=httpx.ConnectError("offline"))
    assert nr.ensure_bgsu_nr_list_available(cache_root=tmp_path) is None


@respx.mock
def test_ensure_download_failure_returns_none_when_no_cache(tmp_path: Path) -> None:
    respx.get(nr.BGSU_NR_CURRENT_RELEASE_URL).mock(
        return_value=httpx.Response(200, text=_RELEASE_PAGE)
    )
    respx.get(nr.bgsu_nr_download_url("4.57")).mock(return_value=httpx.Response(500))
    assert nr.ensure_bgsu_nr_list_available(cache_root=tmp_path) is None
