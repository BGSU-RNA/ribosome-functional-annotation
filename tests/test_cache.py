"""Unit tests for the filesystem cache (spec §17)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from ribosome_state_annotator.cache import DEFAULT_CACHE_ROOT, Cache, CacheInfo

# ---------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------


def test_cache_default_root_points_at_user_cache_dir() -> None:
    cache = Cache()
    assert cache.root == DEFAULT_CACHE_ROOT
    assert str(cache.root).endswith("ribosome-state-annotator")


def test_cache_custom_root_is_honoured(tmp_path: Path) -> None:
    cache = Cache(tmp_path / "custom")
    assert cache.root == tmp_path / "custom"


def test_cache_constructor_does_not_touch_filesystem(tmp_path: Path) -> None:
    """Per design note: the cache is lazy — instantiation must not create
    the root directory, since users may construct one just to query info()."""
    root = tmp_path / "lazy"
    Cache(root)
    assert not root.exists()


# ---------------------------------------------------------------------------
# RCSB
# ---------------------------------------------------------------------------


def test_rcsb_path_uses_lowercase_pdb_id(tmp_path: Path) -> None:
    cache = Cache(tmp_path)
    assert cache.rcsb_path("5J7L") == tmp_path / "rcsb" / "5j7l.json"
    assert cache.rcsb_path("5j7l") == tmp_path / "rcsb" / "5j7l.json"


def test_rcsb_get_returns_none_on_miss(tmp_path: Path) -> None:
    cache = Cache(tmp_path)
    assert cache.get_rcsb_payload("5J7L") is None


def test_rcsb_put_then_get_roundtrips(tmp_path: Path) -> None:
    cache = Cache(tmp_path)
    payload = {"rcsb_id": "5J7L", "exptl": [{"method": "X-RAY DIFFRACTION"}]}
    written = cache.put_rcsb_payload("5J7L", payload)
    assert written.exists()
    assert cache.get_rcsb_payload("5j7l") == payload


def test_rcsb_get_returns_none_on_corrupt_json(tmp_path: Path) -> None:
    """A half-written or hand-edited cache file should not crash the loader."""
    cache = Cache(tmp_path)
    path = cache.rcsb_path("5J7L")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("{not-json")
    assert cache.get_rcsb_payload("5J7L") is None


def test_rcsb_get_returns_none_on_non_object_json(tmp_path: Path) -> None:
    cache = Cache(tmp_path)
    path = cache.rcsb_path("5J7L")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("[1, 2, 3]")
    assert cache.get_rcsb_payload("5J7L") is None


# ---------------------------------------------------------------------------
# BGSU
# ---------------------------------------------------------------------------


def test_bgsu_path_is_sha256_of_key(tmp_path: Path) -> None:
    import hashlib

    cache = Cache(tmp_path)
    key = "units=5J7L|1|AA|G|926|scope=Rfam|resolution=20.0A"
    expected = hashlib.sha256(key.encode()).hexdigest()
    assert cache.bgsu_path(key).name == f"{expected}.json"


def test_bgsu_put_then_get_roundtrips(tmp_path: Path) -> None:
    cache = Cache(tmp_path)
    key = "some-opaque-key"
    payload = {"alignment": [{"reference_unit": "X", "mapped_units": []}]}
    cache.put_bgsu_payload(key, payload)
    assert cache.get_bgsu_payload(key) == payload


def test_bgsu_different_keys_different_files(tmp_path: Path) -> None:
    cache = Cache(tmp_path)
    cache.put_bgsu_payload("k1", {"a": 1})
    cache.put_bgsu_payload("k2", {"b": 2})
    assert cache.get_bgsu_payload("k1") == {"a": 1}
    assert cache.get_bgsu_payload("k2") == {"b": 2}


# ---------------------------------------------------------------------------
# Coordinates
# ---------------------------------------------------------------------------


def test_assembly_coords_path_format(tmp_path: Path) -> None:
    cache = Cache(tmp_path)
    assert cache.assembly_coords_path("5J7L", "1") == tmp_path / "coords" / "5j7l-assembly1.cif.gz"
    assert cache.assembly_coords_path("5FDV", "2") == tmp_path / "coords" / "5fdv-assembly2.cif.gz"


def test_assembly_coords_has_returns_false_on_miss(tmp_path: Path) -> None:
    cache = Cache(tmp_path)
    assert cache.has_assembly_coords("5J7L", "1") is False


def test_assembly_coords_put_then_has_and_read(tmp_path: Path) -> None:
    cache = Cache(tmp_path)
    written = cache.put_assembly_coords("5J7L", "1", b"\x1f\x8bfakegzip")
    assert written.exists()
    assert cache.has_assembly_coords("5J7L", "1")
    assert written.read_bytes() == b"\x1f\x8bfakegzip"


# ---------------------------------------------------------------------------
# Atomic write
# ---------------------------------------------------------------------------


def test_put_leaves_no_tmp_file_after_success(tmp_path: Path) -> None:
    cache = Cache(tmp_path)
    cache.put_rcsb_payload("5J7L", {"x": 1})
    tmp_files = list((tmp_path / "rcsb").glob("*.tmp"))
    assert tmp_files == []


def test_put_overwrites_existing_entry_atomically(tmp_path: Path) -> None:
    cache = Cache(tmp_path)
    cache.put_rcsb_payload("5J7L", {"version": 1})
    cache.put_rcsb_payload("5J7L", {"version": 2})
    assert cache.get_rcsb_payload("5J7L") == {"version": 2}
    assert list((tmp_path / "rcsb").glob("*.tmp")) == []


def test_put_is_deterministic_json_format(tmp_path: Path) -> None:
    """Sorted keys + compact separators — important for hashing/diff stability."""
    cache = Cache(tmp_path)
    path = cache.put_rcsb_payload("5J7L", {"b": 2, "a": 1})
    raw = path.read_text()
    assert raw == '{"a":1,"b":2}'
    assert json.loads(raw) == {"a": 1, "b": 2}


# ---------------------------------------------------------------------------
# info() and clear()
# ---------------------------------------------------------------------------


def test_info_on_missing_root_reports_not_exists(tmp_path: Path) -> None:
    cache = Cache(tmp_path / "missing")
    info = cache.info()
    assert info.exists is False
    assert info.total_entries == 0
    assert info.total_bytes == 0


def test_info_counts_entries_and_bytes(tmp_path: Path) -> None:
    cache = Cache(tmp_path)
    cache.put_rcsb_payload("5J7L", {"a": 1})
    cache.put_rcsb_payload("5TBW", {"b": 2})
    cache.put_bgsu_payload("k1", {"c": 3})
    cache.put_assembly_coords("5J7L", "1", b"some-coordinates-here")
    cache.put_pdbe_payload("5J7L", {"d": 4})
    info = cache.info()
    assert isinstance(info, CacheInfo)
    assert info.exists is True
    assert info.rcsb_entries == 2
    assert info.bgsu_entries == 1
    assert info.coords_entries == 1
    assert info.pdbe_entries == 1
    assert info.total_entries == 5
    assert info.total_bytes > 0


# ---------------------------------------------------------------------------
# PDBe namespace
# ---------------------------------------------------------------------------


def test_pdbe_path_uses_lowercase_pdb_id(tmp_path: Path) -> None:
    cache = Cache(tmp_path)
    assert cache.pdbe_path("5J7L") == tmp_path / "pdbe" / "5j7l.json"


def test_pdbe_get_returns_none_on_miss(tmp_path: Path) -> None:
    assert Cache(tmp_path).get_pdbe_payload("5J7L") is None


def test_pdbe_put_then_get_roundtrips(tmp_path: Path) -> None:
    cache = Cache(tmp_path)
    payload = {"5j7l": {"Rfam": {"RF00177": {"mappings": [{"chain_id": "AA"}]}}}}
    cache.put_pdbe_payload("5J7L", payload)
    assert cache.get_pdbe_payload("5j7l") == payload


def test_clear_removes_cache_root(tmp_path: Path) -> None:
    cache = Cache(tmp_path / "to-clear")
    cache.put_rcsb_payload("5J7L", {"a": 1})
    assert cache.root.exists()
    cache.clear()
    assert not cache.root.exists()


def test_clear_on_missing_root_does_not_raise(tmp_path: Path) -> None:
    cache = Cache(tmp_path / "never-existed")
    cache.clear()  # should not raise
    assert not cache.root.exists()


def test_cacheinfo_total_entries_computed_field() -> None:
    info = CacheInfo(root=Path("/tmp/x"), rcsb_entries=3, bgsu_entries=5, coords_entries=2)
    assert info.total_entries == 10


@pytest.mark.parametrize("pdb_id", ["5J7L", "5j7L", "5J7l", "5j7l"])
def test_rcsb_case_normalisation(tmp_path: Path, pdb_id: str) -> None:
    cache = Cache(tmp_path)
    cache.put_rcsb_payload(pdb_id, {"x": 1})
    # All four spellings should hit the same file.
    assert cache.get_rcsb_payload("5J7L") == {"x": 1}
    assert cache.get_rcsb_payload("5j7l") == {"x": 1}
