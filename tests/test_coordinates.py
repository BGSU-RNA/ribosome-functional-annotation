"""Unit tests for the coordinate retrieval and Gemmi loader (spec §5.3, §10.4)."""

from __future__ import annotations

from pathlib import Path

import gemmi
import httpx
import pytest
import respx

from ribosome_state_annotator import coordinates
from ribosome_state_annotator.cache import Cache
from ribosome_state_annotator.exceptions import (
    CoordinateDownloadError,
    CoordinateParsingError,
)

# ---------------------------------------------------------------------------
# URL builders
# ---------------------------------------------------------------------------


def test_assembly_download_url_lowercase_pdb() -> None:
    assert (
        coordinates.assembly_download_url("5J7L", "1")
        == "https://files.rcsb.org/download/5j7l-assembly1.cif.gz"
    )
    assert (
        coordinates.assembly_download_url("5FDV", "12")
        == "https://files.rcsb.org/download/5fdv-assembly12.cif.gz"
    )


def test_asymmetric_unit_download_url_lowercase_pdb() -> None:
    assert (
        coordinates.asymmetric_unit_download_url("5J7L")
        == "https://files.rcsb.org/download/5j7l.cif.gz"
    )


# ---------------------------------------------------------------------------
# download_assembly_cif (HTTP)
# ---------------------------------------------------------------------------


@respx.mock
def test_download_assembly_returns_bytes_on_200(minimal_cif_bytes: bytes) -> None:
    respx.get("https://files.rcsb.org/download/5j7l-assembly1.cif.gz").mock(
        return_value=httpx.Response(200, content=minimal_cif_bytes)
    )
    result = coordinates.download_assembly_cif("5J7L", "1")
    assert result == minimal_cif_bytes


@respx.mock
def test_download_falls_back_to_asu_on_404(minimal_cif_bytes: bytes) -> None:
    respx.get("https://files.rcsb.org/download/5j7l-assembly1.cif.gz").mock(
        return_value=httpx.Response(404)
    )
    respx.get("https://files.rcsb.org/download/5j7l.cif.gz").mock(
        return_value=httpx.Response(200, content=minimal_cif_bytes)
    )
    result = coordinates.download_assembly_cif("5J7L", "1")
    assert result == minimal_cif_bytes


@respx.mock
def test_download_no_fallback_propagates_404() -> None:
    respx.get("https://files.rcsb.org/download/5j7l-assembly1.cif.gz").mock(
        return_value=httpx.Response(404)
    )
    with pytest.raises(CoordinateDownloadError, match="HTTP 404"):
        coordinates.download_assembly_cif("5J7L", "1", allow_asu_fallback=False)


@respx.mock
def test_download_raises_on_non_200_non_404() -> None:
    respx.get("https://files.rcsb.org/download/5j7l-assembly1.cif.gz").mock(
        return_value=httpx.Response(500)
    )
    with pytest.raises(CoordinateDownloadError, match="HTTP 500"):
        coordinates.download_assembly_cif("5J7L", "1")


@respx.mock
def test_download_raises_on_asu_fallback_failure() -> None:
    respx.get("https://files.rcsb.org/download/5j7l-assembly1.cif.gz").mock(
        return_value=httpx.Response(404)
    )
    respx.get("https://files.rcsb.org/download/5j7l.cif.gz").mock(return_value=httpx.Response(404))
    with pytest.raises(CoordinateDownloadError, match="HTTP 404"):
        coordinates.download_assembly_cif("5J7L", "1")


@respx.mock
def test_download_raises_on_transport_error() -> None:
    respx.get("https://files.rcsb.org/download/5j7l-assembly1.cif.gz").mock(
        side_effect=httpx.ConnectError("connection refused")
    )
    with pytest.raises(CoordinateDownloadError, match="download failed"):
        coordinates.download_assembly_cif("5J7L", "1")


@respx.mock
def test_download_accepts_caller_supplied_client(minimal_cif_bytes: bytes) -> None:
    respx.get("https://files.rcsb.org/download/5j7l-assembly1.cif.gz").mock(
        return_value=httpx.Response(200, content=minimal_cif_bytes)
    )
    with httpx.Client(timeout=5.0, follow_redirects=True) as client:
        result = coordinates.download_assembly_cif("5J7L", "1", client=client)
    assert result == minimal_cif_bytes


# ---------------------------------------------------------------------------
# load_structure / load_structure_from_bytes
# ---------------------------------------------------------------------------


def test_load_structure_reads_uncompressed_cif(minimal_cif_path: Path) -> None:
    structure = coordinates.load_structure(minimal_cif_path)
    assert isinstance(structure, gemmi.Structure)
    chain_names = [c.name for c in structure[0]]
    assert chain_names == ["AA", "DA", "P1"]


def test_load_structure_reads_gzipped_cif(minimal_cif_gz_path: Path) -> None:
    """Spec §10.4: Gemmi auto-decompresses .cif.gz."""
    structure = coordinates.load_structure(minimal_cif_gz_path)
    assert [c.name for c in structure[0]] == ["AA", "DA", "P1"]


def test_load_structure_residues_have_expected_seqids(minimal_cif_path: Path) -> None:
    structure = coordinates.load_structure(minimal_cif_path)
    aa_chain = structure[0]["AA"]
    seqids = [r.seqid.num for r in aa_chain]
    assert seqids == [1, 2, 3, 4, 5]
    # The DA chain uses a 100-based numbering to mimic 23S rRNA mid-sequence
    # numbering and to ensure the loader preserves arbitrary seqid offsets.
    da_chain = structure[0]["DA"]
    assert [r.seqid.num for r in da_chain] == [100, 101, 102, 103, 104]


def test_load_structure_raises_on_missing_file(tmp_path: Path) -> None:
    with pytest.raises(CoordinateParsingError):
        coordinates.load_structure(tmp_path / "nonexistent.cif")


def test_load_structure_raises_on_garbage(tmp_path: Path) -> None:
    bad = tmp_path / "garbage.cif"
    bad.write_text("not a valid cif file at all\n")
    with pytest.raises(CoordinateParsingError):
        coordinates.load_structure(bad)


def test_load_structure_from_bytes_roundtrips(minimal_cif_bytes: bytes) -> None:
    structure = coordinates.load_structure_from_bytes(minimal_cif_bytes)
    assert [c.name for c in structure[0]] == ["AA", "DA", "P1"]


def test_load_structure_from_bytes_cleans_up_tempfile(minimal_cif_bytes: bytes) -> None:
    """Defensive: ensure the helper doesn't leak temp files between calls."""
    import tempfile

    before = set(Path(tempfile.gettempdir()).glob("tmp*.cif.gz"))
    coordinates.load_structure_from_bytes(minimal_cif_bytes)
    after = set(Path(tempfile.gettempdir()).glob("tmp*.cif.gz"))
    # No new lingering tmp files attributable to this call.
    assert not (after - before)


# ---------------------------------------------------------------------------
# get_assembly_structure orchestrator
# ---------------------------------------------------------------------------


def test_orchestrator_local_source_reads_local_path(minimal_cif_path: Path) -> None:
    structure = coordinates.get_assembly_structure(
        "5J7L", "1", source="local", local_path=minimal_cif_path
    )
    assert [c.name for c in structure[0]] == ["AA", "DA", "P1"]


def test_orchestrator_local_source_requires_local_path() -> None:
    with pytest.raises(ValueError, match="local_path"):
        coordinates.get_assembly_structure("5J7L", "1", source="local", local_path=None)


def test_orchestrator_local_source_missing_path_raises(tmp_path: Path) -> None:
    with pytest.raises(CoordinateDownloadError, match="local_path does not exist"):
        coordinates.get_assembly_structure(
            "5J7L", "1", source="local", local_path=tmp_path / "missing.cif"
        )


def test_orchestrator_auto_with_cache_hit_reads_cache(
    tmp_path: Path, minimal_cif_bytes: bytes
) -> None:
    """A cached coords file is read directly without any HTTP call."""
    cache = Cache(tmp_path)
    cache.put_assembly_coords("5J7L", "1", minimal_cif_bytes)
    # respx is NOT active; if download were attempted, this would raise.
    structure = coordinates.get_assembly_structure("5J7L", "1", cache=cache)
    assert [c.name for c in structure[0]] == ["AA", "DA", "P1"]


@respx.mock
def test_orchestrator_auto_with_cache_miss_downloads_and_caches(
    tmp_path: Path, minimal_cif_bytes: bytes
) -> None:
    cache = Cache(tmp_path)
    respx.get("https://files.rcsb.org/download/5j7l-assembly1.cif.gz").mock(
        return_value=httpx.Response(200, content=minimal_cif_bytes)
    )
    structure = coordinates.get_assembly_structure("5J7L", "1", cache=cache)
    assert [c.name for c in structure[0]] == ["AA", "DA", "P1"]
    # After the call, the bytes are in the cache for next time.
    assert cache.has_assembly_coords("5J7L", "1")


@respx.mock
def test_orchestrator_auto_without_cache_returns_structure_via_tempfile(
    minimal_cif_bytes: bytes,
) -> None:
    respx.get("https://files.rcsb.org/download/5j7l-assembly1.cif.gz").mock(
        return_value=httpx.Response(200, content=minimal_cif_bytes)
    )
    structure = coordinates.get_assembly_structure("5J7L", "1")
    assert [c.name for c in structure[0]] == ["AA", "DA", "P1"]
