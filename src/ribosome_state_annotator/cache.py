"""Filesystem cache for RCSB / BGSU / PDBe / coordinate files (spec §17).

The cache is content-addressed and never expires automatically — see spec
§17.2 for the v1 invalidation policy. Four top-level namespaces:

- ``rcsb/<pdb_id_lower>.json``: RCSB GraphQL ``entry`` payloads.
- ``bgsu/<sha256(opaque-key)>.json``: BGSU correspondence responses.
  The opaque key is chosen by the caller; :mod:`bgsu_client` currently
  uses the units list joined with ``,`` plus the scope and resolution.
- ``pdbe/<pdb_id_lower>.json``: PDBe ``/nucleic_mappings/rfam`` responses
  (added in v1 because the spec's RCSB Rfam field is no longer
  populated for rRNA chains; see :mod:`pdbe_client`).
- ``coords/<pdb_id_lower>-assembly<n>.cif.gz``: biological-assembly mmCIFs.

Writes are atomic: each payload is written to a sibling ``.tmp`` file and
renamed into place, so a crash mid-write leaves the previous cached
entry (if any) intact rather than producing a half-written file.
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import shutil
from pathlib import Path
from typing import Any

from pydantic import BaseModel

logger = logging.getLogger(__name__)

DEFAULT_CACHE_ROOT: Path = Path("~/.cache/ribosome-state-annotator").expanduser()
"""Default cache directory per spec §17. Resolved from ``$HOME`` once at import."""

_RCSB_DIR = "rcsb"
_BGSU_DIR = "bgsu"
_COORDS_DIR = "coords"
_PDBE_DIR = "pdbe"
_FR3D_DIR = "fr3d"


class CacheInfo(BaseModel):
    """Summary returned by :meth:`Cache.info` and by ``ribostate cache info``."""

    root: Path
    rcsb_entries: int = 0
    bgsu_entries: int = 0
    coords_entries: int = 0
    pdbe_entries: int = 0
    fr3d_entries: int = 0
    total_bytes: int = 0
    exists: bool = False

    @property
    def total_entries(self) -> int:
        return (
            self.rcsb_entries
            + self.bgsu_entries
            + self.coords_entries
            + self.pdbe_entries
            + self.fr3d_entries
        )


class Cache:
    """On-disk cache for RCSB / BGSU / coordinate files.

    The cache is created lazily — the constructor does not touch the
    filesystem. The first ``put`` call materialises the directory tree.

    Designed to be safe under concurrent reads + interleaved writes from
    independent processes: every put goes through ``write → rename``,
    which is atomic on POSIX and on Python 3.3+ Windows.
    """

    def __init__(self, root: Path | None = None) -> None:
        self.root: Path = Path(root) if root is not None else DEFAULT_CACHE_ROOT

    # ------------------------------------------------------------------ RCSB

    def rcsb_path(self, pdb_id: str) -> Path:
        return self.root / _RCSB_DIR / f"{pdb_id.lower()}.json"

    def get_rcsb_payload(self, pdb_id: str) -> dict[str, Any] | None:
        return _read_json(self.rcsb_path(pdb_id))

    def put_rcsb_payload(self, pdb_id: str, payload: dict[str, Any]) -> Path:
        path = self.rcsb_path(pdb_id)
        _atomic_write_json(path, payload)
        return path

    # ------------------------------------------------------------------ BGSU

    def bgsu_path(self, key: str) -> Path:
        digest = hashlib.sha256(key.encode("utf-8")).hexdigest()
        return self.root / _BGSU_DIR / f"{digest}.json"

    def get_bgsu_payload(self, key: str) -> dict[str, Any] | None:
        return _read_json(self.bgsu_path(key))

    def put_bgsu_payload(self, key: str, payload: dict[str, Any]) -> Path:
        path = self.bgsu_path(key)
        _atomic_write_json(path, payload)
        return path

    # ----------------------------------------------------------- Coordinates

    def assembly_coords_path(self, pdb_id: str, assembly_id: str) -> Path:
        return self.root / _COORDS_DIR / f"{pdb_id.lower()}-assembly{assembly_id}.cif.gz"

    def has_assembly_coords(self, pdb_id: str, assembly_id: str) -> bool:
        return self.assembly_coords_path(pdb_id, assembly_id).is_file()

    def put_assembly_coords(self, pdb_id: str, assembly_id: str, data: bytes) -> Path:
        path = self.assembly_coords_path(pdb_id, assembly_id)
        _atomic_write_bytes(path, data)
        return path

    # ------------------------------------------------------------------ PDBe

    def pdbe_path(self, pdb_id: str) -> Path:
        return self.root / _PDBE_DIR / f"{pdb_id.lower()}.json"

    def get_pdbe_payload(self, pdb_id: str) -> dict[str, Any] | None:
        return _read_json(self.pdbe_path(pdb_id))

    def put_pdbe_payload(self, pdb_id: str, payload: dict[str, Any]) -> Path:
        path = self.pdbe_path(pdb_id)
        _atomic_write_json(path, payload)
        return path

    # ------------------------------------------------------------------ FR3D

    def fr3d_path(self, pdb_id: str) -> Path:
        """Per-PDB FR3D base-pair CSV path."""
        return self.root / _FR3D_DIR / f"{pdb_id.lower()}.csv"

    def get_fr3d_csv(self, pdb_id: str) -> bytes | None:
        """Return the cached FR3D CSV bytes, or ``None`` if absent."""
        path = self.fr3d_path(pdb_id)
        if not path.is_file():
            return None
        try:
            return path.read_bytes()
        except OSError as exc:
            logger.warning("fr3d cache read failed for %s (%s); treating as miss", path, exc)
            return None

    def put_fr3d_csv(self, pdb_id: str, data: bytes) -> Path:
        path = self.fr3d_path(pdb_id)
        _atomic_write_bytes(path, data)
        return path

    # ---------------------------------------------------------- Maintenance

    def info(self) -> CacheInfo:
        """Walk the cache tree and report counts and total disk usage."""
        if not self.root.exists():
            return CacheInfo(root=self.root, exists=False)
        return CacheInfo(
            root=self.root,
            exists=True,
            rcsb_entries=_count_files(self.root / _RCSB_DIR),
            bgsu_entries=_count_files(self.root / _BGSU_DIR),
            coords_entries=_count_files(self.root / _COORDS_DIR),
            pdbe_entries=_count_files(self.root / _PDBE_DIR),
            fr3d_entries=_count_files(self.root / _FR3D_DIR),
            total_bytes=_dir_size(self.root),
        )

    def clear(self) -> None:
        """Delete the cache root and everything below it. No-op if absent."""
        if self.root.exists():
            shutil.rmtree(self.root)


# ---------------------------------------------------------------------------
# Module-level helpers
# ---------------------------------------------------------------------------


def _atomic_write_bytes(path: Path, data: bytes) -> None:
    """Write ``data`` to ``path`` via ``write → rename`` so partial writes
    can't leave a corrupt cache entry."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    try:
        tmp.write_bytes(data)
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            tmp.unlink(missing_ok=True)


def _atomic_write_json(path: Path, payload: Any) -> None:
    data = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    _atomic_write_bytes(path, data)


def _read_json(path: Path) -> dict[str, Any] | None:
    if not path.is_file():
        return None
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        logger.warning("cache read failed for %s (%s); treating as miss", path, exc)
        return None
    if not isinstance(raw, dict):
        logger.warning("cache entry %s is not a JSON object; treating as miss", path)
        return None
    return raw


def _count_files(directory: Path) -> int:
    if not directory.is_dir():
        return 0
    return sum(1 for entry in directory.iterdir() if entry.is_file())


def _dir_size(directory: Path) -> int:
    total = 0
    for current_root, _dirs, files in os.walk(directory):
        for name in files:
            try:
                total += (Path(current_root) / name).stat().st_size
            except OSError:
                continue
    return total


# Public re-export for downstream typing.
__all__ = ["DEFAULT_CACHE_ROOT", "Cache", "CacheInfo"]
