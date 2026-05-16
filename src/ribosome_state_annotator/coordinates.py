"""Biological-assembly coordinate retrieval and Gemmi loading (spec §5.3, §10.4).

This module covers three concerns:

- URL templates for the RCSB ``files.rcsb.org`` download endpoint, plus
  the ASU fallback (§5.1.3).
- :func:`download_assembly_cif` — sync httpx GET that retrieves the
  gzipped mmCIF bytes; falls back to the asymmetric unit when the
  biological-assembly file 404s.
- :func:`load_structure` / :func:`load_structure_from_bytes` — thin Gemmi
  wrappers. Gemmi reads ``.cif`` and ``.cif.gz`` transparently via
  ``gemmi.read_structure(path)`` (spec §10.4).
- :func:`get_assembly_structure` — high-level convenience that combines
  cache lookup, download, and Gemmi parsing in one call.

Per spec §10.4, only the first model of the structure is meaningful —
biological-assembly mmCIFs since the May 2022 wwPDB format change use
chain duplication, not model duplication, to represent operator-expanded
copies. Step 8's neighbour-search layer enforces "read model 0 only";
this module does not need to.
"""

from __future__ import annotations

import logging
import tempfile
from pathlib import Path
from typing import Literal

import gemmi
import httpx

from ribosome_state_annotator.cache import Cache
from ribosome_state_annotator.exceptions import (
    CoordinateDownloadError,
    CoordinateParsingError,
)

logger = logging.getLogger(__name__)

RCSB_ASSEMBLY_DOWNLOAD_TEMPLATE = (
    "https://files.rcsb.org/download/{pdb_id}-assembly{assembly_id}.cif.gz"
)
"""§5.1.3: biological-assembly mmCIF URL template (PDB ID lowercased)."""

RCSB_ASU_DOWNLOAD_TEMPLATE = "https://files.rcsb.org/download/{pdb_id}.cif.gz"
"""§5.1.3: asymmetric-unit fallback URL when the biological-assembly file 404s."""


# ---------------------------------------------------------------------------
# URL builders
# ---------------------------------------------------------------------------


def assembly_download_url(pdb_id: str, assembly_id: str) -> str:
    """Return the biological-assembly mmCIF URL for ``(pdb_id, assembly_id)``."""
    return RCSB_ASSEMBLY_DOWNLOAD_TEMPLATE.format(pdb_id=pdb_id.lower(), assembly_id=assembly_id)


def asymmetric_unit_download_url(pdb_id: str) -> str:
    """Return the asymmetric-unit mmCIF URL for ``pdb_id`` (fallback target)."""
    return RCSB_ASU_DOWNLOAD_TEMPLATE.format(pdb_id=pdb_id.lower())


# ---------------------------------------------------------------------------
# Downloader
# ---------------------------------------------------------------------------


def download_assembly_cif(
    pdb_id: str,
    assembly_id: str,
    *,
    client: httpx.Client | None = None,
    allow_asu_fallback: bool = True,
    timeout: float = 120.0,
) -> bytes:
    """Download the biological-assembly mmCIF (gzipped) from RCSB.

    On HTTP 404 for the assembly file and ``allow_asu_fallback=True``, the
    function retries against the asymmetric unit URL and returns those
    bytes instead. Any other non-200 response (or the fallback also 404ing)
    raises :class:`CoordinateDownloadError`.

    Args:
        pdb_id: PDB accession (case-insensitive; lowercased for the URL).
        assembly_id: Biological assembly identifier (typically ``"1"``).
        client: Optional pre-built :class:`httpx.Client` with
            ``follow_redirects=True``. When omitted, an ephemeral client is
            created and closed inside this function.
        allow_asu_fallback: If ``True`` (default), a 404 on the assembly
            URL triggers a fallback to the asymmetric-unit URL per §5.1.3.
        timeout: Request timeout in seconds. Ignored when ``client`` is supplied.

    Returns:
        The raw gzipped mmCIF bytes — pass to
        :func:`load_structure_from_bytes` or write to disk via
        :meth:`Cache.put_assembly_coords`.
    """
    primary_url = assembly_download_url(pdb_id, assembly_id)
    own_client = client is None
    http: httpx.Client = (
        client if client is not None else httpx.Client(timeout=timeout, follow_redirects=True)
    )
    try:
        response = _safe_get(http, primary_url)
        if response.status_code == 404 and allow_asu_fallback:
            fallback_url = asymmetric_unit_download_url(pdb_id)
            logger.warning(
                "assembly %s for %s not available; falling back to ASU at %s",
                assembly_id,
                pdb_id,
                fallback_url,
            )
            response = _safe_get(http, fallback_url)
        if response.status_code != 200:
            raise CoordinateDownloadError(
                f"download returned HTTP {response.status_code} for {response.url}"
            )
        return response.content
    finally:
        if own_client:
            http.close()


def _safe_get(http: httpx.Client, url: str) -> httpx.Response:
    try:
        return http.get(url)
    except httpx.HTTPError as exc:
        raise CoordinateDownloadError(f"download failed for {url}: {exc}") from exc


# ---------------------------------------------------------------------------
# Gemmi loaders
# ---------------------------------------------------------------------------


def load_structure(path: Path) -> gemmi.Structure:
    """Parse an mmCIF (optionally gzipped) into a :class:`gemmi.Structure`.

    Gemmi auto-decompresses ``.cif.gz`` based on the filename extension
    (spec §10.4). Wraps Gemmi's various parse errors as
    :class:`CoordinateParsingError` so callers don't import gemmi just to
    catch them.
    """
    try:
        return gemmi.read_structure(str(path))
    except (OSError, RuntimeError, ValueError) as exc:
        raise CoordinateParsingError(f"failed to parse {path}: {exc}") from exc


def load_structure_from_bytes(
    data: bytes,
    *,
    suffix: str = ".cif.gz",
) -> gemmi.Structure:
    """Parse mmCIF bytes by writing them to a tempfile and calling Gemmi.

    Gemmi's reader takes a path, not a stream, so this is the only way to
    parse in-memory bytes. The tempfile is deleted before returning even
    on failure — the parsed :class:`gemmi.Structure` doesn't reference
    the file after parsing.
    """
    with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as fh:
        fh.write(data)
        temp_path = Path(fh.name)
    try:
        return load_structure(temp_path)
    finally:
        temp_path.unlink(missing_ok=True)


# ---------------------------------------------------------------------------
# High-level orchestrator
# ---------------------------------------------------------------------------

CoordinateSource = Literal["auto", "local"]


def get_assembly_structure(
    pdb_id: str,
    assembly_id: str,
    *,
    source: CoordinateSource = "auto",
    local_path: Path | None = None,
    cache: Cache | None = None,
    client: httpx.Client | None = None,
) -> gemmi.Structure:
    """Resolve and parse the biological-assembly coordinates for one assembly.

    The ``source`` parameter chooses the resolution strategy:

    - ``"auto"`` (default): check ``cache`` for a hit, otherwise download
      via :func:`download_assembly_cif`. If ``cache`` is provided, the
      downloaded bytes are written into the cache before parsing. If
      ``cache`` is ``None``, the bytes go through a tempfile inside
      :func:`load_structure_from_bytes`.
    - ``"local"``: parse ``local_path`` directly. ``local_path`` must
      exist; otherwise :class:`CoordinateDownloadError` is raised so the
      caller sees the same exception family as a missing remote file.
    """
    if source == "local":
        if local_path is None:
            raise ValueError("source='local' requires local_path")
        if not local_path.exists():
            raise CoordinateDownloadError(f"local_path does not exist: {local_path}")
        return load_structure(local_path)

    # source == "auto"
    if cache is not None and cache.has_assembly_coords(pdb_id, assembly_id):
        cached_path = cache.assembly_coords_path(pdb_id, assembly_id)
        logger.debug("coords cache hit for %s assembly %s at %s", pdb_id, assembly_id, cached_path)
        return load_structure(cached_path)

    logger.info("downloading coordinates for %s assembly %s", pdb_id, assembly_id)
    data = download_assembly_cif(pdb_id, assembly_id, client=client)

    if cache is not None:
        cached_path = cache.put_assembly_coords(pdb_id, assembly_id, data)
        return load_structure(cached_path)

    return load_structure_from_bytes(data)


__all__ = [
    "RCSB_ASSEMBLY_DOWNLOAD_TEMPLATE",
    "RCSB_ASU_DOWNLOAD_TEMPLATE",
    "CoordinateSource",
    "assembly_download_url",
    "asymmetric_unit_download_url",
    "download_assembly_cif",
    "get_assembly_structure",
    "load_structure",
    "load_structure_from_bytes",
]
