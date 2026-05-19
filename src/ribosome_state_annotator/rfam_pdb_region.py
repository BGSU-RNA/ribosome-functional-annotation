"""Rfam ``pdb_full_region`` integration for rRNA Rfam annotations.

EBI publishes a single flat file mapping every PDB chain to its
matching Rfam families with cmsearch alignment scores:

    https://ftp.ebi.ac.uk/pub/databases/Rfam/.preview/pdb_full_region.txt.gz

The file is ~200 KB compressed, contains ~28 k rows across 4.4 k
distinct PDBs, and is updated by EBI on a regular cadence.

For each ``(pdb_id, author_chain_id)`` we pick the **single highest
bit-score** Rfam accession. The reason for single-best selection: rRNA
HMMs from different superkingdoms (bacterial 16S RF00177, archaeal SSU
RF01959, eukaryotic 18S RF01960) hit the same chain in cross-family
scans because they share structural ancestry, but only one is the
biological truth and bit_score reliably picks it. This eliminates the
"MIXED rrna_core" edge case at its source.

Design (mirrors :mod:`.raddb`):

- Cached under ``~/.cache/ribosome-state-annotator/rfam/`` alongside a
  small ``metadata.json`` describing the source URL, the upstream
  ``Last-Modified`` value, and a UTC ``downloaded_at`` timestamp.
- A 7-day staleness window triggers an online check. The check is a
  cheap HEAD request comparing the upstream ``Last-Modified`` header
  with the cached value.
- Every failure mode (no cache, stale cache + offline, malformed file,
  missing columns) returns ``None`` for the lookup so the annotation
  pipeline continues. Rfam integration **must not** crash the pipeline.
"""

from __future__ import annotations

import gzip
import json
import logging
from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any

import httpx

from ribosome_state_annotator.cache import DEFAULT_CACHE_ROOT

logger = logging.getLogger(__name__)

RFAM_PDB_REGION_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/Rfam/.preview/pdb_full_region.txt.gz"
)
"""EBI Rfam preview directory's per-PDB-chain Rfam region file."""

RFAM_PDB_REGION_FILENAME = "pdb_full_region.txt.gz"
RFAM_PDB_REGION_METADATA_FILENAME = "pdb_full_region.metadata.json"

STALENESS_DAYS = 7
"""Local file is considered stale and triggers an online check after this many days."""

DEFAULT_HTTP_TIMEOUT = 60.0


# ---------------------------------------------------------------------------
# Public dataclasses
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class RfamPdbRegionMetadata:
    """Lightweight record describing the cached Rfam pdb_full_region file."""

    source_url: str
    downloaded_at: datetime
    last_modified: str | None = None
    local_file: str = RFAM_PDB_REGION_FILENAME

    def to_dict(self) -> dict[str, str | None]:
        return {
            "source_url": self.source_url,
            "downloaded_at": self.downloaded_at.replace(microsecond=0).isoformat(),
            "last_modified": self.last_modified,
            "local_file": self.local_file,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> RfamPdbRegionMetadata:
        return cls(
            source_url=str(data["source_url"]),
            downloaded_at=_parse_iso(str(data["downloaded_at"])),
            last_modified=data.get("last_modified"),
            local_file=str(data.get("local_file", RFAM_PDB_REGION_FILENAME)),
        )


@dataclass(frozen=True)
class RfamPdbRegionDataset:
    """In-memory Rfam pdb_full_region dataset: per-chain best-score lookup."""

    metadata: RfamPdbRegionMetadata
    lookup: dict[tuple[str, str], str] = field(default_factory=dict)
    """Keyed by ``(pdb_id_lower, auth_chain_id)``. Value is the single
    Rfam accession with the highest bit-score for that chain."""


# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------


def get_rfam_cache_dir(cache_root: Path | None = None) -> Path:
    """Resolve the ``rfam/`` sub-directory inside the package cache."""
    root = cache_root if cache_root is not None else DEFAULT_CACHE_ROOT
    return root / "rfam"


def get_local_rfam_file_path(cache_root: Path | None = None) -> Path:
    return get_rfam_cache_dir(cache_root) / RFAM_PDB_REGION_FILENAME


def get_local_rfam_metadata_path(cache_root: Path | None = None) -> Path:
    return get_rfam_cache_dir(cache_root) / RFAM_PDB_REGION_METADATA_FILENAME


# ---------------------------------------------------------------------------
# Metadata I/O
# ---------------------------------------------------------------------------


def load_rfam_metadata(cache_root: Path | None = None) -> RfamPdbRegionMetadata | None:
    """Read cached Rfam metadata, returning ``None`` if absent or malformed."""
    path = get_local_rfam_metadata_path(cache_root)
    if not path.is_file():
        return None
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        logger.warning("Rfam metadata read failed (%s); treating as missing", exc)
        return None
    if not isinstance(data, dict):
        logger.warning("Rfam metadata is not a JSON object; treating as missing")
        return None
    try:
        return RfamPdbRegionMetadata.from_dict(data)
    except (KeyError, ValueError) as exc:
        logger.warning("Rfam metadata schema mismatch (%s); treating as missing", exc)
        return None


def save_rfam_metadata(
    metadata: RfamPdbRegionMetadata, cache_root: Path | None = None
) -> Path:
    path = get_local_rfam_metadata_path(cache_root)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(metadata.to_dict(), indent=2) + "\n", encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------


def probe_remote_last_modified(
    *, client: httpx.Client | None = None
) -> str | None:
    """HEAD-probe the EBI URL and return the ``Last-Modified`` value.

    Returns ``None`` on any network failure. The value is a raw HTTP-date
    string suitable for direct string comparison against the cached
    metadata's ``last_modified`` field.
    """
    http: httpx.Client = client if client is not None else httpx.Client(timeout=DEFAULT_HTTP_TIMEOUT)
    try:
        response = http.head(RFAM_PDB_REGION_URL, follow_redirects=True)
    except httpx.HTTPError as exc:
        logger.debug("Rfam HEAD %s failed: %s", RFAM_PDB_REGION_URL, exc)
        return None
    finally:
        if client is None:
            http.close()
    if response.status_code != 200:
        logger.warning("Rfam HEAD returned HTTP %d", response.status_code)
        return None
    last_modified = response.headers.get("Last-Modified")
    return last_modified if last_modified is None else str(last_modified)


def download_rfam_pdb_region(
    *,
    cache_root: Path | None = None,
    client: httpx.Client | None = None,
) -> tuple[Path, RfamPdbRegionMetadata]:
    """Download the Rfam ``pdb_full_region.txt.gz`` file and persist alongside metadata.

    Raises :class:`httpx.HTTPError` on network failure or non-200
    status. Callers are responsible for catching and degrading
    gracefully.
    """
    http: httpx.Client = client if client is not None else httpx.Client(timeout=DEFAULT_HTTP_TIMEOUT)
    try:
        response = http.get(RFAM_PDB_REGION_URL, follow_redirects=True)
        response.raise_for_status()
        body = response.content
        last_modified = response.headers.get("Last-Modified")
    finally:
        if client is None:
            http.close()

    file_path = get_local_rfam_file_path(cache_root)
    file_path.parent.mkdir(parents=True, exist_ok=True)
    file_path.write_bytes(body)
    metadata = RfamPdbRegionMetadata(
        source_url=RFAM_PDB_REGION_URL,
        downloaded_at=datetime.now(timezone.utc),
        last_modified=last_modified,
    )
    save_rfam_metadata(metadata, cache_root=cache_root)
    logger.info(
        "downloaded Rfam pdb_full_region (%d bytes; Last-Modified=%s) to %s",
        len(body),
        last_modified,
        file_path,
    )
    return file_path, metadata


# ---------------------------------------------------------------------------
# Refresh orchestration
# ---------------------------------------------------------------------------


def ensure_rfam_pdb_region_available(
    *,
    cache_root: Path | None = None,
    force_refresh: bool = False,
    client: httpx.Client | None = None,
    now: datetime | None = None,
) -> RfamPdbRegionMetadata | None:
    """Make sure a local Rfam pdb_full_region file exists, refreshing if stale.

    Returns the metadata describing the (possibly updated) local file,
    or ``None`` if no local file exists and no download succeeded.

    Behaviour mirrors :func:`.raddb.ensure_raddb_available`:

    1. Local file missing → attempt download.
    2. Local file present + ``force_refresh`` → HEAD-check for newer
       ``Last-Modified``; download if changed.
    3. Local file present + age ≥ ``STALENESS_DAYS`` → HEAD-check for
       newer; download if changed.
    4. Local file present + fresh → return current metadata, no network.
    5. Any network failure with a local file present → keep the local
       file and log a warning.
    6. Any network failure without a local file → return ``None``.
    """
    file_path = get_local_rfam_file_path(cache_root)
    metadata = load_rfam_metadata(cache_root)
    now_dt = now or datetime.now(timezone.utc)

    if not file_path.is_file() or metadata is None:
        logger.info("Rfam pdb_full_region file missing; downloading")
        return _try_refresh(cache_root=cache_root, client=client, current=None)

    if force_refresh:
        logger.info("Rfam refresh forced; checking for newer version")
        refreshed = _try_refresh(cache_root=cache_root, client=client, current=metadata)
        return refreshed if refreshed is not None else metadata

    age = now_dt - _aware(metadata.downloaded_at)
    if age >= timedelta(days=STALENESS_DAYS):
        logger.info(
            "Rfam cache is %d days old; checking for newer version", age.days
        )
        refreshed = _try_refresh(cache_root=cache_root, client=client, current=metadata)
        return refreshed if refreshed is not None else metadata

    logger.info(
        "using cached Rfam pdb_full_region (downloaded %s, Last-Modified=%s)",
        metadata.downloaded_at.date(),
        metadata.last_modified,
    )
    return metadata


def _try_refresh(
    *,
    cache_root: Path | None,
    client: httpx.Client | None,
    current: RfamPdbRegionMetadata | None,
) -> RfamPdbRegionMetadata | None:
    """Check upstream ``Last-Modified`` and download if newer than ``current``.

    Returns the new metadata on successful download; returns ``current``
    unchanged if the upstream copy is the same; returns ``None`` if the
    probe or download fails.
    """
    remote_last_modified = probe_remote_last_modified(client=client)
    if remote_last_modified is None:
        if current is not None:
            logger.warning(
                "Rfam Last-Modified probe failed; keeping cached %s",
                current.last_modified,
            )
        else:
            logger.warning("Rfam unavailable: probe failed and no cached file")
        # Still attempt the download — some servers (or proxies) strip
        # the Last-Modified header. A successful GET is enough to refresh.
        if current is None:
            try:
                _, metadata = download_rfam_pdb_region(
                    cache_root=cache_root, client=client
                )
                return metadata
            except httpx.HTTPError as exc:
                logger.warning("Rfam download failed: %s", exc)
                return None
        return None

    if current is not None and remote_last_modified == current.last_modified:
        logger.info(
            "Rfam cache is already at the latest version (Last-Modified=%s)",
            current.last_modified,
        )
        # Touch the metadata so the next staleness check resets the clock.
        refreshed = RfamPdbRegionMetadata(
            source_url=current.source_url,
            downloaded_at=datetime.now(timezone.utc),
            last_modified=current.last_modified,
            local_file=current.local_file,
        )
        save_rfam_metadata(refreshed, cache_root=cache_root)
        return refreshed

    try:
        _, metadata = download_rfam_pdb_region(cache_root=cache_root, client=client)
    except httpx.HTTPError as exc:
        logger.warning("Rfam download failed: %s", exc)
        return None
    return metadata


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------


def load_rfam_pdb_region_table(
    *,
    cache_root: Path | None = None,
    metadata: RfamPdbRegionMetadata | None = None,
) -> list[tuple[str, str, str, float]]:
    """Read the cached Rfam pdb_full_region file and return raw rows.

    Returns a list of ``(rfam_acc, pdb_id_lower, auth_chain_id,
    bit_score)`` tuples. Rows with malformed bit_score are silently
    skipped. The function does not pick a best-score per chain — that's
    done by :func:`build_rfam_pdb_region_lookup`.
    """
    metadata = metadata or load_rfam_metadata(cache_root)
    if metadata is None:
        return []
    file_path = get_rfam_cache_dir(cache_root) / metadata.local_file
    if not file_path.is_file():
        logger.warning("Rfam file missing at %s", file_path)
        return []

    rows: list[tuple[str, str, str, float]] = []
    with gzip.open(file_path, mode="rt", encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            rfam_acc = parts[0]
            pdb_id = parts[1].lower()
            chain_id = parts[2]
            try:
                bit_score = float(parts[5])
            except ValueError:
                continue
            rows.append((rfam_acc, pdb_id, chain_id, bit_score))
    return rows


def build_rfam_pdb_region_lookup(
    rows: list[tuple[str, str, str, float]],
) -> dict[tuple[str, str], str]:
    """Pick the single best-score Rfam per ``(pdb_id_lower, chain_id)``.

    Ties are broken by the input row order (rare in practice).
    """
    best_score: dict[tuple[str, str], float] = {}
    best_acc: dict[tuple[str, str], str] = {}
    for rfam_acc, pdb_id, chain_id, bit_score in rows:
        key = (pdb_id, chain_id)
        if key not in best_score or bit_score > best_score[key]:
            best_score[key] = bit_score
            best_acc[key] = rfam_acc
    return best_acc


def load_rfam_pdb_region_dataset(
    *,
    cache_root: Path | None = None,
    metadata: RfamPdbRegionMetadata | None = None,
) -> RfamPdbRegionDataset | None:
    """Parse the cached Rfam file into a :class:`RfamPdbRegionDataset`.

    Returns ``None`` if no metadata or file is present.
    """
    metadata = metadata or load_rfam_metadata(cache_root)
    if metadata is None:
        return None
    rows = load_rfam_pdb_region_table(cache_root=cache_root, metadata=metadata)
    lookup = build_rfam_pdb_region_lookup(rows)
    return RfamPdbRegionDataset(metadata=metadata, lookup=lookup)


# ---------------------------------------------------------------------------
# Lookup
# ---------------------------------------------------------------------------


def get_rfam_for_chain(
    dataset: RfamPdbRegionDataset | None,
    pdb_id: str,
    chain_id: str,
) -> str | None:
    """Return the best-score Rfam accession for ``(pdb_id, chain_id)``, or ``None``.

    PDB IDs are normalised to lowercase for the lookup; chain IDs are
    compared case-sensitively (matches the EBI file's casing).
    """
    if dataset is None:
        return None
    return dataset.lookup.get((pdb_id.lower(), chain_id))


def get_rfam_mapping_for_pdb(
    dataset: RfamPdbRegionDataset | None,
    pdb_id: str,
) -> dict[str, list[str]]:
    """Return ``{chain_id: [rfam_acc]}`` (single best-score) for a PDB.

    Output shape matches :func:`pdbe_client.parse_rfam_mappings` so the
    api.py merge step (:func:`api._augment_chain_rfam`) is drop-in
    compatible.
    """
    if dataset is None:
        return {}
    pdb_lower = pdb_id.lower()
    out: dict[str, list[str]] = {}
    for (pdb, chain), rfam_acc in dataset.lookup.items():
        if pdb == pdb_lower:
            out[chain] = [rfam_acc]
    return out


# ---------------------------------------------------------------------------
# CLI helpers (`ribostate rfam info` etc.)
# ---------------------------------------------------------------------------


def list_rfam_files(cache_root: Path | None = None) -> int:
    """Return the size in bytes of the cached Rfam file, or 0 if absent."""
    file_path = get_local_rfam_file_path(cache_root)
    if not file_path.is_file():
        return 0
    return file_path.stat().st_size


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _parse_iso(value: str) -> datetime:
    """Parse an ISO-8601 timestamp; ensure UTC tzinfo on naive inputs."""
    dt = datetime.fromisoformat(value)
    if dt.tzinfo is None:
        return dt.replace(tzinfo=timezone.utc)
    return dt


def _aware(dt: datetime) -> datetime:
    """Coerce a naive datetime to UTC; pass through aware datetimes."""
    if dt.tzinfo is None:
        return dt.replace(tzinfo=timezone.utc)
    return dt
