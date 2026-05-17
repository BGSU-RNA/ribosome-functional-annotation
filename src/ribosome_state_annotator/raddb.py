"""RADdb integration for large-scale ribosome movement metrics.

RADdb (Ribosome Angle Decomposition database;
https://radtool.rc.northeastern.edu/) publishes a weekly CSV that
decomposes the inter-subunit and SSU-head motions of every public
ribosome structure into the canonical body / head rotation, tilt, and
translation parameters introduced by Mears et al.

For v1 the package exposes only the two most widely-used metrics —
``intersubunit_rotation`` (``body rot.``) and ``ssu_head_rotation``
(``head rot.``) — but the full row stays in memory so future versions
can surface tilt / translation without changing the cache or download
path.

Design notes:

- The CSV (~500 KB, ~2k rows) is cached under
  ``~/.cache/ribosome-state-annotator/raddb/`` alongside a small
  ``metadata.json`` describing the source URL, the embedded date
  (``rad_date``), and a UTC ``downloaded_at`` timestamp.
- A 7-day staleness window triggers an online check for a newer
  release. The check walks back at most 60 calendar days from today
  via cheap HEAD requests until it lands on the most recent file.
- Every failure mode (no cache, stale cache + offline, malformed CSV,
  duplicate keys, missing columns) returns ``None`` for the metrics so
  the annotation pipeline can continue. RADdb integration **must not**
  crash the pipeline.
"""

from __future__ import annotations

import csv
import io
import json
import logging
import re
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any

import httpx

from ribosome_state_annotator.cache import DEFAULT_CACHE_ROOT

logger = logging.getLogger(__name__)

RADDB_URL_TEMPLATE = (
    "https://radtool.rc.northeastern.edu/public_database/RADdb.{rad_date}.LSUSSU.csv"
)
"""Format string for the RADdb LSU/SSU CSV. ``rad_date`` is ``YYYYMMDD``."""

RADDB_CSV_FILENAME = "RADdb.LSUSSU.csv"
RADDB_METADATA_FILENAME = "RADdb.LSUSSU.metadata.json"

STALENESS_DAYS = 7
"""Local file is considered stale and triggers an online check after this many days."""

LATEST_LOOKBACK_DAYS = 60
"""Maximum number of days walked back when probing for the latest RADdb file."""

DEFAULT_HTTP_TIMEOUT = 30.0

# CSV column names exposed publicly in v1. The CSV stores them with a
# trailing period — preserved verbatim in the lookup so future versions
# can expose the remaining tilt / translation columns without re-parsing.
_BODY_ROT_COLUMN = "body rot."
_HEAD_ROT_COLUMN = "head rot."

_REQUIRED_COLUMNS = ("RCSB", "LSU chain ID", "SSU chain ID", _BODY_ROT_COLUMN, _HEAD_ROT_COLUMN)

_RAD_DATE_RE = re.compile(r"RADdb\.(\d{8})\.LSUSSU\.csv")


# ---------------------------------------------------------------------------
# Public dataclasses
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class RADdbMetadata:
    """Lightweight record of which RADdb file is currently cached."""

    source_url: str
    downloaded_at: datetime
    rad_date: str
    local_file: str = RADDB_CSV_FILENAME

    def to_dict(self) -> dict[str, str]:
        return {
            "source_url": self.source_url,
            "downloaded_at": self.downloaded_at.replace(microsecond=0).isoformat(),
            "rad_date": self.rad_date,
            "local_file": self.local_file,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> RADdbMetadata:
        return cls(
            source_url=str(data["source_url"]),
            downloaded_at=_parse_iso(str(data["downloaded_at"])),
            rad_date=str(data["rad_date"]),
            local_file=str(data.get("local_file", RADDB_CSV_FILENAME)),
        )


@dataclass(frozen=True)
class RADdbDataset:
    """In-memory RADdb dataset: the per-row lookup plus its metadata."""

    metadata: RADdbMetadata
    lookup: dict[tuple[str, str, str], dict[str, Any]]
    """Keyed by ``(pdb_id_upper, lsu_chain_id, ssu_chain_id)``."""

    duplicate_keys: frozenset[tuple[str, str, str]] = frozenset()
    """Keys that appeared in more than one CSV row. Matched lookups for
    these keys return ``None`` so we never silently pick one of the rows."""


# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------


def get_raddb_cache_dir(cache_root: Path | None = None) -> Path:
    """Resolve the RADdb sub-directory inside the package cache."""
    root = cache_root if cache_root is not None else DEFAULT_CACHE_ROOT
    return root / "raddb"


def get_local_raddb_csv_path(cache_root: Path | None = None) -> Path:
    return get_raddb_cache_dir(cache_root) / RADDB_CSV_FILENAME


def get_local_raddb_metadata_path(cache_root: Path | None = None) -> Path:
    return get_raddb_cache_dir(cache_root) / RADDB_METADATA_FILENAME


# ---------------------------------------------------------------------------
# Metadata I/O
# ---------------------------------------------------------------------------


def load_raddb_metadata(cache_root: Path | None = None) -> RADdbMetadata | None:
    """Read cached RADdb metadata, returning ``None`` if absent or malformed."""
    path = get_local_raddb_metadata_path(cache_root)
    if not path.is_file():
        return None
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        logger.warning("RADdb metadata read failed (%s); treating as missing", exc)
        return None
    if not isinstance(data, dict):
        logger.warning("RADdb metadata is not a JSON object; treating as missing")
        return None
    try:
        return RADdbMetadata.from_dict(data)
    except (KeyError, ValueError) as exc:
        logger.warning("RADdb metadata schema mismatch (%s); treating as missing", exc)
        return None


def save_raddb_metadata(metadata: RADdbMetadata, cache_root: Path | None = None) -> Path:
    path = get_local_raddb_metadata_path(cache_root)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(metadata.to_dict(), indent=2) + "\n", encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# Latest-URL discovery + download
# ---------------------------------------------------------------------------


def find_latest_raddb_url(
    *,
    client: httpx.Client | None = None,
    today: datetime | None = None,
    lookback_days: int = LATEST_LOOKBACK_DAYS,
) -> tuple[str, str] | None:
    """HEAD-walk back ``lookback_days`` days from today.

    Returns ``(url, rad_date_YYYYMMDD)`` for the most recent file that
    responds with HTTP 200, or ``None`` if nothing within the window is
    reachable.
    """
    now = today or datetime.now(timezone.utc)
    http: httpx.Client = client if client is not None else httpx.Client(timeout=DEFAULT_HTTP_TIMEOUT)
    try:
        for delta in range(lookback_days + 1):
            candidate = now - timedelta(days=delta)
            rad_date = candidate.strftime("%Y%m%d")
            url = RADDB_URL_TEMPLATE.format(rad_date=rad_date)
            try:
                response = http.head(url, follow_redirects=True)
            except httpx.HTTPError as exc:
                logger.debug("RADdb HEAD %s failed: %s", url, exc)
                continue
            if response.status_code == 200:
                return url, rad_date
        return None
    finally:
        if client is None:
            http.close()


def download_raddb_csv(
    url: str,
    rad_date: str,
    *,
    cache_root: Path | None = None,
    client: httpx.Client | None = None,
) -> tuple[Path, RADdbMetadata]:
    """Download the RADdb CSV at ``url`` and persist alongside metadata.

    Raises :class:`httpx.HTTPError` on network failure or non-200
    status. Callers are responsible for catching and degrading
    gracefully — the annotation pipeline must never crash on this.
    """
    http: httpx.Client = client if client is not None else httpx.Client(timeout=DEFAULT_HTTP_TIMEOUT)
    try:
        response = http.get(url, follow_redirects=True)
        response.raise_for_status()
        body = response.content
    finally:
        if client is None:
            http.close()

    csv_path = get_local_raddb_csv_path(cache_root)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    csv_path.write_bytes(body)
    metadata = RADdbMetadata(
        source_url=url,
        downloaded_at=datetime.now(timezone.utc),
        rad_date=rad_date,
    )
    save_raddb_metadata(metadata, cache_root=cache_root)
    logger.info("downloaded RADdb %s (%d bytes) to %s", rad_date, len(body), csv_path)
    return csv_path, metadata


# ---------------------------------------------------------------------------
# Refresh orchestration
# ---------------------------------------------------------------------------


def ensure_raddb_available(
    *,
    cache_root: Path | None = None,
    force_refresh: bool = False,
    client: httpx.Client | None = None,
    now: datetime | None = None,
) -> RADdbMetadata | None:
    """Make sure a local RADdb CSV exists, refreshing if stale.

    Returns the metadata describing the (possibly updated) local file,
    or ``None`` if no local file exists and no download succeeded.

    Behaviour, mirroring the spec §RADDB STARTUP BEHAVIOUR section:

    1. Local file missing → attempt download of latest.
    2. Local file present + ``force_refresh`` → check online for newer;
       download if newer is available.
    3. Local file present + age ≥ ``STALENESS_DAYS`` → check online for
       newer; download if newer is available.
    4. Local file present + fresh → return current metadata, no network.
    5. Any network failure with a local file present → keep the local
       file and log a warning.
    6. Any network failure without a local file → return ``None``.
    """
    csv_path = get_local_raddb_csv_path(cache_root)
    metadata = load_raddb_metadata(cache_root)
    now_dt = now or datetime.now(timezone.utc)

    if not csv_path.is_file() or metadata is None:
        logger.info("RADdb file missing; checking for latest version")
        return _try_refresh(cache_root=cache_root, client=client, now=now_dt)

    if force_refresh:
        logger.info("RADdb refresh forced; checking for latest version")
        refreshed = _try_refresh(
            cache_root=cache_root, client=client, current=metadata, now=now_dt
        )
        return refreshed if refreshed is not None else metadata

    age = now_dt - _aware(metadata.downloaded_at)
    if age >= timedelta(days=STALENESS_DAYS):
        logger.info(
            "RADdb cache is %d days old; checking for newer version",
            age.days,
        )
        refreshed = _try_refresh(
            cache_root=cache_root, client=client, current=metadata, now=now_dt
        )
        return refreshed if refreshed is not None else metadata

    logger.info("using cached RADdb %s (downloaded %s)", metadata.rad_date, metadata.downloaded_at.date())
    return metadata


def _try_refresh(
    *,
    cache_root: Path | None,
    client: httpx.Client | None,
    current: RADdbMetadata | None = None,
    now: datetime,
) -> RADdbMetadata | None:
    """Resolve the latest URL and download it if newer than ``current``.

    Returns the new metadata on successful download. Returns ``None`` if
    the online check or download failed (the caller decides how to
    degrade). If the online file is not newer than ``current`` it
    returns ``current`` unchanged.
    """
    try:
        latest = find_latest_raddb_url(client=client, today=now)
    except httpx.HTTPError as exc:
        logger.warning("RADdb latest-version probe failed: %s", exc)
        return None

    if latest is None:
        if current is not None:
            logger.warning("no RADdb release found in the last %d days", LATEST_LOOKBACK_DAYS)
        else:
            logger.warning("RADdb unavailable: no release found and no cached file")
        return None

    url, rad_date = latest
    if current is not None and rad_date <= current.rad_date:
        logger.info("RADdb cache is already at the latest version %s", current.rad_date)
        return current

    try:
        _, metadata = download_raddb_csv(url, rad_date, cache_root=cache_root, client=client)
    except httpx.HTTPError as exc:
        logger.warning("RADdb download failed for %s: %s", url, exc)
        return None
    return metadata


# ---------------------------------------------------------------------------
# CSV parsing
# ---------------------------------------------------------------------------


def load_raddb_table(
    *,
    cache_root: Path | None = None,
    csv_path: Path | None = None,
) -> list[dict[str, str]] | None:
    """Read the cached RADdb CSV as a list of row dicts.

    Returns ``None`` if the file is absent or unreadable.
    """
    path = csv_path or get_local_raddb_csv_path(cache_root)
    if not path.is_file():
        return None
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as exc:
        logger.warning("RADdb CSV read failed (%s); treating as missing", exc)
        return None
    reader = csv.DictReader(io.StringIO(text))
    return list(reader)


def build_raddb_lookup(
    rows: list[dict[str, str]],
) -> tuple[dict[tuple[str, str, str], dict[str, Any]], frozenset[tuple[str, str, str]]]:
    """Build the ``(pdb, lsu, ssu) → row`` lookup dict from raw CSV rows.

    Returns a ``(lookup, duplicate_keys)`` pair. Rows whose key appears
    more than once go into the dict (last-write-wins) but the key is
    recorded in ``duplicate_keys`` so :func:`get_motion_metrics` returns
    ``None`` for them rather than silently picking one row.

    Required columns: ``RCSB``, ``LSU chain ID``, ``SSU chain ID``,
    ``body rot.``, ``head rot.``. Rows missing the PDB ID or either
    chain ID are skipped.
    """
    lookup: dict[tuple[str, str, str], dict[str, Any]] = {}
    duplicates: set[tuple[str, str, str]] = set()
    for row in rows:
        pdb_id = (row.get("RCSB") or "").strip().upper()
        lsu_chain = (row.get("LSU chain ID") or "").strip()
        ssu_chain = (row.get("SSU chain ID") or "").strip()
        if not pdb_id or not lsu_chain or not ssu_chain:
            continue
        key = (pdb_id, lsu_chain, ssu_chain)
        if key in lookup:
            duplicates.add(key)
        lookup[key] = {
            "intersubunit_rotation": _coerce_float(row.get(_BODY_ROT_COLUMN)),
            "ssu_head_rotation": _coerce_float(row.get(_HEAD_ROT_COLUMN)),
            # The full row is kept so future versions can surface tilt /
            # translation without changing the cache format.
            "raw": dict(row),
        }
    if duplicates:
        logger.warning(
            "RADdb contains %d duplicate (pdb, lsu, ssu) keys; lookups for those keys return null",
            len(duplicates),
        )
    return lookup, frozenset(duplicates)


def load_raddb_dataset(
    *,
    cache_root: Path | None = None,
    metadata: RADdbMetadata | None = None,
) -> RADdbDataset | None:
    """Convenience: load metadata + CSV + lookup in one call."""
    meta = metadata if metadata is not None else load_raddb_metadata(cache_root)
    if meta is None:
        return None
    rows = load_raddb_table(cache_root=cache_root)
    if rows is None:
        logger.warning("RADdb metadata present but CSV missing at %s", get_local_raddb_csv_path(cache_root))
        return None
    if not rows:
        logger.warning("RADdb CSV at %s is empty", get_local_raddb_csv_path(cache_root))
        return RADdbDataset(metadata=meta, lookup={}, duplicate_keys=frozenset())
    missing = [c for c in _REQUIRED_COLUMNS if c not in rows[0]]
    if missing:
        logger.warning("RADdb CSV missing required columns: %s", missing)
        return None
    lookup, duplicates = build_raddb_lookup(rows)
    return RADdbDataset(metadata=meta, lookup=lookup, duplicate_keys=duplicates)


# ---------------------------------------------------------------------------
# Lookup helper
# ---------------------------------------------------------------------------


def get_motion_metrics(
    dataset: RADdbDataset | None,
    pdb_id: str,
    lsu_chain_id: str,
    ssu_chain_id: str,
) -> dict[str, Any] | None:
    """Look up RADdb motion metrics for one ``(pdb, lsu, ssu)`` triple.

    Returns ``{"intersubunit_rotation": float, "ssu_head_rotation": float}``
    on a unique match, ``None`` otherwise.
    """
    if dataset is None:
        return None
    key = (pdb_id.upper(), lsu_chain_id, ssu_chain_id)
    if key in dataset.duplicate_keys:
        logger.warning("RADdb lookup ambiguous: multiple rows for %s; returning null metrics", key)
        return None
    row = dataset.lookup.get(key)
    if row is None:
        logger.info("RADdb match not found for %s", key)
        return None
    return {
        "intersubunit_rotation": row["intersubunit_rotation"],
        "ssu_head_rotation": row["ssu_head_rotation"],
    }


def list_raddb_files(cache_root: Path | None = None) -> int:
    """Number of files in the RADdb cache directory (for cache info)."""
    raddb_dir = get_raddb_cache_dir(cache_root)
    if not raddb_dir.is_dir():
        return 0
    return sum(1 for entry in raddb_dir.iterdir() if entry.is_file())


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _coerce_float(value: str | None) -> float | None:
    if value is None:
        return None
    stripped = value.strip()
    if not stripped:
        return None
    try:
        return float(stripped)
    except ValueError:
        return None


def _parse_iso(value: str) -> datetime:
    """Parse an ISO-8601 timestamp into a timezone-aware datetime (UTC if naive)."""
    parsed = datetime.fromisoformat(value)
    return _aware(parsed)


def _aware(dt: datetime) -> datetime:
    return dt if dt.tzinfo is not None else dt.replace(tzinfo=timezone.utc)


def parse_rad_date_from_url(url: str) -> str | None:
    """Extract the ``YYYYMMDD`` date from a RADdb URL, or ``None``."""
    match = _RAD_DATE_RE.search(url)
    return match.group(1) if match else None


__all__ = [
    "RADDB_CSV_FILENAME",
    "RADDB_METADATA_FILENAME",
    "RADDB_URL_TEMPLATE",
    "RADdbDataset",
    "RADdbMetadata",
    "build_raddb_lookup",
    "download_raddb_csv",
    "ensure_raddb_available",
    "find_latest_raddb_url",
    "get_local_raddb_csv_path",
    "get_local_raddb_metadata_path",
    "get_motion_metrics",
    "get_raddb_cache_dir",
    "list_raddb_files",
    "load_raddb_dataset",
    "load_raddb_metadata",
    "load_raddb_table",
    "parse_rad_date_from_url",
    "save_raddb_metadata",
]
