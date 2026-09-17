"""BGSU RNA 3D Hub representative-set ("NR list") integration for rRNA Rfam tags.

The BGSU RNA group publishes a weekly *representative set* release of
every RNA chain in the PDB, grouped into equivalence classes. The
"full" CSV export of a release carries one row per IFE (integrated
functional element — one chain, or several ``+``-joined chains) with a
``rfam`` column giving the Rfam family each chain maps to:

    https://rna.bgsu.edu/rna3dhub/nrlist/download/rna/<release>/all/csv/full

Why this source exists alongside :mod:`.rfam_pdb_region`: EBI's own
``pdb_full_region`` scan is regenerated irregularly and can lag PDB
releases by months (it stalled at 2026-05-18 while BGSU kept shipping
weekly), which left every newer ribosome deposit without any rRNA
identity and skipped as ``partial_ribosome_missing_ssu_or_lsu``. BGSU
maps chains by the same Rfam.pdb file when available and otherwise by
running cmsearch itself against Rfam covariance models that already
have PDB chains, so its families agree with EBI's best-score pick on
~99 % of shared chains and are available within a week of release.

Format rules (BGSU documentation):

- ``ife_id`` is ``PDB|model|chain`` or several joined by ``+``.
- ``rfam`` has one family per chain in the IFE, joined by ``+`` in the
  same order. ``NA`` means the chain maps to no family.
- Chains with no family may be absent from the file altogether (mRNAs
  typically are).

Design mirrors :mod:`.rfam_pdb_region` / :mod:`.raddb`:

- Cached under ``~/.cache/ribosome-state-annotator/bgsu_nr/`` as a
  gzipped copy of the CSV plus a ``metadata.json`` sidecar recording
  the release number and date.
- There is no ``current`` alias for the download, so a refresh first
  reads the release number from the ``nrlist/release/current`` page
  and re-downloads only when it differs from the cached release.
- A 7-day staleness window triggers that check. Every failure mode
  returns ``None`` so the annotation pipeline continues.
"""

from __future__ import annotations

import csv
import gzip
import io
import json
import logging
import re
from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any

import httpx

from ribosome_state_annotator.cache import DEFAULT_CACHE_ROOT

logger = logging.getLogger(__name__)

BGSU_NR_CURRENT_RELEASE_URL = "https://rna.bgsu.edu/rna3dhub/nrlist/release/current"
"""HTML page for the most recent representative-set release; its title
names the release number (e.g. ``Representative set 4.57``)."""

BGSU_NR_DOWNLOAD_URL_TEMPLATE = (
    "https://rna.bgsu.edu/rna3dhub/nrlist/download/rna/{release}/all/csv/full"
)
"""Full CSV export (all resolutions, one row per IFE with Rfam column)."""

BGSU_NR_FILENAME = "nr_all_full.csv.gz"
BGSU_NR_METADATA_FILENAME = "nr_all_full.metadata.json"

STALENESS_DAYS = 7
"""Local file is considered stale and triggers an online check after this many days."""

DEFAULT_HTTP_TIMEOUT = 120.0
"""The full CSV is ~10 MB; allow a generous download window."""

NO_FAMILY = "NA"
"""Sentinel BGSU uses in the ``rfam`` column for chains with no Rfam family."""

_RELEASE_PATTERNS = (
    re.compile(r"Representative set\s+(\d+\.\d+)"),
    re.compile(r"Release\s+(\d+\.\d+)"),
    re.compile(r"/nrlist/release/rna/(\d+\.\d+)/"),
)
_RELEASE_DATE_PATTERN = re.compile(r"Release\s+\d+\.\d+,\s*(\d{4}-\d{2}-\d{2})")


# ---------------------------------------------------------------------------
# Public dataclasses
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class BgsuNrListMetadata:
    """Lightweight record describing the cached BGSU representative-set CSV."""

    source_url: str
    downloaded_at: datetime
    release: str
    release_date: str | None = None
    local_file: str = BGSU_NR_FILENAME

    def to_dict(self) -> dict[str, str | None]:
        return {
            "source_url": self.source_url,
            "downloaded_at": self.downloaded_at.replace(microsecond=0).isoformat(),
            "release": self.release,
            "release_date": self.release_date,
            "local_file": self.local_file,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> BgsuNrListMetadata:
        return cls(
            source_url=str(data["source_url"]),
            downloaded_at=_parse_iso(str(data["downloaded_at"])),
            release=str(data["release"]),
            release_date=data.get("release_date"),
            local_file=str(data.get("local_file", BGSU_NR_FILENAME)),
        )


@dataclass(frozen=True)
class BgsuNrListDataset:
    """In-memory BGSU representative-set dataset: per-chain Rfam lookup."""

    metadata: BgsuNrListMetadata
    lookup: dict[tuple[str, str], str] = field(default_factory=dict)
    """Keyed by ``(pdb_id_lower, auth_chain_id)``. Value is the Rfam
    accession BGSU assigned to that chain. Chains BGSU marks ``NA`` are
    not present."""


# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------


def get_bgsu_nr_cache_dir(cache_root: Path | None = None) -> Path:
    """Resolve the ``bgsu_nr/`` sub-directory inside the package cache."""
    root = cache_root if cache_root is not None else DEFAULT_CACHE_ROOT
    return root / "bgsu_nr"


def get_local_bgsu_nr_file_path(cache_root: Path | None = None) -> Path:
    return get_bgsu_nr_cache_dir(cache_root) / BGSU_NR_FILENAME


def get_local_bgsu_nr_metadata_path(cache_root: Path | None = None) -> Path:
    return get_bgsu_nr_cache_dir(cache_root) / BGSU_NR_METADATA_FILENAME


def bgsu_nr_download_url(release: str) -> str:
    return BGSU_NR_DOWNLOAD_URL_TEMPLATE.format(release=release)


# ---------------------------------------------------------------------------
# Metadata I/O
# ---------------------------------------------------------------------------


def load_bgsu_nr_metadata(cache_root: Path | None = None) -> BgsuNrListMetadata | None:
    """Read cached metadata, returning ``None`` if absent or malformed."""
    path = get_local_bgsu_nr_metadata_path(cache_root)
    if not path.is_file():
        return None
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        logger.warning("BGSU NR-list metadata read failed (%s); treating as missing", exc)
        return None
    if not isinstance(data, dict):
        logger.warning("BGSU NR-list metadata is not a JSON object; treating as missing")
        return None
    try:
        return BgsuNrListMetadata.from_dict(data)
    except (KeyError, ValueError) as exc:
        logger.warning("BGSU NR-list metadata schema mismatch (%s); treating as missing", exc)
        return None


def save_bgsu_nr_metadata(metadata: BgsuNrListMetadata, cache_root: Path | None = None) -> Path:
    path = get_local_bgsu_nr_metadata_path(cache_root)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(metadata.to_dict(), indent=2) + "\n", encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# Release discovery + download
# ---------------------------------------------------------------------------


def parse_current_release(html: str) -> tuple[str, str | None] | None:
    """Extract ``(release, release_date)`` from the current-release page HTML.

    Returns ``None`` when no release number can be found. The date is
    optional (``None`` when the page layout omits it).
    """
    release: str | None = None
    for pattern in _RELEASE_PATTERNS:
        match = pattern.search(html)
        if match:
            release = match.group(1)
            break
    if release is None:
        return None
    date_match = _RELEASE_DATE_PATTERN.search(html)
    return release, (date_match.group(1) if date_match else None)


def probe_current_release(*, client: httpx.Client | None = None) -> tuple[str, str | None] | None:
    """Fetch the current-release page and return ``(release, release_date)``.

    Returns ``None`` on any network failure or if the page cannot be parsed.
    """
    http: httpx.Client = (
        client if client is not None else httpx.Client(timeout=DEFAULT_HTTP_TIMEOUT)
    )
    try:
        response = http.get(BGSU_NR_CURRENT_RELEASE_URL, follow_redirects=True)
    except httpx.HTTPError as exc:
        logger.debug("BGSU NR-list release probe %s failed: %s", BGSU_NR_CURRENT_RELEASE_URL, exc)
        return None
    finally:
        if client is None:
            http.close()
    if response.status_code != 200:
        logger.warning("BGSU NR-list release page returned HTTP %d", response.status_code)
        return None
    parsed = parse_current_release(response.text)
    if parsed is None:
        logger.warning("BGSU NR-list release page did not name a release number")
    return parsed


def download_bgsu_nr_list(
    release: str,
    *,
    release_date: str | None = None,
    cache_root: Path | None = None,
    client: httpx.Client | None = None,
) -> tuple[Path, BgsuNrListMetadata]:
    """Download the full CSV for ``release`` and persist it gzipped with metadata.

    Raises :class:`httpx.HTTPError` on network failure or non-200
    status. Callers are responsible for catching and degrading
    gracefully.
    """
    url = bgsu_nr_download_url(release)
    http: httpx.Client = (
        client if client is not None else httpx.Client(timeout=DEFAULT_HTTP_TIMEOUT)
    )
    try:
        response = http.get(url, follow_redirects=True)
        response.raise_for_status()
        body = response.content
    finally:
        if client is None:
            http.close()

    file_path = get_local_bgsu_nr_file_path(cache_root)
    file_path.parent.mkdir(parents=True, exist_ok=True)
    file_path.write_bytes(gzip.compress(body))
    metadata = BgsuNrListMetadata(
        source_url=url,
        downloaded_at=datetime.now(timezone.utc),
        release=release,
        release_date=release_date,
    )
    save_bgsu_nr_metadata(metadata, cache_root=cache_root)
    logger.info(
        "downloaded BGSU representative set %s (%d bytes; released %s) to %s",
        release,
        len(body),
        release_date or "unknown date",
        file_path,
    )
    return file_path, metadata


# ---------------------------------------------------------------------------
# Refresh orchestration
# ---------------------------------------------------------------------------


def ensure_bgsu_nr_list_available(
    *,
    cache_root: Path | None = None,
    force_refresh: bool = False,
    client: httpx.Client | None = None,
    now: datetime | None = None,
) -> BgsuNrListMetadata | None:
    """Make sure a local BGSU representative-set CSV exists, refreshing if stale.

    Returns the metadata describing the (possibly updated) local file,
    or ``None`` if no local file exists and no download succeeded.

    1. Local file missing → resolve current release, download.
    2. Local file present + ``force_refresh`` → resolve current release;
       download if it differs from the cached release.
    3. Local file present + age ≥ ``STALENESS_DAYS`` → same as 2.
    4. Local file present + fresh → return current metadata, no network.
    5. Any network failure with a local file present → keep the local
       file and log a warning.
    6. Any network failure without a local file → return ``None``.
    """
    file_path = get_local_bgsu_nr_file_path(cache_root)
    metadata = load_bgsu_nr_metadata(cache_root)
    now_dt = now or datetime.now(timezone.utc)

    if not file_path.is_file() or metadata is None:
        logger.info("BGSU representative-set file missing; downloading")
        return _try_refresh(cache_root=cache_root, client=client, current=None)

    if force_refresh:
        logger.info("BGSU representative-set refresh forced; checking for a newer release")
        refreshed = _try_refresh(cache_root=cache_root, client=client, current=metadata)
        return refreshed if refreshed is not None else metadata

    age = now_dt - _aware(metadata.downloaded_at)
    if age >= timedelta(days=STALENESS_DAYS):
        logger.info(
            "BGSU representative-set cache is %d days old; checking for a newer release",
            age.days,
        )
        refreshed = _try_refresh(cache_root=cache_root, client=client, current=metadata)
        return refreshed if refreshed is not None else metadata

    logger.info(
        "using cached BGSU representative set %s (downloaded %s)",
        metadata.release,
        metadata.downloaded_at.date(),
    )
    return metadata


def _try_refresh(
    *,
    cache_root: Path | None,
    client: httpx.Client | None,
    current: BgsuNrListMetadata | None,
) -> BgsuNrListMetadata | None:
    """Resolve the current release and download if it differs from ``current``.

    Returns the new metadata on successful download; returns refreshed
    ``current`` metadata if the release is unchanged; returns ``None``
    if the probe or download fails.
    """
    probed = probe_current_release(client=client)
    if probed is None:
        if current is not None:
            logger.warning(
                "BGSU representative-set release probe failed; keeping cached release %s",
                current.release,
            )
        else:
            logger.warning("BGSU representative set unavailable: probe failed and no cached file")
        return None

    release, release_date = probed
    if current is not None and release == current.release:
        logger.info("BGSU representative set is already at release %s", release)
        # Touch the metadata so the next staleness check resets the clock.
        refreshed = BgsuNrListMetadata(
            source_url=current.source_url,
            downloaded_at=datetime.now(timezone.utc),
            release=current.release,
            release_date=current.release_date or release_date,
            local_file=current.local_file,
        )
        save_bgsu_nr_metadata(refreshed, cache_root=cache_root)
        return refreshed

    try:
        _, metadata = download_bgsu_nr_list(
            release, release_date=release_date, cache_root=cache_root, client=client
        )
    except httpx.HTTPError as exc:
        logger.warning("BGSU representative-set download failed: %s", exc)
        return None
    return metadata


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------


def parse_bgsu_nr_rows(text: str) -> dict[tuple[str, str], str]:
    """Build the ``(pdb_id_lower, chain) -> rfam_acc`` lookup from CSV text.

    Composite IFEs (``A+B``) are split and zipped positionally with
    their ``+``-joined Rfam families. Rows whose chain and family counts
    differ are skipped with a warning; ``NA`` families are dropped.
    """
    reader = csv.DictReader(io.StringIO(text))
    if (
        reader.fieldnames is None
        or "ife_id" not in reader.fieldnames
        or "rfam" not in reader.fieldnames
    ):
        logger.warning("BGSU representative-set CSV lacks ife_id/rfam columns; ignoring file")
        return {}
    lookup: dict[tuple[str, str], str] = {}
    skipped = 0
    for row in reader:
        ife_field = (row.get("ife_id") or "").strip()
        rfam_field = (row.get("rfam") or "").strip()
        if not ife_field:
            continue
        ifes = ife_field.split("+")
        families = rfam_field.split("+") if rfam_field else []
        if len(ifes) != len(families):
            skipped += 1
            continue
        for ife, family in zip(ifes, families, strict=True):
            parts = ife.split("|")
            if len(parts) < 3:
                skipped += 1
                continue
            family = family.strip().rstrip("*")
            if not family or family == NO_FAMILY:
                continue
            lookup[(parts[0].lower(), parts[2])] = family
    if skipped:
        logger.warning(
            "BGSU representative-set CSV: skipped %d malformed IFE/Rfam entries", skipped
        )
    return lookup


def load_bgsu_nr_list_dataset(
    *,
    cache_root: Path | None = None,
    metadata: BgsuNrListMetadata | None = None,
) -> BgsuNrListDataset | None:
    """Parse the cached CSV into a :class:`BgsuNrListDataset`.

    Returns ``None`` if no metadata or file is present.
    """
    metadata = metadata or load_bgsu_nr_metadata(cache_root)
    if metadata is None:
        return None
    file_path = get_bgsu_nr_cache_dir(cache_root) / metadata.local_file
    if not file_path.is_file():
        logger.warning("BGSU representative-set file missing at %s", file_path)
        return None
    with gzip.open(file_path, mode="rt", encoding="utf-8", newline="") as fh:
        text = fh.read()
    return BgsuNrListDataset(metadata=metadata, lookup=parse_bgsu_nr_rows(text))


# ---------------------------------------------------------------------------
# Lookup
# ---------------------------------------------------------------------------


def get_bgsu_rfam_for_chain(
    dataset: BgsuNrListDataset | None,
    pdb_id: str,
    chain_id: str,
) -> str | None:
    """Return BGSU's Rfam accession for ``(pdb_id, chain_id)``, or ``None``."""
    if dataset is None:
        return None
    return dataset.lookup.get((pdb_id.lower(), chain_id))


def get_bgsu_rfam_mapping_for_pdb(
    dataset: BgsuNrListDataset | None,
    pdb_id: str,
) -> dict[str, list[str]]:
    """Return ``{chain_id: [rfam_acc]}`` for a PDB.

    Same shape as :func:`rfam_pdb_region.get_rfam_mapping_for_pdb` so
    the api.py merge step can combine the two sources.
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
# CLI helpers (`ribostate nrlist info` etc.)
# ---------------------------------------------------------------------------


def list_bgsu_nr_files(cache_root: Path | None = None) -> int:
    """Return the size in bytes of the cached CSV, or 0 if absent."""
    file_path = get_local_bgsu_nr_file_path(cache_root)
    if not file_path.is_file():
        return 0
    return file_path.stat().st_size


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _parse_iso(value: str) -> datetime:
    dt = datetime.fromisoformat(value)
    if dt.tzinfo is None:
        return dt.replace(tzinfo=timezone.utc)
    return dt


def _aware(dt: datetime) -> datetime:
    if dt.tzinfo is None:
        return dt.replace(tzinfo=timezone.utc)
    return dt
