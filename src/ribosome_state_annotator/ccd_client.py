"""Per-component Chemical Component Dictionary (CCD) client.

Fetches and parses individual residue definitions from the PDB's
per-component CCD files at

    https://files.rcsb.org/ligands/view/<comp_id>.cif

The full CCD entry carries fields that the structure mmCIF doesn't
ship — most importantly ``_chem_comp.mon_nstd_parent_comp_id``, which
identifies the canonical parent monomer for non-standard residues
(e.g. ``U8U`` → ``U``, ``OMG`` → ``G``, ``PSU`` → ``U``). This is the
authoritative source for parent-base resolution, and the ~2 KB cost
per lookup is small enough to fetch lazily when Gemmi's built-in
tabulated dictionary returns a blank or missing entry.

Every failure path (network error, non-200 status, parse error)
returns ``None`` so the caller can degrade gracefully — the codon /
anticodon extraction must never crash because one residue's CCD
definition is unreachable.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import gemmi
import httpx

from ribosome_state_annotator.cache import Cache

logger = logging.getLogger(__name__)

CCD_URL_TEMPLATE = "https://files.rcsb.org/ligands/view/{comp_id}.cif"
"""Per-component CCD CIF endpoint. ``comp_id`` is upper-cased."""

DEFAULT_HTTP_TIMEOUT = 30.0


@dataclass(frozen=True)
class ChemCompInfo:
    """Subset of the CCD entry that the package needs.

    All fields except ``id`` may be ``None`` — the PDB occasionally
    ships partial entries (e.g. one_letter_code blank for some
    phosphate-monomer forms).
    """

    id: str
    parent_comp_id: str | None
    one_letter_code: str | None
    type: str | None
    name: str | None


def fetch_chem_comp(
    comp_id: str,
    *,
    cache: Cache | None = None,
    client: httpx.Client | None = None,
) -> ChemCompInfo | None:
    """Return the parsed CCD info for ``comp_id``, or ``None`` on failure.

    Caches the raw CIF bytes per ``comp_id`` under the ``ccd/``
    namespace; the second call for the same comp_id is served from
    disk.
    """
    comp_upper = comp_id.upper().strip()
    if not comp_upper:
        return None
    if cache is not None:
        cached = cache.get_ccd_cif(comp_upper)
        if cached is not None:
            return _parse_ccd_cif(comp_upper, cached)

    url = CCD_URL_TEMPLATE.format(comp_id=comp_upper)
    http = client if client is not None else httpx.Client(timeout=DEFAULT_HTTP_TIMEOUT)
    try:
        try:
            response = http.get(url, follow_redirects=True)
        except httpx.HTTPError as exc:
            logger.warning("CCD fetch failed for %s: %s", comp_upper, exc)
            return None
        if response.status_code != 200:
            logger.warning("CCD returned HTTP %d for %s", response.status_code, comp_upper)
            return None
        body = response.content
    finally:
        if client is None:
            http.close()

    if cache is not None:
        cache.put_ccd_cif(comp_upper, body)
    return _parse_ccd_cif(comp_upper, body)


def _parse_ccd_cif(comp_id: str, data: bytes) -> ChemCompInfo | None:
    """Extract the relevant ``_chem_comp.*`` fields from CCD CIF bytes."""
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as exc:
        logger.warning("CCD CIF is not valid UTF-8 for %s (%s)", comp_id, exc)
        return None
    try:
        doc = gemmi.cif.read_string(text)
        block = doc.sole_block()
    except (RuntimeError, ValueError) as exc:
        logger.warning("CCD CIF parse failed for %s: %s", comp_id, exc)
        return None

    return ChemCompInfo(
        id=comp_id,
        parent_comp_id=_pluck(block, "_chem_comp.mon_nstd_parent_comp_id"),
        one_letter_code=_pluck(block, "_chem_comp.one_letter_code"),
        type=_pluck(block, "_chem_comp.type"),
        name=_pluck(block, "_chem_comp.name"),
    )


def _pluck(block: gemmi.cif.Block, key: str) -> str | None:
    """Return a clean string for the given CIF tag, or ``None`` when
    the entry is missing, blank, or the CCD's "?" / "." sentinel.

    Note: ``Block.find_value`` may return ``None`` at runtime even
    though Gemmi's Python stubs declare ``str`` — we tolerate both.
    """
    raw: str | None = block.find_value(key)
    if not raw:
        return None
    value = raw.strip().strip('"').strip("'")
    if not value or value in {"?", "."}:
        return None
    return value


__all__ = [
    "CCD_URL_TEMPLATE",
    "ChemCompInfo",
    "fetch_chem_comp",
]
