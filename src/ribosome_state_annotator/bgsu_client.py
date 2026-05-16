"""BGSU correspondence-mapping API client (spec §5.2).

This module owns the HTTP layer for the BGSU RNA 3D Hub
`map_across_chains` endpoint and the parser that converts its JSON
response shape into a flat ``{reference_unit: [mapped_unit, ...]}`` dict.
The downstream filtering (§5.2.2) and :class:`CorrespondenceResult`
assembly live in :mod:`.correspondence`.

The live BGSU endpoint returns JSON only when the request includes
``format=json``; without it the response is tab-delimited ``text/plain``
regardless of the ``Accept`` header. :func:`fetch_correspondence`
always sends ``format=json``.

The spec acknowledges that BGSU's exact JSON field names may evolve
(spec §5.2.1 last paragraph). :func:`parse_alignment_response` is
therefore tolerant of:

- additional top-level fields (e.g. ``query`` echo);
- additional per-row fields;
- ``aligned_units`` as a synonym for ``mapped_units``;
- the same ``reference_unit`` appearing across multiple rows (merged).

Anything fundamentally malformed (no top-level ``alignment`` list,
non-object payload) raises :class:`~.exceptions.CorrespondenceMappingError`.
An empty alignment list is **not** an error — it is a valid response
that simply yields no mappings, and the §5.2.3 missing-anchor warning
machinery handles it downstream.
"""

from __future__ import annotations

import json
import logging
from collections.abc import Iterable
from typing import Any

import httpx

from ribosome_state_annotator.exceptions import (
    ApiRequestError,
    CorrespondenceMappingError,
)

logger = logging.getLogger(__name__)

BGSU_CORRESPONDENCE_URL = "https://rna.bgsu.edu/correspondence/map_across_chains"
"""BGSU RNA 3D Hub `map_across_chains` endpoint (spec §5.2)."""

DEFAULT_BGSU_SCOPE = "Rfam"
"""``scope=Rfam`` makes the alignment go across species via Infernal /
covariance models — the cross-species transfer we want for §3.3."""

DEFAULT_BGSU_RESOLUTION = "20.0A"
"""``resolution=20.0A`` permits low-resolution cryo-EM structures. Most
ribosome structures are ≤4 Å; this loose threshold is conservatively safe."""

DEFAULT_BGSU_DEPTH = 700
"""``depth=N`` widens BGSU's correspondence search beyond the default
top-N equivalence-class representatives. Without ``depth``, the
endpoint returns ~120 NR representatives per query and never includes
common ribosome PDBs like 5J7L or 4V5Q themselves. With ``depth=700``
the response grows to ~1260 rows and reliably includes the major
bacterial / eukaryotic ribosome deposits we care about.

Still returns at most ONE chain per PDB even at high depth — multi-chain
PDBs (e.g. 5FDV's ``1a`` vs ``2a``, or 5J7L's ``AA`` vs ``BA``) are
handled by the chain-substitution fallback in
:func:`correspondence.build_correspondence_result`."""


# ---------------------------------------------------------------------------
# HTTP
# ---------------------------------------------------------------------------


def fetch_correspondence(
    reference_units: Iterable[str],
    *,
    scope: str = DEFAULT_BGSU_SCOPE,
    resolution: str = DEFAULT_BGSU_RESOLUTION,
    depth: int = DEFAULT_BGSU_DEPTH,
    client: httpx.Client | None = None,
    timeout: float = 60.0,
) -> dict[str, list[str]]:
    """GET BGSU correspondence mappings for a batch of reference units.

    Returns a ``{reference_unit: [mapped_unit, ...]}`` dict (after passing
    through :func:`parse_alignment_response`). The returned mapped-unit
    lists span every aligned PDB in the same Rfam family — caller is
    expected to apply the §5.2.2 PDB-prefix + assembly-chain filter via
    :mod:`.correspondence`.

    Single batch query containing all the unit IDs (spec §5.2 default).

    Args:
        reference_units: Non-empty iterable of BGSU unit IDs to query.
        scope: BGSU scope parameter; default ``"Rfam"``.
        resolution: BGSU resolution parameter; default ``"20.0A"``.
        depth: BGSU depth parameter; default :data:`DEFAULT_BGSU_DEPTH`.
        client: Optional pre-built :class:`httpx.Client`. Used by tests and
            by the cache layer (step 6). When omitted, an ephemeral client
            is created and closed inside this function.
        timeout: Request timeout in seconds. Ignored when ``client`` is supplied.

    Raises:
        ValueError: if ``reference_units`` is empty.
        ApiRequestError: HTTP failure, non-200 status, non-JSON body, or
            non-object JSON.
        CorrespondenceMappingError: if the parsed JSON does not contain a
            top-level ``alignment`` or ``mappings`` list (see
            :func:`parse_alignment_response`).
    """
    units_list = list(reference_units)
    if not units_list:
        raise ValueError("reference_units must contain at least one unit ID")

    own_client = client is None
    http: httpx.Client = client if client is not None else httpx.Client(timeout=timeout)
    try:
        return _query_bgsu(
            http, units_list, scope=scope, resolution=resolution, depth=depth
        )
    finally:
        if own_client:
            http.close()


def _query_bgsu(
    http: httpx.Client,
    units: list[str],
    *,
    scope: str,
    resolution: str,
    depth: int,
) -> dict[str, list[str]]:
    """Single HTTP call to BGSU; raises the typed exceptions on failure."""
    params = {
        "id": ",".join(units),
        "scope": scope,
        "resolution": resolution,
        # See DEFAULT_BGSU_DEPTH; without this the response is limited to
        # ~120 NR equivalence-class representatives and excludes common
        # ribosome PDBs like 5J7L and 4V5Q.
        "depth": str(depth),
        # Live BGSU ignores the ``Accept`` header — the JSON response is
        # gated by an explicit ``format=json`` query parameter. Without
        # it the endpoint returns tab-delimited ``text/plain``.
        "format": "json",
    }
    headers = {"Accept": "application/json"}
    try:
        response = http.get(BGSU_CORRESPONDENCE_URL, params=params, headers=headers)
    except httpx.HTTPError as exc:
        raise ApiRequestError(f"BGSU correspondence request failed: {exc}") from exc
    if response.status_code != 200:
        raise ApiRequestError(f"BGSU correspondence returned HTTP {response.status_code}")
    try:
        payload: Any = response.json()
    except (json.JSONDecodeError, ValueError) as exc:
        raise ApiRequestError(f"BGSU correspondence returned non-JSON body: {exc}") from exc
    if not isinstance(payload, dict):
        raise ApiRequestError("BGSU correspondence returned non-object JSON")
    return parse_alignment_response(payload, query_units=units)


# ---------------------------------------------------------------------------
# Pure parser
# ---------------------------------------------------------------------------


def parse_alignment_response(
    payload: dict[str, Any], *, query_units: list[str] | None = None
) -> dict[str, list[str]]:
    """Parse a BGSU correspondence JSON response into ``{ref_unit: [mapped, …]}``.

    Live BGSU shape (May 2026, with ``format=json``)::

        {
          "query": [...],
          "mappings": [
            {
              "unit_id_list": ["9PJ7|1|A|G|910", "9PJ7|1|A|4OC|1389", ...],
              "pdb": "9PJ7",
              "rfam_EC_chain": [["RF00177", "NR_20.0_86966.1", "9PJ7|1|A"]],
              "sequence": "..."
            },
            ...
          ]
        }

    The mapped units inside each row's ``unit_id_list`` are positionally
    aligned with ``query_units`` (or with ``payload["query"]`` if that
    argument isn't supplied).

    Also accepts the spec §5.2.1 idealised shape::

        {
          "alignment": [
            {"reference_unit": "...", "mapped_units": ["...", ...]},
            ...
          ]
        }

    Tolerances applied to the idealised shape:

    - ``query`` is optional.
    - Per-row fields beyond ``reference_unit`` / ``mapped_units`` are ignored.
    - ``aligned_units`` is accepted as a synonym for ``mapped_units``.
    - A given ``reference_unit`` appearing across multiple rows has its
      ``mapped_units`` concatenated.
    - Malformed individual rows are silently skipped.

    Raises :class:`CorrespondenceMappingError` only for fundamental shape
    errors: non-dict payload, missing both ``mappings`` and ``alignment``,
    or ``mappings`` shape without a way to determine query units.
    """
    if not isinstance(payload, dict):
        raise CorrespondenceMappingError("BGSU response is not a JSON object")

    if isinstance(payload.get("mappings"), list):
        return _parse_mappings_shape(payload, query_units=query_units)

    if isinstance(payload.get("alignment"), list):
        return _parse_alignment_shape(payload)

    raise CorrespondenceMappingError(
        "BGSU response missing required 'mappings' or 'alignment' list"
    )


def _parse_mappings_shape(
    payload: dict[str, Any], *, query_units: list[str] | None
) -> dict[str, list[str]]:
    """Parse the live BGSU ``mappings`` shape (May 2026).

    Each ``query_unit`` is also injected into its own mapped list, so the
    annotating-the-reference-PDB case works: BGSU returns mappings to
    OTHER ribosomes, but when the target PDB IS the reference PDB (e.g.
    annotating 5J7L itself), the query unit needs to pass through the
    §5.2.2 filter on its own merits. For non-reference targets the
    self-injection is harmless — the PDB-prefix filter discards it.
    """
    if query_units is None:
        echoed_query = payload.get("query")
        if not isinstance(echoed_query, list) or not all(isinstance(q, str) for q in echoed_query):
            raise CorrespondenceMappingError(
                "BGSU `mappings` response provided without query_units and "
                "without a top-level `query` echo; cannot align positions"
            )
        query_units = [q for q in echoed_query if isinstance(q, str)]

    # Seed each query unit's mapped list with itself so the
    # annotating-the-reference-PDB case works.
    out: dict[str, list[str]] = {unit: [unit] for unit in query_units}
    mappings = payload["mappings"]
    if not isinstance(mappings, list):
        return out
    for row in mappings:
        if not isinstance(row, dict):
            continue
        unit_list = row.get("unit_id_list")
        if not isinstance(unit_list, list):
            continue
        for offset, query_unit in enumerate(query_units):
            if offset >= len(unit_list):
                break
            candidate = unit_list[offset]
            if isinstance(candidate, str) and candidate:
                out[query_unit].append(candidate)
    return out


def _parse_alignment_shape(payload: dict[str, Any]) -> dict[str, list[str]]:
    """Parse the spec §5.2.1 idealised ``alignment`` shape (also used by the
    cache layer's stored payloads)."""
    alignment = payload["alignment"]
    out: dict[str, list[str]] = {}
    if not isinstance(alignment, list):
        return out
    for row in alignment:
        if not isinstance(row, dict):
            continue
        ref = row.get("reference_unit")
        if not isinstance(ref, str) or not ref:
            continue
        mapped_raw = row.get("mapped_units")
        if mapped_raw is None:
            mapped_raw = row.get("aligned_units")
        if not isinstance(mapped_raw, list):
            continue
        valid = [u for u in mapped_raw if isinstance(u, str) and u]
        existing = out.get(ref)
        if existing is not None:
            existing.extend(valid)
        else:
            out[ref] = list(valid)
    return out
