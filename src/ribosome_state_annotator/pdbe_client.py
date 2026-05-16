"""PDBe REST API client for Rfam chain mappings.

The spec (§5.1.1) describes Rfam annotations as a field on RCSB's
GraphQL ``polymer_entity``. The live RCSB schema no longer supplies
Rfam for rRNA chains — 5J7L's 16S chain in 2026 carries only a
GenBank cross-reference. We therefore query PDBe's REST endpoint
``/nucleic_mappings/rfam/<pdb_id>`` for the authoritative chain-to-Rfam
mapping and augment ``ChainRef.rfam_accessions`` in :mod:`api` before
classification runs.

The PDBe response shape (one PDB ID per call):

::

    {
      "<pdb_id_lower>": {
        "Rfam": {
          "RF00177": {
            "identifier": "Bacterial small subunit ribosomal RNA",
            "family":     "SSU_rRNA_bacteria",
            "mappings": [
              {"chain_id": "AA", "entity_id": 1, ...},
              {"chain_id": "BA", "entity_id": 1, ...}
            ],
            ...
          },
          "RF02541": {...},
          ...
        }
      }
    }
"""

from __future__ import annotations

import json
import logging
from typing import Any

import httpx

from ribosome_state_annotator.exceptions import ApiRequestError

logger = logging.getLogger(__name__)

PDBE_RFAM_URL_TEMPLATE = "https://www.ebi.ac.uk/pdbe/api/v2/nucleic_mappings/rfam/{pdb_id}"
"""Endpoint template; substitute ``pdb_id`` lower-cased."""


# ---------------------------------------------------------------------------
# HTTP
# ---------------------------------------------------------------------------


def fetch_rfam_mappings(
    pdb_id: str,
    *,
    client: httpx.Client | None = None,
    timeout: float = 30.0,
) -> dict[str, list[str]]:
    """Return ``{chain_id: [Rfam_accession, ...]}`` for one PDB entry.

    Args:
        pdb_id: PDB accession (case-insensitive; lower-cased for the URL).
        client: Optional pre-built :class:`httpx.Client` (used by tests
            and by the cache layer). When omitted, an ephemeral client
            is created and closed inside this function.
        timeout: Request timeout in seconds. Ignored when ``client`` is supplied.

    Returns:
        Empty dict when PDBe has no Rfam mappings for the entry (the
        common case for protein-only structures). Otherwise a chain →
        list-of-accessions mapping suitable for augmenting
        :class:`~.models.ChainRef`.

    Raises:
        ApiRequestError: HTTP failure, non-200/404 status, non-JSON body.
            A 404 is **not** an error — PDBe returns 404 for entries it
            doesn't recognise, and the caller should treat that as "no
            Rfam mappings" rather than a hard failure.
    """
    url = PDBE_RFAM_URL_TEMPLATE.format(pdb_id=pdb_id.lower())
    own_client = client is None
    http: httpx.Client = client if client is not None else httpx.Client(timeout=timeout)
    try:
        try:
            response = http.get(url)
        except httpx.HTTPError as exc:
            raise ApiRequestError(f"PDBe Rfam request failed for {pdb_id}: {exc}") from exc
        if response.status_code == 404:
            # PDBe answers 404 when it has no nucleic-acid mappings for the
            # entry — treat as "no Rfam data" rather than failing the whole
            # annotation pipeline.
            logger.debug("PDBe returned 404 (no Rfam mappings) for %s", pdb_id)
            return {}
        if response.status_code != 200:
            raise ApiRequestError(f"PDBe Rfam returned HTTP {response.status_code} for {pdb_id}")
        try:
            payload: Any = response.json()
        except (json.JSONDecodeError, ValueError) as exc:
            raise ApiRequestError(f"PDBe Rfam returned non-JSON body for {pdb_id}: {exc}") from exc
        if not isinstance(payload, dict):
            raise ApiRequestError(f"PDBe Rfam returned non-object JSON for {pdb_id}")
        return parse_rfam_mappings(payload, pdb_id=pdb_id)
    finally:
        if own_client:
            http.close()


# ---------------------------------------------------------------------------
# Pure parser
# ---------------------------------------------------------------------------


def parse_rfam_mappings(payload: dict[str, Any], *, pdb_id: str) -> dict[str, list[str]]:
    """Flatten the PDBe response into ``{chain_id: [Rfam_accession, ...]}``.

    The response keys at the top level by ``pdb_id_lower``; we look up
    via ``pdb_id.lower()``. Tolerant of:

    - missing top-level ``pdb_id`` key → empty mapping;
    - missing ``"Rfam"`` sub-key → empty mapping;
    - malformed individual entries → silently skipped.

    Multiple Rfam accessions on the same chain (rare; e.g. tRNA + 5S
    clan) are collected in input order, deduplicated.
    """
    pdb_block = payload.get(pdb_id.lower())
    if not isinstance(pdb_block, dict):
        return {}
    rfam_block = pdb_block.get("Rfam")
    if not isinstance(rfam_block, dict):
        return {}

    out: dict[str, list[str]] = {}
    for accession, body in rfam_block.items():
        if not isinstance(accession, str) or not isinstance(body, dict):
            continue
        mappings = body.get("mappings")
        if not isinstance(mappings, list):
            continue
        for mapping in mappings:
            if not isinstance(mapping, dict):
                continue
            chain_id = mapping.get("chain_id")
            if not isinstance(chain_id, str) or not chain_id:
                continue
            accessions_for_chain = out.setdefault(chain_id, [])
            if accession not in accessions_for_chain:
                accessions_for_chain.append(accession)
    return out


__all__ = [
    "PDBE_RFAM_URL_TEMPLATE",
    "fetch_rfam_mappings",
    "parse_rfam_mappings",
]
