"""RCSB GraphQL Data API client (spec §5.1).

This module owns:

- the GraphQL endpoint URL and the exact ``RibosomeEntry`` query from
  §5.1.1 (do not edit fields out of the query without updating §5.1);
- :func:`fetch_entry_payload`, a thin HTTP wrapper that returns the raw
  ``entry`` dict from RCSB and raises :class:`~.exceptions.ApiRequestError`
  on transport, JSON, or GraphQL errors;
- :func:`parse_assemblies`, a pure function that turns one RCSB entry
  payload into a list of :class:`~.models.AssemblyContext` objects.

The split is deliberate (spec §23.2): the HTTP call has no biological
knowledge, and the parser has no network state. Tests mock the HTTP layer
with :mod:`respx` and exercise the parser against the §5.1.2 example
response verbatim.
"""

from __future__ import annotations

import json
import logging
import re
from typing import Any

import httpx

from ribosome_state_annotator.classify import matches_ribosomal_protein_narrow
from ribosome_state_annotator.exceptions import ApiRequestError
from ribosome_state_annotator.models import AssemblyContext, ChainRef, LigandRef

logger = logging.getLogger(__name__)

RCSB_GRAPHQL_URL = "https://data.rcsb.org/graphql"
"""RCSB GraphQL Data API endpoint (spec §5.1)."""

# Spec §5.1.1 — verbatim. Do not change field paths without updating §5.1.
RIBOSOME_ENTRY_QUERY = """\
query RibosomeEntry($entry_id: String!) {
  entry(entry_id: $entry_id) {
    rcsb_id
    exptl { method }
    assemblies {
      rcsb_id
      polymer_entity_instances {
        rcsb_polymer_entity_instance_container_identifiers {
          auth_asym_id
          asym_id
          entity_id
        }
        polymer_entity {
          rcsb_id
          entity_poly { rcsb_entity_polymer_type }
          rcsb_polymer_entity { pdbx_description }
          rcsb_polymer_entity_annotation {
            annotation_id
            type
            name
            description
          }
          rcsb_entity_source_organism {
            ncbi_taxonomy_id
            ncbi_scientific_name
            ncbi_parent_scientific_name
          }
          uniprots {
            rcsb_uniprot_protein { name { value } }
          }
        }
      }
      nonpolymer_entity_instances {
        rcsb_nonpolymer_entity_instance_container_identifiers {
          auth_asym_id
          comp_id
        }
        nonpolymer_entity {
          pdbx_entity_nonpoly { name }
          rcsb_nonpolymer_entity_annotation {
            type
            annotation_id
            description
          }
        }
      }
    }
  }
}
"""

# Assembly rcsb_id format: "<PDB>-<assembly_n>", e.g. "5J7L-1".
_ASSEMBLY_RCSB_ID_RE = re.compile(r"^[A-Za-z0-9]+-(\w+)$")


# ---------------------------------------------------------------------------
# HTTP
# ---------------------------------------------------------------------------


def fetch_entry_payload(
    pdb_id: str,
    *,
    client: httpx.Client | None = None,
    timeout: float = 30.0,
) -> dict[str, Any]:
    """Fetch the raw ``entry`` payload for one PDB ID from RCSB.

    The package's caching layer (step 6) wraps this function; until then
    callers always pay the network round-trip.

    Args:
        pdb_id: PDB accession. Case-insensitive; uppercased before sending.
        client: Optional pre-built :class:`httpx.Client`. Used by tests so a
            single :mod:`respx` mock can intercept the call. When omitted, a
            client is created and closed inside this function.
        timeout: Request timeout in seconds. Ignored when ``client`` is supplied.

    Returns:
        The ``data.entry`` dict from the GraphQL response.

    Raises:
        ApiRequestError: HTTP failure, non-200 status, non-JSON body,
            GraphQL ``errors`` array, or a ``null`` ``entry`` field.
    """
    body = {
        "query": RIBOSOME_ENTRY_QUERY,
        "variables": {"entry_id": pdb_id.upper()},
    }
    own_client = client is None
    http: httpx.Client = client if client is not None else httpx.Client(timeout=timeout)
    try:
        try:
            response = http.post(RCSB_GRAPHQL_URL, json=body)
        except httpx.HTTPError as exc:
            raise ApiRequestError(f"RCSB GraphQL request failed for {pdb_id}: {exc}") from exc
        if response.status_code != 200:
            raise ApiRequestError(f"RCSB GraphQL returned HTTP {response.status_code} for {pdb_id}")
        try:
            payload: Any = response.json()
        except (json.JSONDecodeError, ValueError) as exc:
            raise ApiRequestError(
                f"RCSB GraphQL returned non-JSON body for {pdb_id}: {exc}"
            ) from exc
        if not isinstance(payload, dict):
            raise ApiRequestError(f"RCSB GraphQL returned non-object JSON for {pdb_id}")
        errors = payload.get("errors")
        if errors:
            raise ApiRequestError(f"RCSB GraphQL returned errors for {pdb_id}: {errors}")
        data = payload.get("data")
        if not isinstance(data, dict):
            raise ApiRequestError(f"RCSB GraphQL response for {pdb_id} has no `data` object")
        entry = data.get("entry")
        if entry is None:
            raise ApiRequestError(
                f"RCSB GraphQL returned no entry for {pdb_id} (PDB ID may not exist)"
            )
        if not isinstance(entry, dict):
            raise ApiRequestError(f"RCSB GraphQL `entry` for {pdb_id} is not an object")
        return entry
    finally:
        if own_client:
            http.close()


# ---------------------------------------------------------------------------
# Pure parser
# ---------------------------------------------------------------------------


def parse_assemblies(entry_payload: dict[str, Any]) -> list[AssemblyContext]:
    """Convert one raw RCSB ``entry`` dict into a list of :class:`AssemblyContext`.

    Pure function — no network, no filesystem. The classification layer
    (step 4) decides which assemblies survive the NMR / archaeal / partial
    filters; this function returns every assembly the entry declares.
    """
    pdb_id_raw = entry_payload.get("rcsb_id")
    if not isinstance(pdb_id_raw, str) or not pdb_id_raw:
        raise ApiRequestError("entry payload missing rcsb_id")
    pdb_id = pdb_id_raw.upper()

    methods: list[str] = []
    for exptl_entry in entry_payload.get("exptl") or []:
        if isinstance(exptl_entry, dict):
            method = exptl_entry.get("method")
            if isinstance(method, str) and method:
                methods.append(method)

    assemblies_payload = entry_payload.get("assemblies") or []
    return [_parse_assembly(pdb_id, methods, a) for a in assemblies_payload if isinstance(a, dict)]


def _parse_assembly(
    pdb_id: str,
    methods: list[str],
    assembly_payload: dict[str, Any],
) -> AssemblyContext:
    rcsb_id = assembly_payload.get("rcsb_id")
    rcsb_id_str = rcsb_id if isinstance(rcsb_id, str) else ""
    assembly_id = _parse_assembly_id(rcsb_id_str, fallback="1")

    rna_chains: list[ChainRef] = []
    protein_chains: list[ChainRef] = []

    for instance in assembly_payload.get("polymer_entity_instances") or []:
        if not isinstance(instance, dict):
            continue
        chain = _parse_polymer_instance(pdb_id, assembly_id, instance)
        if chain is None:
            continue
        polymer_type = (chain.polymer_type or "").upper()
        if polymer_type == "RNA":
            rna_chains.append(chain)
        elif polymer_type == "PROTEIN":
            protein_chains.append(chain)
        else:
            # DNA / NA-hybrid / Other — kept off both lists in v1. Downstream
            # logic (mRNA assignment, factor labelling) does not look at DNA
            # chains; if a real ribosome ever needs DNA-chain handling, route
            # it into a dedicated bucket here.
            logger.debug(
                "skipping non-RNA / non-Protein chain %s (polymer_type=%r)",
                chain.ife,
                chain.polymer_type,
            )

    ligands: list[LigandRef] = []
    seen_comp_ids: set[str] = set()
    for instance in assembly_payload.get("nonpolymer_entity_instances") or []:
        if not isinstance(instance, dict):
            continue
        lig = _parse_nonpolymer_instance(instance)
        if lig is None or lig.comp_id in seen_comp_ids:
            continue
        ligands.append(lig)
        seen_comp_ids.add(lig.comp_id)

    return AssemblyContext(
        pdb_id=pdb_id,
        assembly_id=assembly_id,
        experimental_methods=list(methods),
        rna_chains=rna_chains,
        protein_chains=protein_chains,
        ligands=ligands,
    )


def _parse_polymer_instance(
    pdb_id: str,
    assembly_id: str,
    instance: dict[str, Any],
) -> ChainRef | None:
    ids = instance.get("rcsb_polymer_entity_instance_container_identifiers") or {}
    if not isinstance(ids, dict):
        return None
    auth_asym_id = _str_or_none(ids.get("auth_asym_id"))
    if not auth_asym_id:
        logger.warning("polymer instance in %s/%s missing auth_asym_id", pdb_id, assembly_id)
        return None

    entity = instance.get("polymer_entity")
    if not isinstance(entity, dict):
        entity = {}

    polymer_type = _str_or_none(_get_nested(entity, "entity_poly", "rcsb_entity_polymer_type"))
    description = _str_or_none(_get_nested(entity, "rcsb_polymer_entity", "pdbx_description"))

    rfam_accessions: list[str] = []
    # Live RCSB schema renamed the field from the spec §5.1.1 plural
    # ``rcsb_polymer_entity_annotations`` to the singular
    # ``rcsb_polymer_entity_annotation``; accept either for forwards/
    # backwards compatibility.
    annotations = (
        entity.get("rcsb_polymer_entity_annotation")
        or entity.get("rcsb_polymer_entity_annotations")
        or []
    )
    if isinstance(annotations, list):
        for ann in annotations:
            if isinstance(ann, dict) and ann.get("type") == "Rfam":
                accession = _str_or_none(ann.get("annotation_id"))
                if accession:
                    rfam_accessions.append(accession)
    # Note: the live RCSB schema no longer supplies Rfam annotations for
    # rRNA chains, so ``rfam_accessions`` here is often empty. The api
    # layer augments these via PDBe's `/nucleic_mappings/rfam` endpoint
    # (see :mod:`pdbe_client`) before classification runs.

    sources = entity.get("rcsb_entity_source_organism") or []
    source_list = sources if isinstance(sources, list) else []
    tax_id = _first_non_null_int(source_list, "ncbi_taxonomy_id")
    scientific_name = _first_non_null_str(source_list, "ncbi_scientific_name")
    superkingdom = _first_non_null_str(source_list, "ncbi_parent_scientific_name")

    uniprot_name = _extract_first_uniprot_name(entity.get("uniprots"))
    is_ribosomal = matches_ribosomal_protein_narrow(description, uniprot_name)

    return ChainRef(
        pdb_id=pdb_id,
        assembly_id=assembly_id,
        auth_asym_id=auth_asym_id,
        label_asym_id=_str_or_none(ids.get("asym_id")),
        entity_id=_str_or_none(ids.get("entity_id")),
        polymer_type=polymer_type,
        description=description,
        rfam_accessions=rfam_accessions,
        tax_id=tax_id,
        scientific_name=scientific_name,
        superkingdom=superkingdom,
        uniprot_name=uniprot_name,
        is_ribosomal_protein=is_ribosomal,
    )


def _parse_nonpolymer_instance(instance: dict[str, Any]) -> LigandRef | None:
    ids = instance.get("rcsb_nonpolymer_entity_instance_container_identifiers") or {}
    if not isinstance(ids, dict):
        return None
    comp_id = _str_or_none(ids.get("comp_id"))
    if not comp_id:
        return None

    entity = instance.get("nonpolymer_entity")
    if not isinstance(entity, dict):
        entity = {}

    name = _str_or_none(_get_nested(entity, "pdbx_entity_nonpoly", "name"))

    drugbank_id: str | None = None
    drugbank_description: str | None = None
    annotations = entity.get("rcsb_nonpolymer_entity_annotation") or []
    if isinstance(annotations, list):
        for ann in annotations:
            if isinstance(ann, dict) and ann.get("type") == "DrugBank":
                drugbank_id = _str_or_none(ann.get("annotation_id"))
                drugbank_description = _str_or_none(ann.get("description"))
                break

    return LigandRef(
        comp_id=comp_id,
        name=name,
        auth_asym_id=_str_or_none(ids.get("auth_asym_id")),
        drugbank_id=drugbank_id,
        drugbank_description=drugbank_description,
    )


def _parse_assembly_id(rcsb_id: str, *, fallback: str) -> str:
    """Parse an assembly rcsb_id like ``"5J7L-1"`` into ``"1"``.

    Returns ``fallback`` if the format is unrecognised (e.g. RCSB starts
    issuing differently-shaped IDs in the future).
    """
    match = _ASSEMBLY_RCSB_ID_RE.match(rcsb_id)
    if match is None:
        if rcsb_id:
            logger.warning("unrecognised assembly rcsb_id %r; using fallback %r", rcsb_id, fallback)
        return fallback
    return match.group(1)


# ---------------------------------------------------------------------------
# Local helpers
# ---------------------------------------------------------------------------


def _get_nested(payload: dict[str, Any] | None, *path: str) -> Any:
    """Safely walk a nested dict, returning ``None`` at the first missing key."""
    current: Any = payload
    for key in path:
        if not isinstance(current, dict):
            return None
        current = current.get(key)
    return current


def _str_or_none(value: Any) -> str | None:
    return value if isinstance(value, str) and value else None


def _first_non_null_str(items: list[Any], key: str) -> str | None:
    for item in items:
        if isinstance(item, dict):
            val = item.get(key)
            if isinstance(val, str) and val:
                return val
    return None


def _first_non_null_int(items: list[Any], key: str) -> int | None:
    for item in items:
        if isinstance(item, dict):
            val = item.get(key)
            # ``isinstance(True, int)`` is True; reject bools explicitly so a
            # JSON ``true`` cannot land in tax_id.
            if isinstance(val, int) and not isinstance(val, bool):
                return val
    return None


def _extract_first_uniprot_name(uniprots: Any) -> str | None:
    if not isinstance(uniprots, list):
        return None
    for entry in uniprots:
        name = _str_or_none(_get_nested(entry, "rcsb_uniprot_protein", "name", "value"))
        if name:
            return name
    return None
