"""JSON / JSONL / CSV output emitters (spec §15).

This module is pure formatting. It consumes :class:`RibosomeAnnotation`
objects (the canonical in-memory shape) and emits the four output
artifacts v1 promises:

- **JSON**: one assembly per object; single PDB → array of objects;
  batch → array of objects across PDBs.
- **JSONL**: one JSON object per line. Suitable for streaming and
  for large batches.
- **Chain-level CSV** (``ribosome_chain_annotation.csv``): the legacy
  13-column shape from the prototype, derived from the canonical
  role-based JSON fields per the §15.3 mapping. One row per **annotated**
  assembly; skipped/failed assemblies are omitted.
- **Assembly-level CSV** (``ribosome_assembly_annotation.csv``): one row
  per ``(pdb_id, assembly_id, chain, property, value)`` tuple. Includes
  ``species_name``, ``non_ribosomal_proteins``, ``bound_ligands``,
  ``unmapped_rna_chains``, plus the v1-only ``ribosome_classification``,
  ``dominant_superkingdom``, and ``warning`` rows appended at the end so
  the byte-stable prefix matches the prototype.

CSV line endings use Python's csv default ``\\r\\n`` (matches the
prototype and spec §15.3). Quoting uses ``csv.QUOTE_MINIMAL`` (only
values containing a delimiter, quote, or newline get quoted) which
again matches the prototype.
"""

from __future__ import annotations

import csv
import io
import json
import logging
from collections.abc import Iterable
from pathlib import Path

from ribosome_state_annotator import constants as C
from ribosome_state_annotator.models import (
    ChainRef,
    LigandRef,
    RibosomeAnnotation,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# JSON / JSONL
# ---------------------------------------------------------------------------


def render_json(annotations: list[RibosomeAnnotation], *, indent: int | None = 2) -> str:
    """Render a list of annotations as a single JSON array string."""
    payload = [a.model_dump(mode="json") for a in annotations]
    return json.dumps(payload, indent=indent)


def write_json(
    annotations: list[RibosomeAnnotation], path: Path, *, indent: int | None = 2
) -> None:
    """Write a list of annotations as a JSON array file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(render_json(annotations, indent=indent) + "\n", encoding="utf-8")


def render_jsonl(annotations: Iterable[RibosomeAnnotation]) -> str:
    """Render annotations as JSON-Lines (one object per line)."""
    lines = [json.dumps(a.model_dump(mode="json")) for a in annotations]
    return "\n".join(lines) + "\n" if lines else ""


def write_jsonl(annotations: Iterable[RibosomeAnnotation], path: Path) -> None:
    """Write annotations as a JSON-Lines file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(render_jsonl(annotations), encoding="utf-8")


# ---------------------------------------------------------------------------
# Chain-level CSV (spec §15.3)
# ---------------------------------------------------------------------------

CHAIN_CSV_FIELDS: tuple[str, ...] = (
    "pdb_id",
    "assembly_id",
    "ssu_chain",
    "lsu_large_chain",
    "lsu_medium_chain",
    "lsu_small_chain",
    "mrna",
    "aminoacyl_trna",
    "peptidyl_trna",
    "exit_trna",
    "aminoacyl_trna_state",
    "peptidyl_trna_state",
    "exit_trna_state",
)
"""Legacy chain-level CSV columns in the prototype's exact order. Do not
reorder — :file:`tests/fixtures/legacy_csv/ribosome_chain_annotation.csv`
is the v1 contract."""


def _ife_or_empty(chain: ChainRef | None) -> str:
    return chain.ife if chain is not None else ""


def _singleton_or_empty(chains: list[ChainRef]) -> str:
    """Spec §15.3: legacy single-chain CSV columns are populated only when
    exactly one rRNA chain fills that role; otherwise empty, with a warning
    recorded upstream."""
    return chains[0].ife if len(chains) == 1 else ""


def _associated_chain_by_rfam(chains: list[ChainRef], rfam_accession: str) -> str:
    """Pick the chain whose ``rfam_accessions`` contains ``rfam_accession``,
    or empty if none. Used to fill ``lsu_medium_chain`` (5.8S, ``RF00002``)
    and ``lsu_small_chain`` (5S, ``RF00001``) from
    ``lsu_associated_rrna_chains``."""
    for chain in chains:
        if rfam_accession in chain.rfam_accessions:
            return chain.ife
    return ""


def chain_csv_row(annotation: RibosomeAnnotation) -> dict[str, str]:
    """Build the chain-level CSV row for one annotation.

    Per §15.3 only ``status="annotated"`` annotations should be rendered;
    this helper does NOT enforce that (the caller filters). Skipped/failed
    annotations passed in produce a row with mostly-empty values.
    """
    return {
        "pdb_id": annotation.pdb_id,
        "assembly_id": annotation.assembly_id or "",
        "ssu_chain": _singleton_or_empty(annotation.ssu_main_rrna_chains),
        "lsu_large_chain": _singleton_or_empty(annotation.lsu_main_rrna_chains),
        "lsu_medium_chain": _associated_chain_by_rfam(
            annotation.lsu_associated_rrna_chains, C.LSU_MEDIUM_RFAM
        ),
        "lsu_small_chain": _associated_chain_by_rfam(
            annotation.lsu_associated_rrna_chains, C.LSU_SMALL_RFAM
        ),
        "mrna": _ife_or_empty(annotation.mrna_chain),
        "aminoacyl_trna": _ife_or_empty(annotation.aminoacyl_trna_chain),
        "peptidyl_trna": _ife_or_empty(annotation.peptidyl_trna_chain),
        "exit_trna": _ife_or_empty(annotation.exit_trna_chain),
        "aminoacyl_trna_state": annotation.aminoacyl_trna_state or "",
        "peptidyl_trna_state": annotation.peptidyl_trna_state or "",
        "exit_trna_state": annotation.exit_trna_state or "",
    }


def render_chain_csv(annotations: Iterable[RibosomeAnnotation]) -> str:
    """Render the chain-level CSV (skipped/failed annotations omitted)."""
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=list(CHAIN_CSV_FIELDS))
    writer.writeheader()
    for annotation in annotations:
        if annotation.status != "annotated":
            continue
        writer.writerow(chain_csv_row(annotation))
    return buffer.getvalue()


def write_chain_csv(annotations: Iterable[RibosomeAnnotation], path: Path) -> None:
    """Write the chain-level CSV file (skipped/failed annotations omitted)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(render_chain_csv(annotations), encoding="utf-8", newline="")


# ---------------------------------------------------------------------------
# Assembly-level CSV (spec §15.3)
# ---------------------------------------------------------------------------

ASSEMBLY_CSV_FIELDS: tuple[str, ...] = (
    "pdb_id",
    "assembly_id",
    "chain",
    "property",
    "value",
)
"""Assembly-level CSV columns. Always five fields per the prototype layout."""


def _species_name_value(annotation: RibosomeAnnotation) -> str:
    """Per spec §15.3.1: ``scientific_name`` of the **first SSU main rRNA chain**."""
    if not annotation.ssu_main_rrna_chains:
        return ""
    return annotation.ssu_main_rrna_chains[0].scientific_name or ""


def assembly_csv_rows(annotation: RibosomeAnnotation) -> list[dict[str, str]]:
    """Build all assembly-level CSV rows for one annotation.

    Row order matches the prototype: ``species_name`` first, then per-chain
    ``non_ribosomal_proteins`` rows, then ``bound_ligands``, then
    ``unmapped_rna_chains``. The v1-only ``ribosome_classification``,
    ``dominant_superkingdom``, and ``warning`` rows are appended at the
    end so prefix bytes match the legacy fixture.
    """
    pdb_id = annotation.pdb_id
    assembly_id = annotation.assembly_id or ""

    def _row(chain: str, prop: str, value: str) -> dict[str, str]:
        return {
            "pdb_id": pdb_id,
            "assembly_id": assembly_id,
            "chain": chain,
            "property": prop,
            "value": value,
        }

    rows: list[dict[str, str]] = []
    species = _species_name_value(annotation)
    if species:
        rows.append(_row("", "species_name", species))

    for protein in annotation.non_ribosomal_proteins:
        rows.append(_row(protein.ife, "non_ribosomal_proteins", protein.description or ""))

    # Bound ligands: one row per unique ligand name; chain stays empty per
    # the prototype example.
    seen_ligand_names: set[str] = set()
    for ligand in annotation.bound_ligands:
        ligand_name = _ligand_display_name(ligand)
        if ligand_name and ligand_name not in seen_ligand_names:
            seen_ligand_names.add(ligand_name)
            rows.append(_row("", "bound_ligands", ligand_name))

    for chain in annotation.other_rna_chains:
        rows.append(_row(chain.ife, "unmapped_rna_chains", chain.description or ""))

    # v1 extension rows (spec §15.3 explicitly permits these at the end).
    if annotation.ribosome_classification is not None:
        rows.append(_row("", "ribosome_classification", annotation.ribosome_classification))
    dominant = annotation.classification_evidence.get("dominant_ribosomal_protein_superkingdom")
    if isinstance(dominant, str) and dominant:
        rows.append(_row("", "dominant_superkingdom", dominant))
    for warning in annotation.warnings:
        rows.append(_row("", "warning", warning))

    return rows


def _ligand_display_name(ligand: LigandRef) -> str:
    """Pick a non-empty display value: prefer ``name``, fall back to ``comp_id``."""
    return ligand.name or ligand.comp_id


def render_assembly_csv(annotations: Iterable[RibosomeAnnotation]) -> str:
    """Render the assembly-level CSV (skipped/failed annotations omitted)."""
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=list(ASSEMBLY_CSV_FIELDS))
    writer.writeheader()
    for annotation in annotations:
        if annotation.status != "annotated":
            continue
        for row in assembly_csv_rows(annotation):
            writer.writerow(row)
    return buffer.getvalue()


def write_assembly_csv(annotations: Iterable[RibosomeAnnotation], path: Path) -> None:
    """Write the assembly-level CSV file (skipped/failed annotations omitted)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(render_assembly_csv(annotations), encoding="utf-8", newline="")


__all__ = [
    "ASSEMBLY_CSV_FIELDS",
    "CHAIN_CSV_FIELDS",
    "assembly_csv_rows",
    "chain_csv_row",
    "render_assembly_csv",
    "render_chain_csv",
    "render_json",
    "render_jsonl",
    "write_assembly_csv",
    "write_chain_csv",
    "write_json",
    "write_jsonl",
]
