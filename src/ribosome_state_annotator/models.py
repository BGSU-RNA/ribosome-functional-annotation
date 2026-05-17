"""Pydantic v2 data models (spec §9.1, §13.3).

These types form the in-memory contract between the package's stages
(RCSB parsing → classification → correspondence → contact detection →
inference → output). Field names and structure deliberately match the
spec; do not rename or reorder fields without updating §9.1.

Design notes:

- The canonical rRNA chain shape is **list-based** (``ssu_main_rrna_chains``,
  ``lsu_main_rrna_chains``, ``lsu_associated_rrna_chains``) because the
  number of rRNA chains is not fixed across bacterial / eukaryotic
  cytoplasmic / organellar systems (spec §6.2). The single-chain
  ``ssu_chain`` / ``lsu_chain`` aliases on :class:`RibosomeAnnotation` are
  derived convenience properties that resolve to ``None`` for anything
  other than a singleton list.
- ``ChainRef`` and ``LigandRef`` are intentionally mutable (no
  ``model_config = ConfigDict(frozen=True)``) because the RCSB parsing
  layer builds them incrementally — e.g. ``is_ribosomal_protein`` is set
  after the polymer entity has been parsed. ``RibosomeAnnotation`` is
  similarly mutable so the inference layer can append warnings as it goes.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, Field, computed_field

# Status / classification literal types are re-exported so callers don't
# have to import the verbatim string list from §8 / §9 themselves.
RibosomeStatus = Literal["annotated", "skipped", "failed"]
RibosomeClassification = Literal[
    "bacterial_ribosome",
    "eukaryotic_ribosome",
    "eukaryotic_organellar_ribosome",
]


class ChainRef(BaseModel):
    """A single polymer chain in a biological assembly.

    Chain identifiers are author-style (``auth_asym_id``) because BGSU
    unit IDs and downstream consumers use author-style chain labels. The
    ``ife`` computed field returns the BGSU-style IFE string
    ``<pdb_id>|1|<auth_asym_id>`` used in CSV outputs (spec §15.3). The
    model number is hard-coded to ``1`` because v1 reads only the primary
    model of biological assembly mmCIF files (spec §10.4).
    """

    pdb_id: str
    assembly_id: str
    auth_asym_id: str
    label_asym_id: str | None = None
    entity_id: str | None = None
    polymer_type: str | None = None  # "RNA" | "Protein" | "DNA" | ...
    description: str | None = None
    rfam_accessions: list[str] = Field(default_factory=list)
    tax_id: int | None = None
    scientific_name: str | None = None
    superkingdom: str | None = None
    # First non-null UniProt protein name (when available). Required for the
    # §8.4 voting-eligibility rule, which checks both ``description`` and
    # ``uniprot_name`` for the "ribosomal protein" substring.
    uniprot_name: str | None = None
    is_ribosomal_protein: bool = False

    @computed_field  # type: ignore[prop-decorator]
    @property
    def ife(self) -> str:
        """BGSU-style IFE chain identifier (e.g. ``"5J7L|1|AA"``)."""
        return f"{self.pdb_id}|1|{self.auth_asym_id}"


class LigandRef(BaseModel):
    """A bound non-polymer ligand in a biological assembly (spec §13.3).

    Deduplicated by ``comp_id``: one entry per unique chemical component,
    keeping any of its observed ``auth_asym_id`` occurrences.
    """

    comp_id: str
    name: str | None = None
    auth_asym_id: str | None = None
    drugbank_id: str | None = None
    drugbank_description: str | None = None


class AssemblyContext(BaseModel):
    """All chains and ligands belonging to one biological assembly.

    This is the working context the inference layer operates on. RNA and
    protein chains are kept separate to make downstream "candidate RNA
    pool" / "candidate protein factor pool" logic explicit (spec §11, §12.4).
    """

    pdb_id: str
    assembly_id: str
    experimental_methods: list[str] = Field(default_factory=list)
    rna_chains: list[ChainRef] = Field(default_factory=list)
    protein_chains: list[ChainRef] = Field(default_factory=list)
    ligands: list[LigandRef] = Field(default_factory=list)
    coordinate_path: Path | None = None


class LargeScaleMovements(BaseModel):
    """RADdb-derived large-scale ribosome motion metrics for one assembly.

    Always emitted (with ``rad_date=None`` and null metrics) when RADdb
    is unavailable, so the JSON schema stays stable across runs.
    """

    source: Literal["RADdb"] = "RADdb"
    rad_date: str | None = None
    intersubunit_rotation: float | None = None
    ssu_head_rotation: float | None = None


class CorrespondenceResult(BaseModel):
    """Output of the BGSU correspondence layer for one functional-site key
    (e.g. ``"ssu_atrna"``).

    ``mapped_units`` is the flat list of query-PDB unit IDs after PDB-prefix
    and assembly-chain filtering (spec §5.2.2). ``mapped_units_by_chain``
    groups those units by their author chain ID so the Gemmi neighbour-
    search layer can iterate by chain. ``warnings`` collects per-anchor
    issues (e.g. ``"correspondence_missing_for_..."``) for surfacing in
    the final annotation result.
    """

    reference_key: str
    reference_units: list[str] = Field(default_factory=list)
    mapped_units: list[str] = Field(default_factory=list)
    mapped_units_by_chain: dict[str, list[str]] = Field(default_factory=dict)
    warnings: list[str] = Field(default_factory=list)


class RibosomeAnnotation(BaseModel):
    """Final annotation result for one biological assembly (spec §9.1).

    The rRNA chain fields are canonical role-based **lists** (§6.2). The
    single-chain ``ssu_chain`` / ``lsu_chain`` aliases are convenience
    computed properties for the common one-chain-per-role case; they
    return ``None`` if the corresponding list has zero or 2+ entries.
    Consumers that need to handle organellar or fragmented rRNA
    assemblies must read the canonical list fields directly.

    ``assembly_id`` is ``None`` only for entry-level skip results that
    apply before any assembly is considered (e.g. NMR exclusion, §7.3).
    Assembly-level skips still carry the assembly_id.
    """

    pdb_id: str
    assembly_id: str | None = None
    status: RibosomeStatus
    skip_reason: str | None = None
    ribosome_classification: RibosomeClassification | None = None

    # Canonical role-based rRNA outputs (§6.2). Always lists; may be empty.
    ssu_main_rrna_chains: list[ChainRef] = Field(default_factory=list)
    lsu_main_rrna_chains: list[ChainRef] = Field(default_factory=list)
    lsu_associated_rrna_chains: list[ChainRef] = Field(default_factory=list)  # 5S, 5.8S
    other_rna_chains: list[ChainRef] = Field(default_factory=list)

    # Functional chains assigned by contact-transfer (§11).
    mrna_chain: ChainRef | None = None
    aminoacyl_trna_chain: ChainRef | None = None
    peptidyl_trna_chain: ChainRef | None = None
    exit_trna_chain: ChainRef | None = None

    # tRNA functional states (§12). Format: "<SSU>/<LSU>", e.g. "A/A", "ap/AP",
    # "A/Elongation factor Tu", "A/*", "A/**".
    aminoacyl_trna_state: str | None = None
    peptidyl_trna_state: str | None = None
    exit_trna_state: str | None = None

    non_ribosomal_proteins: list[ChainRef] = Field(default_factory=list)
    bound_ligands: list[LigandRef] = Field(default_factory=list)
    # RADdb-derived inter-subunit rotation + SSU head rotation (raddb spec).
    # Always emitted (with null metrics) for annotated assemblies so the
    # output schema stays stable whether RADdb is reachable or not.
    large_scale_movements: LargeScaleMovements | None = None
    classification_evidence: dict[str, Any] = Field(default_factory=dict)
    warnings: list[str] = Field(default_factory=list)

    @computed_field  # type: ignore[prop-decorator]
    @property
    def ssu_chain(self) -> ChainRef | None:
        """The single SSU main rRNA chain if and only if exactly one is present."""
        return self.ssu_main_rrna_chains[0] if len(self.ssu_main_rrna_chains) == 1 else None

    @computed_field  # type: ignore[prop-decorator]
    @property
    def lsu_chain(self) -> ChainRef | None:
        """The single LSU main rRNA chain if and only if exactly one is present."""
        return self.lsu_main_rrna_chains[0] if len(self.lsu_main_rrna_chains) == 1 else None
