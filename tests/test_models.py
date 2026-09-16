"""Unit tests for Pydantic v2 data models (spec §9.1)."""

from __future__ import annotations

from pathlib import Path

import pytest
from pydantic import ValidationError

from ribosome_state_annotator.models import (
    AssemblyContext,
    ChainRef,
    CorrespondenceResult,
    LigandRef,
    RibosomeAnnotation,
)

# ---------------------------------------------------------------------------
# ChainRef
# ---------------------------------------------------------------------------


def _chain(pdb: str = "5J7L", asym: str = "AA", **kw: object) -> ChainRef:
    return ChainRef(pdb_id=pdb, assembly_id="1", auth_asym_id=asym, **kw)  # type: ignore[arg-type]


def test_chainref_ife_format() -> None:
    chain = _chain("5J7L", "AA")
    assert chain.ife == "5J7L|1|AA"


def test_chainref_ife_uses_auth_asym_not_label_asym() -> None:
    chain = _chain("5J7L", "AA", label_asym_id="AB")
    assert chain.ife == "5J7L|1|AA"


def test_chainref_defaults() -> None:
    chain = _chain()
    assert chain.rfam_accessions == []
    assert chain.is_ribosomal_protein is False
    assert chain.label_asym_id is None
    assert chain.entity_id is None
    assert chain.tax_id is None


def test_chainref_default_lists_are_per_instance() -> None:
    """Pydantic v2 must deep-copy mutable defaults — mutating one instance's
    list must not leak into others."""
    a = _chain()
    b = _chain()
    a.rfam_accessions.append("RF00177")
    assert b.rfam_accessions == []


def test_chainref_dump_includes_ife() -> None:
    chain = _chain(rfam_accessions=["RF00177"], superkingdom="Bacteria")
    dump = chain.model_dump()
    assert dump["ife"] == "5J7L|1|AA"
    assert dump["rfam_accessions"] == ["RF00177"]
    assert dump["superkingdom"] == "Bacteria"


def test_chainref_required_fields_validated() -> None:
    with pytest.raises(ValidationError):
        ChainRef.model_validate({"pdb_id": "5J7L", "assembly_id": "1"})  # no auth_asym_id


# ---------------------------------------------------------------------------
# LigandRef
# ---------------------------------------------------------------------------


def test_ligandref_minimal() -> None:
    lig = LigandRef(comp_id="MG")
    assert lig.comp_id == "MG"
    assert lig.name is None


def test_ligandref_full() -> None:
    lig = LigandRef(
        comp_id="STR",
        name="STREPTOMYCIN",
        auth_asym_id="Z",
        drugbank_id="DB01082",
        drugbank_description="Aminoglycoside antibiotic",
    )
    assert lig.drugbank_id == "DB01082"


# ---------------------------------------------------------------------------
# AssemblyContext
# ---------------------------------------------------------------------------


def test_assembly_context_minimal() -> None:
    ctx = AssemblyContext(pdb_id="5J7L", assembly_id="1")
    assert ctx.rna_chains == []
    assert ctx.protein_chains == []
    assert ctx.coordinate_path is None


def test_assembly_context_with_chains_and_path(tmp_path: Path) -> None:
    rna = _chain(asym="AA")
    prot = _chain(asym="P1")
    ctx = AssemblyContext(
        pdb_id="5J7L",
        assembly_id="1",
        rna_chains=[rna],
        protein_chains=[prot],
        coordinate_path=tmp_path / "5j7l-assembly1.cif.gz",
    )
    assert ctx.rna_chains[0].ife == "5J7L|1|AA"
    assert ctx.coordinate_path is not None
    assert ctx.coordinate_path.name == "5j7l-assembly1.cif.gz"


# ---------------------------------------------------------------------------
# CorrespondenceResult
# ---------------------------------------------------------------------------


def test_correspondence_result_minimal() -> None:
    cr = CorrespondenceResult(reference_key="ssu_atrna")
    assert cr.warnings == []
    assert cr.mapped_units == []
    assert cr.mapped_units_by_chain == {}


def test_correspondence_result_round_trip() -> None:
    cr = CorrespondenceResult(
        reference_key="ssu_atrna",
        reference_units=["5J7L|1|AA|G|530"],
        mapped_units=["7K00|1|a|G|530"],
        mapped_units_by_chain={"a": ["7K00|1|a|G|530"]},
        warnings=["correspondence_missing_for_ssu_atrna_5J7L|1|AA|A|1492"],
    )
    again = CorrespondenceResult.model_validate(cr.model_dump())
    assert again == cr


# ---------------------------------------------------------------------------
# RibosomeAnnotation
# ---------------------------------------------------------------------------


def test_annotation_minimal_skip() -> None:
    """The NMR / partial entry-level skip case: assembly_id is None."""
    ann = RibosomeAnnotation(
        pdb_id="1ABC",
        assembly_id=None,
        status="skipped",
        skip_reason="nmr_structure_not_supported",
    )
    assert ann.assembly_id is None
    assert ann.status == "skipped"
    assert ann.ribosome_classification is None
    assert ann.warnings == []


def test_annotation_rejects_unknown_status() -> None:
    with pytest.raises(ValidationError):
        RibosomeAnnotation(pdb_id="5J7L", assembly_id="1", status="ok")  # type: ignore[arg-type]


def test_annotation_rejects_unknown_classification() -> None:
    with pytest.raises(ValidationError):
        RibosomeAnnotation(
            pdb_id="5J7L",
            assembly_id="1",
            status="annotated",
            ribosome_classification="archaeal_ribosome",  # type: ignore[arg-type]
        )


def test_ssu_chain_alias_none_when_empty() -> None:
    ann = RibosomeAnnotation(pdb_id="5J7L", assembly_id="1", status="annotated")
    assert ann.ssu_chain is None
    assert ann.lsu_chain is None


def test_ssu_chain_alias_returns_singleton() -> None:
    ssu = _chain(asym="AA")
    lsu = _chain(asym="DA")
    ann = RibosomeAnnotation(
        pdb_id="5J7L",
        assembly_id="1",
        status="annotated",
        ssu_main_rrna_chains=[ssu],
        lsu_main_rrna_chains=[lsu],
    )
    assert ann.ssu_chain == ssu
    assert ann.lsu_chain == lsu
    assert ann.ssu_chain is not None
    assert ann.ssu_chain.ife == "5J7L|1|AA"


def test_ssu_chain_alias_none_when_multiple() -> None:
    """The list-based form remains canonical when multiple SSU main rRNA chains
    are present; the convenience alias must be ``None`` (§9.1)."""
    ssu_a = _chain(asym="AA")
    ssu_b = _chain(asym="BB")
    ann = RibosomeAnnotation(
        pdb_id="5J7L",
        assembly_id="1",
        status="annotated",
        ssu_main_rrna_chains=[ssu_a, ssu_b],
    )
    assert ann.ssu_chain is None
    assert len(ann.ssu_main_rrna_chains) == 2


def test_annotation_round_trip_preserves_all_fields() -> None:
    ssu = _chain(asym="AA", rfam_accessions=["RF00177"], polymer_type="RNA")
    lsu = _chain(asym="DA", rfam_accessions=["RF02541"], polymer_type="RNA")
    fivev = _chain(asym="BA", rfam_accessions=["RF00001"], polymer_type="RNA")
    mrna = _chain(asym="X", polymer_type="RNA", description="mRNA")
    atrna = _chain(asym="V", polymer_type="RNA", description="tRNA-Phe")
    ann = RibosomeAnnotation(
        pdb_id="5J7L",
        assembly_id="1",
        status="annotated",
        ribosome_classification="bacterial_ribosome",
        ssu_main_rrna_chains=[ssu],
        lsu_main_rrna_chains=[lsu],
        lsu_associated_rrna_chains=[fivev],
        mrna_chain=mrna,
        aminoacyl_trna_chain=atrna,
        aminoacyl_trna_state="A/A",
        bound_ligands=[LigandRef(comp_id="MG", name="MAGNESIUM ION")],
        classification_evidence={
            "ssu_rfam": ["RF00177"],
            "lsu_rfam": ["RF02541"],
            "rrna_core": "bacterial_like",
            "rule": "bacterial_like_rfam_plus_bacterial_proteins",
        },
        warnings=["something_minor_to_report"],
    )
    again = RibosomeAnnotation.model_validate(ann.model_dump())
    assert again == ann
    assert again.ssu_chain is not None
    assert again.ssu_chain.ife == "5J7L|1|AA"


def test_annotation_dump_includes_aliases() -> None:
    ssu = _chain(asym="AA")
    lsu = _chain(asym="DA")
    ann = RibosomeAnnotation(
        pdb_id="5J7L",
        assembly_id="1",
        status="annotated",
        ssu_main_rrna_chains=[ssu],
        lsu_main_rrna_chains=[lsu],
    )
    dump = ann.model_dump()
    assert "ssu_chain" in dump
    assert "lsu_chain" in dump
    assert dump["ssu_chain"]["ife"] == "5J7L|1|AA"
    assert dump["lsu_chain"]["ife"] == "5J7L|1|DA"


# ---------------------------------------------------------------------------
# RibosomeAnnotation.summary()
# ---------------------------------------------------------------------------


def _full_annotation() -> RibosomeAnnotation:
    from ribosome_state_annotator.models import (
        Anticodon,
        AnticodonResidue,
        AssemblyTaxonomy,
        BasePair,
        Codon,
        LargeScaleMovements,
        TaxonNode,
        TRNAmRNAInteraction,
    )

    def chain(asym: str, **kw: object) -> ChainRef:
        return ChainRef(pdb_id="5UYM", assembly_id="1", auth_asym_id=asym, **kw)  # type: ignore[arg-type]

    pair = BasePair(
        codon_position=3,
        trna_position=34,
        codon_unit_id="5UYM|1|V|C|21",
        trna_unit_id="5UYM|1|Y|G|34",
        codon_base="C",
        trna_parent_base="G",
        trna_chem_comp_id="G",
        trna_is_modified=False,
        fr3d_interaction="cWW",
        basepair="C-G",
        is_wobble_position=True,
        assignment_status="assigned",
    )
    interaction = TRNAmRNAInteraction(
        site="A",
        mrna_chain_id="V",
        trna_chain_id="Y",
        codon=Codon(sequence="UUC", assignment_status="complete"),
        anticodon=Anticodon(
            sequence_parent="GAA",
            residues=[
                AnticodonResidue(
                    trna_position=34,
                    unit_id="5UYM|1|Y|G|34",
                    parent_base="G",
                    trna_chem_comp_id="G",
                    is_modified=False,
                )
            ],
        ),
        pairs=[pair],
    )
    return RibosomeAnnotation(
        pdb_id="5UYM",
        assembly_id="1",
        status="annotated",
        ribosome_classification="bacterial_ribosome",
        ssu_main_rrna_chains=[chain("A")],
        lsu_main_rrna_chains=[chain("01")],
        lsu_associated_rrna_chains=[chain("02")],
        mrna_chain=chain("V"),
        aminoacyl_trna_chain=chain("Y"),
        peptidyl_trna_chain=chain("W"),
        aminoacyl_trna_state="A/Elongation factor Tu 1",
        peptidyl_trna_state="P/P",
        non_ribosomal_proteins=[
            chain("Z", description="Elongation factor Tu 2", uniprot_name="Elongation factor Tu 1")
        ],
        assembly_taxonomy=AssemblyTaxonomy(
            lineage=(TaxonNode(tax_id=2, name="Bacteria", depth=1),),
            domain="Bacteria",
            species="Escherichia coli",
        ),
        large_scale_movements=LargeScaleMovements(
            rad_date="20260508", intersubunit_rotation=0.9, ssu_head_rotation=2.5
        ),
        trna_mrna_interactions=[interaction],
        warnings=["something"],
    )


def test_summary_annotated_lists_sites_states_and_evidence() -> None:
    text = _full_annotation().summary()
    lines = text.splitlines()
    assert lines[0] == "5UYM assembly 1: annotated — bacterial_ribosome (complete)"
    assert "organism:  Escherichia coli" in text
    assert "SSU rRNA:  5UYM|1|A" in text
    assert "LSU rRNA:  5UYM|1|01" in text
    assert "5S/5.8S:   5UYM|1|02" in text
    assert "mRNA:      5UYM|1|V" in text
    assert "A-site tRNA: 5UYM|1|Y   state A/Elongation factor Tu 1   codon UUC / anticodon GAA (1 FR3D pair(s))" in text
    assert "P-site tRNA: 5UYM|1|W   state P/P" in text
    assert "E-site tRNA: -" in text
    # UniProt name wins over the depositor description, matching the
    # state label ("A/Elongation factor Tu 1") so one chain never shows
    # under two names in the same block.
    assert "factors:   Elongation factor Tu 1 [Z]" in text
    assert "Tu 2" not in text
    assert "rotation:  intersubunit 0.9°, SSU head 2.5°" in text
    assert "warnings:  1 (see .warnings)" in text
    # No pydantic repr noise (the tester's "TaxonNode" flood).
    assert "TaxonNode" not in text
    assert "ChainRef" not in text


def test_summary_failed_is_one_line_with_reason() -> None:
    ann = RibosomeAnnotation(
        pdb_id="5J7L",
        assembly_id="1",
        status="failed",
        skip_reason="correspondence_failure (lsu): BGSU correspondence request timed out",
    )
    assert ann.summary() == (
        "5J7L assembly 1: failed — correspondence_failure (lsu): "
        "BGSU correspondence request timed out"
    )


def test_summary_skipped_entry_level_uses_dash_for_missing_assembly() -> None:
    ann = RibosomeAnnotation(pdb_id="2N0L", assembly_id=None, status="skipped", skip_reason="nmr")
    assert ann.summary() == "2N0L assembly -: skipped — nmr"


def test_summary_annotated_without_optional_blocks() -> None:
    ann = RibosomeAnnotation(pdb_id="5J7L", assembly_id="1", status="annotated")
    text = ann.summary()
    assert "unclassified" in text
    assert "SSU rRNA:  -" in text
    assert "mRNA:      -" in text
    assert "organism" not in text
    assert "rotation" not in text
    assert "warnings" not in text


def test_summary_factor_falls_back_to_description_without_uniprot_name() -> None:
    ann = RibosomeAnnotation(
        pdb_id="5UYM",
        assembly_id="1",
        status="annotated",
        non_ribosomal_proteins=[
            ChainRef(pdb_id="5UYM", assembly_id="1", auth_asym_id="Z", description="Some factor")
        ],
    )
    assert "factors:   Some factor [Z]" in ann.summary()
