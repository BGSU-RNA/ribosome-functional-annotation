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
