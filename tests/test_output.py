"""Unit tests for the JSON / JSONL / CSV emitters (spec §15)."""

from __future__ import annotations

import json
from pathlib import Path

from ribosome_state_annotator import output
from ribosome_state_annotator.models import (
    ChainRef,
    LigandRef,
    RibosomeAnnotation,
)

# ---------------------------------------------------------------------------
# Fixture builders
# ---------------------------------------------------------------------------


def _chain(
    pdb_id: str,
    auth: str,
    *,
    rfam: tuple[str, ...] = (),
    description: str | None = None,
    scientific_name: str | None = None,
    is_ribosomal_protein: bool = False,
) -> ChainRef:
    return ChainRef(
        pdb_id=pdb_id,
        assembly_id="1",
        auth_asym_id=auth,
        rfam_accessions=list(rfam),
        description=description,
        scientific_name=scientific_name,
        is_ribosomal_protein=is_ribosomal_protein,
    )


def _bacterial_annotated(pdb_id: str = "5J7L") -> RibosomeAnnotation:
    """A fully-populated bacterial annotation for chain-CSV tests."""
    ssu = _chain(pdb_id, "AA", rfam=("RF00177",), scientific_name="Escherichia coli")
    lsu = _chain(pdb_id, "DA", rfam=("RF02541",))
    five_s = _chain(pdb_id, "BA", rfam=("RF00001",))
    mrna = _chain(pdb_id, "X", description="mRNA")
    atrna = _chain(pdb_id, "V", description="tRNA-Phe")
    ptrna = _chain(pdb_id, "W", description="tRNA-Lys")
    etrna = _chain(pdb_id, "Y", description="tRNA-Ala")
    return RibosomeAnnotation(
        pdb_id=pdb_id,
        assembly_id="1",
        status="annotated",
        ribosome_classification="bacterial_ribosome",
        ssu_main_rrna_chains=[ssu],
        lsu_main_rrna_chains=[lsu],
        lsu_associated_rrna_chains=[five_s],
        mrna_chain=mrna,
        aminoacyl_trna_chain=atrna,
        peptidyl_trna_chain=ptrna,
        exit_trna_chain=etrna,
        aminoacyl_trna_state="A/A",
        peptidyl_trna_state="P/P",
        exit_trna_state="E/E",
        non_ribosomal_proteins=[
            _chain(pdb_id, "T", description="Elongation factor Tu"),
        ],
        bound_ligands=[
            LigandRef(comp_id="MG", name="MAGNESIUM ION"),
            LigandRef(comp_id="STR", name="STREPTOMYCIN"),
        ],
        classification_evidence={
            "rrna_core": "bacterial_like",
            "dominant_ribosomal_protein_superkingdom": "Bacteria",
            "rule": "bacterial_like_rfam_plus_bacterial_proteins",
        },
        warnings=["organellar_classification_without_keyword_support"],
    )


def _skipped(
    pdb_id: str = "XXXX", skip_reason: str = "nmr_structure_not_supported"
) -> RibosomeAnnotation:
    return RibosomeAnnotation(
        pdb_id=pdb_id,
        assembly_id=None,
        status="skipped",
        skip_reason=skip_reason,
    )


# ---------------------------------------------------------------------------
# JSON / JSONL
# ---------------------------------------------------------------------------


def test_render_json_returns_array_of_objects() -> None:
    annotations = [_bacterial_annotated(), _skipped()]
    rendered = output.render_json(annotations)
    parsed = json.loads(rendered)
    assert isinstance(parsed, list)
    assert len(parsed) == 2
    assert parsed[0]["pdb_id"] == "5J7L"
    assert parsed[0]["status"] == "annotated"
    assert parsed[1]["status"] == "skipped"
    assert parsed[1]["assembly_id"] is None  # null preserved


def test_render_json_includes_computed_ssu_chain_alias() -> None:
    ann = _bacterial_annotated()
    parsed = json.loads(output.render_json([ann]))[0]
    assert parsed["ssu_chain"]["ife"] == "5J7L|1|AA"
    assert parsed["lsu_chain"]["ife"] == "5J7L|1|DA"


def test_render_json_skipped_has_null_chain_aliases() -> None:
    parsed = json.loads(output.render_json([_skipped()]))[0]
    assert parsed["ssu_chain"] is None
    assert parsed["lsu_chain"] is None
    assert parsed["aminoacyl_trna_chain"] is None


def test_write_json_creates_parent_dirs(tmp_path: Path) -> None:
    out = tmp_path / "nested" / "deep" / "out.json"
    output.write_json([_bacterial_annotated()], out)
    assert out.exists()
    parsed = json.loads(out.read_text())
    assert parsed[0]["pdb_id"] == "5J7L"


def test_render_jsonl_one_object_per_line() -> None:
    rendered = output.render_jsonl([_bacterial_annotated(), _skipped()])
    lines = [line for line in rendered.splitlines() if line]
    assert len(lines) == 2
    for line in lines:
        parsed = json.loads(line)
        assert isinstance(parsed, dict)


def test_render_jsonl_empty_returns_empty_string() -> None:
    assert output.render_jsonl([]) == ""


# ---------------------------------------------------------------------------
# Chain-level CSV
# ---------------------------------------------------------------------------


def test_chain_csv_header_matches_spec_15_3_order() -> None:
    rendered = output.render_chain_csv([_bacterial_annotated()])
    header_line = rendered.splitlines()[0]
    expected = "pdb_id,assembly_id,ssu_chain,lsu_large_chain,lsu_medium_chain,lsu_small_chain,mrna,aminoacyl_trna,peptidyl_trna,exit_trna,aminoacyl_trna_state,peptidyl_trna_state,exit_trna_state"
    assert header_line == expected


def test_chain_csv_row_has_ife_strings_for_bacterial() -> None:
    row = output.chain_csv_row(_bacterial_annotated())
    assert row["pdb_id"] == "5J7L"
    assert row["assembly_id"] == "1"
    assert row["ssu_chain"] == "5J7L|1|AA"
    assert row["lsu_large_chain"] == "5J7L|1|DA"
    assert row["lsu_medium_chain"] == ""  # bacterial: no 5.8S
    assert row["lsu_small_chain"] == "5J7L|1|BA"  # 5S
    assert row["mrna"] == "5J7L|1|X"
    assert row["aminoacyl_trna"] == "5J7L|1|V"
    assert row["peptidyl_trna"] == "5J7L|1|W"
    assert row["exit_trna"] == "5J7L|1|Y"
    assert row["aminoacyl_trna_state"] == "A/A"
    assert row["peptidyl_trna_state"] == "P/P"
    assert row["exit_trna_state"] == "E/E"


def test_chain_csv_row_eukaryotic_has_5_8s_in_lsu_medium() -> None:
    """Eukaryotic cytoplasmic ribosomes have both 5.8S (RF00002) and 5S (RF00001)."""
    pdb_id = "5TBW"
    ssu = _chain(pdb_id, "A", rfam=("RF01960",))
    lsu = _chain(pdb_id, "1", rfam=("RF02543",))
    five_eight_s = _chain(pdb_id, "B", rfam=("RF00002",))
    five_s = _chain(pdb_id, "C", rfam=("RF00001",))
    ann = RibosomeAnnotation(
        pdb_id=pdb_id,
        assembly_id="1",
        status="annotated",
        ribosome_classification="eukaryotic_ribosome",
        ssu_main_rrna_chains=[ssu],
        lsu_main_rrna_chains=[lsu],
        lsu_associated_rrna_chains=[five_eight_s, five_s],
    )
    row = output.chain_csv_row(ann)
    assert row["lsu_medium_chain"] == "5TBW|1|B"
    assert row["lsu_small_chain"] == "5TBW|1|C"


def test_chain_csv_row_multiple_ssu_leaves_singleton_columns_empty() -> None:
    pdb_id = "FAKE"
    ssu1 = _chain(pdb_id, "AA", rfam=("RF00177",))
    ssu2 = _chain(pdb_id, "AB", rfam=("RF00177",))
    lsu = _chain(pdb_id, "DA", rfam=("RF02541",))
    ann = RibosomeAnnotation(
        pdb_id=pdb_id,
        assembly_id="1",
        status="annotated",
        ssu_main_rrna_chains=[ssu1, ssu2],
        lsu_main_rrna_chains=[lsu],
    )
    row = output.chain_csv_row(ann)
    assert row["ssu_chain"] == ""  # multiple SSU → empty
    assert row["lsu_large_chain"] == "FAKE|1|DA"  # single LSU still fills


def test_chain_csv_uses_crlf_line_endings() -> None:
    rendered = output.render_chain_csv([_bacterial_annotated()])
    assert "\r\n" in rendered


def test_chain_csv_excludes_skipped_annotations() -> None:
    rendered = output.render_chain_csv([_bacterial_annotated(), _skipped()])
    lines = rendered.splitlines()
    # Header + 1 data row only.
    assert len(lines) == 2
    assert "XXXX" not in rendered


def test_chain_csv_includes_assemblies_without_trna_assignments() -> None:
    """Spec §15.3: SSU+LSU-only annotations (no mRNA/tRNAs) are still rows."""
    pdb_id = "PART"
    ssu = _chain(pdb_id, "AA", rfam=("RF00177",))
    lsu = _chain(pdb_id, "DA", rfam=("RF02541",))
    ann = RibosomeAnnotation(
        pdb_id=pdb_id,
        assembly_id="1",
        status="annotated",
        ssu_main_rrna_chains=[ssu],
        lsu_main_rrna_chains=[lsu],
    )
    rendered = output.render_chain_csv([ann])
    assert "PART,1,PART|1|AA,PART|1|DA," in rendered


def test_write_chain_csv_writes_to_disk(tmp_path: Path) -> None:
    out = tmp_path / "chain.csv"
    output.write_chain_csv([_bacterial_annotated()], out)
    assert out.exists()
    text = out.read_text()
    assert "5J7L|1|AA" in text


# ---------------------------------------------------------------------------
# Assembly-level CSV
# ---------------------------------------------------------------------------


def test_assembly_csv_header() -> None:
    rendered = output.render_assembly_csv([_bacterial_annotated()])
    assert rendered.splitlines()[0] == "pdb_id,assembly_id,chain,property,value"


def test_assembly_csv_emits_species_name() -> None:
    rows = output.assembly_csv_rows(_bacterial_annotated())
    species_rows = [r for r in rows if r["property"] == "species_name"]
    assert len(species_rows) == 1
    assert species_rows[0]["chain"] == ""
    assert species_rows[0]["value"] == "Escherichia coli"


def test_assembly_csv_emits_non_ribosomal_proteins_per_chain() -> None:
    rows = output.assembly_csv_rows(_bacterial_annotated())
    factor_rows = [r for r in rows if r["property"] == "non_ribosomal_proteins"]
    assert len(factor_rows) == 1
    assert factor_rows[0]["chain"] == "5J7L|1|T"
    assert factor_rows[0]["value"] == "Elongation factor Tu"


def test_assembly_csv_emits_bound_ligands_unique() -> None:
    """One row per unique ligand name; chain stays empty per spec example."""
    rows = output.assembly_csv_rows(_bacterial_annotated())
    ligand_rows = [r for r in rows if r["property"] == "bound_ligands"]
    values = {r["value"] for r in ligand_rows}
    assert values == {"MAGNESIUM ION", "STREPTOMYCIN"}
    for r in ligand_rows:
        assert r["chain"] == ""


def test_assembly_csv_emits_unmapped_rna_chains() -> None:
    pdb_id = "FAKE"
    ssu = _chain(pdb_id, "AA", rfam=("RF00177",), scientific_name="Escherichia coli")
    lsu = _chain(pdb_id, "DA", rfam=("RF02541",))
    leftover = _chain(pdb_id, "Z", description="tRNA-Phe")
    ann = RibosomeAnnotation(
        pdb_id=pdb_id,
        assembly_id="1",
        status="annotated",
        ssu_main_rrna_chains=[ssu],
        lsu_main_rrna_chains=[lsu],
        other_rna_chains=[leftover],
    )
    rows = output.assembly_csv_rows(ann)
    unmapped_rows = [r for r in rows if r["property"] == "unmapped_rna_chains"]
    assert len(unmapped_rows) == 1
    assert unmapped_rows[0]["chain"] == "FAKE|1|Z"
    assert unmapped_rows[0]["value"] == "tRNA-Phe"


def test_assembly_csv_v1_extension_rows_appended() -> None:
    """v1 adds ribosome_classification, dominant_superkingdom, and warning
    rows at the end (per §15.3 to preserve byte-stable prefix vs prototype)."""
    rows = output.assembly_csv_rows(_bacterial_annotated())
    properties = [r["property"] for r in rows]
    # Last few rows are the v1-extension rows in this order.
    assert "ribosome_classification" in properties
    assert "dominant_superkingdom" in properties
    assert "warning" in properties
    # Verify they come after the legacy properties.
    last_legacy = max(
        i
        for i, prop in enumerate(properties)
        if prop
        in {"species_name", "non_ribosomal_proteins", "bound_ligands", "unmapped_rna_chains"}
    )
    first_extension = min(
        i
        for i, prop in enumerate(properties)
        if prop in {"ribosome_classification", "dominant_superkingdom", "warning"}
    )
    assert last_legacy < first_extension


def test_assembly_csv_excludes_skipped_annotations() -> None:
    rendered = output.render_assembly_csv([_bacterial_annotated(), _skipped()])
    assert "XXXX" not in rendered


def test_assembly_csv_skips_species_when_no_ssu_taxonomy() -> None:
    pdb_id = "FAKE"
    ssu = _chain(pdb_id, "AA", rfam=("RF00177",))  # no scientific_name
    lsu = _chain(pdb_id, "DA", rfam=("RF02541",))
    ann = RibosomeAnnotation(
        pdb_id=pdb_id,
        assembly_id="1",
        status="annotated",
        ssu_main_rrna_chains=[ssu],
        lsu_main_rrna_chains=[lsu],
    )
    rows = output.assembly_csv_rows(ann)
    assert all(r["property"] != "species_name" for r in rows)


def test_write_assembly_csv_writes_to_disk(tmp_path: Path) -> None:
    out = tmp_path / "assembly.csv"
    output.write_assembly_csv([_bacterial_annotated()], out)
    assert out.exists()
    text = out.read_text()
    assert "Escherichia coli" in text
