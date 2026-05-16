"""Unit tests for chain assignment and tRNA state inference (spec §11, §12)."""

from __future__ import annotations

import gemmi
import pytest

from ribosome_state_annotator import infer
from ribosome_state_annotator.models import (
    AssemblyContext,
    ChainRef,
    CorrespondenceResult,
    LigandRef,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _chain(
    auth: str,
    *,
    polymer_type: str = "RNA",
    rfam: tuple[str, ...] = (),
    description: str | None = None,
    is_ribosomal_protein: bool = False,
) -> ChainRef:
    return ChainRef(
        pdb_id="RIBOFIXTURE",
        assembly_id="1",
        auth_asym_id=auth,
        polymer_type=polymer_type,
        rfam_accessions=list(rfam),
        description=description,
        is_ribosomal_protein=is_ribosomal_protein,
    )


def _assembly(rna_chains: list[ChainRef], protein_chains: list[ChainRef]) -> AssemblyContext:
    return AssemblyContext(
        pdb_id="RIBOFIXTURE",
        assembly_id="1",
        experimental_methods=["X-RAY DIFFRACTION"],
        rna_chains=rna_chains,
        protein_chains=protein_chains,
        ligands=[LigandRef(comp_id="MG", name="MAGNESIUM ION")],
    )


def _correspondence(site_key: str, mapped_units: list[str]) -> CorrespondenceResult:
    return CorrespondenceResult(
        reference_key=site_key,
        reference_units=mapped_units,  # placeholder
        mapped_units=mapped_units,
        warnings=[],
    )


def _build_ribosome_inputs() -> tuple[
    AssemblyContext,
    dict[str, list[ChainRef]],
    dict[str, CorrespondenceResult],
]:
    """Standard inputs for the ribosome_fixture: assembly, by_role partition,
    correspondence_by_site dict with anchors pointing at the fixture's
    pre-computed anchor residues.

    Anchor residue choices (encoded in the fixture):
        ssu_mrna     → S/10
        ssu_atrna    → S/50
        ssu_ptrna    → S/30
        ssu_etrna    → S/70
        lsu_atrna    → L/20
        lsu_ptrna    → L/40
        lsu_etrna    → L/60
    """
    s_chain = _chain("S", rfam=("RF00177",))  # bacterial SSU
    l_chain = _chain("L", rfam=("RF02541",))  # bacterial LSU
    mrna = _chain("M", description="mRNA")
    atrna = _chain("TA", description="tRNA-Phe")
    ptrna = _chain("TP", description="tRNA-Lys")
    etrna = _chain("TE", description="tRNA-Ala")
    asl = _chain("ASL", description="tRNA-Phe anticodon stem-loop")

    eftu = _chain(
        "EFTU",
        polymer_type="Protein",
        description="Elongation factor Tu",
        is_ribosomal_protein=False,
    )
    l1 = _chain(
        "L1",
        polymer_type="Protein",
        description="50S ribosomal protein L1",
        is_ribosomal_protein=True,
    )

    assembly = _assembly(
        rna_chains=[s_chain, l_chain, mrna, atrna, ptrna, etrna, asl],
        protein_chains=[eftu, l1],
    )
    by_role: dict[str, list[ChainRef]] = {
        "ssu_main_rrna": [s_chain],
        "lsu_main_rrna": [l_chain],
        "lsu_associated_rrna": [],
        "trna": [],
        "unmapped": [mrna, atrna, ptrna, etrna, asl],
    }
    correspondence_by_site = {
        "ssu_mrna": _correspondence("ssu_mrna", ["RIBOFIXTURE|1|S|U|10"]),
        "ssu_atrna": _correspondence("ssu_atrna", ["RIBOFIXTURE|1|S|U|50"]),
        "ssu_ptrna": _correspondence("ssu_ptrna", ["RIBOFIXTURE|1|S|U|30"]),
        "ssu_etrna": _correspondence("ssu_etrna", ["RIBOFIXTURE|1|S|U|70"]),
        "lsu_atrna": _correspondence("lsu_atrna", ["RIBOFIXTURE|1|L|U|20"]),
        "lsu_ptrna": _correspondence("lsu_ptrna", ["RIBOFIXTURE|1|L|U|40"]),
        "lsu_etrna": _correspondence("lsu_etrna", ["RIBOFIXTURE|1|L|U|60"]),
    }
    return assembly, by_role, correspondence_by_site


# ---------------------------------------------------------------------------
# assign_functional_chains — §11 happy path
# ---------------------------------------------------------------------------


def test_assignment_happy_path_picks_all_four_chains(
    ribosome_fixture: gemmi.Structure,
) -> None:
    assembly, by_role, correspondence = _build_ribosome_inputs()
    result = infer.assign_functional_chains(ribosome_fixture, assembly, by_role, correspondence)
    assert result.mrna_chain is not None
    assert result.mrna_chain.auth_asym_id == "M"
    assert result.aminoacyl_trna_chain is not None
    assert result.aminoacyl_trna_chain.auth_asym_id == "TA"
    assert result.peptidyl_trna_chain is not None
    assert result.peptidyl_trna_chain.auth_asym_id == "TP"
    assert result.exit_trna_chain is not None
    assert result.exit_trna_chain.auth_asym_id == "TE"


def test_assignment_excludes_rrna_chains_from_pool(
    ribosome_fixture: gemmi.Structure,
) -> None:
    """Even though anchors are on S/L rRNA chains, mRNA/tRNA assignments
    must pull from the non-rRNA candidate pool only."""
    assembly, by_role, correspondence = _build_ribosome_inputs()
    result = infer.assign_functional_chains(ribosome_fixture, assembly, by_role, correspondence)
    assigned_ifes = {
        c.ife
        for c in (
            result.mrna_chain,
            result.aminoacyl_trna_chain,
            result.peptidyl_trna_chain,
            result.exit_trna_chain,
        )
        if c is not None
    }
    rrna_ifes = {by_role["ssu_main_rrna"][0].ife, by_role["lsu_main_rrna"][0].ife}
    assert assigned_ifes.isdisjoint(rrna_ifes)


def test_assignment_no_chain_used_twice(ribosome_fixture: gemmi.Structure) -> None:
    """Pool removal must enforce that the same chain can't fill two roles."""
    assembly, by_role, correspondence = _build_ribosome_inputs()
    result = infer.assign_functional_chains(ribosome_fixture, assembly, by_role, correspondence)
    assigned_ifes = [
        c.ife
        for c in (
            result.mrna_chain,
            result.aminoacyl_trna_chain,
            result.peptidyl_trna_chain,
            result.exit_trna_chain,
        )
        if c is not None
    ]
    assert len(assigned_ifes) == len(set(assigned_ifes))


def test_assignment_no_candidates_returns_all_none(
    ribosome_fixture: gemmi.Structure,
) -> None:
    """When the candidate pool contains only rRNA, nothing should be assigned."""
    s_chain = _chain("S", rfam=("RF00177",))
    l_chain = _chain("L", rfam=("RF02541",))
    assembly = _assembly(rna_chains=[s_chain, l_chain], protein_chains=[])
    by_role: dict[str, list[ChainRef]] = {
        "ssu_main_rrna": [s_chain],
        "lsu_main_rrna": [l_chain],
        "lsu_associated_rrna": [],
        "trna": [],
        "unmapped": [],
    }
    correspondence = {
        "ssu_mrna": _correspondence("ssu_mrna", ["RIBOFIXTURE|1|S|U|10"]),
    }
    result = infer.assign_functional_chains(ribosome_fixture, assembly, by_role, correspondence)
    assert result.mrna_chain is None
    assert result.aminoacyl_trna_chain is None
    assert result.peptidyl_trna_chain is None
    assert result.exit_trna_chain is None


def test_assignment_missing_anchors_skip_their_site(
    ribosome_fixture: gemmi.Structure,
) -> None:
    """Sites absent from correspondence_by_site silently produce no assignment."""
    assembly, by_role, correspondence = _build_ribosome_inputs()
    correspondence = {k: v for k, v in correspondence.items() if k != "ssu_mrna"}
    result = infer.assign_functional_chains(ribosome_fixture, assembly, by_role, correspondence)
    assert result.mrna_chain is None
    # Other assignments still succeed because their anchors are present.
    assert result.aminoacyl_trna_chain is not None
    assert result.peptidyl_trna_chain is not None
    assert result.exit_trna_chain is not None


def test_assignment_missing_residue_in_coords_emits_warning(
    ribosome_fixture: gemmi.Structure,
) -> None:
    """A mapped unit ID the coordinates don't have should surface as a warning."""
    assembly, by_role, correspondence = _build_ribosome_inputs()
    correspondence["ssu_mrna"] = _correspondence("ssu_mrna", ["RIBOFIXTURE|1|S|U|9999"])
    result = infer.assign_functional_chains(ribosome_fixture, assembly, by_role, correspondence)
    assert any("missing_residue_in_coords_ssu_mrna" in w for w in result.warnings)


# ---------------------------------------------------------------------------
# compute_trna_states — §12 classic states (A/A, P/P, E/E)
# ---------------------------------------------------------------------------


def test_classic_atrna_state_is_A_over_A(
    ribosome_fixture: gemmi.Structure,
) -> None:
    """A-tRNA in classic A/A state: contacts SSU ssu_atrna anchor (not the
    SSU ptrna anchor), contacts LSU atrna anchor (not LSU ptrna anchor)."""
    assembly, by_role, correspondence = _build_ribosome_inputs()
    assignments = infer.assign_functional_chains(
        ribosome_fixture, assembly, by_role, correspondence
    )
    states = infer.compute_trna_states(ribosome_fixture, assembly, assignments, correspondence)
    assert states.aminoacyl_trna_state == "A/A"


def test_classic_ptrna_state_is_P_over_P(
    ribosome_fixture: gemmi.Structure,
) -> None:
    assembly, by_role, correspondence = _build_ribosome_inputs()
    assignments = infer.assign_functional_chains(
        ribosome_fixture, assembly, by_role, correspondence
    )
    states = infer.compute_trna_states(ribosome_fixture, assembly, assignments, correspondence)
    assert states.peptidyl_trna_state == "P/P"


def test_classic_etrna_state_is_E_over_E(
    ribosome_fixture: gemmi.Structure,
) -> None:
    """If E-tRNA is assigned at all, the state is always E/E (§12.3)."""
    assembly, by_role, correspondence = _build_ribosome_inputs()
    assignments = infer.assign_functional_chains(
        ribosome_fixture, assembly, by_role, correspondence
    )
    states = infer.compute_trna_states(ribosome_fixture, assembly, assignments, correspondence)
    assert states.exit_trna_state == "E/E"


# ---------------------------------------------------------------------------
# compute_trna_states — protein-factor LSU label (§12.4)
# ---------------------------------------------------------------------------


def test_atrna_factor_label_when_no_LSU_contact(
    ribosome_fixture: gemmi.Structure,
) -> None:
    """Drop the LSU atrna + LSU ptrna anchors so the A-tRNA contacts neither
    LSU site; §12.4 should label LSU with the factor description."""
    assembly, by_role, correspondence = _build_ribosome_inputs()
    correspondence["lsu_atrna"] = _correspondence("lsu_atrna", [])
    correspondence["lsu_ptrna"] = _correspondence("lsu_ptrna", [])
    assignments = infer.assign_functional_chains(
        ribosome_fixture, assembly, by_role, correspondence
    )
    states = infer.compute_trna_states(ribosome_fixture, assembly, assignments, correspondence)
    assert states.aminoacyl_trna_state == "A/Elongation factor Tu"
    # Evidence dict carries chain IFE, raw description, and distance.
    assert states.trna_state_evidence["aminoacyl_trna_factor_chain"] == ("RIBOFIXTURE|1|EFTU")
    assert states.trna_state_evidence["aminoacyl_trna_factor_description"] == "Elongation factor Tu"
    assert states.trna_state_evidence["aminoacyl_trna_factor_distance"] > 0.0


def test_atrna_factor_search_excludes_ribosomal_proteins(
    ribosome_fixture: gemmi.Structure,
) -> None:
    """L1 sits CLOSER to the CCA end than EFTU (0.7 vs 1.0 Å) but is
    flagged is_ribosomal_protein=True. §12.4 must exclude it and pick EFTU."""
    assembly, by_role, correspondence = _build_ribosome_inputs()
    correspondence["lsu_atrna"] = _correspondence("lsu_atrna", [])
    correspondence["lsu_ptrna"] = _correspondence("lsu_ptrna", [])
    assignments = infer.assign_functional_chains(
        ribosome_fixture, assembly, by_role, correspondence
    )
    states = infer.compute_trna_states(ribosome_fixture, assembly, assignments, correspondence)
    # If L1 had been picked the state would be "A/50S ribosomal protein L1".
    assert states.aminoacyl_trna_state == "A/Elongation factor Tu"


def test_atrna_factor_star_when_no_qualifying_protein(
    ribosome_fixture: gemmi.Structure,
) -> None:
    """No LSU contact AND no non-ribosomal protein nearby → "A/*"."""
    assembly, by_role, correspondence = _build_ribosome_inputs()
    correspondence["lsu_atrna"] = _correspondence("lsu_atrna", [])
    correspondence["lsu_ptrna"] = _correspondence("lsu_ptrna", [])
    # Drop EFTU and L1 from the protein pool so §12.4 finds no candidate.
    assembly_no_factor = AssemblyContext(
        pdb_id=assembly.pdb_id,
        assembly_id=assembly.assembly_id,
        experimental_methods=list(assembly.experimental_methods),
        rna_chains=list(assembly.rna_chains),
        protein_chains=[],
        ligands=list(assembly.ligands),
    )
    assignments = infer.assign_functional_chains(
        ribosome_fixture, assembly_no_factor, by_role, correspondence
    )
    states = infer.compute_trna_states(
        ribosome_fixture, assembly_no_factor, assignments, correspondence
    )
    assert states.aminoacyl_trna_state == "A/*"


def test_atrna_double_star_for_asl_fragment(
    ribosome_fixture: gemmi.Structure,
) -> None:
    """An A-tRNA chain with fewer than 30 residues triggers the ASL fragment
    branch and produces "<SSU>/**" regardless of nearby protein factors."""
    s_chain = _chain("S", rfam=("RF00177",))
    l_chain = _chain("L", rfam=("RF02541",))
    asl_chain = _chain("ASL", description="tRNA anticodon stem-loop")
    eftu = _chain(
        "EFTU",
        polymer_type="Protein",
        description="Elongation factor Tu",
    )
    assembly = _assembly(rna_chains=[s_chain, l_chain, asl_chain], protein_chains=[eftu])
    # by_role partition would be unused — the assignment step is bypassed below.
    _ = {
        "ssu_main_rrna": [s_chain],
        "lsu_main_rrna": [l_chain],
        "lsu_associated_rrna": [],
        "trna": [],
        "unmapped": [asl_chain],
    }
    # ASL is far from every anchor in the fixture, so the orchestrator
    # won't assign it. To exercise the **-fragment branch we hand-craft
    # the assignments.
    assignments = infer.ChainAssignments(aminoacyl_trna_chain=asl_chain)
    correspondence = {
        "lsu_atrna": _correspondence("lsu_atrna", []),
        "lsu_ptrna": _correspondence("lsu_ptrna", []),
        "ssu_ptrna": _correspondence("ssu_ptrna", []),
    }
    states = infer.compute_trna_states(ribosome_fixture, assembly, assignments, correspondence)
    assert states.aminoacyl_trna_state == "A/**"


# ---------------------------------------------------------------------------
# compute_trna_states — partial states (hybrid + chimeric)
# ---------------------------------------------------------------------------


def test_atrna_hybrid_state_A_over_P(
    ribosome_fixture: gemmi.Structure,
) -> None:
    """A-tRNA contacts LSU ptrna anchor but NOT LSU atrna anchor → SSU=A,
    LSU=P (the A/P hybrid case)."""
    assembly, by_role, correspondence = _build_ribosome_inputs()
    # Drop only the lsu_atrna anchor — keep lsu_ptrna in play. The
    # fixture's A-tRNA CCA end sits next to L/20 (lsu_atrna) only; to
    # exercise the hybrid we redirect lsu_ptrna onto L/20 instead.
    correspondence["lsu_atrna"] = _correspondence("lsu_atrna", [])
    correspondence["lsu_ptrna"] = _correspondence("lsu_ptrna", ["RIBOFIXTURE|1|L|U|20"])
    assignments = infer.assign_functional_chains(
        ribosome_fixture, assembly, by_role, correspondence
    )
    states = infer.compute_trna_states(ribosome_fixture, assembly, assignments, correspondence)
    assert states.aminoacyl_trna_state == "A/P"


def test_atrna_chimeric_ap_over_AP(
    ribosome_fixture: gemmi.Structure,
) -> None:
    """A-tRNA contacts BOTH ssu_atrna AND ssu_ptrna anchors → SSU=ap.
    Same residue covers both LSU sites → LSU=AP.

    Bypasses :func:`assign_functional_chains` because the chimeric
    geometry (one tRNA contacting both SSU sites) would otherwise pull
    the same chain into both P-tRNA and A-tRNA roles via the pool
    removal — the assignment side is correctly defensive about that. The
    state-computation logic is what's under test here.
    """
    assembly, _, correspondence = _build_ribosome_inputs()
    # Anchors that A-tRNA's chain TA contacts:
    #   ssu_atrna stays at S/50; move ssu_ptrna onto S/50 too.
    #   lsu_atrna stays at L/20; move lsu_ptrna onto L/20 too.
    correspondence["ssu_ptrna"] = _correspondence("ssu_ptrna", ["RIBOFIXTURE|1|S|U|50"])
    correspondence["lsu_ptrna"] = _correspondence("lsu_ptrna", ["RIBOFIXTURE|1|L|U|20"])
    ta_chain = next(c for c in assembly.rna_chains if c.auth_asym_id == "TA")
    assignments = infer.ChainAssignments(aminoacyl_trna_chain=ta_chain)
    states = infer.compute_trna_states(ribosome_fixture, assembly, assignments, correspondence)
    assert states.aminoacyl_trna_state == "ap/AP"


def test_ptrna_hybrid_state_P_over_E(
    ribosome_fixture: gemmi.Structure,
) -> None:
    assembly, by_role, correspondence = _build_ribosome_inputs()
    # P-tRNA contacts lsu_etrna anchor but not lsu_ptrna anchor.
    correspondence["lsu_ptrna"] = _correspondence("lsu_ptrna", [])
    correspondence["lsu_etrna"] = _correspondence(
        "lsu_etrna",
        ["RIBOFIXTURE|1|L|U|40"],  # the P-tRNA's CCA-end residue
    )
    assignments = infer.assign_functional_chains(
        ribosome_fixture, assembly, by_role, correspondence
    )
    states = infer.compute_trna_states(ribosome_fixture, assembly, assignments, correspondence)
    assert states.peptidyl_trna_state == "P/E"


# ---------------------------------------------------------------------------
# compute_trna_states — None states when no tRNA assigned
# ---------------------------------------------------------------------------


def test_states_are_none_when_no_chains_assigned(
    ribosome_fixture: gemmi.Structure,
) -> None:
    assembly, _, correspondence = _build_ribosome_inputs()
    empty_assignments = infer.ChainAssignments()
    states = infer.compute_trna_states(
        ribosome_fixture, assembly, empty_assignments, correspondence
    )
    assert states.aminoacyl_trna_state is None
    assert states.peptidyl_trna_state is None
    assert states.exit_trna_state is None


def test_etrna_state_is_none_when_only_etrna_missing(
    ribosome_fixture: gemmi.Structure,
) -> None:
    assembly, by_role, correspondence = _build_ribosome_inputs()
    assignments = infer.assign_functional_chains(
        ribosome_fixture, assembly, by_role, correspondence
    )
    # Manually clear just the E-tRNA assignment.
    assignments = assignments.model_copy(update={"exit_trna_chain": None})
    states = infer.compute_trna_states(ribosome_fixture, assembly, assignments, correspondence)
    assert states.exit_trna_state is None
    assert states.aminoacyl_trna_state == "A/A"


# ---------------------------------------------------------------------------
# Parameterised: SSU + LSU truth table for A-tRNA
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("ssu_ptrna_anchor", "lsu_atrna_anchor", "lsu_ptrna_anchor", "expected"),
    [
        # SSU=ap (touches both ssu_atrna+ssu_ptrna), LSU=AP (touches both LSU sites)
        ("RIBOFIXTURE|1|S|U|50", "RIBOFIXTURE|1|L|U|20", "RIBOFIXTURE|1|L|U|20", "ap/AP"),
        # SSU=A, LSU=AP
        ("", "RIBOFIXTURE|1|L|U|20", "RIBOFIXTURE|1|L|U|20", "A/AP"),
        # SSU=A, LSU=A
        ("", "RIBOFIXTURE|1|L|U|20", "", "A/A"),
        # SSU=A, LSU=P (hybrid)
        ("", "", "RIBOFIXTURE|1|L|U|20", "A/P"),
    ],
)
def test_atrna_truth_table(
    ribosome_fixture: gemmi.Structure,
    ssu_ptrna_anchor: str,
    lsu_atrna_anchor: str,
    lsu_ptrna_anchor: str,
    expected: str,
) -> None:
    """State-computation truth table. Bypasses chain assignment because the
    ``ap/AP`` case has TA touching anchors normally owned by P-tRNA — the
    correct defensive behaviour of :func:`assign_functional_chains` (no
    chain in two roles) would block it."""
    assembly, _, correspondence = _build_ribosome_inputs()
    correspondence["ssu_ptrna"] = _correspondence(
        "ssu_ptrna", [ssu_ptrna_anchor] if ssu_ptrna_anchor else []
    )
    correspondence["lsu_atrna"] = _correspondence(
        "lsu_atrna", [lsu_atrna_anchor] if lsu_atrna_anchor else []
    )
    correspondence["lsu_ptrna"] = _correspondence(
        "lsu_ptrna", [lsu_ptrna_anchor] if lsu_ptrna_anchor else []
    )
    ta_chain = next(c for c in assembly.rna_chains if c.auth_asym_id == "TA")
    assignments = infer.ChainAssignments(aminoacyl_trna_chain=ta_chain)
    states = infer.compute_trna_states(ribosome_fixture, assembly, assignments, correspondence)
    assert states.aminoacyl_trna_state == expected
