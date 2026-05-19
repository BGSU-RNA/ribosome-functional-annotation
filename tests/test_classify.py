"""Unit tests for the complete-ribosome filter and classifier (spec §7, §8)."""

from __future__ import annotations

import pytest

from ribosome_state_annotator import classify
from ribosome_state_annotator import constants as C
from ribosome_state_annotator.config import CompletenessThresholds
from ribosome_state_annotator.models import AssemblyContext, ChainRef, LigandRef

# ---------------------------------------------------------------------------
# Fixture builders
# ---------------------------------------------------------------------------


def _chain(
    asym: str,
    *,
    polymer_type: str = "Protein",
    description: str | None = None,
    uniprot_name: str | None = None,
    superkingdom: str | None = None,
    rfam: tuple[str, ...] = (),
    pdb_id: str = "5J7L",
    assembly_id: str = "1",
) -> ChainRef:
    return ChainRef(
        pdb_id=pdb_id,
        assembly_id=assembly_id,
        auth_asym_id=asym,
        polymer_type=polymer_type,
        description=description,
        uniprot_name=uniprot_name,
        superkingdom=superkingdom,
        rfam_accessions=list(rfam),
    )


def _bacterial_ssu(asym: str = "AA", **kw: object) -> ChainRef:
    return _chain(asym, polymer_type="RNA", rfam=("RF00177",), **kw)  # type: ignore[arg-type]


def _bacterial_lsu(asym: str = "DA", **kw: object) -> ChainRef:
    return _chain(asym, polymer_type="RNA", rfam=("RF02541",), **kw)  # type: ignore[arg-type]


def _eukaryotic_ssu(asym: str = "A", **kw: object) -> ChainRef:
    return _chain(asym, polymer_type="RNA", rfam=("RF01960",), **kw)  # type: ignore[arg-type]


def _eukaryotic_lsu(asym: str = "1", **kw: object) -> ChainRef:
    return _chain(asym, polymer_type="RNA", rfam=("RF02543",), **kw)  # type: ignore[arg-type]


def _ribosomal_protein(
    asym: str,
    superkingdom: str = "Bacteria",
    description: str = "50S ribosomal protein L1",
) -> ChainRef:
    return _chain(
        asym,
        polymer_type="Protein",
        description=description,
        superkingdom=superkingdom,
    )


def _assembly(
    *,
    rna_chains: list[ChainRef] | None = None,
    protein_chains: list[ChainRef] | None = None,
    ligands: list[LigandRef] | None = None,
    methods: list[str] | None = None,
    pdb_id: str = "5J7L",
    assembly_id: str = "1",
) -> AssemblyContext:
    return AssemblyContext(
        pdb_id=pdb_id,
        assembly_id=assembly_id,
        experimental_methods=methods if methods is not None else ["X-RAY DIFFRACTION"],
        rna_chains=rna_chains or [],
        protein_chains=protein_chains or [],
        ligands=ligands or [],
    )


def _n_proteins(n: int, superkingdom: str = "Bacteria") -> list[ChainRef]:
    return [_ribosomal_protein(f"P{i}", superkingdom=superkingdom) for i in range(n)]


# ---------------------------------------------------------------------------
# matches_ribosomal_protein_narrow (§13.1)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("description", "uniprot", "expected"),
    [
        ("50S ribosomal protein L1", None, True),
        ("30S RIBOSOMAL PROTEIN S1", None, True),
        (None, "ribosomal protein L11", True),
        # Ban-nomenclature pattern used by UniProt for modern ribosomal
        # protein names (cytoplasmic and mitoribosomal). The narrow rule
        # must accept "ribosomal subunit protein" as a synonym of
        # "ribosomal protein", otherwise every mitoribosomal protein gets
        # mis-classified as a non-ribosomal factor candidate (cf. 3J9M
        # human 55S mitoribosome).
        ("uL2m", "Large ribosomal subunit protein uL2m", True),
        ("uS5m", "Small ribosomal subunit protein uS5m", True),
        (None, "Large ribosomal subunit protein eL40", True),
        # Substring is case-insensitive.
        (None, "LARGE RIBOSOMAL SUBUNIT PROTEIN UL2M", True),
        # Ban-nomenclature short-form pattern: when RCSB stores only the
        # bare systematic name as pdbx_description and provides no
        # UniProt cross-reference, the narrow rule must still recognise
        # it as a ribosomal protein. See 3J9M chain X (bL28m) for the
        # in-the-wild example, and REFERENCES.md for the Ban 2014
        # nomenclature paper.
        ("bL28m", None, True),
        ("uL2m", None, True),
        ("mS35", None, True),
        ("eL40", None, True),
        ("bS6m", None, True),
        # Looks Ban-like but isn't — should NOT match.
        ("mLpRetraction", None, False),
        ("uLtimately", None, False),
        # Non-Ban prefixes used by older deposits — not matched by the
        # narrow rule (the broader §8.4 rule handles them for the
        # superkingdom vote).
        ("RPS5", None, False),
        ("MRPS6", None, False),
        ("Elongation factor Tu", None, False),
        ("16S ribosomal RNA", None, False),  # "ribosomal" but not "ribosomal protein"
        ("S1", None, False),  # §13.1 known false negative
        (None, None, False),
        ("", "", False),
    ],
)
def test_matches_ribosomal_protein_narrow(
    description: str | None, uniprot: str | None, expected: bool
) -> None:
    assert classify.matches_ribosomal_protein_narrow(description, uniprot) is expected


# ---------------------------------------------------------------------------
# matches_ribosomal_protein_broad (§8.4)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("description", "uniprot", "expected", "label"),
    [
        ("50S ribosomal protein L1", None, True, "rule-1-substring"),
        ("ribosomal subunit protein S3", None, True, "rule-1-subunit-substring"),
        ("S1", None, True, "rule-2-regex-bare-S1"),
        ("L11", None, True, "rule-2-regex-bare-L11"),
        ("uS3", None, True, "rule-2-regex-modern-uS3"),
        ("mL11", None, True, "rule-2-regex-modern-mL11"),
        ("S12A", None, True, "rule-2-regex-S12-with-suffix"),
        ("RPS", None, True, "rule-2-regex-RPS-bare"),
        ("RPL", None, True, "rule-2-regex-RPL-bare"),
        ("MRPL", None, True, "rule-2-regex-MRPL-bare"),
        ("30S ribosomal subunit", None, True, "rule-4-30S+subunit"),
        ("40S protein component", None, True, "rule-4-40S+protein+component"),
        ("60S subunit component", None, True, "rule-4-60S+subunit+component"),
        (None, "ribosomal protein L11", True, "rule-3-uniprot"),
        # NEGATIVES
        ("Elongation factor Tu", None, False, "factor-not-ribosomal"),
        ("16S ribosomal RNA", None, False, "rRNA-not-ribosomal-protein"),
        ("RPS3 isoform 1", None, False, "rule-2-RPS3-has-trailing-digit-no-match"),
        ("Some 30S binding study", None, False, "30S-without-subunit-word"),
        ("30S", None, False, "30S-bare-no-subunit-word"),
        (None, None, False, "both-null"),
    ],
)
def test_matches_ribosomal_protein_broad(
    description: str | None,
    uniprot: str | None,
    expected: bool,
    label: str,
) -> None:
    chain = ChainRef(
        pdb_id="X",
        assembly_id="1",
        auth_asym_id="X",
        description=description,
        uniprot_name=uniprot,
    )
    assert classify.matches_ribosomal_protein_broad(chain) is expected, label


def test_broad_is_superset_of_narrow() -> None:
    """Every chain narrow-classifies as ribosomal must also broad-classify."""
    for desc, up in (
        ("50S ribosomal protein L1", None),
        ("RIBOSOMAL PROTEIN L11", None),
        (None, "ribosomal protein S3"),
    ):
        chain = ChainRef(
            pdb_id="X",
            assembly_id="1",
            auth_asym_id="X",
            description=desc,
            uniprot_name=up,
        )
        assert classify.matches_ribosomal_protein_narrow(desc, up)
        assert classify.matches_ribosomal_protein_broad(chain)


# ---------------------------------------------------------------------------
# role_for_chain / partition_rna_chains_by_role
# ---------------------------------------------------------------------------


def test_role_for_chain_maps_to_each_role() -> None:
    assert classify.role_for_chain(_bacterial_ssu()) == "ssu_main_rrna"
    assert classify.role_for_chain(_bacterial_lsu()) == "lsu_main_rrna"
    assert (
        classify.role_for_chain(_chain("B", rfam=("RF00001",), polymer_type="RNA"))
        == "lsu_associated_rrna"
    )
    assert classify.role_for_chain(_chain("T", rfam=("RF00005",), polymer_type="RNA")) == "trna"


def test_role_for_chain_returns_none_when_no_match() -> None:
    chain = _chain("X", polymer_type="RNA", rfam=("RF99999",))
    assert classify.role_for_chain(chain) is None


def test_role_for_chain_picks_first_when_multiple_match() -> None:
    """When a chain has Rfam accessions in multiple ribosomal roles, the
    function logs a warning and returns the alphabetically-first role for
    determinism."""
    chain = _chain("X", polymer_type="RNA", rfam=("RF00177", "RF00001"))
    # roles are {"ssu_main_rrna", "lsu_associated_rrna"}; alphabetical first = "lsu_associated_rrna"
    assert classify.role_for_chain(chain) == "lsu_associated_rrna"


def test_partition_rna_chains_buckets_all_chains() -> None:
    chains = [
        _bacterial_ssu(),
        _bacterial_lsu(),
        _chain("5S", polymer_type="RNA", rfam=("RF00001",)),
        _chain("T1", polymer_type="RNA", rfam=("RF00005",)),
        _chain("M", polymer_type="RNA"),  # unmapped (no Rfam)
    ]
    buckets = classify.partition_rna_chains_by_role(chains)
    assert {b.auth_asym_id for b in buckets["ssu_main_rrna"]} == {"AA"}
    assert {b.auth_asym_id for b in buckets["lsu_main_rrna"]} == {"DA"}
    assert {b.auth_asym_id for b in buckets["lsu_associated_rrna"]} == {"5S"}
    assert {b.auth_asym_id for b in buckets["trna"]} == {"T1"}
    assert {b.auth_asym_id for b in buckets["unmapped"]} == {"M"}


# ---------------------------------------------------------------------------
# determine_rrna_core (§8.2)
# ---------------------------------------------------------------------------


def test_rrna_core_bacterial() -> None:
    assert classify.determine_rrna_core([_bacterial_ssu()], [_bacterial_lsu()]) == "bacterial_like"


def test_rrna_core_eukaryotic() -> None:
    assert (
        classify.determine_rrna_core([_eukaryotic_ssu()], [_eukaryotic_lsu()]) == "eukaryotic_like"
    )


def test_rrna_core_mixed() -> None:
    """Same chain hits both bacterial and eukaryotic Rfam accessions — rare but
    happens with mis-annotated entries."""
    ssu = _chain("AA", polymer_type="RNA", rfam=("RF00177", "RF01960"))
    lsu = _chain("DA", polymer_type="RNA", rfam=("RF02541", "RF02543"))
    assert classify.determine_rrna_core([ssu], [lsu]) == "mixed"


def test_rrna_core_ambiguous_when_only_ssu() -> None:
    """Bacterial SSU but no matching LSU Rfam → ambiguous."""
    assert classify.determine_rrna_core([_bacterial_ssu()], []) == "ambiguous"


def test_rrna_core_ambiguous_when_mismatched() -> None:
    """Bacterial SSU but eukaryotic LSU → ambiguous (no matching core)."""
    assert classify.determine_rrna_core([_bacterial_ssu()], [_eukaryotic_lsu()]) == "ambiguous"


# ---------------------------------------------------------------------------
# compute_dominant_superkingdom (§8.3)
# ---------------------------------------------------------------------------


def test_dominant_superkingdom_clean_majority() -> None:
    chains = _n_proteins(5, "Bacteria")
    sk, counts = classify.compute_dominant_superkingdom(chains)
    assert sk == "Bacteria"
    assert counts == {"Bacteria": 5}


def test_dominant_superkingdom_eukaryote_clean() -> None:
    sk, _ = classify.compute_dominant_superkingdom(_n_proteins(5, "Eukaryota"))
    assert sk == "Eukaryota"


def test_dominant_superkingdom_bacteria_wins_tie_over_eukaryota() -> None:
    chains = _n_proteins(3, "Bacteria") + _n_proteins(3, "Eukaryota")
    sk, _ = classify.compute_dominant_superkingdom(chains)
    assert sk == "Bacteria"


def test_dominant_superkingdom_eukaryota_wins_tie_over_archaea() -> None:
    chains = _n_proteins(3, "Eukaryota") + _n_proteins(3, "Archaea")
    sk, _ = classify.compute_dominant_superkingdom(chains)
    assert sk == "Eukaryota"


def test_dominant_superkingdom_unknown_when_below_min_voters() -> None:
    """Fewer than 3 voters → 'unknown' regardless of their distribution."""
    sk, counts = classify.compute_dominant_superkingdom(_n_proteins(2, "Bacteria"))
    assert sk == "unknown"
    assert counts == {"Bacteria": 2}


def test_dominant_superkingdom_unknown_when_no_taxonomy() -> None:
    """All chains have null superkingdom → 'unknown' with empty counts."""
    chains = [_ribosomal_protein(f"P{i}", superkingdom="") for i in range(5)]
    # superkingdom="" gets filtered as falsy — set to None instead
    chains = [
        ChainRef(pdb_id="X", assembly_id="1", auth_asym_id=f"P{i}", superkingdom=None)
        for i in range(5)
    ]
    sk, counts = classify.compute_dominant_superkingdom(chains)
    assert sk == "unknown"
    assert counts == {}


def test_dominant_superkingdom_min_voters_threshold_configurable() -> None:
    """The caller can lower the threshold for testing."""
    sk, _ = classify.compute_dominant_superkingdom(_n_proteins(2, "Bacteria"), min_voters=1)
    assert sk == "Bacteria"


def test_dominant_superkingdom_empty_voters_with_zero_threshold() -> None:
    """Edge case: min_voters=0 and no voters at all → still unknown, not a crash."""
    sk, counts = classify.compute_dominant_superkingdom([], min_voters=0)
    assert sk == "unknown"
    assert counts == {}


def test_dominant_superkingdom_falls_back_when_leader_outside_canonical_set() -> None:
    """A non-canonical superkingdom (e.g. 'Viruses') leads → returned as-is so
    the §8.5 truth table can route it to insufficient_taxonomy."""
    chains = [_ribosomal_protein(f"P{i}", superkingdom="Viruses") for i in range(5)]
    sk, counts = classify.compute_dominant_superkingdom(chains)
    assert sk == "Viruses"
    assert counts == {"Viruses": 5}


def test_classification_result_is_skip_property() -> None:
    skip_result = classify.ClassificationResult(skip_reason=C.SKIP_NMR)
    assert skip_result.is_skip is True
    ok_result = classify.ClassificationResult(classification="bacterial_ribosome")
    assert ok_result.is_skip is False


def test_supporting_terms_skips_chains_with_no_description() -> None:
    chains = [
        _ribosomal_protein("P1", description="Mitochondrial ribosomal protein"),
        ChainRef(pdb_id="X", assembly_id="1", auth_asym_id="P2"),  # no description
    ]
    assert classify.extract_supporting_terms(chains) == ["mitochondrial"]


# ---------------------------------------------------------------------------
# extract_supporting_terms (§8.6)
# ---------------------------------------------------------------------------


def test_supporting_terms_finds_keyword() -> None:
    chains = [_ribosomal_protein("P1", description="Mitochondrial ribosomal protein L11")]
    assert classify.extract_supporting_terms(chains) == ["mitochondrial"]


def test_supporting_terms_word_boundary_not_substring() -> None:
    """'mitochondrially' must NOT count for 'mitochondrial' (§8.6 says exact word match)."""
    chains = [_ribosomal_protein("P1", description="Encoded mitochondrially in S. cerevisiae")]
    assert classify.extract_supporting_terms(chains) == []


def test_supporting_terms_case_insensitive() -> None:
    chains = [_ribosomal_protein("P1", description="CHLOROPLAST ribosomal protein")]
    assert classify.extract_supporting_terms(chains) == ["chloroplast"]


def test_supporting_terms_dedup_across_chains() -> None:
    chains = [
        _ribosomal_protein("P1", description="Mitochondrial ribosomal protein L11"),
        _ribosomal_protein("P2", description="Mitochondrial ribosomal protein S2"),
    ]
    assert classify.extract_supporting_terms(chains) == ["mitochondrial"]


def test_supporting_terms_returns_sorted() -> None:
    chains = [
        _ribosomal_protein("P1", description="Plastid ribosomal protein"),
        _ribosomal_protein("P2", description="Mitochondrial ribosomal protein"),
    ]
    assert classify.extract_supporting_terms(chains) == ["mitochondrial", "plastid"]


def test_supporting_terms_empty_when_no_keywords() -> None:
    chains = [_ribosomal_protein("P1", description="50S ribosomal protein L1")]
    assert classify.extract_supporting_terms(chains) == []


# ---------------------------------------------------------------------------
# classify_assembly — §8.5 truth table (one test per rule)
# ---------------------------------------------------------------------------


def test_bacterial_ribosome_classified() -> None:
    a = _assembly(
        rna_chains=[_bacterial_ssu(), _bacterial_lsu()],
        protein_chains=_n_proteins(20, "Bacteria"),
    )
    result = classify.classify_assembly(a)
    assert result.classification == "bacterial_ribosome"
    assert result.skip_reason is None
    assert result.evidence["rule"] == C.CLASSIFICATION_RULE_BACTERIAL
    assert result.evidence["rrna_core"] == "bacterial_like"
    assert result.evidence["dominant_ribosomal_protein_superkingdom"] == "Bacteria"
    assert result.warnings == []


def test_eukaryotic_ribosome_classified() -> None:
    a = _assembly(
        rna_chains=[_eukaryotic_ssu(), _eukaryotic_lsu()],
        protein_chains=_n_proteins(35, "Eukaryota"),
    )
    result = classify.classify_assembly(a)
    assert result.classification == "eukaryotic_ribosome"
    assert result.evidence["rule"] == C.CLASSIFICATION_RULE_EUKARYOTIC


def test_eukaryotic_organellar_classified() -> None:
    """Bacterial-like rRNA + Eukaryota proteins → organellar."""
    a = _assembly(
        rna_chains=[_bacterial_ssu(), _bacterial_lsu()],
        protein_chains=[
            _ribosomal_protein(f"P{i}", "Eukaryota", "Mitochondrial ribosomal protein L11")
            for i in range(25)
        ],
    )
    result = classify.classify_assembly(a)
    assert result.classification == "eukaryotic_organellar_ribosome"
    assert result.evidence["rule"] == C.CLASSIFICATION_RULE_ORGANELLAR
    assert "mitochondrial" in result.evidence["supporting_terms"]


def test_organellar_without_keyword_emits_warning() -> None:
    """Bacterial-like rRNA + Eukaryota proteins with no organellar keywords →
    classify as organellar but emit warning."""
    a = _assembly(
        rna_chains=[_bacterial_ssu(), _bacterial_lsu()],
        protein_chains=[
            _ribosomal_protein(f"P{i}", "Eukaryota", f"50S ribosomal protein L{i}")
            for i in range(25)
        ],
    )
    result = classify.classify_assembly(a)
    assert result.classification == "eukaryotic_organellar_ribosome"
    assert classify.WARNING_ORGANELLAR_NO_KEYWORDS in result.warnings


def test_eukaryotic_like_rrna_plus_bacteria_proteins_is_ambiguous() -> None:
    a = _assembly(
        rna_chains=[_eukaryotic_ssu(), _eukaryotic_lsu()],
        protein_chains=_n_proteins(30, "Bacteria"),
    )
    result = classify.classify_assembly(a)
    assert result.classification is None
    assert result.skip_reason == C.SKIP_AMBIGUOUS_CLASSIFICATION


# ---------------------------------------------------------------------------
# classify_assembly — skip paths
# ---------------------------------------------------------------------------


def test_nmr_short_circuits_with_skip_reason() -> None:
    a = _assembly(
        methods=["SOLUTION NMR"],
        rna_chains=[_bacterial_ssu(), _bacterial_lsu()],
        protein_chains=_n_proteins(20, "Bacteria"),
    )
    result = classify.classify_assembly(a)
    assert result.classification is None
    assert result.skip_reason == C.SKIP_NMR


def test_missing_ssu_short_circuits() -> None:
    a = _assembly(
        rna_chains=[_bacterial_lsu()],
        protein_chains=_n_proteins(20, "Bacteria"),
    )
    result = classify.classify_assembly(a)
    assert result.skip_reason == C.SKIP_PARTIAL_MISSING_SSU_OR_LSU


def test_missing_lsu_short_circuits() -> None:
    a = _assembly(
        rna_chains=[_bacterial_ssu()],
        protein_chains=_n_proteins(20, "Bacteria"),
    )
    result = classify.classify_assembly(a)
    assert result.skip_reason == C.SKIP_PARTIAL_MISSING_SSU_OR_LSU


def test_ambiguous_rrna_core_short_circuits() -> None:
    """Bacterial SSU + eukaryotic LSU → ambiguous_rrna_core."""
    a = _assembly(
        rna_chains=[_bacterial_ssu(), _eukaryotic_lsu()],
        protein_chains=_n_proteins(20, "Bacteria"),
    )
    result = classify.classify_assembly(a)
    assert result.skip_reason == C.SKIP_AMBIGUOUS_RRNA_CORE


def test_mixed_rrna_core_resolved_to_bacterial_when_proteins_bacterial() -> None:
    """MIXED rrna core resolves to bacterial when the protein vote is Bacteria."""
    ssu = _chain("AA", polymer_type="RNA", rfam=("RF00177", "RF01960"))
    lsu = _chain("DA", polymer_type="RNA", rfam=("RF02541", "RF02543"))
    a = _assembly(rna_chains=[ssu, lsu], protein_chains=_n_proteins(20, "Bacteria"))
    result = classify.classify_assembly(a)
    assert result.classification == "bacterial_ribosome"
    assert classify.WARNING_MIXED_RRNA_CORE in result.warnings
    assert result.evidence["rrna_core_resolved"] == classify.RRNA_CORE_BACTERIAL


def test_mixed_rrna_core_resolved_to_eukaryotic_when_proteins_eukaryotic() -> None:
    """9B0S-style PDBe over-annotation: MIXED rrna core + Eukaryota proteins
    should resolve to eukaryotic_ribosome, not eukaryotic_organellar_ribosome."""
    ssu = _chain("AA", polymer_type="RNA", rfam=("RF00177", "RF01959", "RF01960"))
    lsu = _chain("DA", polymer_type="RNA", rfam=("RF02540", "RF02541", "RF02543"))
    a = _assembly(rna_chains=[ssu, lsu], protein_chains=_n_proteins(20, "Eukaryota"))
    result = classify.classify_assembly(a)
    assert result.classification == "eukaryotic_ribosome"
    assert classify.WARNING_MIXED_RRNA_CORE in result.warnings
    assert result.evidence["rrna_core_resolved"] == classify.RRNA_CORE_EUKARYOTIC


def test_archaeal_short_circuits() -> None:
    """Dominant superkingdom = Archaea → archaeal_ribosome_not_supported."""
    a = _assembly(
        rna_chains=[_bacterial_ssu(), _bacterial_lsu()],
        protein_chains=_n_proteins(20, "Archaea"),
    )
    result = classify.classify_assembly(a)
    assert result.skip_reason == C.SKIP_ARCHAEAL


def test_insufficient_taxonomy_short_circuits() -> None:
    """Fewer than 3 voters with non-null superkingdom → insufficient_taxonomy."""
    a = _assembly(
        rna_chains=[_bacterial_ssu(), _bacterial_lsu()],
        protein_chains=_n_proteins(2, "Bacteria"),
    )
    result = classify.classify_assembly(a)
    assert result.skip_reason == C.SKIP_INSUFFICIENT_TAXONOMY


# ---------------------------------------------------------------------------
# classify_assembly — completeness threshold (§7.4)
# ---------------------------------------------------------------------------


def test_low_protein_count_warns_in_default_mode() -> None:
    """Default (non-strict): low protein count emits warning, classification proceeds."""
    a = _assembly(
        rna_chains=[_bacterial_ssu(), _bacterial_lsu()],
        protein_chains=_n_proteins(5, "Bacteria"),  # below bacterial threshold of 15
    )
    result = classify.classify_assembly(a, strict_complete_check=False)
    assert result.classification == "bacterial_ribosome"
    assert classify.WARNING_LOW_PROTEIN_COUNT in result.warnings


def test_low_protein_count_skips_in_strict_mode() -> None:
    a = _assembly(
        rna_chains=[_bacterial_ssu(), _bacterial_lsu()],
        protein_chains=_n_proteins(5, "Bacteria"),
    )
    result = classify.classify_assembly(a, strict_complete_check=True)
    assert result.classification is None
    assert result.skip_reason == C.SKIP_LOW_PROTEIN_COUNT


def test_organellar_threshold_applied() -> None:
    """Organellar threshold is 20, not 15. 18 Eukaryota proteins should warn (organellar)."""
    a = _assembly(
        rna_chains=[_bacterial_ssu(), _bacterial_lsu()],
        protein_chains=[
            _ribosomal_protein(f"P{i}", "Eukaryota", "Mitochondrial ribosomal protein L11")
            for i in range(18)
        ],
    )
    result = classify.classify_assembly(a)
    assert result.classification == "eukaryotic_organellar_ribosome"
    assert classify.WARNING_LOW_PROTEIN_COUNT in result.warnings


def test_custom_thresholds_honoured() -> None:
    a = _assembly(
        rna_chains=[_bacterial_ssu(), _bacterial_lsu()],
        protein_chains=_n_proteins(5, "Bacteria"),
    )
    result = classify.classify_assembly(
        a,
        thresholds=CompletenessThresholds(bacterial=3),
    )
    assert result.classification == "bacterial_ribosome"
    assert classify.WARNING_LOW_PROTEIN_COUNT not in result.warnings


# ---------------------------------------------------------------------------
# classify_assembly — evidence dict fields (§8.7)
# ---------------------------------------------------------------------------


def test_evidence_dict_has_required_keys_on_success() -> None:
    a = _assembly(
        rna_chains=[_bacterial_ssu(), _bacterial_lsu()],
        protein_chains=_n_proteins(20, "Bacteria"),
    )
    ev = classify.classify_assembly(a).evidence
    expected = {
        "ssu_rfam",
        "lsu_rfam",
        "rrna_core",
        "dominant_ribosomal_protein_superkingdom",
        "ribosomal_protein_superkingdom_votes",
        "ribosomal_protein_chains_voting",
        "supporting_terms",
        "rule",
    }
    assert expected.issubset(ev.keys())
    assert ev["ssu_rfam"] == ["RF00177"]
    assert ev["lsu_rfam"] == ["RF02541"]
    assert ev["ribosomal_protein_chains_voting"] == 20
    assert ev["ribosomal_protein_superkingdom_votes"] == {"Bacteria": 20}


def test_evidence_dict_populated_even_on_skip() -> None:
    """A partial assembly is skipped but ssu_rfam/lsu_rfam are still reported."""
    a = _assembly(
        rna_chains=[_bacterial_ssu()],  # SSU only
        protein_chains=_n_proteins(20, "Bacteria"),
    )
    result = classify.classify_assembly(a)
    assert result.skip_reason == C.SKIP_PARTIAL_MISSING_SSU_OR_LSU
    assert result.evidence["ssu_rfam"] == ["RF00177"]
    assert result.evidence["lsu_rfam"] == []
