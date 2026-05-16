"""Unit tests for the curated reference constants (spec §6.1, §6.3, §8.6)."""

from __future__ import annotations

import re

import pytest

from ribosome_state_annotator import constants as C

# Format: PDB|model|chain|residue_name|residue_number
# residue_name can include modified-nt codes like "4OC", "OMG"; residue_number is digits.
BGSU_UNIT_ID_RE = re.compile(r"^[A-Za-z0-9]{4}\|\d+\|[A-Za-z0-9]+\|[A-Za-z0-9]+\|-?\d+$")


# ---------------------------------------------------------------------------
# Rfam role tables
# ---------------------------------------------------------------------------


def test_role_sets_partition_role_map() -> None:
    """Every key in RFAM_ROLE_MAP appears in exactly one role set, and every
    role-set member appears in RFAM_ROLE_MAP."""
    role_to_set = {
        "ssu_main_rrna": C.SSU_MAIN_RRNA,
        "lsu_main_rrna": C.LSU_MAIN_RRNA,
        "lsu_associated_rrna": C.LSU_ASSOCIATED_RRNA,
        "trna": C.TRNA,
    }
    for accession, role in C.RFAM_ROLE_MAP.items():
        assert role in role_to_set, f"{accession} has unknown role {role!r}"
        assert accession in role_to_set[role], (
            f"{accession} mapped to role {role!r} but is not in the corresponding set"
        )
    flat_union = set().union(*role_to_set.values())
    assert flat_union == set(C.RFAM_ROLE_MAP.keys())


def test_role_sets_are_disjoint() -> None:
    sets = [C.SSU_MAIN_RRNA, C.LSU_MAIN_RRNA, C.LSU_ASSOCIATED_RRNA, C.TRNA]
    for i, a in enumerate(sets):
        for b in sets[i + 1 :]:
            assert a.isdisjoint(b), f"role sets overlap: {a & b}"


def test_core_constants_match_role_sets() -> None:
    assert C.BACTERIAL_RRNA_CORE_SSU in C.SSU_MAIN_RRNA
    assert C.BACTERIAL_RRNA_CORE_LSU in C.LSU_MAIN_RRNA
    assert C.EUKARYOTIC_RRNA_CORE_SSU in C.SSU_MAIN_RRNA
    assert C.EUKARYOTIC_RRNA_CORE_LSU in C.LSU_MAIN_RRNA
    assert C.LSU_MEDIUM_RFAM in C.LSU_ASSOCIATED_RRNA
    assert C.LSU_SMALL_RFAM in C.LSU_ASSOCIATED_RRNA
    assert C.LSU_MEDIUM_RFAM != C.LSU_SMALL_RFAM


# ---------------------------------------------------------------------------
# Functional-site anchor sets
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("reference_units", "expected_pdb"),
    [
        (C.ECOLI_REFERENCE_UNITS, "5J7L"),
        (C.YEAST_REFERENCE_UNITS, "7ZW0"),
    ],
)
def test_reference_units_have_expected_site_keys(
    reference_units: dict[str, tuple[str, ...]], expected_pdb: str
) -> None:
    assert set(reference_units.keys()) == set(C.REFERENCE_SITE_KEYS)
    for site, units in reference_units.items():
        assert len(units) >= 1, f"{site}: at least one anchor required"
        for unit_id in units:
            assert BGSU_UNIT_ID_RE.match(unit_id), (
                f"{site}: {unit_id!r} is not a valid BGSU unit ID"
            )
            assert unit_id.startswith(f"{expected_pdb}|1|"), (
                f"{site}: {unit_id!r} should start with {expected_pdb}|1|"
            )


def test_reference_unit_ssu_lsu_chain_distinct() -> None:
    """Within each reference dict, SSU sites should live on a different chain
    from LSU sites — both reference structures place 16S/18S on chain A/AA
    and 23S/25S on a different chain (DA / chain ``1``)."""
    for ref in (C.ECOLI_REFERENCE_UNITS, C.YEAST_REFERENCE_UNITS):
        ssu_chains = {u.split("|")[2] for site in ref if site.startswith("ssu_") for u in ref[site]}
        lsu_chains = {u.split("|")[2] for site in ref if site.startswith("lsu_") for u in ref[site]}
        assert ssu_chains.isdisjoint(lsu_chains), (
            f"SSU/LSU anchors share a chain: {ssu_chains & lsu_chains}"
        )


def test_reference_units_contain_modified_nucleotides() -> None:
    """Sanity: at least one modified-nucleotide unit is present in each reference,
    so the unit-ID parser (step 7) is exercised against real inputs."""
    flat_ecoli = [u for site in C.ECOLI_REFERENCE_UNITS.values() for u in site]
    assert any("|4OC|" in u for u in flat_ecoli), "expected 4OC modified nt in E. coli reference"
    assert any("|OMG|" in u for u in flat_ecoli), "expected OMG modified nt in E. coli reference"


# ---------------------------------------------------------------------------
# Organellar / classification
# ---------------------------------------------------------------------------


def test_organellar_keywords_are_lowercase_distinct() -> None:
    assert all(k == k.lower() for k in C.ORGANELLAR_KEYWORDS)
    assert len(set(C.ORGANELLAR_KEYWORDS)) == len(C.ORGANELLAR_KEYWORDS)


def test_superkingdom_tiebreak_order_starts_with_bacteria() -> None:
    """Spec §8.3: Bacteria > Eukaryota > Archaea on ties."""
    assert C.SUPERKINGDOM_TIEBREAK_ORDER == ("Bacteria", "Eukaryota", "Archaea")


# ---------------------------------------------------------------------------
# Misc defaults
# ---------------------------------------------------------------------------


def test_default_contact_cutoff_is_5_angstrom() -> None:
    """Spec §10.2 says 5.0 Å."""
    assert C.DEFAULT_CONTACT_CUTOFF_ANGSTROM == 5.0


def test_asl_fragment_threshold_is_30() -> None:
    """Spec §12.1: chain length < 30 → ``**`` LSU state."""
    assert C.ASL_FRAGMENT_MAX_LENGTH == 30


def test_water_is_excluded_by_default() -> None:
    assert "HOH" in C.DEFAULT_LIGAND_EXCLUSIONS


def test_skip_reason_constants_are_unique() -> None:
    reasons = {
        C.SKIP_NMR,
        C.SKIP_ARCHAEAL,
        C.SKIP_PARTIAL_MISSING_SSU_OR_LSU,
        C.SKIP_LOW_PROTEIN_COUNT,
        C.SKIP_AMBIGUOUS_CLASSIFICATION,
        C.SKIP_AMBIGUOUS_RRNA_CORE,
        C.SKIP_INSUFFICIENT_TAXONOMY,
    }
    assert len(reasons) == 7
