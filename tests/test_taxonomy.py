"""Unit tests for :mod:`ribosome_state_annotator.taxonomy` (spec §34).

Synthetic ChainRef fixtures with hand-crafted lineages exercise the
strict-LCA → majority-mode-constrained-to-branch aggregator without any
RCSB I/O.
"""

from __future__ import annotations

import pytest

from ribosome_state_annotator.models import (
    AssemblyTaxonomy,
    ChainRef,
    TaxonNode,
)
from ribosome_state_annotator.taxonomy import aggregate_assembly_lineage


# ---------------------------------------------------------------------------
# Lineage fixtures — pinned to NCBI tax_ids returned by the live RCSB
# GraphQL schema (verified manually before adding here, see commit
# message). These need not match exactly across NCBI revisions, but the
# tax_ids picked are stable canonical ones.
# ---------------------------------------------------------------------------

_ECOLI_LINEAGE = (
    TaxonNode(tax_id=131567, name="cellular organisms", depth=1),
    TaxonNode(tax_id=2, name="Bacteria", depth=2),
    TaxonNode(tax_id=3379134, name="Pseudomonadati", depth=3),
    TaxonNode(tax_id=1224, name="Pseudomonadota", depth=4),
    TaxonNode(tax_id=1236, name="Gammaproteobacteria", depth=5),
    TaxonNode(tax_id=91347, name="Enterobacterales", depth=6),
    TaxonNode(tax_id=543, name="Enterobacteriaceae", depth=7),
    TaxonNode(tax_id=561, name="Escherichia", depth=8),
    TaxonNode(tax_id=562, name="Escherichia coli", depth=9),
    TaxonNode(tax_id=83333, name="Escherichia coli K-12", depth=10),
)

_SALMONELLA_LINEAGE = (
    TaxonNode(tax_id=131567, name="cellular organisms", depth=1),
    TaxonNode(tax_id=2, name="Bacteria", depth=2),
    TaxonNode(tax_id=3379134, name="Pseudomonadati", depth=3),
    TaxonNode(tax_id=1224, name="Pseudomonadota", depth=4),
    TaxonNode(tax_id=1236, name="Gammaproteobacteria", depth=5),
    TaxonNode(tax_id=91347, name="Enterobacterales", depth=6),
    TaxonNode(tax_id=543, name="Enterobacteriaceae", depth=7),
    TaxonNode(tax_id=590, name="Salmonella", depth=8),
    TaxonNode(tax_id=28901, name="Salmonella enterica", depth=9),
)

_BACILLUS_LINEAGE = (
    TaxonNode(tax_id=131567, name="cellular organisms", depth=1),
    TaxonNode(tax_id=2, name="Bacteria", depth=2),
    TaxonNode(tax_id=1783272, name="Bacillati", depth=3),
    TaxonNode(tax_id=1239, name="Bacillota", depth=4),
    TaxonNode(tax_id=91061, name="Bacilli", depth=5),
    TaxonNode(tax_id=1385, name="Bacillales", depth=6),
    TaxonNode(tax_id=186817, name="Bacillaceae", depth=7),
    TaxonNode(tax_id=1386, name="Bacillus", depth=8),
    TaxonNode(tax_id=1423, name="Bacillus subtilis", depth=9),
)

_HUMAN_LINEAGE = (
    TaxonNode(tax_id=131567, name="cellular organisms", depth=1),
    TaxonNode(tax_id=2759, name="Eukaryota", depth=2),
    TaxonNode(tax_id=33154, name="Opisthokonta", depth=3),
    TaxonNode(tax_id=33208, name="Metazoa", depth=4),
    TaxonNode(tax_id=9606, name="Homo sapiens", depth=31),
)


def _chain(
    pdb_id: str,
    auth: str,
    lineage: tuple[TaxonNode, ...],
    *,
    tax_id: int | None = None,
    scientific_name: str | None = None,
    rfam: tuple[str, ...] = (),
) -> ChainRef:
    """Build a synthetic ChainRef with the given lineage."""
    return ChainRef(
        pdb_id=pdb_id,
        assembly_id="1",
        auth_asym_id=auth,
        polymer_type="RNA",
        rfam_accessions=list(rfam),
        tax_id=tax_id if tax_id is not None else (lineage[-1].tax_id if lineage else None),
        scientific_name=scientific_name
        if scientific_name is not None
        else (lineage[-1].name if lineage else None),
        taxonomy_lineage=lineage,
    )


# ---------------------------------------------------------------------------
# Single-organism happy path
# ---------------------------------------------------------------------------


def test_single_organism_resolves_full_lineage() -> None:
    rrna = [
        _chain("5J7L", "AA", _ECOLI_LINEAGE),
        _chain("5J7L", "DA", _ECOLI_LINEAGE),
        _chain("5J7L", "BA", _ECOLI_LINEAGE),
    ]
    result = aggregate_assembly_lineage(rrna)
    assert result is not None
    assert result.domain == "Bacteria"
    assert result.species == "Escherichia coli K-12"
    assert result.is_mixed is False
    assert result.strict_lca_depth == _ECOLI_LINEAGE[-1].depth
    assert len(result.lineage) == len(_ECOLI_LINEAGE)
    assert tuple(n.tax_id for n in result.lineage) == tuple(n.tax_id for n in _ECOLI_LINEAGE)
    assert result.n_voting_chains == 3
    assert result.voting_chains == ("5J7L|1|AA", "5J7L|1|DA", "5J7L|1|BA")


def test_eukaryote_returns_eukaryota_domain() -> None:
    rrna = [_chain("9B0S", "S2", _HUMAN_LINEAGE)]
    result = aggregate_assembly_lineage(rrna)
    assert result is not None
    assert result.domain == "Eukaryota"
    assert result.species == "Homo sapiens"


# ---------------------------------------------------------------------------
# No-lineage / all-synthetic
# ---------------------------------------------------------------------------


def test_all_chains_lack_lineage_returns_none() -> None:
    rrna = [
        _chain("FAKE", "A", ()),
        _chain("FAKE", "B", ()),
    ]
    assert aggregate_assembly_lineage(rrna) is None


def test_empty_voting_set_returns_none() -> None:
    assert aggregate_assembly_lineage([]) is None


# ---------------------------------------------------------------------------
# Strict-LCA + majority-mode hybrid (the chimera case)
# ---------------------------------------------------------------------------


def test_strain_disagreement_strict_to_species_then_majority() -> None:
    """5 chains *E. coli K-12* + 4 chains *E. coli* (no K-12 strain row).

    Strict-LCA reaches species (depth 9), majority extends to strain.
    """
    ecoli_no_strain = _ECOLI_LINEAGE[:-1]  # drop K-12
    rrna = [
        *(_chain("X", f"A{i}", _ECOLI_LINEAGE) for i in range(5)),
        *(_chain("X", f"B{i}", ecoli_no_strain) for i in range(4)),
    ]
    result = aggregate_assembly_lineage(rrna)
    assert result is not None
    assert result.domain == "Bacteria"
    assert result.species == "Escherichia coli K-12"
    assert result.is_mixed is True
    # Strict reaches Escherichia coli (depth 9); strain depth 10 is majority.
    assert result.strict_lca_depth == 9
    assert len(result.lineage) == 10


def test_two_species_same_family_strict_to_family_majority_to_dominant_species() -> None:
    """5 E. coli + 4 Salmonella — both in Enterobacteriaceae.

    Strict-LCA stops at Enterobacteriaceae (depth 7). Majority continues
    on Escherichia branch since 5 > 4.
    """
    rrna = [
        *(_chain("X", f"E{i}", _ECOLI_LINEAGE) for i in range(5)),
        *(_chain("X", f"S{i}", _SALMONELLA_LINEAGE) for i in range(4)),
    ]
    result = aggregate_assembly_lineage(rrna)
    assert result is not None
    assert result.is_mixed is True
    assert result.strict_lca_depth == 7  # Enterobacteriaceae
    assert result.species == "Escherichia coli K-12"
    # Lineage must end in K-12 (the majority branch), NOT mix Salmonella nodes.
    leaf_ids = [n.tax_id for n in result.lineage]
    assert 590 not in leaf_ids  # Salmonella genus excluded
    assert 561 in leaf_ids  # Escherichia included


def test_three_branch_chimera_constrains_to_strict_subtree() -> None:
    """5 E. coli + 4 Salmonella + 4 Bacillus — three branches.

    Strict-LCA stops at Bacteria (depth 2). Majority picks
    Pseudomonadati (9/13) at depth 3. Subsequent depths constrained to
    that subtree, so Bacillus chains (4) drop out and the depth-8
    contest becomes 5 Escherichia vs 4 Salmonella → Escherichia wins.
    The returned lineage must be a valid root-to-leaf path through the
    Pseudomonadati subtree, NOT a mix of Pseudomonadati and Bacillati
    descendants.
    """
    rrna = [
        *(_chain("X", f"E{i}", _ECOLI_LINEAGE) for i in range(5)),
        *(_chain("X", f"S{i}", _SALMONELLA_LINEAGE) for i in range(4)),
        *(_chain("X", f"B{i}", _BACILLUS_LINEAGE) for i in range(4)),
    ]
    result = aggregate_assembly_lineage(rrna)
    assert result is not None
    assert result.is_mixed is True
    assert result.strict_lca_depth == 2  # Bacteria
    # The lineage should follow the Pseudomonadati subtree (Bacillati
    # had only 4 chains vs Pseudomonadati's 9).
    names = [n.name for n in result.lineage]
    assert "Pseudomonadati" in names
    assert "Bacillati" not in names
    assert "Bacillus" not in names
    # And E. coli wins the depth-8 contest within Pseudomonadati.
    assert result.species == "Escherichia coli K-12"
    assert result.domain == "Bacteria"


# ---------------------------------------------------------------------------
# Source organism collection
# ---------------------------------------------------------------------------


def test_source_organisms_dedupes_across_voting_and_non_voting_chains() -> None:
    rrna = [_chain("5J7L", "AA", _ECOLI_LINEAGE)]
    non_voting = [
        _chain("5J7L", "V", (), tax_id=4932, scientific_name="Saccharomyces cerevisiae"),
        _chain("5J7L", "W", (), tax_id=4932, scientific_name="Saccharomyces cerevisiae"),
        _chain("5J7L", "X", ()),  # synthetic: tax_id=None, scientific_name=None
    ]
    result = aggregate_assembly_lineage(rrna, all_chains=rrna + non_voting)
    assert result is not None
    organisms = {(o.tax_id, o.scientific_name) for o in result.source_organisms}
    assert (83333, "Escherichia coli K-12") in organisms
    assert (4932, "Saccharomyces cerevisiae") in organisms
    assert (None, None) in organisms
    # The yeast chain appears twice but the organism is reported once.
    assert sum(1 for o in result.source_organisms if o.tax_id == 4932) == 1


# ---------------------------------------------------------------------------
# AssemblyTaxonomy as serialisable model
# ---------------------------------------------------------------------------


def test_assembly_taxonomy_round_trips_through_model_dump() -> None:
    rrna = [_chain("9B0S", "S2", _HUMAN_LINEAGE)]
    result = aggregate_assembly_lineage(rrna)
    assert result is not None
    dumped = result.model_dump()
    rebuilt = AssemblyTaxonomy.model_validate(dumped)
    assert rebuilt == result
