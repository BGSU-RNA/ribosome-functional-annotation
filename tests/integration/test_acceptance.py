"""Acceptance tests (spec §25.1) — eight PDB entries via the live APIs.

Marked ``network`` via :mod:`tests.integration.conftest` and excluded
from the default ``pytest`` run. To execute:

::

    pytest -m network tests/integration/test_acceptance.py

Each test asserts the **shape** of the result (status, classification,
chain-assignment count) rather than specific chain letter codes, per
spec §25.2 — the BGSU correspondence and RCSB metadata evolve and
exact IDs may shift. Tightening any assertion to a specific chain or
state string is a one-line change once the live behaviour is verified.
"""

from __future__ import annotations

import pytest

from ribosome_state_annotator import api
from ribosome_state_annotator import constants as C
from ribosome_state_annotator.models import RibosomeAnnotation


def _only(results: list[RibosomeAnnotation]) -> RibosomeAnnotation:
    assert len(results) == 1, f"expected exactly one result, got {len(results)}: " + ", ".join(
        f"{r.pdb_id}/{r.assembly_id}={r.status}" for r in results
    )
    return results[0]


def _by_assembly_id(
    results: list[RibosomeAnnotation], assembly_id: str
) -> RibosomeAnnotation:
    """Pick out one specific assembly from a multi-assembly entry's results."""
    matches = [r for r in results if r.assembly_id == assembly_id]
    assert matches, f"no assembly {assembly_id} in: " + ", ".join(
        f"{r.pdb_id}/{r.assembly_id}={r.status}" for r in results
    )
    return matches[0]


# ---------------------------------------------------------------------------
# Cases 1-4: annotated assemblies
# ---------------------------------------------------------------------------


def test_case_1_5j7l_bacterial_classic() -> None:
    """Case 1: *E. coli* 70S — bacterial_ribosome classification.

    Spec §25.1 describes 5J7L as carrying mRNA + A/P/E tRNAs, but the
    actual deposited 5J7L is the apo (no-substrate) high-resolution
    structure; the prototype's legacy CSV correctly shows no functional
    chains assigned. We assert only the classification + rRNA shape,
    which is the part of the spec's claim that the data actually supports.
    """
    results = api.annotate_pdb("5J7L", no_cache=True)
    assert len(results) >= 1
    assembly_1 = _by_assembly_id(results, "1")
    assert assembly_1.status == "annotated"
    assert assembly_1.ribosome_classification == "bacterial_ribosome"
    assert len(assembly_1.ssu_main_rrna_chains) == 1
    assert len(assembly_1.lsu_main_rrna_chains) == 1


def test_case_2_5tbw_yeast_80s() -> None:
    """Case 2: Yeast cytoplasmic 80S, eukaryotic_ribosome."""
    results = api.annotate_pdb("5TBW", no_cache=True)
    assert len(results) >= 1
    assembly_1 = _by_assembly_id(results, "1")
    assert assembly_1.status == "annotated"
    assert assembly_1.ribosome_classification == "eukaryotic_ribosome"
    # Yeast 80S has 5.8S — eukaryotic ribosomes populate lsu_associated_rrna_chains.
    assert len(assembly_1.lsu_associated_rrna_chains) >= 1


def test_case_3_6zmi_human_80s() -> None:
    """Case 3: Human cytoplasmic 80S — multi-character chain IDs (S2, L5)."""
    results = api.annotate_pdb("6ZMI", no_cache=True)
    assert len(results) >= 1
    assembly_1 = _by_assembly_id(results, "1")
    assert assembly_1.status == "annotated"
    assert assembly_1.ribosome_classification == "eukaryotic_ribosome"


def test_case_4_6zm6_mitoribosome_organellar() -> None:
    """Case 4: Human mitoribosome — bacterial-like rRNA + Eukaryota proteins
    classifies as `eukaryotic_organellar_ribosome` (§8.5)."""
    results = api.annotate_pdb("6ZM6", no_cache=True)
    assert len(results) >= 1
    assembly_1 = _by_assembly_id(results, "1")
    assert assembly_1.status == "annotated"
    assert assembly_1.ribosome_classification == "eukaryotic_organellar_ribosome"


# ---------------------------------------------------------------------------
# Cases 5-7: structured skip results
# ---------------------------------------------------------------------------


def test_case_5_isolated_ssu_skipped() -> None:
    """Case 5: a bacterial 30S-only deposit. See tests/fixtures/acceptance_set.md
    for the PDB-pick rationale."""
    result = _only(api.annotate_pdb("1J5E", no_cache=True))
    assert result.status == "skipped"
    assert result.skip_reason == C.SKIP_PARTIAL_MISSING_SSU_OR_LSU


def test_case_6_isolated_lsu_skipped() -> None:
    """Case 6: bacterial 50S-only deposit."""
    result = _only(api.annotate_pdb("1NKW", no_cache=True))
    assert result.status == "skipped"
    assert result.skip_reason == C.SKIP_PARTIAL_MISSING_SSU_OR_LSU


def test_case_7_archaeal_skipped() -> None:
    """Case 7: an archaeal 70S — `archaeal_ribosome_not_supported`.

    Uses *Pyrococcus furiosus* 70S (4V6U), which carries both archaeal
    SSU rRNA (Rfam RF01959) and archaeal LSU rRNA (Rfam RF02540). These
    Rfam families are in ``RFAM_ROLE_MAP`` so the SSU/LSU detection
    succeeds; the §8.3 protein-superkingdom vote then identifies the
    Archaea taxonomy and the §8.5 truth table emits the archaeal skip.
    """
    results = api.annotate_pdb("4V6U", no_cache=True)
    assert len(results) >= 1
    for r in results:
        assert r.status == "skipped"
        assert r.skip_reason == C.SKIP_ARCHAEAL


# ---------------------------------------------------------------------------
# Case 8: multi-assembly entry
# ---------------------------------------------------------------------------


def test_case_8_5fdv_multi_assembly() -> None:
    """Case 8: an entry whose biological assemblies each carry one ribosome.
    Should produce **two** annotated `RibosomeAnnotation` objects."""
    results = api.annotate_pdb("5FDV", no_cache=True)
    assert len(results) == 2
    assert {r.assembly_id for r in results} == {"1", "2"}
    for r in results:
        # Both assemblies should annotate, not skip.
        assert r.status == "annotated", f"5FDV assembly {r.assembly_id} skipped: {r.skip_reason}"


# ---------------------------------------------------------------------------
# Smoke test — quickest live check
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("pdb_id", ["5J7L"])
def test_smoke_annotate_pdb_returns_at_least_one_result(pdb_id: str) -> None:
    """Bare-minimum live check: annotate_pdb returns something non-empty."""
    results = api.annotate_pdb(pdb_id, no_cache=True)
    assert len(results) >= 1
    assert all(isinstance(r, RibosomeAnnotation) for r in results)
