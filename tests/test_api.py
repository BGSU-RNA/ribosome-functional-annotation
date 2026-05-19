"""Integration tests for the high-level api orchestrator.

Every external dependency is mocked: RCSB GraphQL (respx), BGSU
correspondence (respx), RCSB file download (respx). The Gemmi structure
served back to the loader is the ``ribosome_fixture`` from
``conftest.py`` — same geometry the step-9 inference tests rely on, so
the chain assignments and tRNA states are deterministic.
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import Any

import gemmi
import httpx
import pytest
import respx

from ribosome_state_annotator import api
from ribosome_state_annotator import constants as C
from ribosome_state_annotator.bgsu_client import BGSU_CORRESPONDENCE_URL
from ribosome_state_annotator.cache import Cache
from ribosome_state_annotator.coordinates import RCSB_ASSEMBLY_DOWNLOAD_TEMPLATE
from ribosome_state_annotator.rcsb_client import RCSB_GRAPHQL_URL
from ribosome_state_annotator.rfam_pdb_region import RFAM_PDB_REGION_URL

FIXTURE_PDB_ID = "RIBOFIXTURE"

# Anchor unit IDs in the synthetic ribosome_fixture, picked to map each
# E. coli §6.3 reference site onto a residue the fixture actually has.
_ANCHOR_BY_SITE: dict[str, str] = {
    "ssu_mrna": f"{FIXTURE_PDB_ID}|1|S|U|10",
    "ssu_atrna": f"{FIXTURE_PDB_ID}|1|S|U|50",
    "ssu_ptrna": f"{FIXTURE_PDB_ID}|1|S|U|30",
    "ssu_etrna": f"{FIXTURE_PDB_ID}|1|S|U|70",
    "lsu_atrna": f"{FIXTURE_PDB_ID}|1|L|U|20",
    "lsu_ptrna": f"{FIXTURE_PDB_ID}|1|L|U|40",
    "lsu_etrna": f"{FIXTURE_PDB_ID}|1|L|U|60",
}


# ---------------------------------------------------------------------------
# Fixture builders
# ---------------------------------------------------------------------------


def _polymer_instance(
    auth: str,
    polymer_type: str,
    description: str,
    *,
    rfam: tuple[str, ...] = (),
    superkingdom: str | None = "Bacteria",
) -> dict[str, Any]:
    return {
        "rcsb_polymer_entity_instance_container_identifiers": {
            "auth_asym_id": auth,
            "asym_id": auth,
            "entity_id": auth,
        },
        "polymer_entity": {
            "entity_poly": {"rcsb_entity_polymer_type": polymer_type},
            "rcsb_polymer_entity": {"pdbx_description": description},
            "rcsb_polymer_entity_annotations": [
                {"annotation_id": acc, "type": "Rfam", "name": "rRNA"} for acc in rfam
            ],
            "rcsb_entity_source_organism": (
                [
                    {
                        "ncbi_taxonomy_id": 562,
                        "ncbi_scientific_name": "Escherichia coli",
                        "ncbi_parent_scientific_name": superkingdom,
                    }
                ]
                if superkingdom is not None
                else []
            ),
            "uniprots": None,
        },
    }


def _bacterial_entry_payload() -> dict[str, Any]:
    """Synthetic RCSB entry matching the chain layout of ``ribosome_fixture``.

    - S, L: rRNA chains (RF00177 / RF02541)
    - M, TA, TP, TE: unmapped RNA candidates (mRNA + tRNAs)
    - EFTU: non-ribosomal protein (elongation factor)
    - L1: ribosomal protein
    - RP1..RP14: extra ribosomal proteins so the §7.4 threshold (15) is met
      and the §8.3 vote has ≥3 voters with Bacteria taxonomy.
    """
    instances = [
        _polymer_instance("S", "RNA", "16S ribosomal RNA", rfam=("RF00177",)),
        _polymer_instance("L", "RNA", "23S ribosomal RNA", rfam=("RF02541",)),
        _polymer_instance("M", "RNA", "mRNA"),
        _polymer_instance("TA", "RNA", "tRNA-Phe"),
        _polymer_instance("TP", "RNA", "tRNA-Lys"),
        _polymer_instance("TE", "RNA", "tRNA-Ala"),
        _polymer_instance("EFTU", "Protein", "Elongation factor Tu"),
        _polymer_instance("L1", "Protein", "50S ribosomal protein L1"),
    ]
    instances.extend(
        _polymer_instance(f"RP{i}", "Protein", f"50S ribosomal protein L{i}") for i in range(2, 16)
    )
    return {
        "rcsb_id": FIXTURE_PDB_ID,
        "exptl": [{"method": "X-RAY DIFFRACTION"}],
        "assemblies": [
            {
                "rcsb_id": f"{FIXTURE_PDB_ID}-1",
                "polymer_entity_instances": instances,
                "nonpolymer_entity_instances": [],
            }
        ],
    }


def _combined_bgsu_alignment(anchor_by_site: dict[str, str]) -> dict[str, Any]:
    """Build one BGSU response containing alignments for every E. coli
    site at once. Each reference unit maps to itself plus the fixture's
    anchor residue, so the §5.2.2 filter keeps the anchor for our
    target PDB."""
    alignment = []
    for site_key, anchor in anchor_by_site.items():
        for ref_unit in C.ECOLI_REFERENCE_UNITS[site_key]:
            alignment.append({"reference_unit": ref_unit, "mapped_units": [ref_unit, anchor]})
    return {"alignment": alignment}


def _ribosome_cif_bytes(structure: gemmi.Structure, tmp_path: Path) -> bytes:
    """Write a Gemmi structure to disk and return the gzipped CIF bytes."""
    cif_path = tmp_path / "fixture.cif"
    structure.make_mmcif_document().write_file(str(cif_path))
    return gzip.compress(cif_path.read_bytes())


def _pdbe_payload(pdb_id: str = FIXTURE_PDB_ID) -> dict[str, Any]:
    """PDBe-shape Rfam mapping for the synthetic fixture's S and L chains.

    Duplicates the RCSB-supplied accessions on those chains; api.augment
    dedupes so the duplication is harmless and exercises the merge path.
    """
    return {
        pdb_id.lower(): {
            "Rfam": {
                "RF00177": {"mappings": [{"chain_id": "S"}]},
                "RF02541": {"mappings": [{"chain_id": "L"}]},
            }
        }
    }


def _install_mocks(
    *,
    entry_payload: dict[str, Any] | None = None,
    bgsu_payload: dict[str, Any] | None = None,
    pdbe_payload: dict[str, Any] | None = None,
    cif_bytes: bytes | None = None,
    pdb_id: str = FIXTURE_PDB_ID,
    assembly_id: str = "1",
) -> dict[str, respx.Route]:
    """Install all respx routes and return them by name.

    The ``pdbe_payload`` argument is preserved for backward compatibility
    with older test signatures but is unused — the package no longer
    queries PDBe's REST endpoint for Rfam annotations (it uses the EBI
    pdb_full_region file instead, mocked here as unavailable so tests
    don't try to fetch it).
    """
    del pdbe_payload  # historical compatibility, no longer used
    entry = entry_payload if entry_payload is not None else _bacterial_entry_payload()
    bgsu = bgsu_payload if bgsu_payload is not None else _combined_bgsu_alignment(_ANCHOR_BY_SITE)
    cif = cif_bytes if cif_bytes is not None else b""

    rcsb_route = respx.post(RCSB_GRAPHQL_URL).mock(
        return_value=httpx.Response(200, json={"data": {"entry": entry}})
    )
    bgsu_route = respx.get(BGSU_CORRESPONDENCE_URL).mock(
        return_value=httpx.Response(200, json=bgsu)
    )
    download_url = RCSB_ASSEMBLY_DOWNLOAD_TEMPLATE.format(
        pdb_id=pdb_id.lower(), assembly_id=assembly_id
    )
    coord_route = respx.get(download_url).mock(return_value=httpx.Response(200, content=cif))
    # Mock the EBI Rfam pdb_full_region URL as a 404 — tests don't
    # exercise the Rfam-file augmentation path; RCSB-supplied Rfam tags
    # on the entry payload are sufficient for classification.
    rfam_route = respx.route(url=RFAM_PDB_REGION_URL).mock(
        return_value=httpx.Response(404)
    )
    return {
        "rcsb": rcsb_route,
        "bgsu": bgsu_route,
        "coord": coord_route,
        "rfam": rfam_route,
    }


# ---------------------------------------------------------------------------
# annotate_pdb — happy path
# ---------------------------------------------------------------------------


@respx.mock
def test_annotate_pdb_happy_path(ribosome_fixture: gemmi.Structure, tmp_path: Path) -> None:
    cif_bytes = _ribosome_cif_bytes(ribosome_fixture, tmp_path)
    _install_mocks(cif_bytes=cif_bytes)

    results = api.annotate_pdb(FIXTURE_PDB_ID, no_cache=True)

    assert len(results) == 1
    ann = results[0]
    assert ann.status == "annotated"
    assert ann.pdb_id == FIXTURE_PDB_ID
    assert ann.assembly_id == "1"
    assert ann.ribosome_classification == "bacterial_ribosome"

    # rRNA partition surfaced canonically.
    assert [c.auth_asym_id for c in ann.ssu_main_rrna_chains] == ["S"]
    assert [c.auth_asym_id for c in ann.lsu_main_rrna_chains] == ["L"]

    # Functional chain assignment matches the fixture's anchor geometry.
    assert ann.mrna_chain is not None and ann.mrna_chain.auth_asym_id == "M"
    assert ann.aminoacyl_trna_chain is not None and ann.aminoacyl_trna_chain.auth_asym_id == "TA"
    assert ann.peptidyl_trna_chain is not None and ann.peptidyl_trna_chain.auth_asym_id == "TP"
    assert ann.exit_trna_chain is not None and ann.exit_trna_chain.auth_asym_id == "TE"
    assert ann.aminoacyl_trna_state == "A/A"
    assert ann.peptidyl_trna_state == "P/P"
    assert ann.exit_trna_state == "E/E"

    # Non-ribosomal proteins exclude any chain flagged is_ribosomal_protein.
    assert [c.auth_asym_id for c in ann.non_ribosomal_proteins] == ["EFTU"]


@respx.mock
def test_annotate_pdb_emits_classification_evidence(
    ribosome_fixture: gemmi.Structure, tmp_path: Path
) -> None:
    cif_bytes = _ribosome_cif_bytes(ribosome_fixture, tmp_path)
    _install_mocks(cif_bytes=cif_bytes)
    ann = api.annotate_pdb(FIXTURE_PDB_ID, no_cache=True)[0]
    ev = ann.classification_evidence
    assert ev["ssu_rfam"] == ["RF00177"]
    assert ev["lsu_rfam"] == ["RF02541"]
    assert ev["rrna_core"] == "bacterial_like"
    assert ev["dominant_ribosomal_protein_superkingdom"] == "Bacteria"
    assert ev["rule"] == C.CLASSIFICATION_RULE_BACTERIAL


# ---------------------------------------------------------------------------
# annotate_pdb — skip paths
# ---------------------------------------------------------------------------


@respx.mock
def test_annotate_pdb_nmr_entry_skips_all_assemblies() -> None:
    entry = _bacterial_entry_payload()
    entry["exptl"] = [{"method": "SOLUTION NMR"}]
    _install_mocks(entry_payload=entry)

    results = api.annotate_pdb(FIXTURE_PDB_ID, no_cache=True)
    assert len(results) == 1
    assert results[0].status == "skipped"
    assert results[0].skip_reason == C.SKIP_NMR
    assert results[0].assembly_id is None


@respx.mock
def test_annotate_pdb_partial_assembly_skips() -> None:
    """An assembly with SSU but no LSU rRNA → partial_ribosome_missing_ssu_or_lsu."""
    entry = _bacterial_entry_payload()
    # Drop the LSU rRNA instance.
    entry["assemblies"][0]["polymer_entity_instances"] = [
        inst
        for inst in entry["assemblies"][0]["polymer_entity_instances"]
        if inst["rcsb_polymer_entity_instance_container_identifiers"]["auth_asym_id"] != "L"
    ]
    _install_mocks(entry_payload=entry)

    results = api.annotate_pdb(FIXTURE_PDB_ID, no_cache=True)
    assert len(results) == 1
    assert results[0].status == "skipped"
    assert results[0].skip_reason == C.SKIP_PARTIAL_MISSING_SSU_OR_LSU


@respx.mock
def test_annotate_pdb_rcsb_failure_returns_failed_annotation() -> None:
    respx.post(RCSB_GRAPHQL_URL).mock(return_value=httpx.Response(500))
    results = api.annotate_pdb(FIXTURE_PDB_ID, no_cache=True)
    assert len(results) == 1
    assert results[0].status == "failed"
    assert "rcsb_fetch_failed" in (results[0].skip_reason or "")


# ---------------------------------------------------------------------------
# annotate_assembly (single-assembly wrapper)
# ---------------------------------------------------------------------------


@respx.mock
def test_annotate_assembly_returns_single_annotation(
    ribosome_fixture: gemmi.Structure, tmp_path: Path
) -> None:
    cif_bytes = _ribosome_cif_bytes(ribosome_fixture, tmp_path)
    _install_mocks(cif_bytes=cif_bytes)
    ann = api.annotate_assembly(FIXTURE_PDB_ID, "1", no_cache=True)
    assert ann.assembly_id == "1"
    assert ann.status == "annotated"


@respx.mock
def test_annotate_assembly_missing_assembly_returns_failed() -> None:
    _install_mocks()
    ann = api.annotate_assembly(FIXTURE_PDB_ID, "99", no_cache=True)
    assert ann.status == "failed"
    assert "assembly_not_found" in (ann.skip_reason or "")


def test_demote_no_contact_fragments_clears_double_star_states() -> None:
    """When a tRNA state is `**/**` (fragment with no anchor contact on
    either subunit) the safeguard demotes the chain to unmapped so it
    doesn't claim a tRNA role."""
    from ribosome_state_annotator.infer import ChainAssignments, TRNAStates
    from ribosome_state_annotator.models import ChainRef

    fragment = ChainRef(
        pdb_id="TEST",
        assembly_id="1",
        auth_asym_id="X",
        polymer_type="RNA",
    )
    assignments = ChainAssignments(aminoacyl_trna_chain=fragment)
    states = TRNAStates(aminoacyl_trna_state="**/**")
    new_assignments, new_states = api._demote_no_contact_fragments(assignments, states)
    assert new_assignments.aminoacyl_trna_chain is None
    assert new_states.aminoacyl_trna_state is None


def test_demote_no_contact_fragments_preserves_real_assignments() -> None:
    """A canonical `A/A` state is untouched by the safeguard."""
    from ribosome_state_annotator.infer import ChainAssignments, TRNAStates
    from ribosome_state_annotator.models import ChainRef

    real = ChainRef(
        pdb_id="TEST",
        assembly_id="1",
        auth_asym_id="TA",
        polymer_type="RNA",
    )
    assignments = ChainAssignments(aminoacyl_trna_chain=real)
    states = TRNAStates(aminoacyl_trna_state="A/A")
    new_assignments, new_states = api._demote_no_contact_fragments(assignments, states)
    assert new_assignments.aminoacyl_trna_chain is real
    assert new_states.aminoacyl_trna_state == "A/A"


# ---------------------------------------------------------------------------
# annotate_many
# ---------------------------------------------------------------------------


@respx.mock
def test_annotate_many_aggregates(ribosome_fixture: gemmi.Structure, tmp_path: Path) -> None:
    cif_bytes = _ribosome_cif_bytes(ribosome_fixture, tmp_path)
    _install_mocks(cif_bytes=cif_bytes)
    # Single PDB called twice → 2 annotations.
    results = api.annotate_many(
        [FIXTURE_PDB_ID, FIXTURE_PDB_ID],
        no_cache=True,
    )
    assert len(results) == 2
    for ann in results:
        assert ann.status == "annotated"


@respx.mock
def test_annotate_pdb_attaches_large_scale_movements_when_dataset_matches(
    ribosome_fixture: gemmi.Structure, tmp_path: Path
) -> None:
    """A RADdb dataset whose key matches the fixture's (pdb, lsu, ssu)
    should populate ``large_scale_movements`` with non-null metrics."""
    from datetime import datetime, timezone

    from ribosome_state_annotator.models import LargeScaleMovements
    from ribosome_state_annotator.raddb import (
        RADdbDataset,
        RADdbMetadata,
        build_raddb_lookup,
    )

    cif_bytes = _ribosome_cif_bytes(ribosome_fixture, tmp_path)
    _install_mocks(cif_bytes=cif_bytes)

    lookup, dupes = build_raddb_lookup(
        [
            {
                "RCSB": FIXTURE_PDB_ID,
                "LSU chain ID": "L",
                "SSU chain ID": "S",
                "body rot.": "5.8",
                "head rot.": "9.4",
            }
        ]
    )
    dataset = RADdbDataset(
        metadata=RADdbMetadata(
            source_url="https://example.test/RADdb.20260508.LSUSSU.csv",
            downloaded_at=datetime(2026, 5, 8, tzinfo=timezone.utc),
            rad_date="20260508",
        ),
        lookup=lookup,
        duplicate_keys=dupes,
    )

    ann = api.annotate_pdb(FIXTURE_PDB_ID, no_cache=True, raddb_dataset=dataset)[0]
    assert isinstance(ann.large_scale_movements, LargeScaleMovements)
    assert ann.large_scale_movements.rad_date == "20260508"
    assert ann.large_scale_movements.intersubunit_rotation == 5.8
    assert ann.large_scale_movements.ssu_head_rotation == 9.4


@respx.mock
def test_annotate_pdb_no_raddb_emits_null_metrics(
    ribosome_fixture: gemmi.Structure, tmp_path: Path
) -> None:
    """``no_raddb=True`` skips integration entirely; the JSON block still
    appears with ``rad_date=None`` for schema stability."""
    cif_bytes = _ribosome_cif_bytes(ribosome_fixture, tmp_path)
    _install_mocks(cif_bytes=cif_bytes)

    ann = api.annotate_pdb(FIXTURE_PDB_ID, no_cache=True, no_raddb=True)[0]
    assert ann.large_scale_movements is not None
    assert ann.large_scale_movements.rad_date is None
    assert ann.large_scale_movements.intersubunit_rotation is None
    assert ann.large_scale_movements.ssu_head_rotation is None


@respx.mock
def test_annotate_many_continue_on_error_catches_per_entry_failures() -> None:
    """When ``continue_on_error=True``, a per-PDB error becomes a failed
    annotation rather than aborting the whole batch."""
    # RCSB fails for any request.
    respx.post(RCSB_GRAPHQL_URL).mock(side_effect=httpx.ConnectError("down"))
    results = api.annotate_many(
        ["AAAA", "BBBB"],
        continue_on_error=True,
        no_cache=True,
    )
    assert len(results) == 2
    for ann in results:
        assert ann.status == "failed"


# ---------------------------------------------------------------------------
# Cache wiring
# ---------------------------------------------------------------------------


@respx.mock
def test_rcsb_cache_avoids_second_http_call(
    ribosome_fixture: gemmi.Structure, tmp_path: Path
) -> None:
    cif_bytes = _ribosome_cif_bytes(ribosome_fixture, tmp_path)
    routes = _install_mocks(cif_bytes=cif_bytes)
    cache = Cache(tmp_path / "shared-cache")

    api.annotate_pdb(FIXTURE_PDB_ID, cache=cache)
    assert routes["rcsb"].call_count == 1
    # Per-subunit BGSU batching: 1 call for the four SSU sites concatenated,
    # 1 call for the three LSU sites concatenated.
    assert routes["bgsu"].call_count == 2
    assert routes["coord"].call_count == 1

    api.annotate_pdb(FIXTURE_PDB_ID, cache=cache)
    # All caches should hit on the second call — no new HTTP traffic.
    assert routes["rcsb"].call_count == 1
    assert routes["bgsu"].call_count == 2
    assert routes["coord"].call_count == 1


def test_resolve_cache_no_cache_returns_none() -> None:
    assert api._resolve_cache(None, None, no_cache=True) is None


def test_resolve_cache_explicit_takes_precedence(tmp_path: Path) -> None:
    explicit = Cache(tmp_path / "explicit")
    result = api._resolve_cache(explicit, tmp_path / "dir", no_cache=False)
    assert result is explicit


def test_resolve_cache_dir_builds_cache(tmp_path: Path) -> None:
    result = api._resolve_cache(None, tmp_path / "custom-dir", no_cache=False)
    assert result is not None
    assert result.root == tmp_path / "custom-dir"


# ---------------------------------------------------------------------------
# Chain-substitution map builder (§5.2.2 multi-assembly fallback)
# ---------------------------------------------------------------------------


def test_build_chain_substitution_maps_ssu_and_lsu() -> None:
    """For E. coli references → another bacterial assembly, the map should
    translate reference chains (AA for SSU, DA for LSU) into the target's
    SSU/LSU chains via the by_role partition."""
    from ribosome_state_annotator.models import ChainRef

    by_role: dict[str, list[ChainRef]] = {
        "ssu_main_rrna": [
            ChainRef(pdb_id="5J7L", assembly_id="2", auth_asym_id="BA", rfam_accessions=["RF00177"])
        ],
        "lsu_main_rrna": [
            ChainRef(pdb_id="5J7L", assembly_id="2", auth_asym_id="CA", rfam_accessions=["RF02541"])
        ],
    }
    sub = api._build_chain_substitution(C.ECOLI_REFERENCE_UNITS, by_role)
    # E. coli reference uses AA (SSU) and DA (LSU); target assembly uses BA / CA.
    assert sub["AA"] == "BA"
    assert sub["DA"] == "CA"


def test_build_chain_substitution_skips_role_with_multiple_chains() -> None:
    """If the target has multiple SSU rRNA chains, substitution is ambiguous
    — skip the SSU entry rather than picking arbitrarily."""
    from ribosome_state_annotator.models import ChainRef

    by_role: dict[str, list[ChainRef]] = {
        "ssu_main_rrna": [
            ChainRef(pdb_id="X", assembly_id="1", auth_asym_id="A", rfam_accessions=["RF00177"]),
            ChainRef(pdb_id="X", assembly_id="1", auth_asym_id="B", rfam_accessions=["RF00177"]),
        ],
        "lsu_main_rrna": [
            ChainRef(pdb_id="X", assembly_id="1", auth_asym_id="C", rfam_accessions=["RF02541"]),
        ],
    }
    sub = api._build_chain_substitution(C.ECOLI_REFERENCE_UNITS, by_role)
    assert "AA" not in sub  # SSU skipped due to ambiguity
    assert sub["DA"] == "C"


def test_build_chain_substitution_empty_when_no_rrna() -> None:
    sub = api._build_chain_substitution(C.ECOLI_REFERENCE_UNITS, {})
    assert sub == {}


def test_bgsu_cache_key_includes_depth() -> None:
    """Cache key must reflect the depth parameter so a wider-depth request
    doesn't incorrectly hit a narrower-depth cache entry."""
    k1 = api._bgsu_cache_key(["5J7L|1|AA|G|926"], "Rfam", "20.0A", 700)
    k2 = api._bgsu_cache_key(["5J7L|1|AA|G|926"], "Rfam", "20.0A", 1500)
    assert k1 != k2
    assert "depth=700" in k1
    assert "depth=1500" in k2


# ---------------------------------------------------------------------------
# Reference selection (§6.4)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("classification", "expected_pdb"),
    [
        ("bacterial_ribosome", "5J7L"),
        ("eukaryotic_organellar_ribosome", "5J7L"),  # filtered E. coli anchors
        ("eukaryotic_ribosome", "7ZW0"),
    ],
)
def test_reference_selection(classification: str, expected_pdb: str) -> None:
    refs = api._select_reference_units(classification)
    # Every reference unit ID should start with the expected reference PDB.
    # (Sites may be empty for organellar — ssu_etrna has no mt-friendly anchor.)
    for site_units in refs.values():
        for unit_id in site_units:
            assert unit_id.startswith(f"{expected_pdb}|"), (
                f"{classification} site referenced {unit_id}, expected {expected_pdb}"
            )


def test_reference_selection_unknown_classification_returns_empty() -> None:
    assert api._select_reference_units("archaeal_ribosome") == {}
