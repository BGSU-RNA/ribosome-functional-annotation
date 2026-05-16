"""Unit tests for the RCSB client (spec §5.1).

The HTTP layer is mocked with :mod:`respx`; no test in this file makes a
live network call. Integration tests against the real RCSB endpoint live
in ``tests/integration/`` and are marked with the ``network`` marker.
"""

from __future__ import annotations

import json
from typing import Any

import httpx
import pytest
import respx

from ribosome_state_annotator import rcsb_client
from ribosome_state_annotator.exceptions import ApiRequestError

# ---------------------------------------------------------------------------
# Fixture builders
# ---------------------------------------------------------------------------


def _polymer_instance(
    auth_asym: str,
    *,
    asym: str | None = None,
    entity_id: str | None = "1",
    polymer_type: str = "RNA",
    description: str | None = None,
    rfam: tuple[str, ...] = (),
    sources: tuple[dict[str, Any], ...] = (),
    uniprot_name: str | None = None,
) -> dict[str, Any]:
    """Build a polymer_entity_instance payload matching the RCSB GraphQL shape."""
    rfam_annotations = [
        {
            "annotation_id": accession,
            "type": "Rfam",
            "name": f"Rfam {accession}",
            "description": None,
        }
        for accession in rfam
    ]
    uniprots: list[dict[str, Any]] | None
    if uniprot_name is None:
        uniprots = None
    else:
        uniprots = [{"rcsb_uniprot_protein": {"name": {"value": uniprot_name}}}]
    return {
        "rcsb_polymer_entity_instance_container_identifiers": {
            "auth_asym_id": auth_asym,
            "asym_id": asym if asym is not None else auth_asym,
            "entity_id": entity_id,
        },
        "polymer_entity": {
            "rcsb_id": f"5J7L_{entity_id}",
            "entity_poly": {"rcsb_entity_polymer_type": polymer_type},
            "rcsb_polymer_entity": {"pdbx_description": description},
            "rcsb_polymer_entity_annotations": rfam_annotations,
            "rcsb_entity_source_organism": list(sources),
            "uniprots": uniprots,
        },
    }


def _nonpolymer_instance(
    auth_asym: str,
    comp_id: str,
    *,
    name: str | None = None,
    drugbank: tuple[str, str] | None = None,
) -> dict[str, Any]:
    annotations: list[dict[str, Any]] = []
    if drugbank is not None:
        annotations.append(
            {
                "type": "DrugBank",
                "annotation_id": drugbank[0],
                "description": drugbank[1],
            }
        )
    return {
        "rcsb_nonpolymer_entity_instance_container_identifiers": {
            "auth_asym_id": auth_asym,
            "comp_id": comp_id,
        },
        "nonpolymer_entity": {
            "pdbx_entity_nonpoly": {"name": name},
            "rcsb_nonpolymer_entity_annotation": annotations,
        },
    }


def _entry(
    *,
    pdb_id: str = "5J7L",
    methods: tuple[str, ...] = ("X-RAY DIFFRACTION",),
    assemblies: tuple[dict[str, Any], ...] = (),
) -> dict[str, Any]:
    return {
        "rcsb_id": pdb_id,
        "exptl": [{"method": m} for m in methods],
        "assemblies": list(assemblies),
    }


def _assembly(
    n: int,
    *,
    polymer_instances: tuple[dict[str, Any], ...] = (),
    nonpolymer_instances: tuple[dict[str, Any], ...] = (),
    pdb_id: str = "5J7L",
) -> dict[str, Any]:
    return {
        "rcsb_id": f"{pdb_id}-{n}",
        "polymer_entity_instances": list(polymer_instances),
        "nonpolymer_entity_instances": list(nonpolymer_instances),
    }


# The §5.1.2 example response, verbatim.
SPEC_5_1_2_ENTRY: dict[str, Any] = {
    "rcsb_id": "5J7L",
    "exptl": [{"method": "X-RAY DIFFRACTION"}],
    "assemblies": [
        {
            "rcsb_id": "5J7L-1",
            "polymer_entity_instances": [
                {
                    "rcsb_polymer_entity_instance_container_identifiers": {
                        "auth_asym_id": "AA",
                        "asym_id": "AA",
                        "entity_id": "1",
                    },
                    "polymer_entity": {
                        "entity_poly": {"rcsb_entity_polymer_type": "RNA"},
                        "rcsb_polymer_entity": {"pdbx_description": "16S ribosomal RNA"},
                        "rcsb_polymer_entity_annotations": [
                            {
                                "annotation_id": "RF00177",
                                "type": "Rfam",
                                "name": "Bacterial small subunit ribosomal RNA",
                            }
                        ],
                        "rcsb_entity_source_organism": [
                            {
                                "ncbi_taxonomy_id": 562,
                                "ncbi_scientific_name": "Escherichia coli",
                                "ncbi_parent_scientific_name": "Bacteria",
                            }
                        ],
                        "uniprots": None,
                    },
                }
            ],
            "nonpolymer_entity_instances": [
                {
                    "rcsb_nonpolymer_entity_instance_container_identifiers": {
                        "auth_asym_id": "AA",
                        "comp_id": "MG",
                    },
                    "nonpolymer_entity": {
                        "pdbx_entity_nonpoly": {"name": "MAGNESIUM ION"},
                        "rcsb_nonpolymer_entity_annotation": [],
                    },
                }
            ],
        }
    ],
}


# ---------------------------------------------------------------------------
# parse_assemblies — verbatim §5.1.2 example
# ---------------------------------------------------------------------------


def test_parse_spec_5_1_2_response_yields_one_assembly() -> None:
    assemblies = rcsb_client.parse_assemblies(SPEC_5_1_2_ENTRY)
    assert len(assemblies) == 1
    a = assemblies[0]
    assert a.pdb_id == "5J7L"
    assert a.assembly_id == "1"
    assert a.experimental_methods == ["X-RAY DIFFRACTION"]
    assert len(a.rna_chains) == 1
    assert a.protein_chains == []
    assert len(a.ligands) == 1


def test_parse_spec_5_1_2_response_rna_chain_fields() -> None:
    chain = rcsb_client.parse_assemblies(SPEC_5_1_2_ENTRY)[0].rna_chains[0]
    assert chain.pdb_id == "5J7L"
    assert chain.assembly_id == "1"
    assert chain.auth_asym_id == "AA"
    assert chain.label_asym_id == "AA"
    assert chain.entity_id == "1"
    assert chain.polymer_type == "RNA"
    assert chain.description == "16S ribosomal RNA"
    assert chain.rfam_accessions == ["RF00177"]
    assert chain.tax_id == 562
    assert chain.scientific_name == "Escherichia coli"
    assert chain.superkingdom == "Bacteria"
    # 16S rRNA description contains "ribosomal" but not "ribosomal protein".
    assert chain.is_ribosomal_protein is False
    assert chain.ife == "5J7L|1|AA"


def test_parse_spec_5_1_2_response_ligand_fields() -> None:
    lig = rcsb_client.parse_assemblies(SPEC_5_1_2_ENTRY)[0].ligands[0]
    assert lig.comp_id == "MG"
    assert lig.name == "MAGNESIUM ION"
    assert lig.auth_asym_id == "AA"
    assert lig.drugbank_id is None
    assert lig.drugbank_description is None


# ---------------------------------------------------------------------------
# parse_assemblies — synthetic edge cases
# ---------------------------------------------------------------------------


def test_polymer_type_routing_splits_rna_and_protein() -> None:
    entry = _entry(
        assemblies=(
            _assembly(
                1,
                polymer_instances=(
                    _polymer_instance("AA", polymer_type="RNA", description="16S rRNA"),
                    _polymer_instance(
                        "L1",
                        entity_id="2",
                        polymer_type="Protein",
                        description="50S ribosomal protein L1",
                    ),
                    _polymer_instance(
                        "Z1",
                        entity_id="3",
                        polymer_type="DNA",
                        description="DNA probe",
                    ),
                ),
            ),
        ),
    )
    a = rcsb_client.parse_assemblies(entry)[0]
    assert {c.auth_asym_id for c in a.rna_chains} == {"AA"}
    assert {c.auth_asym_id for c in a.protein_chains} == {"L1"}


def test_ribosomal_protein_flag_from_description() -> None:
    entry = _entry(
        assemblies=(
            _assembly(
                1,
                polymer_instances=(
                    _polymer_instance(
                        "L1",
                        polymer_type="Protein",
                        description="50S ribosomal protein L1",
                    ),
                ),
            ),
        ),
    )
    chain = rcsb_client.parse_assemblies(entry)[0].protein_chains[0]
    assert chain.is_ribosomal_protein is True


def test_ribosomal_protein_flag_from_uniprot_when_description_empty() -> None:
    entry = _entry(
        assemblies=(
            _assembly(
                1,
                polymer_instances=(
                    _polymer_instance(
                        "P1",
                        polymer_type="Protein",
                        description=None,
                        uniprot_name="30S ribosomal protein S1",
                    ),
                ),
            ),
        ),
    )
    chain = rcsb_client.parse_assemblies(entry)[0].protein_chains[0]
    assert chain.is_ribosomal_protein is True


def test_ribosomal_protein_flag_false_for_factor_chains() -> None:
    """A translation factor like EF-Tu must NOT be flagged as ribosomal —
    the §13.1 substring rule needs the exact phrase "ribosomal protein"."""
    entry = _entry(
        assemblies=(
            _assembly(
                1,
                polymer_instances=(
                    _polymer_instance(
                        "T",
                        polymer_type="Protein",
                        description="Elongation factor Tu",
                    ),
                ),
            ),
        ),
    )
    chain = rcsb_client.parse_assemblies(entry)[0].protein_chains[0]
    assert chain.is_ribosomal_protein is False


def test_ribosomal_protein_flag_case_insensitive() -> None:
    entry = _entry(
        assemblies=(
            _assembly(
                1,
                polymer_instances=(
                    _polymer_instance(
                        "L1",
                        polymer_type="Protein",
                        description="50S RIBOSOMAL PROTEIN L1",
                    ),
                ),
            ),
        ),
    )
    chain = rcsb_client.parse_assemblies(entry)[0].protein_chains[0]
    assert chain.is_ribosomal_protein is True


def test_chimeric_source_takes_first_non_null_for_each_field() -> None:
    """When the source list has nulls in different fields, the parser takes
    each field's first non-null value independently."""
    entry = _entry(
        assemblies=(
            _assembly(
                1,
                polymer_instances=(
                    _polymer_instance(
                        "AA",
                        polymer_type="RNA",
                        description="16S rRNA",
                        sources=(
                            {
                                "ncbi_taxonomy_id": None,
                                "ncbi_scientific_name": None,
                                "ncbi_parent_scientific_name": "Eukaryota",
                            },
                            {
                                "ncbi_taxonomy_id": 9606,
                                "ncbi_scientific_name": "Homo sapiens",
                                "ncbi_parent_scientific_name": "Eukaryota",
                            },
                        ),
                    ),
                ),
            ),
        ),
    )
    chain = rcsb_client.parse_assemblies(entry)[0].rna_chains[0]
    assert chain.superkingdom == "Eukaryota"
    assert chain.tax_id == 9606
    assert chain.scientific_name == "Homo sapiens"


def test_uniprots_null_yields_none_uniprot_name() -> None:
    entry = _entry(
        assemblies=(
            _assembly(
                1,
                polymer_instances=(
                    _polymer_instance("X", polymer_type="Protein", description="something"),
                ),
            ),
        ),
    )
    chain = rcsb_client.parse_assemblies(entry)[0].protein_chains[0]
    assert chain.is_ribosomal_protein is False


def test_multiple_rfam_annotations_all_captured() -> None:
    entry = _entry(
        assemblies=(
            _assembly(
                1,
                polymer_instances=(_polymer_instance("AA", rfam=("RF00177", "RF00001")),),
            ),
        ),
    )
    chain = rcsb_client.parse_assemblies(entry)[0].rna_chains[0]
    assert chain.rfam_accessions == ["RF00177", "RF00001"]


def test_non_rfam_annotations_are_filtered_out() -> None:
    """Annotations with type != "Rfam" must not leak into rfam_accessions."""
    instance = _polymer_instance("AA")
    # Inject a non-Rfam annotation alongside one Rfam annotation.
    instance["polymer_entity"]["rcsb_polymer_entity_annotations"] = [
        {"annotation_id": "GO:0006412", "type": "GO", "name": "translation"},
        {"annotation_id": "RF00177", "type": "Rfam", "name": "16S"},
    ]
    entry = _entry(assemblies=(_assembly(1, polymer_instances=(instance,)),))
    chain = rcsb_client.parse_assemblies(entry)[0].rna_chains[0]
    assert chain.rfam_accessions == ["RF00177"]


def test_ligands_deduped_by_comp_id() -> None:
    """Multiple MG occurrences across different auth_asym_ids collapse to one entry."""
    entry = _entry(
        assemblies=(
            _assembly(
                1,
                nonpolymer_instances=(
                    _nonpolymer_instance("Z", "MG", name="MAGNESIUM ION"),
                    _nonpolymer_instance("Y", "MG", name="MAGNESIUM ION"),
                    _nonpolymer_instance("X", "ATP", name="ADENOSINE TRIPHOSPHATE"),
                ),
            ),
        ),
    )
    ligs = rcsb_client.parse_assemblies(entry)[0].ligands
    assert {lig.comp_id for lig in ligs} == {"MG", "ATP"}


def test_drugbank_annotation_extracted_for_ligand() -> None:
    entry = _entry(
        assemblies=(
            _assembly(
                1,
                nonpolymer_instances=(
                    _nonpolymer_instance(
                        "Z",
                        "STR",
                        name="STREPTOMYCIN",
                        drugbank=("DB01082", "Aminoglycoside antibiotic"),
                    ),
                ),
            ),
        ),
    )
    lig = rcsb_client.parse_assemblies(entry)[0].ligands[0]
    assert lig.drugbank_id == "DB01082"
    assert lig.drugbank_description == "Aminoglycoside antibiotic"


def test_entry_with_multiple_assemblies_returns_multiple_contexts() -> None:
    entry = _entry(
        assemblies=(
            _assembly(1, polymer_instances=(_polymer_instance("A1"),)),
            _assembly(2, polymer_instances=(_polymer_instance("A2"),)),
        ),
    )
    contexts = rcsb_client.parse_assemblies(entry)
    assert [a.assembly_id for a in contexts] == ["1", "2"]
    assert contexts[0].rna_chains[0].auth_asym_id == "A1"
    assert contexts[1].rna_chains[0].auth_asym_id == "A2"


def test_polymer_instance_missing_auth_asym_id_is_skipped() -> None:
    instance = _polymer_instance("A")
    instance["rcsb_polymer_entity_instance_container_identifiers"]["auth_asym_id"] = None
    entry = _entry(assemblies=(_assembly(1, polymer_instances=(instance,)),))
    a = rcsb_client.parse_assemblies(entry)[0]
    assert a.rna_chains == []
    assert a.protein_chains == []


def test_methods_propagated_to_every_assembly() -> None:
    entry = _entry(
        methods=("ELECTRON MICROSCOPY", "X-RAY DIFFRACTION"),
        assemblies=(
            _assembly(1, polymer_instances=(_polymer_instance("A"),)),
            _assembly(2, polymer_instances=(_polymer_instance("B"),)),
        ),
    )
    for a in rcsb_client.parse_assemblies(entry):
        assert a.experimental_methods == ["ELECTRON MICROSCOPY", "X-RAY DIFFRACTION"]


def test_entry_payload_missing_rcsb_id_raises() -> None:
    with pytest.raises(ApiRequestError, match="missing rcsb_id"):
        rcsb_client.parse_assemblies({})


def test_pdb_id_is_uppercased_from_entry_payload() -> None:
    """RCSB returns rcsb_id uppercased, but if a payload arrives lower-cased
    the parser must normalise so downstream IFE strings are stable."""
    entry = _entry(
        pdb_id="5j7l",
        assemblies=(_assembly(1, polymer_instances=(_polymer_instance("AA"),)),),
    )
    a = rcsb_client.parse_assemblies(entry)[0]
    assert a.pdb_id == "5J7L"
    assert a.rna_chains[0].ife == "5J7L|1|AA"


def test_int_field_rejects_json_boolean() -> None:
    """JSON true/false must not slip into tax_id (booleans are ints in Python)."""
    instance = _polymer_instance(
        "AA",
        sources=({"ncbi_taxonomy_id": True, "ncbi_scientific_name": "X"},),
    )
    entry = _entry(assemblies=(_assembly(1, polymer_instances=(instance,)),))
    chain = rcsb_client.parse_assemblies(entry)[0].rna_chains[0]
    assert chain.tax_id is None
    assert chain.scientific_name == "X"


# ---------------------------------------------------------------------------
# _parse_assembly_id
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("rcsb_id", "expected"),
    [
        ("5J7L-1", "1"),
        ("5J7L-2", "2"),
        ("5FDV-12", "12"),
        ("6ZMI-A", "A"),  # alphabetic assembly suffixes are accepted
    ],
)
def test_parse_assembly_id_standard_cases(rcsb_id: str, expected: str) -> None:
    assert rcsb_client._parse_assembly_id(rcsb_id, fallback="X") == expected


def test_parse_assembly_id_falls_back_when_malformed() -> None:
    assert rcsb_client._parse_assembly_id("not-a-pdb-id", fallback="1") == "1"
    assert rcsb_client._parse_assembly_id("", fallback="1") == "1"


# ---------------------------------------------------------------------------
# fetch_entry_payload (HTTP, mocked)
# ---------------------------------------------------------------------------


@respx.mock
def test_fetch_entry_payload_posts_to_graphql_with_uppercased_id() -> None:
    route = respx.post(rcsb_client.RCSB_GRAPHQL_URL).mock(
        return_value=httpx.Response(200, json={"data": {"entry": SPEC_5_1_2_ENTRY}})
    )
    result = rcsb_client.fetch_entry_payload("5j7l")
    assert route.called
    request = route.calls.last.request
    body = json.loads(request.content)
    assert "RibosomeEntry" in body["query"]
    assert body["variables"] == {"entry_id": "5J7L"}
    assert result["rcsb_id"] == "5J7L"


@respx.mock
def test_fetch_entry_payload_raises_on_http_500() -> None:
    respx.post(rcsb_client.RCSB_GRAPHQL_URL).mock(return_value=httpx.Response(500))
    with pytest.raises(ApiRequestError, match="HTTP 500"):
        rcsb_client.fetch_entry_payload("5J7L")


@respx.mock
def test_fetch_entry_payload_raises_on_graphql_errors() -> None:
    respx.post(rcsb_client.RCSB_GRAPHQL_URL).mock(
        return_value=httpx.Response(
            200,
            json={"errors": [{"message": "Field 'foo' is missing"}], "data": None},
        )
    )
    with pytest.raises(ApiRequestError, match="returned errors"):
        rcsb_client.fetch_entry_payload("5J7L")


@respx.mock
def test_fetch_entry_payload_raises_on_null_entry() -> None:
    respx.post(rcsb_client.RCSB_GRAPHQL_URL).mock(
        return_value=httpx.Response(200, json={"data": {"entry": None}})
    )
    with pytest.raises(ApiRequestError, match="no entry"):
        rcsb_client.fetch_entry_payload("XXXX")


@respx.mock
def test_fetch_entry_payload_raises_on_non_json_body() -> None:
    respx.post(rcsb_client.RCSB_GRAPHQL_URL).mock(
        return_value=httpx.Response(200, text="<html>not json</html>")
    )
    with pytest.raises(ApiRequestError, match="non-JSON"):
        rcsb_client.fetch_entry_payload("5J7L")


@respx.mock
def test_fetch_entry_payload_raises_on_network_error() -> None:
    respx.post(rcsb_client.RCSB_GRAPHQL_URL).mock(
        side_effect=httpx.ConnectError("connection refused")
    )
    with pytest.raises(ApiRequestError, match="request failed"):
        rcsb_client.fetch_entry_payload("5J7L")


@respx.mock
def test_fetch_entry_payload_accepts_caller_supplied_client() -> None:
    """Tests / future cache layer can pass their own client and reuse the connection pool."""
    respx.post(rcsb_client.RCSB_GRAPHQL_URL).mock(
        return_value=httpx.Response(200, json={"data": {"entry": SPEC_5_1_2_ENTRY}})
    )
    with httpx.Client(timeout=5.0) as client:
        entry = rcsb_client.fetch_entry_payload("5J7L", client=client)
    assert entry["rcsb_id"] == "5J7L"


@respx.mock
def test_fetch_then_parse_round_trip() -> None:
    """End-to-end mocked flow: fetch_entry_payload → parse_assemblies."""
    respx.post(rcsb_client.RCSB_GRAPHQL_URL).mock(
        return_value=httpx.Response(200, json={"data": {"entry": SPEC_5_1_2_ENTRY}})
    )
    entry = rcsb_client.fetch_entry_payload("5J7L")
    assemblies = rcsb_client.parse_assemblies(entry)
    assert len(assemblies) == 1
    assert assemblies[0].rna_chains[0].ife == "5J7L|1|AA"
