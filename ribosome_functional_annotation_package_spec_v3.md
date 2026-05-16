# Specification: Production Python Package for Ribosome Functional State Annotation

## 1. Goal

Build a production-quality Python package that annotates functional states of complete bacterial and eukaryotic ribosome assemblies. The package should be usable both as:

1. A command-line tool.
2. A Python library callable from another program, pipeline, notebook, or workflow manager.

The implementation should modernise the existing prototype logic while removing the dependency on the local database. External data should be retrieved via public APIs and structural contact detection should be performed locally using Gemmi.

## 2. Package name

Suggested package name:

```text
ribosome-state-annotator
```

Suggested import name:

```python
import ribosome_state_annotator
```

Alternative shorter import name:

```python
import ribostate
```

## 3. Scope for v1

### 3.1 In scope

The package should annotate whole bacterial and eukaryotic ribosome assemblies from PDB entries.

For each valid assembly, it should determine:

- Whether the assembly is a complete ribosome containing both SSU and LSU rRNA components.
- Whether the assembly should be classified as:
  - `bacterial_ribosome`
  - `eukaryotic_ribosome`
  - `eukaryotic_organellar_ribosome`
- SSU and LSU rRNA chains.
- mRNA chain, if present.
- A-site / aminoacyl tRNA chain, if present.
- P-site / peptidyl tRNA chain, if present.
- E-site / exit tRNA chain, if present.
- tRNA positional states inferred from SSU/LSU contact patterns.
- Non-ribosomal protein factors present in the assembly.
- Bound ligands present in the assembly.
- Structured skip/failure reasons for unsupported entries or assemblies.

### 3.2 Out of scope for v1

The package should not attempt to support:

- Archaeal ribosomes.
- Ribosome fragments.
- Partial ribosome assemblies missing either SSU or LSU.
- Entries solved only by NMR.
- Multiple complete ribosomes in the same assembly.
- Complex symmetry/operator expansion beyond the biological assembly coordinates downloaded from PDB/RCSB.
- Advanced tRNA state labels beyond the existing contact-based scheme.
- Local database-backed annotation.

## 3.3 Central biological concept: contact-transfer annotation

The most important conceptual idea in this package is that annotation is performed using **conserved ribosomal functional contacts**, not chain names or purely sequence-based heuristics.

The pipeline uses curated sets of reference ribosomal nucleotides that are known to physically contact:

- mRNA
- A-site tRNA
- P-site tRNA
- E-site tRNA

in experimentally determined ribosome structures.

These reference nucleotides act as:

```text
functional-site anchors
```

The workflow is therefore:

```text
reference functional contacts
    ->
evolutionary correspondence mapping
    ->
local structural neighbourhood detection
    ->
functional chain assignment
    ->
functional-state inference
```

This is fundamentally different from:

- chain-name-based annotation;
- taxonomy-only annotation;
- simple sequence alignment;
- static residue numbering transfer.

Instead, the algorithm transfers conserved ribosomal interaction sites across homologous ribosomes using BGSU nucleotide correspondence mapping.

For example:

```python
"ssu_mrna": [
    "5J7L|1|AA|G|926",
    "5J7L|1|AA|4OC|1402",
    "5J7L|1|AA|C|1403",
]
```

does NOT simply represent arbitrary conserved residues.

These are:

```text
SSU rRNA nucleotides known to contact mRNA
```

in the bacterial reference ribosome.

The pipeline therefore performs:

```text
Find equivalent SSU nucleotides in the query ribosome
    ->
Find which RNA chain contacts those nucleotides
    ->
Infer that the chain is mRNA
```

The same concept applies to:

- `ssu_atrna`
- `ssu_ptrna`
- `ssu_etrna`
- `lsu_atrna`
- `lsu_ptrna`
- `lsu_etrna`

These define conserved ribosomal contact regions associated with tRNA positioning and movement.

This contact-transfer strategy is the core scientific abstraction of the package and should be clearly reflected in:

- documentation;
- naming conventions;
- API design;
- code structure;
- testing strategy.

In particular:

- the reference nucleotide dictionaries should be treated as curated biological reference data;
- the correspondence layer should be treated as functional-site transfer;
- the Gemmi neighbour-search layer should be treated as local interaction detection.

The implementation should therefore avoid terminology suggesting these are merely generic "reference units" or generic residue mappings.

Preferred terminology:

```text
functional-site anchors
reference contact nucleotides
conserved ribosomal interaction sites
contact-transfer annotation
```

## 4. High-level workflow

For each input PDB ID:

1. Retrieve entry-level metadata from RCSB/PDB APIs.
2. Exclude entries solved by NMR.
3. Retrieve assembly-level polymer, non-polymer, Rfam, taxonomy, and method metadata.
4. Iterate over biological assemblies in the entry.
5. For each assembly:
   1. Identify RNA and protein chains belonging to that assembly.
   2. Identify Rfam mappings for assembly RNA chains.
   3. Check whether the assembly contains both SSU and LSU rRNA components.
   4. Exclude archaeal, fragment, partial, or ambiguous assemblies.
   5. Determine ribosome classification.
   6. Select reference nucleotide set:
      - Use *E. coli* reference nucleotides for bacterial-like ribosomes.
      - Use yeast reference nucleotides for cytoplasmic eukaryotic ribosomes.
      - For eukaryotic organellar ribosomes, use *E. coli* reference nucleotides initially because the rRNA core is bacterial-like.
   7. Query BGSU correspondence API to map reference nucleotides to the query PDB entry.
   8. Filter mapped units to chains present in the current biological assembly.
   9. Download or load the biological assembly coordinate file.
   10. Use Gemmi to find neighbouring chains within 5 Å of mapped reference nucleotides.
   11. Infer mRNA, A-tRNA, P-tRNA, and E-tRNA chains.
   12. Infer tRNA states from SSU/LSU contact patterns.
   13. Emit structured annotation output.

## 5. Required external data sources

### 5.1 RCSB/PDB APIs

The package uses two RCSB endpoints:

- **GraphQL Data API:** `https://data.rcsb.org/graphql` (POST). This is the
  primary endpoint and returns everything the package needs in 1–2 queries
  per PDB entry.
- **File download:** `https://files.rcsb.org/download/{pdb_id_lower}-assembly{n}.cif.gz`
  for biological assembly mmCIF files (e.g. `5j7l-assembly1.cif.gz`).
  Fall back to `https://files.rcsb.org/download/{pdb_id_lower}.cif.gz` for
  the asymmetric unit if the assembly file 404s.

The package must retrieve, for each PDB entry:

- Entry experimental method (`exptl.method`).
- Biological assemblies (`assemblies`).
- Polymer entity instances per assembly, including author and label asym IDs.
- Polymer type per entity (RNA / protein / DNA).
- Entity description (`rcsb_polymer_entity.pdbx_description`).
- Rfam annotations per RNA chain (`rcsb_polymer_entity_annotations` where
  `type == "Rfam"`).
- UniProt protein name per protein chain, where available.
- Taxonomy per entity: `ncbi_taxonomy_id`, `ncbi_scientific_name`, and
  `ncbi_parent_scientific_name`. The latter is the **superkingdom** for
  cellular organisms — this is the key signal for classification (§8) and
  is documented in the RCSB taxonomy schema as resolving to `Bacteria`,
  `Eukaryota`, or `Archaea`.
- Non-polymer entity instances per assembly, including `comp_id`,
  ligand `name`, optional DrugBank ID and description.

#### 5.1.1 Reference GraphQL query

The following single query returns everything required for one entry. The
implementation must use this (or a close variant) rather than inventing
field paths from scratch:

```graphql
query RibosomeEntry($entry_id: String!) {
  entry(entry_id: $entry_id) {
    rcsb_id
    exptl { method }
    assemblies {
      rcsb_id
      polymer_entity_instances {
        rcsb_polymer_entity_instance_container_identifiers {
          auth_asym_id
          asym_id
          entity_id
        }
        polymer_entity {
          rcsb_id
          entity_poly { rcsb_entity_polymer_type }
          rcsb_polymer_entity { pdbx_description }
          rcsb_polymer_entity_annotation {
            annotation_id
            type
            name
            description
          }
          rcsb_entity_source_organism {
            ncbi_taxonomy_id
            ncbi_scientific_name
            ncbi_parent_scientific_name
          }
          uniprots {
            rcsb_uniprot_protein { name { value } }
          }
        }
      }
      nonpolymer_entity_instances {
        rcsb_nonpolymer_entity_instance_container_identifiers {
          auth_asym_id
          comp_id
        }
        nonpolymer_entity {
          pdbx_entity_nonpoly { name }
          rcsb_nonpolymer_entity_annotation {
            type
            annotation_id
            description
          }
        }
      }
    }
  }
}
```

Notes for the implementer:

- `rcsb_polymer_entity_annotation` (singular — RCSB renamed this field
  from the original `rcsb_polymer_entity_annotations` plural in the live
  schema; the implementation accepts both names defensively) is a list;
  filter to entries where `type == "Rfam"` and read `annotation_id` (the
  accession, e.g. `"RF00177"`).
- **v3.1 live-API note (rRNA Rfam):** the live RCSB schema currently
  returns no `Rfam` annotations for rRNA polymer entities — only for
  small ncRNAs. Because this package's rRNA-typing rules (§6 and §8)
  depend on those Rfam accessions, the implementation augments the
  per-entity annotation list with a call to PDBe's REST endpoint
  `https://www.ebi.ac.uk/pdbe/api/v2/nucleic_mappings/rfam/<pdb_id>`,
  whose response maps each nucleic-acid chain ID to a list of Rfam
  accessions. See §5.3 (PDBe Rfam fallback) and the §28 addendum.
- `rcsb_entity_source_organism` is also a list (chimeric entities can have
  multiple sources); when computing the "dominant superkingdom" vote in §8,
  treat each occurrence as one vote weighted by 1.
- `ncbi_parent_scientific_name` is documented by RCSB as the superkingdom
  for cellular organisms, taking values `"Bacteria"`, `"Eukaryota"`, or
  `"Archaea"` (and clade names for viruses, which are out of scope for v1).
- DrugBank annotations appear under `rcsb_nonpolymer_entity_annotation` with
  `type == "DrugBank"`.

#### 5.1.2 Example minimal response shape

```json
{
  "data": {
    "entry": {
      "rcsb_id": "5J7L",
      "exptl": [{"method": "X-RAY DIFFRACTION"}],
      "assemblies": [
        {
          "rcsb_id": "5J7L-1",
          "polymer_entity_instances": [
            {
              "rcsb_polymer_entity_instance_container_identifiers": {
                "auth_asym_id": "AA", "asym_id": "AA", "entity_id": "1"
              },
              "polymer_entity": {
                "entity_poly": {"rcsb_entity_polymer_type": "RNA"},
                "rcsb_polymer_entity": {"pdbx_description": "16S ribosomal RNA"},
                "rcsb_polymer_entity_annotations": [
                  {"annotation_id": "RF00177", "type": "Rfam",
                   "name": "Bacterial small subunit ribosomal RNA"}
                ],
                "rcsb_entity_source_organism": [
                  {"ncbi_taxonomy_id": 562,
                   "ncbi_scientific_name": "Escherichia coli",
                   "ncbi_parent_scientific_name": "Bacteria"}
                ],
                "uniprots": null
              }
            }
          ],
          "nonpolymer_entity_instances": [
            {
              "rcsb_nonpolymer_entity_instance_container_identifiers": {
                "auth_asym_id": "AA", "comp_id": "MG"
              },
              "nonpolymer_entity": {
                "pdbx_entity_nonpoly": {"name": "MAGNESIUM ION"},
                "rcsb_nonpolymer_entity_annotation": []
              }
            }
          ]
        }
      ]
    }
  }
}
```

#### 5.1.3 Assembly file download

```text
https://files.rcsb.org/download/{pdb_id_lower}-assembly{assembly_n}.cif.gz
```

For example, for PDB `5J7L` assembly `1`:

```text
https://files.rcsb.org/download/5j7l-assembly1.cif.gz
```

The file is gzipped mmCIF. Gemmi can read it directly via
`gemmi.read_structure(path)` (gemmi auto-decompresses `.gz`).

### 5.2 BGSU RNA 3D Hub correspondence API

The BGSU correspondence API is used to map curated reference nucleotides
(§6.3) from the bacterial or yeast reference ribosome onto homologous
positions in the query PDB. This is the "evolutionary correspondence
mapping" layer of the contact-transfer workflow (§3.3).

Endpoint:

```text
GET https://rna.bgsu.edu/correspondence/map_across_chains
    ?id=<comma-separated-unit-ids>
    &scope=Rfam
    &resolution=20.0A
    &depth=700
    &format=json
```

- `scope=Rfam` makes the alignment go across species via Infernal /
  covariance models, which is what we want for cross-species transfer.
- `resolution=20.0A` permits low-resolution cryo-EM structures (most
  ribosome structures are ≤4 Å, but the loose threshold is safer).
- **`depth=700` is mandatory in practice (v3.1).** Without an explicit
  depth the live endpoint caps the response at ~120 NR
  equivalence-class representatives and silently excludes common
  ribosome PDBs such as 5J7L and 4V5Q. With `depth=700` the response
  grows to ~1260 rows and reliably includes the major bacterial /
  eukaryotic / archaeal ribosome deposits we care about. The depth is
  exposed via the `DEFAULT_BGSU_DEPTH` constant and is overridable per
  call.
- **`format=json` is mandatory (v3.1).** The live BGSU endpoint ignores
  the `Accept: application/json` header and returns tab-delimited
  `text/plain` unless the `format=json` query parameter is set.

Example request:

```text
https://rna.bgsu.edu/correspondence/map_across_chains?id=5J7L|1|AA|G|926,5J7L|1|AA|4OC|1402,5J7L|1|AA|C|1403&scope=Rfam&resolution=20.0A&depth=700&format=json
```

#### 5.2.1 Response shape

The endpoint requires `format=json` (see §5.2). Each row of the response
corresponds to a single reference unit and lists its mapped equivalent
in every aligned chain across the Rfam family.

**Spec idealised shape** (kept here as the canonical model; the
implementation accepts it as a fallback because cached / unit-test
payloads use it):

```json
{
  "query": ["5J7L|1|AA|G|926", "5J7L|1|AA|4OC|1402", "5J7L|1|AA|C|1403"],
  "alignment": [
    {
      "reference_unit": "5J7L|1|AA|G|926",
      "mapped_units": [
        "5J7L|1|AA|G|926",
        "7K00|1|a|G|926",
        "4V51|1|AA|G|926",
        "6ZMI|1|S2|G|1207",
        "5TBW|1|A|A|1150"
      ]
    },
    {
      "reference_unit": "5J7L|1|AA|4OC|1402",
      "mapped_units": [
        "5J7L|1|AA|4OC|1402",
        "7K00|1|a|C|1402",
        "6ZMI|1|S2|C|1703"
      ]
    }
  ]
}
```

**Live BGSU shape (v3.1, May 2026)** — the actual response is keyed by
target-PDB row, not by reference unit, and the per-row `unit_id_list`
is positionally aligned with the original `query`:

```json
{
  "query": ["5J7L|1|AA|G|926", "5J7L|1|AA|4OC|1402", "5J7L|1|AA|C|1403"],
  "mappings": [
    {
      "pdb": "9PJ7",
      "unit_id_list": ["9PJ7|1|A|G|910", "9PJ7|1|A|4OC|1389", "9PJ7|1|A|C|1390"],
      "rfam_EC_chain": [["RF00177", "NR_20.0_86966.1", "9PJ7|1|A"]],
      "sequence": "..."
    }
  ]
}
```

The BGSU client converts the live `mappings` shape into the spec
`{reference_unit: [mapped_unit, ...]}` map by zipping each row's
`unit_id_list` with the top-level `query`. The reference PDB's own
units are self-injected into their own mapped list so that annotating
the reference PDB itself still works.

**Two load-bearing semantic constraints (v3.1) implementers must respect:**

- **NR-representative chains only.** Even at high `depth`, BGSU
  returns at most ONE chain per PDB — specifically the chain BGSU
  picked as the non-redundant equivalence-class representative.
  Common PDB entries with multiple identical SSU/LSU copies in the
  asymmetric unit (e.g. 5J7L's `AA` vs `BA`, 5FDV's `1a` vs `2a`,
  5TBW's `A` vs `sR`) only have ONE of their chains visible in any
  given response. This means correspondence rows referencing the
  "missing" chain are simply absent from the response, not present
  but filtered out. The implementation handles this with a
  **chain-substitution fallback** (see `correspondence.py`
  `_try_chain_substitution`) that translates the assembly's
  equivalent chain ID to the BGSU-NR-representative chain before
  applying the §5.2.2 filter.
- **Batches are intersected.** A single batch call returns rows ONLY
  for target PDBs that match at EVERY query position. If even one
  query position lies in a poorly-conserved region of the Rfam
  alignment for a particular target organism, that whole PDB is
  silently dropped from the response. The implementation cannot
  detect this from the response alone — it manifests as a missing
  target. The mitigation is to choose reference anchors that lie in
  well-conserved regions across the Rfam family (§6.3 yeast
  `lsu_etrna` was trimmed from 3 anchors to 2 for exactly this
  reason — see §6.3 and §28).

The exact field names may evolve; the package's BGSU client should
treat the response as a list of (reference_unit, mapped_units[]) pairs
and be tolerant of additional fields. It accepts:

- the `mappings` / `unit_id_list` live shape;
- the `alignment` / `mapped_units` idealised shape (with
  `aligned_units` as a synonym for `mapped_units`);
- repeated `reference_unit` rows (merged);
- additional top-level and per-row fields.

#### 5.2.2 Required filtering

After parsing, the package must apply two filters in this order:

1. **Filter by PDB ID:** Keep only `mapped_units` whose unit-ID
   begins with the query PDB ID (case-insensitive). Discard rows from
   any other PDB entry.
2. **Filter by assembly chain set:** From the chain segment of each
   surviving unit ID, keep only those chains that appear in the
   biological assembly currently being processed (per the RCSB
   `polymer_entity_instances` listing for that assembly). This is
   essential because a PDB entry can contain multiple assemblies and
   correspondence rows cross all of them.

**v3.1 chain-substitution fallback:** because BGSU returns only the
NR-representative chain per PDB (see §5.2.1), step 2 can leave the
result empty for assemblies whose chain naming differs from the NR
representative. The canonical example is multi-assembly PDBs: 4V5Q's
16S rRNA copies are deposited under author chains `AA` and `CA`
(across the two biological assemblies), but BGSU collapses them to a
single NR representative — say `AA`. When the package is annotating
the assembly that contains `CA` (and not `AA`), every BGSU row
references the wrong chain and the filter discards them all. In that
case the implementation runs a one-shot substitution: it looks up the
reference's BGSU-NR-representative chain, the assembly's equivalent
chain via `polymer_entity_instances`, and rewrites the chain segment
of every mapped unit ID (`AA` → `CA`) before reapplying step 2. The
result keeps the spec's two-filter contract while accommodating
BGSU's NR redundancy collapse. See `correspondence._try_chain_substitution`.

#### 5.2.3 Missing or partial mappings

If a given reference unit has no mapping in the query PDB, record a
warning (`f"correspondence_missing_for_{reference_key}_{unit_id}"`)
and continue. A site (e.g. `ssu_atrna`) with at least one mapped
nucleotide in the assembly is sufficient for that site's contact
search; if zero of a site's reference units map, skip just that site
(do not skip the whole assembly).

### 5.3 PDBe REST API — rRNA Rfam fallback (v3.1)

The live RCSB GraphQL schema returns no `Rfam` annotation entries for
rRNA polymer entities (it still returns Rfam for small ncRNAs such as
tRNA / SRP RNA). The package therefore augments its Rfam-by-chain map
with a call to PDBe's REST endpoint:

```text
GET https://www.ebi.ac.uk/pdbe/api/v2/nucleic_mappings/rfam/<pdb_id>
```

Response shape:

```json
{
  "5j7l": {
    "Rfam": {
      "RF00177": {
        "identifier": "Bacterial small subunit ribosomal RNA",
        "mappings": [
          {"entity_id": 1, "chain_id": "AA", "start": {"residue_number": 1}, "end": {"residue_number": 1542}, ...},
          ...
        ]
      },
      "RF02541": { ... }
    }
  }
}
```

Treatment:

- Parsed into `{auth_chain_id: [Rfam_accession, ...]}`.
- Merged into the chain-level Rfam set produced from the RCSB
  `rcsb_polymer_entity_annotation` field (de-duplicated).
- The merged map is what §6.1 (Rfam role table) and §8.2 (rRNA-core
  determination) consult.
- An HTTP 404 from the endpoint is treated as "no PDBe Rfam mappings
  for this entry" — not a fatal error — and the RCSB-only map is used
  unchanged.
- Responses are cached under the `pdbe/` namespace alongside `rcsb/`
  and `bgsu/`.

### 5.4 PDB/mmCIF coordinate files

Use biological assembly mmCIF files where possible. The coordinate source should be configurable:

- `auto`: download from RCSB/PDB if not already cached.
- `local`: read from a local file or local directory.
- `url`: read from a user-provided coordinate URL.

Gemmi should be used for parsing coordinates and performing neighbour searches.

## 6. Reference data

### 6.1 Rfam accessions

The package should define Rfam accessions centrally in a constants/config module.

The implementation must avoid assuming that every ribosome has the same number of LSU rRNA chains. Different ribosome systems can split the LSU rRNA into different numbers of chains/fragments, especially across bacterial, eukaryotic cytoplasmic, mitochondrial, and plastid ribosomes. Therefore, the classification should be **role-based and accession-based**, not count-based.

Recommended role model:

```python
SSU_MAIN_RRNA = {"RF00177", "RF01960", "RF01959"}
LSU_MAIN_RRNA = {"RF02541", "RF02543", "RF02540"}
LSU_ASSOCIATED_RRNA = {"RF00001", "RF00002"}  # 5S and 5.8S-like components
TRNA = {"RF00005"}
```

More explicit mapping:

```python
RFAM_ROLE_MAP = {
    "RF00177": "ssu_main_rrna",        # bacterial-like 16S SSU rRNA
    "RF01960": "ssu_main_rrna",        # eukaryotic-like 18S SSU rRNA
    "RF01959": "ssu_main_rrna",        # archaeal SSU rRNA
    "RF02541": "lsu_main_rrna",        # bacterial-like 23S LSU rRNA
    "RF02543": "lsu_main_rrna",        # eukaryotic-like 25S/28S LSU rRNA
    "RF02540": "lsu_main_rrna",        # archaeal LSU rRNA
    "RF00002": "lsu_associated_rrna",  # 5.8S rRNA-like component
    "RF00001": "lsu_associated_rrna",  # 5S rRNA-like component
    "RF00005": "trna",
}
```

> **v3.1 note — archaeal Rfam families added.** Spec §3.2 keeps
> archaeal ribosomes out-of-scope, but the SSU/LSU detection must
> still **recognise** them so they reach the §8.5 archaeal-protein
> short-circuit (`archaeal_ribosome_not_supported`) rather than being
> dropped at the §7.4 partial-ribosome check with
> `partial_ribosome_missing_ssu_or_lsu`. Without RF01959/RF02540 in
> the role map, an archaeal 70S like 4V6U skips with the wrong reason.

The key requirement for a complete supported ribosome in v1 should be:

```text
at least one SSU main rRNA + at least one LSU main rRNA
```

Associated LSU rRNAs such as 5S and 5.8S should be captured as optional components, not used as mandatory completeness criteria.

### 6.2 Number-agnostic rRNA component model

The output model should not hard-code only `LSU_large_chain`, `LSU_medium_chain`, and `LSU_small_chain` as if every ribosome has exactly those components.

Use this instead:

```python
ssu_main_rrna_chains: list[ChainRef]
lsu_main_rrna_chains: list[ChainRef]
lsu_associated_rrna_chains: list[ChainRef]
trna_chains: list[ChainRef]
other_rna_chains: list[ChainRef]
```

For convenience, if there is exactly one SSU main rRNA and one LSU main rRNA, expose aliases:

```python
ssu_chain: ChainRef | None
lsu_chain: ChainRef | None
```

But the canonical output should remain list-based.

This makes the package robust to:

- bacterial ribosomes with 16S + 23S + 5S;
- eukaryotic cytoplasmic ribosomes with 18S + 25S/28S + 5.8S + 5S;
- organellar ribosomes where LSU rRNA organisation and protein content can differ;
- assemblies where some associated rRNA chains are unresolved or absent;
- future extension to fragmented or split rRNA systems.

### 6.3 Reference nucleotide sets for contact mapping

The package must include the reference nucleotide sets explicitly. These are the unit IDs passed to the BGSU correspondence API. The API returns mapped units in the query PDB, which are then filtered to the biological assembly being processed.

#### Bacterial reference: *Escherichia coli* / 5J7L

```python
ECOLI_REFERENCE_UNITS = {
    "ssu_mrna": [
        "5J7L|1|AA|G|926",
        "5J7L|1|AA|4OC|1402",
        "5J7L|1|AA|C|1403",
    ],
    "ssu_ptrna": [
        "5J7L|1|AA|G|1338",
        "5J7L|1|AA|A|1339",
        "5J7L|1|AA|C|1400",
    ],
    "ssu_atrna": [
        "5J7L|1|AA|G|530",
        "5J7L|1|AA|A|1492",
        "5J7L|1|AA|A|1493",
    ],
    "ssu_etrna": [
        "5J7L|1|AA|G|693",
        "5J7L|1|AA|A|694",
    ],
    "lsu_atrna": [
        "5J7L|1|DA|G|2553",
        "5J7L|1|DA|G|2583",
        "5J7L|1|DA|U|2585",
    ],
    "lsu_ptrna": [
        "5J7L|1|DA|OMG|2251",
        "5J7L|1|DA|G|2252",
        "5J7L|1|DA|G|2253",
    ],
    "lsu_etrna": [
        "5J7L|1|DA|G|2112",
        "5J7L|1|DA|G|2421",
        "5J7L|1|DA|C|2422",
    ],
}
```

#### Eukaryotic cytoplasmic reference: *Saccharomyces cerevisiae* / 7ZW0

```python
YEAST_REFERENCE_UNITS = {
    "ssu_mrna": [
        "7ZW0|1|2|G|1150",
        "7ZW0|1|2|C|1639",
        "7ZW0|1|2|C|1640",
    ],
    "ssu_atrna": [
        "7ZW0|1|2|G|577",
        "7ZW0|1|2|A|1755",
        "7ZW0|1|2|A|1756",
    ],
    "ssu_ptrna": [
        "7ZW0|1|2|G|1575",
        "7ZW0|1|2|A|1576",
        "7ZW0|1|2|C|1637",
    ],
    "ssu_etrna": [
        "7ZW0|1|2|G|904",
        "7ZW0|1|2|A|905",
    ],
    "lsu_atrna": [
        "7ZW0|1|LA|G|2922",
        "7ZW0|1|LA|G|2952",
        "7ZW0|1|LA|U|2954",
    ],
    "lsu_ptrna": [
        "7ZW0|1|LA|G|2619",
        "7ZW0|1|LA|G|2620",
        "7ZW0|1|LA|G|2621",
    ],
    "lsu_etrna": [
        "7ZW0|1|LA|G|2793",
        "7ZW0|1|LA|G|2794",
    ],
}
```

> **v3.1 note — yeast reference changed from 5TBW to 7ZW0.** Three
> reasons documented in `memory/project_bgsu_yeast_ssu_chain.md`:
> 1. 5TBW assembly-1 chain `A` isn't BGSU's NR representative for the
>    yeast 18S (BGSU picked assembly-2 chain `sR`), so queries with
>    `5TBW|1|A|...` returned zero mappings.
> 2. 5TBW chain `1` is missing residue `2454` (one of the originally
>    listed lsu_etrna anchors), and BGSU's batch endpoint is
>    all-or-nothing — one bad anchor poisons the whole query.
> 3. 7ZW0 is an active 80S (FAP-80S, rotated state) with bound mRNA +
>    two tRNAs, which is more representative for the contact-transfer
>    workflow than 5TBW's empty 80S.
>
> Same S. cerevisiae residue numbering as 5TBW, so individual anchor
> positions match the v3 values one-to-one. `lsu_etrna` was also
> trimmed from 3 anchors to 2 — position `2454` is in a less-conserved
> region of LSU 25S rRNA that BGSU's Rfam alignment doesn't map across
> all eukaryotic representatives (specifically, human 80S targets like
> 6Y57 lack a 2454 equivalent, so including it in a batched query
> causes BGSU's intersection-semantic to silently drop the target).
> §5.2.3's "one anchor per site is sufficient" rule applies.

### 6.4 Reference selection rule

Use this rule in v1:

```text
if ribosome_classification == bacterial_ribosome:
    reference_units = ECOLI_REFERENCE_UNITS
elif ribosome_classification == eukaryotic_organellar_ribosome:
    reference_units = ECOLI_REFERENCE_UNITS
elif ribosome_classification == eukaryotic_ribosome:
    reference_units = YEAST_REFERENCE_UNITS
else:
    skip as unsupported
```

Rationale:

- bacterial ribosomes use bacterial-like rRNA site references;
- eukaryotic cytoplasmic ribosomes use yeast-like rRNA site references;
- eukaryotic organellar ribosomes are initially treated as bacterial-like for reference-site mapping because their rRNA core is generally bacterial-derived, even though their protein taxonomy is eukaryotic.

## 7. Complete ribosome validation

An assembly should be considered a whole supported ribosome only if all of the following are true:

1. The entry is not solved by NMR.
2. The assembly contains at least one recognised SSU rRNA chain.
3. The assembly contains at least one recognised LSU rRNA chain.
4. The assembly does not contain multiple recognised SSU/LSU copies that would indicate multiple ribosomes or ambiguity.
5. The assembly is not only an isolated SSU or isolated LSU.
6. The assembly contains enough ribosomal protein content to plausibly represent a whole ribosome.
7. The assembly is not archaeal according to rRNA Rfam classification and taxonomy signals.

### 7.1 SSU/LSU detection

Use assembly-level Rfam mappings.

Supported bacterial-like core:

```text
SSU: RF00177
LSU: RF02541
```

Supported eukaryotic cytoplasmic core:

```text
SSU: RF01960
LSU: RF02543
```

Eukaryotic organellar ribosomes are expected to look bacterial-like at the rRNA Rfam level but eukaryotic at the dominant organism/superkingdom level.

### 7.2 Excluding archaeal ribosomes

If the dominant rRNA/protein taxonomy is Archaea, skip the assembly with:

```json
{
  "status": "skipped",
  "skip_reason": "archaeal_ribosome_not_supported"
}
```

### 7.3 Excluding NMR structures

If the entry experimental method contains `NMR`, skip all assemblies with:

```json
{
  "status": "skipped",
  "skip_reason": "nmr_structure_not_supported"
}
```

### 7.4 Handling partial/fragments

If only SSU or only LSU is detected, skip with:

```json
{
  "status": "skipped",
  "skip_reason": "partial_ribosome_missing_ssu_or_lsu"
}
```

If both SSU and LSU are present but ribosomal protein content is too low, emit a warning or skip depending on strictness mode:

```json
{
  "status": "skipped",
  "skip_reason": "likely_partial_ribosome_low_ribosomal_protein_count"
}
```

#### Minimum ribosomal protein counts (v1 thresholds)

Counts use the §8.4 voting definition of "ribosomal protein" (the broader
one), not the §13.1 definition. Counts apply per assembly. Provisional
classification (after the §8.5 rule has run) determines which threshold
applies:

| Classification                  | Minimum ribosomal protein chains |
|--------------------------------|-----------------------------------|
| bacterial_ribosome              | 15                                |
| eukaryotic_ribosome             | 30                                |
| eukaryotic_organellar_ribosome  | 20                                |

These thresholds are deliberately conservative — real complete ribosomes
typically have 50+ protein chains. The aim is to exclude obviously
incomplete deposits, not to police biological completeness.

#### Strictness mode

Strictness is configurable per call and via the CLI `--strict` flag:

```python
strict_complete_ribosome_check: bool = False  # default: warn only
```

- `strict=False` (default): below-threshold protein count emits warning
  `low_ribosomal_protein_count` and **continues** with annotation.
- `strict=True`: below-threshold protein count skips the assembly with
  `skip_reason="likely_partial_ribosome_low_ribosomal_protein_count"`.

The threshold values themselves are exposed as a config object so users
can override them without monkey-patching:

```python
class CompletenessThresholds(BaseModel):
    bacterial: int = 15
    eukaryotic: int = 30
    eukaryotic_organellar: int = 20
```

## 8. Ribosome classification

The output field `ribosome_classification` must have one of exactly three values:

```text
bacterial_ribosome
eukaryotic_ribosome
eukaryotic_organellar_ribosome
```

If the assembly cannot be assigned to one of these confidently, the assembly should be skipped with a structured reason (`skip_reason="ambiguous_ribosome_classification"` with the evidence in `classification_evidence`).

### 8.1 Classification signals

Use the following signals, in this order of trust:

1. Rfam accessions of SSU/LSU rRNAs — determines whether the rRNA core is
   **bacterial-like** (`RF00177` SSU and/or `RF02541` LSU) or
   **eukaryotic-like** (`RF01960` SSU and/or `RF02543` LSU).
2. Dominant superkingdom among **ribosomal protein chains** (see §8.4 for
   the operational definition).
3. Protein description supporting terms (`mitochondrial`, `chloroplast`,
   `plastid`, `organellar`) as **tiebreakers only**.

The package must not classify based on free-text species names or
`rcsb_entity_source_organism.ncbi_scientific_name` alone. Use NCBI
superkingdom (`ncbi_parent_scientific_name`, per §5.1) which RCSB
documents as resolving to `"Bacteria"`, `"Eukaryota"`, or `"Archaea"` for
cellular organisms.

### 8.2 rRNA-core determination (precise rule)

Look across **all** rRNA chains in the assembly (SSU main and LSU main
roles, per §6.1) and collect their Rfam accessions:

```text
rrna_core = "bacterial_like"  if any SSU chain has RF00177 AND any LSU chain has RF02541
rrna_core = "eukaryotic_like" if any SSU chain has RF01960 AND any LSU chain has RF02543
rrna_core = "bacterial_like"  if any SSU chain has RF01959 AND any LSU chain has RF02540
                              (archaeal — bucketed with bacterial_like so the §8.5
                              dispatch can fire the dominant-superkingdom check
                              and emit `archaeal_ribosome_not_supported`)
rrna_core = "mixed"           if bacterial_like AND eukaryotic_like are both true
                              (rare; flag warning)
rrna_core = "ambiguous"       otherwise
```

`"mixed"` → treat as `"bacterial_like"` for §8.3 purposes and emit warning
`mixed_rrna_core_treated_as_bacterial_like`.

`"ambiguous"` → skip the assembly with
`skip_reason="ambiguous_rrna_core"`.

> **v3.1 note — archaeal bucketing.** The archaeal-as-`bacterial_like`
> rule above is the v3.1 addition. Without it, an archaeal 70S whose
> rRNA chains carry RF01959/RF02540 (but neither RF00177 nor
> RF01960/RF02541/RF02543) would fall into the `ambiguous` bucket and
> skip with `ambiguous_rrna_core` instead of the spec-intended
> `archaeal_ribosome_not_supported`. Biologically defensible: archaeal
> ribosomes are evolutionarily closer to bacterial than eukaryotic, so
> bucketing them with bacterial_like during rRNA-core determination is
> the right v1 choice. The §8.5 truth table then routes them to the
> archaeal skip via the dominant-protein-superkingdom check.

### 8.3 Dominant ribosomal protein superkingdom (precise rule)

Compute a vote over **ribosomal protein chains in the assembly**:

```python
def dominant_protein_superkingdom(ribosomal_protein_chains: list[ChainRef]) -> tuple[str, dict]:
    """Return (dominant_superkingdom, vote_counts).

    Each ribosomal protein chain contributes 1 vote for the superkingdom in
    its first non-null `ncbi_parent_scientific_name`. Chains with no taxonomy
    do not vote. Ties are broken in this order:
        Bacteria > Eukaryota > Archaea > unknown
    (Eukaryota loses to Bacteria on a tie because organellar ribosomes
    typically have a clean Eukaryota plurality; a tie usually means thin or
    chimeric data and bacterial-like is the safer default.)
    """
    counts = Counter()
    for chain in ribosomal_protein_chains:
        if chain.superkingdom:
            counts[chain.superkingdom] += 1

    if not counts:
        return ("unknown", {})

    max_votes = max(counts.values())
    leaders = [k for k, v in counts.items() if v == max_votes]
    tiebreak_order = ["Bacteria", "Eukaryota", "Archaea"]
    for sk in tiebreak_order:
        if sk in leaders:
            return (sk, dict(counts))
    return (leaders[0], dict(counts))
```

The vote requires **at least 3 ribosomal protein chains with non-null
taxonomy**; otherwise return `"unknown"` and proceed to §8.6 (skip with
`insufficient_taxonomy_evidence`).

### 8.4 What counts as a "ribosomal protein" for this vote

This is the operational definition used **only** for the §8.3 vote. It is
deliberately broader than the §13.1 definition to avoid the circularity of
relying on `"ribosomal protein"` substring matching, which fails for
organellar proteins named `"mitochondrial ribosomal protein L11"` (those
correctly contain `"ribosomal protein"`) but also for legitimately-named
factors like `"S1"` (no `"ribosomal"` substring).

A protein chain qualifies for the §8.3 vote if **any** of:

1. Its entity description (`pdbx_description`, case-insensitive) contains
   `"ribosomal protein"`, `"ribosomal subunit protein"`, or matches the
   regex `^(MRP|MRPS|MRPL|RPS|RPL|S\d{1,2}[A-Za-z]?|L\d{1,2}[A-Za-z]?|uS\d{1,2}|uL\d{1,2}|bS\d{1,2}|bL\d{1,2}|eS\d{1,2}|eL\d{1,2}|mS\d{1,2}|mL\d{1,2})\b`.
2. Its UniProt protein name (when available) contains
   `"ribosomal protein"`.
3. Its description contains `"30S"`, `"40S"`, `"50S"`, or `"60S"` AND
   another word from `{"subunit", "protein", "component"}`.

The §13.1 definition (used for `is_ribosomal_protein` in the output) is
the simpler substring rule and is unchanged.

### 8.5 Primary classification rule (final)

```text
rrna_core = "bacterial_like"  AND  dominant_protein_sk = "Bacteria"
    -> bacterial_ribosome

rrna_core = "eukaryotic_like" AND  dominant_protein_sk = "Eukaryota"
    -> eukaryotic_ribosome

rrna_core = "bacterial_like"  AND  dominant_protein_sk = "Eukaryota"
    -> eukaryotic_organellar_ribosome

rrna_core = "eukaryotic_like" AND  dominant_protein_sk = "Bacteria"
    -> ambiguous_ribosome_classification  (skip; warn)

dominant_protein_sk = "Archaea"
    -> archaeal_ribosome_not_supported  (skip)

dominant_protein_sk = "unknown"
    -> insufficient_taxonomy_evidence  (skip)
```

### 8.6 Organellar supporting evidence

After classification, regardless of outcome, populate
`classification_evidence.supporting_terms` by scanning ribosomal protein
descriptions for these keywords (case-insensitive, exact word match):

```python
ORGANELLAR_KEYWORDS = ["mitochondrial", "chloroplast", "plastid", "organellar"]
```

If `dominant_protein_sk == "Eukaryota"` and `rrna_core == "bacterial_like"`
and none of these keywords appears, still classify as
`eukaryotic_organellar_ribosome` but emit warning
`organellar_classification_without_keyword_support`.

### 8.7 Classification evidence output

The `classification_evidence` dict must always contain these keys when a
classification is emitted:

```json
{
  "ssu_rfam": ["RF00177"],
  "lsu_rfam": ["RF02541"],
  "rrna_core": "bacterial_like",
  "dominant_ribosomal_protein_superkingdom": "Eukaryota",
  "ribosomal_protein_superkingdom_votes": {"Eukaryota": 78, "Bacteria": 2},
  "ribosomal_protein_chains_voting": 80,
  "supporting_terms": ["mitochondrial"],
  "rule": "bacterial_like_rfam_plus_eukaryotic_proteins"
}
```

`rule` is one of:

- `bacterial_like_rfam_plus_bacterial_proteins`
- `eukaryotic_like_rfam_plus_eukaryotic_proteins`
- `bacterial_like_rfam_plus_eukaryotic_proteins`

## 9. Assembly-level data model

### 9.1 Internal dataclasses / Pydantic models

Use Pydantic v2 models or dataclasses with type hints. Pydantic is preferred for validation and JSON serialisation.

Suggested models:

```python
class ChainRef(BaseModel):
    pdb_id: str
    assembly_id: str
    auth_asym_id: str
    label_asym_id: str | None = None
    entity_id: str | None = None
    polymer_type: str | None = None
    description: str | None = None
    rfam_accessions: list[str] = []
    tax_id: int | None = None
    scientific_name: str | None = None
    superkingdom: str | None = None
    is_ribosomal_protein: bool = False

    @computed_field
    @property
    def ife(self) -> str:
        """BGSU-style IFE chain identifier: ``<pdb>|<model>|<auth_asym_id>``.

        Model is hard-coded to ``1`` because v1 only reads model 0 of the
        biological assembly file (§10.4). The IFE string is the canonical
        chain identifier in CSV outputs and matches the format the original
        prototype emitted (e.g. ``5J7L|1|AA``).
        """
        return f"{self.pdb_id}|1|{self.auth_asym_id}"

class AssemblyContext(BaseModel):
    pdb_id: str
    assembly_id: str
    experimental_methods: list[str]
    rna_chains: list[ChainRef]
    protein_chains: list[ChainRef]
    ligands: list[LigandRef]
    coordinate_path: Path | None = None

class CorrespondenceResult(BaseModel):
    reference_key: str
    reference_units: list[str]
    mapped_units: list[str]
    mapped_units_by_chain: dict[str, list[str]]
    warnings: list[str] = []

class RibosomeAnnotation(BaseModel):
    """Final annotation result for one biological assembly.

    The rRNA chain fields are canonical role-based **lists**, consistent with
    §6.2. The single-chain aliases (`ssu_chain`, `lsu_chain`) are convenience
    computed properties for the common case of exactly one SSU main rRNA and
    one LSU main rRNA. Do NOT add separate `LSU_large_chain` /
    `LSU_medium_chain` / `LSU_small_chain` fields — that pattern is explicitly
    deprecated in §6.2 because it does not generalise to organellar ribosomes,
    fragmented rRNA systems, or assemblies with unresolved associated rRNAs.
    """

    pdb_id: str
    assembly_id: str | None = None  # None for entry-level skip results (e.g. NMR)
    status: Literal["annotated", "skipped", "failed"]
    skip_reason: str | None = None
    ribosome_classification: Literal[
        "bacterial_ribosome",
        "eukaryotic_ribosome",
        "eukaryotic_organellar_ribosome",
    ] | None = None

    # Canonical role-based rRNA outputs (§6.2). Always lists; may be empty.
    ssu_main_rrna_chains: list[ChainRef] = []
    lsu_main_rrna_chains: list[ChainRef] = []
    lsu_associated_rrna_chains: list[ChainRef] = []  # e.g. 5S, 5.8S
    other_rna_chains: list[ChainRef] = []  # RNAs not assigned to any role

    # Functional chains assigned by contact-transfer
    mrna_chain: ChainRef | None = None
    aminoacyl_trna_chain: ChainRef | None = None
    peptidyl_trna_chain: ChainRef | None = None
    exit_trna_chain: ChainRef | None = None

    # tRNA functional states (see §12)
    aminoacyl_trna_state: str | None = None
    peptidyl_trna_state: str | None = None
    exit_trna_state: str | None = None

    non_ribosomal_proteins: list[ChainRef] = []
    bound_ligands: list[LigandRef] = []
    classification_evidence: dict = {}
    warnings: list[str] = []

    @computed_field
    @property
    def ssu_chain(self) -> ChainRef | None:
        """Convenience alias: the single SSU main rRNA chain, if exactly one."""
        return self.ssu_main_rrna_chains[0] if len(self.ssu_main_rrna_chains) == 1 else None

    @computed_field
    @property
    def lsu_chain(self) -> ChainRef | None:
        """Convenience alias: the single LSU main rRNA chain, if exactly one."""
        return self.lsu_main_rrna_chains[0] if len(self.lsu_main_rrna_chains) == 1 else None
```

The `LigandRef` model (referenced above) is:

```python
class LigandRef(BaseModel):
    comp_id: str                              # e.g. "MG", "ATP", "STR"
    name: str | None = None                   # human-readable name
    auth_asym_id: str | None = None
    drugbank_id: str | None = None
    drugbank_description: str | None = None
```

## 10. Contact detection using Gemmi

### 10.1 Goal

Replace database lookups of unit-pair interactions with local geometry-based neighbour detection using Gemmi.

### 10.2 Contact definition

A chain is considered neighbouring a target nucleotide if any atom in the candidate chain is within 5.0 Å of any atom in the target nucleotide residue.

Default cutoff:

```python
contact_cutoff_angstrom = 5.0
```

This should be configurable.

### 10.3 Unit ID parsing

The package must parse BGSU-style unit IDs:

```text
PDB|model|chain|residue_name|residue_number
```

Example:

```text
5J7L|1|AA|G|926
```

The parser should return:

```python
UnitId(
    pdb_id="5J7L",
    model="1",
    chain="AA",
    residue_name="G",
    residue_number="926"
)
```

The parser must handle modified nucleotides such as `4OC`.

### 10.4 Mapping unit IDs to Gemmi residues

Use author chain IDs and author residue numbers where possible because BGSU unit IDs are chain/residue identifiers in author-style notation.

The implementation should:

1. Parse the mmCIF biological assembly with Gemmi via
   `gemmi.read_structure(path)` (which auto-decompresses `.cif.gz`).
2. **Iterate over `structure[0]`** (the first / primary model). Biological
   assembly mmCIF files since the May 2022 wwPDB format change use the
   asymmetric unit's chain IDs and place operator-expanded copies as
   additional chains within the same model, not as additional models.
   Therefore: read only model 0 (index) / 1 (label), and do not iterate
   over multiple models. If the file unexpectedly contains >1 model, log a
   warning and use only model 0.
3. Locate the chain using the chain ID from the unit ID (this is an author
   asym ID).
4. Locate the residue using author sequence ID / residue number.
5. Validate residue name if possible, but do not fail solely because the residue name differs due to modification or normalisation.

If a mapped unit cannot be found in the coordinate file, record a warning and continue.

### 10.5 Restrict neighbour search to current assembly chains

When processing assembly `n`, only candidate chains that belong to assembly `n` should be considered. This is essential because a PDB entry can contain multiple biological assemblies.

### 10.6 Neighbour search output

For each reference site, return all neighbouring chains grouped by polymer type:

```json
{
  "site": "ssu_atrna",
  "target_units": ["7XXX|1|A|G|530", "7XXX|1|A|A|1492"],
  "neighbouring_rna_chains": ["mRNA", "tRNA-A"],
  "neighbouring_protein_chains": ["EF-Tu"],
  "min_distances": {
    "X": 3.2,
    "Y": 4.8
  }
}
```

## 11. Inference logic

### 11.1 Component assignment

Use Rfam mappings (per the `RFAM_ROLE_MAP` in §6.1) to populate the role-based
rRNA fields:

- `ssu_main_rrna_chains` — every chain whose Rfam role is `ssu_main_rrna`.
- `lsu_main_rrna_chains` — every chain whose Rfam role is `lsu_main_rrna`.
- `lsu_associated_rrna_chains` — every chain whose Rfam role is
  `lsu_associated_rrna` (typically 5S, 5.8S).
- `other_rna_chains` — any RNA chain in the assembly that is not assigned
  one of the roles above AND is not subsequently assigned as mRNA / A-tRNA /
  P-tRNA / E-tRNA.

Use Gemmi contact detection to assign the functional chains:

- `mrna_chain`
- `aminoacyl_trna_chain`
- `peptidyl_trna_chain`
- `exit_trna_chain`

A single chain must not appear in more than one role. As each chain is
assigned, remove it from the pool of candidate RNA chains before the next
assignment step (this preserves the behaviour of the original prototype's
`remove_ife` logic).

The legacy CSV columns `lsu_medium_chain` and `lsu_small_chain` (§15.3)
are derived from `lsu_associated_rrna_chains` by Rfam role:

- `lsu_medium_chain` ← the chain with Rfam accession `RF00002` (5.8S-like)
- `lsu_small_chain` ← the chain with Rfam accession `RF00001` (5S-like)

For bacterial ribosomes this means `lsu_medium_chain` is typically empty
and `lsu_small_chain` carries the 5S. For eukaryotic cytoplasmic
ribosomes both are usually populated (5.8S and 5S respectively).
Organellar ribosomes often have neither.

### 11.2 mRNA assignment

Find RNA chains neighbouring the mapped `ssu_mrna` reference site. Exclude SSU/LSU rRNA chains and known tRNA chains where possible. If multiple candidates remain, rank by:

1. Chain length compatible with mRNA fragment.
2. Description containing `mRNA`, `messenger`, `codon`, or similar.
3. Stronger/minimum contact distance to the mRNA reference site.

> **v3.1 note — the "exclude known tRNA chains" requirement is
> load-bearing**, not optional. The SSU mRNA-anchor nucleotides (the
> decoding-loop residues that interact with the mRNA codon) ALSO
> contact tRNA anticodons. With overlapping contact geometries, an
> incoming tRNA chain can win the closest-distance race against the
> actual mRNA chain. Concrete operational rule:
>
> **Any chain mapped to Rfam `RF00005` (tRNA) is excluded from the
> mRNA candidate pool, regardless of geometric contact.** This is a
> per-pool filter that runs BEFORE the §11.2 ranking, in addition to
> the rRNA-role exclusion. Such chains remain eligible for §11.3 / §11.4
> / §11.5 (P/A/E-tRNA) assignment.
>
> Discovered empirically on 6Y57: without this filter, the P/E hybrid
> tRNA `B4` was selected as mRNA because it contacted the SSU mRNA
> anchor before the real mRNA chain `A4`.

### 11.3 P-site tRNA assignment

Find RNA chains neighbouring the mapped `ssu_ptrna` reference site. Exclude core rRNAs and assigned mRNA.

### 11.4 A-site tRNA assignment

Find RNA chains neighbouring the mapped `ssu_atrna` reference site. Exclude core rRNAs and assigned mRNA. If both mRNA and tRNA neighbour A1492/A1493-like positions, prefer the non-mRNA chain as A-site tRNA.

### 11.5 E-site tRNA assignment

Find RNA chains neighbouring the mapped `lsu_etrna` reference site. Exclude core rRNAs and assigned P-site tRNA. If the only candidate is the P-site tRNA, do not assign an independent E-site tRNA.

## 12. tRNA state inference

Keep the existing state logic but implement it with Gemmi-derived contacts rather than database unit-pair interactions.

### 12.1 A-site / aminoacyl tRNA

For assigned A-tRNA:

- Check whether it contacts `lsu_atrna` site.
- Check whether it contacts neighbouring P-site reference regions:
  - `ssu_ptrna`
  - `lsu_ptrna`

Rules:

```text
SSU state:
- if A-tRNA contacts SSU P-site reference nts: ap
- else: A

LSU state:
- if contacts LSU A-site and LSU P-site: AP
- if contacts LSU A-site only: A
- if contacts LSU P-site only: P
- if contacts neither:
    - if chain length < 30: **
    - else compute the protein-factor LSU label (see §12.4);
      if none found, use *
```

Final state:

```text
<SSU state>/<LSU state>
```

Examples:

```text
A/A
A/P
ap/AP
A/*
A/**
```

### 12.2 P-site / peptidyl tRNA

For assigned P-tRNA:

- Check whether it contacts `lsu_ptrna` site.
- Check whether it contacts neighbouring E-site reference regions:
  - `ssu_etrna`
  - `lsu_etrna`

Rules:

```text
SSU state:
- if P-tRNA contacts SSU E-site reference nts: pe
- else: P

LSU state:
- if contacts LSU P-site and LSU E-site: PE
- if contacts LSU P-site only: P
- if contacts LSU E-site only: E
- if contacts neither:
    - if chain length < 30: **
    - else compute the protein-factor LSU label (see §12.4);
      if none found, use *
```

Final state:

```text
<SSU state>/<LSU state>
```

Examples:

```text
P/P
P/E
pe/PE
P/*
P/**
```

### 12.3 E-site / exit tRNA

For assigned E-tRNA:

```text
E/E
```

If no E-tRNA is assigned:

```python
exit_trna_state = None
```

### 12.4 Protein-factor LSU label (used by §12.1 and §12.2)

When the LSU state cannot be derived from rRNA contacts (the "contacts
neither" branch) and the tRNA is full-length (≥30 nt), the LSU state is
labelled by the nearest non-ribosomal protein factor at the tRNA's
acceptor end. The rule is:

1. Identify the **last three residues** of the tRNA chain (the
   approximate CCA end) using Gemmi's per-chain residue iterator and
   taking the three highest author-sequence-numbered nucleotide residues
   on the chain. **Filter to `gemmi.EntityType.Polymer`** before sorting
   by seqid — mmCIF chains routinely include HETATM cofactors (Mg ions,
   waters, modified ligands) under the same chain ID as the polymer,
   and those typically have the highest seqids (e.g. tRNA residues
   1..76 followed by `MG/101..` and `HOH/201..`). Without the polymer
   filter the "last three residues" selection picks waters/ions whose
   distance to any nearby factor is meaningless, and §12.4 silently
   collapses to `*`. Modified nucleotides (`4SU`, `PSU`, `MIA`, `7MG`,
   `5MU`, `H2U`) keep `EntityType.Polymer` and are correctly retained.
2. Run a Gemmi neighbour search around all atoms of those three residues
   with the configured `contact_cutoff_angstrom` (default 5.0 Å).
3. Collect protein chains in the same biological assembly with at least
   one atom inside the cutoff.
4. Exclude any chain classified as a ribosomal protein by the §13.1 rule.
5. Among the remaining non-ribosomal protein chains, **rank by minimum
   atom–atom distance to the CCA-end atoms** and pick the closest.
6. The LSU label is the chosen chain's entity description
   (`pdbx_description`) **as-is** — no canonical-token mapping, no
   truncation. This matches the original prototype's behaviour and the
   expectations of downstream consumers, who treat this field as a free
   text label.

Sample LSU labels that appear in real outputs from the original
prototype (one per line; this is descriptive, not a closed enumeration):

```text
Elongation factor Tu
Elongation factor G
Elongation factor 1-alpha 1
Elongation factor 4
Elongation factor Tu 1
Elongation factor Tu-A
Selenocysteine-specific elongation factor
Peptide chain release factor 1
Ribosome-recycling factor
Translation initiation factor IF-2
Eukaryotic translation initiation factor 5B
GTP pyrophosphokinase
RPS3 isoform 1
Macrolide efflux protein
dipeptide
ermDL
```

The "non-factor" entries at the bottom (`dipeptide`, `ermDL`) are real:
they arise when the nearest protein near the CCA end is a nascent
peptide or a peptide product, not a translation factor. The package
must not filter these out — the label is information about what is
near the acceptor end, not specifically a "factor name".

Resulting tRNA state strings then look like:

```text
A/A                                  # classic A/A
A/P                                  # A/P hybrid
ap/AP                                # chimeric
A/Elongation factor Tu               # A-tRNA with EF-Tu bound at LSU
P/Ribosome-recycling factor          # P-tRNA with RRF bound
P/Eukaryotic translation initiation factor 5B
A/dipeptide                          # nascent dipeptide near CCA end
A/*                                  # full-length, no factor found
A/**                                 # ASL fragment <30 nt
```

The implementation must record, in `classification_evidence` or a
dedicated `trna_state_evidence` block:

- the chain ID of the protein factor used;
- its raw description (the same string used as the LSU label);
- the minimum distance to the CCA-end atoms.

#### Note on the prior token-table proposal

An earlier draft of this spec proposed mapping descriptions to canonical
short tokens (`T`, `R`, `G`, etc.). That proposal is **explicitly
rejected**: real outputs use the full description, downstream consumers
expect that, and a closed token table cannot accommodate the long tail of
real labels (peptide products, organellar factors, mutant variants,
species-specific paralogues). A future v2 may add an optional canonical
token alongside the raw description, but the raw description remains the
primary output.

## 13. Non-ribosomal proteins and ligands

### 13.1 Ribosomal protein detection

This is the rule used for the `is_ribosomal_protein` flag on `ChainRef`
and for excluding ribosomal proteins from §12.4's protein-factor search.
It is **intentionally narrower** than the §8.4 voting definition.

A protein chain is considered ribosomal if **any** of:

- UniProt protein name (when available) contains the substring
  `"ribosomal protein"` (case-insensitive).
- PDB entity description (`pdbx_description`) contains the substring
  `"ribosomal protein"` (case-insensitive).

Non-ribosomal proteins are all protein chains not classified as ribosomal
by this rule.

#### Known limitations

This simple substring rule has documented false positives and false
negatives that the v1 package accepts:

- **False positive**: `"ribosomal recycling factor"` contains the
  substring and will be classified as ribosomal. RRF is technically a
  factor, not a structural ribosomal protein, but its proximity behaviour
  is similar enough that excluding it from §12.4 is acceptable in v1.
- **False negative**: chains named simply `"S1"`, `"L11"`, etc. without
  the `"ribosomal"` substring will be missed. These are rare in
  modern PDB depositions but appear in older entries. A future v2 may
  adopt the broader §8.4 definition here as well.

Implementations must not silently "fix" these by adopting different
rules; document the limitation in the README and let users override via
a configurable allow-list / deny-list of entity descriptions:

```python
extra_ribosomal_descriptions: list[str] = []   # treat as ribosomal
non_ribosomal_overrides: list[str] = []        # force non-ribosomal
```

### 13.2 Protein factors near tRNA CCA end

See §12.4 for the canonical rule. §13.2 is a cross-reference only; the
implementation must not duplicate the logic.

### 13.3 Ligands

Report unique bound ligands per assembly. Each ligand entry contains:

```python
class LigandRef(BaseModel):
    comp_id: str                       # e.g. "MG", "ATP", "STR"
    name: str | None = None            # pdbx_entity_nonpoly.name
    auth_asym_id: str | None = None    # ligand asym ID
    drugbank_id: str | None = None
    drugbank_description: str | None = None
```

Deduplicate by `comp_id`: report one `LigandRef` per unique `comp_id`
present in the assembly, choosing any of its `auth_asym_id` occurrences.
Common buffer/solvent comp_ids may be filtered via a configurable
exclusion list (default: `["HOH"]` for water; magnesium and other
biologically meaningful cofactors are kept).

## 14. API design

### 14.1 Main Python API

```python
from ribosome_state_annotator import annotate_pdb, annotate_assembly

annotations = annotate_pdb("5J7L")

annotation = annotate_assembly(
    pdb_id="5J7L",
    assembly_id="1",
    contact_cutoff_angstrom=5.0,
    cache_dir="./.cache/ribosome_state_annotator"
)
```

### 14.2 Batch API

```python
from ribosome_state_annotator import annotate_many

results = annotate_many(
    pdb_ids=["5J7L", "5TBW", "6ZMI"],
    output_format="jsonl",
    continue_on_error=True,
)
```

### 14.3 CLI

The CLI must support both single-PDB and batch modes.

#### Single PDB mode

Annotate all valid biological assemblies for one PDB entry:

```bash
ribostate annotate 5J7L --output 5j7l_annotation.json
```

Annotate a specific biological assembly:

```bash
ribostate annotate 5J7L --assembly-id 1 --output 5j7l_assembly1_annotation.json
```

The `--output` flag is required unless the user explicitly requests stdout output using `--stdout`.

#### Batch mode

Batch mode must accept a plain text file containing one PDB ID per line:

```text
5J7L
5TBW
6ZMI
```

Run batch annotation:

```bash
ribostate annotate-batch pdb_ids.txt --output ribosome_annotations.json
```

The batch output should be a single JSON file containing a list of annotation objects. Each object should correspond to one processed assembly or one skipped/failed PDB/assembly result.

Example:

```json
[
  {
    "pdb_id": "5J7L",
    "assembly_id": "1",
    "status": "annotated"
  },
  {
    "pdb_id": "XXXX",
    "assembly_id": null,
    "status": "skipped",
    "skip_reason": "nmr_structure_not_supported"
  }
]
```

The CLI should also support a continue-on-error mode for batch processing:

```bash
ribostate annotate-batch pdb_ids.txt --output ribosome_annotations.json --continue-on-error
```

#### Output flag

The CLI must provide an output flag that determines where the JSON file is written:

```bash
--output /path/to/output.json
```

Short alias:

```bash
-o /path/to/output.json
```

Examples:

```bash
ribostate annotate 5J7L -o results/5j7l.json
ribostate annotate-batch pdb_ids.txt -o results/batch_annotations.json
```

If the parent output directory does not exist, the CLI should create it automatically unless `--no-create-dirs` is passed.

#### Additional useful CLI flags

```bash
--assembly-id <id>              # optional, single assembly only
--cutoff 5.0                    # neighbour-search cutoff in Å
--cache-dir <path>              # cache directory
--no-cache                      # disable cache
--strict                        # strict completeness filtering
--continue-on-error             # batch mode only
--stdout                        # write JSON to stdout instead of file
--verbose                       # INFO logging
--debug                         # DEBUG logging
```

## 15. Output formats

### 15.1 JSON output

Each assembly should produce one JSON object.

Example:

```json
{
  "pdb_id": "5J7L",
  "assembly_id": "1",
  "status": "annotated",
  "ribosome_classification": "bacterial_ribosome",
  "ssu_main_rrna_chains": [
    {"auth_asym_id": "AA", "ife": "5J7L|1|AA",
     "rfam_accessions": ["RF00177"], "polymer_type": "RNA"}
  ],
  "lsu_main_rrna_chains": [
    {"auth_asym_id": "DA", "ife": "5J7L|1|DA",
     "rfam_accessions": ["RF02541"], "polymer_type": "RNA"}
  ],
  "lsu_associated_rrna_chains": [
    {"auth_asym_id": "BA", "ife": "5J7L|1|BA",
     "rfam_accessions": ["RF00001"], "polymer_type": "RNA"}
  ],
  "other_rna_chains": [],
  "ssu_chain": {"auth_asym_id": "AA", "ife": "5J7L|1|AA"},
  "lsu_chain": {"auth_asym_id": "DA", "ife": "5J7L|1|DA"},
  "mrna_chain": {"auth_asym_id": "X", "ife": "5J7L|1|X"},
  "aminoacyl_trna_chain": {"auth_asym_id": "V", "ife": "5J7L|1|V"},
  "peptidyl_trna_chain": {"auth_asym_id": "W", "ife": "5J7L|1|W"},
  "exit_trna_chain": {"auth_asym_id": "Y", "ife": "5J7L|1|Y"},
  "aminoacyl_trna_state": "A/A",
  "peptidyl_trna_state": "P/P",
  "exit_trna_state": "E/E",
  "non_ribosomal_proteins": [],
  "bound_ligands": [],
  "classification_evidence": {
    "ssu_rfam": ["RF00177"],
    "lsu_rfam": ["RF02541"],
    "rrna_core": "bacterial_like",
    "dominant_ribosomal_protein_superkingdom": "Bacteria",
    "ribosomal_protein_superkingdom_votes": {"Bacteria": 52},
    "ribosomal_protein_chains_voting": 52,
    "supporting_terms": [],
    "rule": "bacterial_like_rfam_plus_bacterial_proteins"
  },
  "warnings": []
}
```

Note: `ssu_chain` and `lsu_chain` are derived convenience fields. Consumers
that need to handle organellar or fragmented rRNA assemblies must read the
canonical list fields (`ssu_main_rrna_chains`, `lsu_main_rrna_chains`,
`lsu_associated_rrna_chains`) directly. The `ife` field on each `ChainRef`
is the IFE string used in the CSV outputs (§15.3).

### 15.2 JSONL output

For batch mode, emit one JSON object per line.

### 15.3 CSV output

Two CSV files are emitted, with column layouts that **match the original
prototype's output files exactly** so existing downstream consumers
continue to work. The original CSVs were the de-facto contract before
this rewrite, and reproducing them is a hard requirement.

Both files use standard CSV with `\r\n` line endings (`csv.writer`
default) and quote any value containing commas, double quotes, or
whitespace runs. Empty cells are empty strings, not `null` or `NA`.

#### Chain-level annotation (`ribosome_chain_annotation.csv`)

One row per assembly. Chain identifiers are emitted as **IFE strings**
(`<pdb_id>|1|<auth_asym_id>`) via the `ChainRef.ife` computed field
(§9.1). Empty when the chain is not assigned.

```text
pdb_id,assembly_id,ssu_chain,lsu_large_chain,lsu_medium_chain,lsu_small_chain,mrna,aminoacyl_trna,peptidyl_trna,exit_trna,aminoacyl_trna_state,peptidyl_trna_state,exit_trna_state
```

Example rows reproduced from the prototype's real output:

```text
4V48,1,4V48|1|BA,4V48|1|A0,,4V48|1|A9,,,,,,,
5HCR,1,5HCR|1|1a,5HCR|1|1A,,5HCR|1|1B,5HCR|1|1v,,5HCR|1|1x,,,P/P,
4U6F,1,4U6F|1|2,4U6F|1|1,4U6F|1|4,4U6F|1|3,,,,,,,
8QBT,1,8QBT|1|i,8QBT|1|A,,8QBT|1|B,8QBT|1|4,8QBT|1|5,8QBT|1|6,8QBT|1|7,A/A,P/P,E/E
```

The four legacy single-chain columns are populated from the canonical
role-based fields as follows:

| CSV column            | Source                                                          |
|-----------------------|-----------------------------------------------------------------|
| `ssu_chain`           | `ssu_main_rrna_chains[0].ife` if exactly one; else empty        |
| `lsu_large_chain`     | `lsu_main_rrna_chains[0].ife` if exactly one; else empty        |
| `lsu_medium_chain`    | The 5.8S-role chain (Rfam `RF00002`) in `lsu_associated_rrna_chains`, if present |
| `lsu_small_chain`     | The 5S-role chain (Rfam `RF00001`) in `lsu_associated_rrna_chains`, if present |
| `mrna`                | `mrna_chain.ife` if assigned                                    |
| `aminoacyl_trna`      | `aminoacyl_trna_chain.ife` if assigned                          |
| `peptidyl_trna`       | `peptidyl_trna_chain.ife` if assigned                           |
| `exit_trna`           | `exit_trna_chain.ife` if assigned                               |

If an assembly has multiple SSU main rRNA chains or multiple LSU main
rRNA chains, the corresponding legacy column is left empty and a warning
(`multiple_chains_for_legacy_csv_column`) is recorded — the canonical
JSON list fields still carry the full information.

The chain CSV must include rows for assemblies regardless of whether
mRNA/tRNAs were assigned. Looking at the prototype's real output of
1,104 rows, only ~34% have a full mRNA + tRNA assignment; ~27% have no
tRNAs at all (just SSU/LSU). Both are valid annotations, not skips.

A `--no-csv` flag suppresses CSV emission entirely (JSON-only mode).

#### Assembly-level annotation (`ribosome_assembly_annotation.csv`)

One row per (property, value) pair, matching the prototype layout
exactly:

```text
pdb_id,assembly_id,chain,property,value
```

`chain` is either an IFE string or empty (some properties are
assembly-level and have no chain). `value` is a free-text string.

Properties emitted, with semantics matching the prototype:

| `property`              | `chain`           | `value`                                                                |
|-------------------------|-------------------|------------------------------------------------------------------------|
| `species_name`          | empty             | Scientific name of the dominant rRNA source organism (§15.3.1 below)   |
| `non_ribosomal_proteins`| IFE of the chain  | The protein's `pdbx_description`                                       |
| `bound_ligands`         | empty             | Ligand name (`name` from `LigandRef`), one row per unique ligand       |
| `unmapped_rna_chains`   | IFE of the chain  | The RNA chain's `pdbx_description` (e.g. `tRNA-Phe`, `mRNA`)           |

Real example rows from the prototype:

```text
4V48,1,,species_name,Escherichia coli
4V49,1,4V49|1|BT,non_ribosomal_proteins,general stress protein Ctc
4V49,1,4V49|1|AV,unmapped_rna_chains,mRNA
5HCR,1,5HCR|1|1z,non_ribosomal_proteins,Oncocin 10wt
5HCR,1,,bound_ligands,Iron/sulfur cluster
```

Additional properties beyond the prototype set are permitted and should
be appended at the end of the row stream rather than interleaved, to
preserve byte-equivalent prefixes when comparing against legacy outputs:

- `ribosome_classification` (assembly-level; `chain` empty; `value` is
  one of the three classification strings)
- `dominant_superkingdom` (assembly-level; `value` is `Bacteria`,
  `Eukaryota`, or `Archaea`)
- `warning` (one row per warning string)

#### 15.3.1 `species_name` derivation

`species_name` is taken from the `rcsb_entity_source_organism.ncbi_scientific_name`
of the **first SSU main rRNA chain** in the assembly. This deliberately
re-uses the rRNA source as the assembly's species label (rather than
voting across all chains) because that is what the prototype did and
what consumers expect. If the SSU has multiple source organisms
(chimeric), use the first non-null one and emit a warning
`chimeric_ssu_taxonomy`. Unlike the prototype, v1 does **not** restrict
the value to a hard-coded list of four species — any scientific name
from RCSB is emitted as-is.

## 16. Error handling

The package must not use `sys.exit()` except in the CLI entrypoint.

Library functions should return structured results or raise typed exceptions only for unexpected failures.

Suggested exception hierarchy:

```python
class RibosomeAnnotatorError(Exception): ...
class ApiRequestError(RibosomeAnnotatorError): ...
class CoordinateDownloadError(RibosomeAnnotatorError): ...
class CoordinateParsingError(RibosomeAnnotatorError): ...
class CorrespondenceMappingError(RibosomeAnnotatorError): ...
class UnsupportedRibosomeError(RibosomeAnnotatorError): ...
```

Expected unsupported cases should become `status="skipped"`, not hard failures.

## 17. Caching

Implement a cache layer for:

- RCSB GraphQL API responses.
- BGSU correspondence API responses.
- Downloaded coordinate (assembly mmCIF) files.

Default cache directory:

```text
~/.cache/ribosome-state-annotator/
```

The cache must be optional and configurable:

```bash
ribostate annotate 5J7L --cache-dir /path/to/cache
ribostate annotate 5J7L --no-cache
```

### 17.1 Cache key format

Use deterministic file-system keys:

```text
rcsb/<pdb_id_lower>.json
bgsu/<sha256(query_url)>.json
coords/<pdb_id_lower>-assembly<n>.cif.gz
```

### 17.2 Invalidation policy (v1)

The v1 cache is **content-addressed and never expires automatically**.
PDB entries change rarely after deposition; BGSU correspondence is
recomputed weekly but the changes are small. Users wanting fresh data
must either delete the cache directory or pass `--no-cache`.

The package must not honour HTTP `ETag` / `Last-Modified` in v1 — these
add complexity without meaningful benefit at this scale. A future v2
may add an explicit `--refresh` flag.

### 17.3 Cache size

No automatic eviction in v1. The CLI must support `ribostate cache info`
(report total bytes and entry count) and `ribostate cache clear` (delete
the cache directory) as housekeeping commands.

## 18. Logging

Use the standard Python `logging` module.

The package must include `logger.info()` statements at major workflow steps so users can follow progress during single-entry and batch runs.

Do not use `print()` except in the CLI rendering layer.

### 18.1 Logging requirements

Each major stage should emit an informative `logger.info()` message.

Recommended workflow logging points:

```python
logger.info("Starting ribosome annotation for PDB entry %s", pdb_id)
logger.info("Fetching entry metadata for %s", pdb_id)
logger.info("Checking experimental method for %s", pdb_id)
logger.info("Found %d biological assemblies for %s", len(assemblies), pdb_id)
logger.info("Processing assembly %s for PDB entry %s", assembly_id, pdb_id)
logger.info("Identifying RNA and protein chains for %s assembly %s", pdb_id, assembly_id)
logger.info("Checking whether assembly %s is a complete supported ribosome", assembly_id)
logger.info("Classified %s assembly %s as %s", pdb_id, assembly_id, ribosome_classification)
logger.info("Selected %s functional-site anchor set for %s assembly %s", reference_name, pdb_id, assembly_id)
logger.info("Fetching BGSU correspondence mappings for %s assembly %s", pdb_id, assembly_id)
logger.info("Downloading/loading biological assembly coordinates for %s assembly %s", pdb_id, assembly_id)
logger.info("Running Gemmi neighbour search with %.1f Å cutoff", contact_cutoff_angstrom)
logger.info("Assigning mRNA and tRNA chains for %s assembly %s", pdb_id, assembly_id)
logger.info("Inferring tRNA functional states for %s assembly %s", pdb_id, assembly_id)
logger.info("Finished annotation for %s assembly %s with status %s", pdb_id, assembly_id, status)
```

Batch mode should also log:

```python
logger.info("Starting batch ribosome annotation for %d PDB entries", len(pdb_ids))
logger.info("Processing batch item %d/%d: %s", index, total, pdb_id)
logger.info("Writing JSON output to %s", output_path)
logger.info("Finished batch annotation: %d annotated, %d skipped, %d failed", n_annotated, n_skipped, n_failed)
```

### 18.2 Log levels

Use the following conventions:

- `INFO`: high-level workflow progress.
- `WARNING`: missing correspondence, missing coordinate residue, ambiguous chain assignment, weak classification evidence.
- `ERROR`: failed API request, failed coordinate parsing, unexpected annotation failure.
- `DEBUG`: detailed neighbour search, chain scoring, raw correspondence parsing, cache hits/misses.

### 18.3 Logger setup

Each module should define a module-level logger:

```python
import logging

logger = logging.getLogger(__name__)
```

The library must not configure global logging by default. Logging configuration should be handled by the CLI or calling application.

The CLI should configure logging based on flags:

```bash
--verbose   # INFO
--debug     # DEBUG
```

## 19. Docstrings and code documentation

All public functions, classes, and methods must have clear docstrings.

Docstrings should explain:

- what the function does;
- important biological assumptions;
- input parameters;
- return values;
- possible exceptions or skip behaviour;
- whether the function performs network or file-system operations.

Use Google-style or NumPy-style docstrings consistently. Google-style is recommended.

Example:

```python
def annotate_pdb(
    pdb_id: str,
    *,
    assembly_id: str | None = None,
    contact_cutoff_angstrom: float = 5.0,
    cache_dir: Path | None = None,
    strict: bool = False,
) -> list[RibosomeAnnotation]:
    """Annotate ribosome functional states for one PDB entry.

    This is the main library entry point. It retrieves entry and assembly
    metadata, filters unsupported assemblies, transfers conserved ribosomal
    functional-site anchors using BGSU correspondence, detects local contacts
    with Gemmi, and infers mRNA/tRNA assignments and tRNA functional states.

    Args:
        pdb_id: Four-character PDB accession.
        assembly_id: Optional biological assembly identifier. If omitted, all
            biological assemblies in the entry are processed.
        contact_cutoff_angstrom: Distance cutoff used for Gemmi neighbour
            detection around mapped functional-site anchor nucleotides.
        cache_dir: Optional cache directory for API responses and coordinate
            files.
        strict: If True, apply stricter completeness filtering to likely
            partial ribosome assemblies.

    Returns:
        A list of ribosome annotation results, one per processed assembly or
        structured skipped/failed result.

    Raises:
        ApiRequestError: If required external API metadata cannot be retrieved
            and the error is not recoverable.
        CoordinateParsingError: If the coordinate file cannot be parsed.
    """
```

Private helper functions should also have short docstrings when their behaviour is non-trivial, especially for biological inference logic.

## 20. Suggested package structure

```text
ribosome-state-annotator/
├── pyproject.toml
├── README.md
├── src/
│   └── ribosome_state_annotator/
│       ├── __init__.py
│       ├── api.py
│       ├── cli.py
│       ├── config.py
│       ├── constants.py
│       ├── models.py
│       ├── exceptions.py
│       ├── rcsb_client.py
│       ├── bgsu_client.py
│       ├── coordinates.py
│       ├── gemmi_contacts.py
│       ├── classify.py
│       ├── correspondence.py
│       ├── infer.py
│       ├── output.py
│       └── cache.py
├── tests/
│   ├── test_unit_id_parser.py
│   ├── test_classification.py
│   ├── test_complete_ribosome_filter.py
│   ├── test_correspondence_parser.py
│   ├── test_gemmi_contacts.py
│   ├── test_inference.py
│   ├── test_legacy_csv_regression.py
│   └── fixtures/
│       ├── golden/                              # expected JSON outputs
│       ├── legacy_csv/
│       │   ├── ribosome_chain_annotation.csv    # prototype output, 1104 rows
│       │   └── ribosome_assembly_annotation.csv # prototype output, 6043 rows
│       └── acceptance_pdbs.txt
└── docs/
    ├── workflow.md
    ├── output_schema.md
    └── examples.md
```

## 21. Dependencies

Required:

```text
gemmi
requests or httpx
pydantic>=2
rich
typer
pandas optional for CSV/dataframe convenience
```

Recommended:

```text
pytest
pytest-cov
respx or requests-mock
ruff
mypy
pre-commit
```

## 22. Testing strategy

### 22.1 Unit tests

Test:

- Unit ID parser (including modified nucleotides like `4OC`, `OMG`).
- Rfam-based SSU/LSU detection.
- Ribosome classification rules (each of the §8.5 branches).
- §8.3 superkingdom voting (including ties).
- NMR exclusion.
- Partial ribosome exclusion (both `missing_ssu_or_lsu` and low protein count).
- BGSU correspondence response parsing.
- Assembly-chain filtering.
- Gemmi residue lookup.
- Neighbouring-chain detection.
- tRNA state inference rules (each branch of §12.1, §12.2).
- §12.4 protein-factor labelling (CCA-end residue selection, ribosomal-protein exclusion, minimum-distance ranking, raw description as label).

### 22.2 Integration tests

Use the acceptance-test PDB set from §25 below. All eight cases must pass
before v1 is considered ready.

### 22.3 Golden outputs

Two kinds of golden corpus:

**JSON golden outputs.** For selected acceptance-test entries (§25),
store expected JSON outputs in `tests/fixtures/golden/`. Tests should
compare stable fields and allow dynamic fields such as API timestamps to
differ.

**CSV regression corpus.** The original prototype's
`ribosome_chain_annotation.csv` (1,104 rows) and
`ribosome_assembly_annotation.csv` (6,043 rows) are checked in to
`tests/fixtures/legacy_csv/` as the de-facto contract for the
prototype's behaviour on a broad PDB sample. A regression test must:

1. Pick a sample of ~20 entries from the legacy chain CSV covering each
   of these row shapes:
   - bacterial complete (e.g. `4V48`, `4V49`) — SSU+LSU+5S, no mRNA/tRNAs
   - bacterial with full mRNA + A/P/E (e.g. `8QBT`, `7K50`, `4V4Y`)
   - bacterial with classic states `A/A`, `P/P`, `E/E`
   - bacterial with hybrid states (`P/E`, `A/P`)
   - bacterial with chimeric states (`pe/E`, `ap/*`)
   - bacterial with protein-factor LSU labels (`A/Elongation factor Tu`,
     `P/Ribosome-recycling factor`, `A/ermDL`)
   - eukaryotic cytoplasmic (e.g. `4U6F`, `4U4R`) — populated
     `lsu_medium_chain` (5.8S)
   - multi-assembly (e.g. `4V9P`, `4V9O`, `5E81`)
2. Re-annotate each entry through the new package.
3. Compare the new chain CSV and assembly CSV rows against the legacy
   rows.

Exact chain-ID equality is **not** required (the prototype's
chain-assignment was sensitive to ordering and could disagree with v1
on ambiguous cases). Instead the assertions are:

- `ssu_chain`, `lsu_large_chain`: same chain or empty in both.
- `mrna`, `aminoacyl_trna`, `peptidyl_trna`, `exit_trna`: when both are
  non-empty, they refer to the same chain; when one is empty and the
  other is not, log as a divergence (not a failure unless >10% of the
  sample diverges).
- tRNA states: same first character on each side of the `/` (so
  `A/Elongation factor Tu` and `A/Elongation factor Tu 2` match;
  `A/A` and `ap/A` do not — chimeric/classic must agree).
- Assembly CSV: each prototype `species_name` / `bound_ligands` /
  `non_ribosomal_proteins` row has at least one corresponding new row
  for the same `(pdb_id, assembly_id)`. Exact value-string equality
  not required (DrugBank names may differ).

This regression suite is the strongest signal that the rewrite
preserves the prototype's scientific behaviour at scale.

## 23. Implementation notes for Claude Code

### 23.1 Engineering priorities

Implement in this order:

1. Create package skeleton with `pyproject.toml`, `src/`, tests, and CLI.
2. Implement models and constants.
3. Implement RCSB client and assembly parsing.
4. Implement complete-ribosome filter and classification logic.
5. Implement BGSU correspondence client and parser.
6. Implement coordinate download/cache and Gemmi parsing.
7. Implement Unit ID parser and mapped-residue lookup.
8. Implement neighbour-chain detection.
9. Implement chain assignment and tRNA-state inference.
10. Implement JSON/JSONL/CSV output.
11. Add tests and documentation.

### 23.2 Design principles

- Keep pure logic separate from API calls.
- Keep API response parsing separate from biological inference.
- Do not hard-code network calls deep inside inference functions.
- Make all thresholds configurable.
- Return structured results for skipped assemblies.
- Avoid global mutable state.
- Avoid `sys.exit()` outside the CLI.
- Use type hints throughout.
- Use Pydantic models or dataclasses for every major data object.
- Use dependency injection where practical, especially for clients and cache directories.

## 24. Open questions / decisions confirmed for v1

The original spec listed open questions; v1 decisions are now baked in:

| Question | v1 decision |
|---|---|
| Reference set for organellar ribosomes | *E. coli* reference (§6.4). Separate mito/chloroplast references deferred to v2. |
| Minimum ribosomal protein counts | Bacterial 15, Eukaryotic 30, Organellar 20 (§7.4). |
| Low protein count: skip or warn? | Warn by default, skip with `--strict` (§7.4). |
| Local mmCIF files | First-class input via `--input-file` from v1. |
| Author vs label asym IDs | Both stored in `ChainRef`; JSON output uses `auth_asym_id` as the human-facing chain identifier and includes `label_asym_id` as a sibling field. |
| Classification ties | Tiebreak Bacteria > Eukaryota > Archaea (§8.3). |
| tRNA state factor labels | Use the full `pdbx_description` of the nearest non-ribosomal protein at the CCA end, verbatim. No canonical token table in v1 — this matches the original prototype and downstream-consumer expectations (§12.4). |

## 25. Acceptance tests

These are concrete PDB entries that the v1 implementation must handle
correctly. The expected outcomes are derived from the original prototype's
known-good outputs and from manual inspection of each entry. The
implementation passes acceptance when **all eight** produce the expected
status, classification, and (where applicable) the expected functional
chain assignments.

The acceptance tests must be run as integration tests (network-enabled)
and are separate from the cached golden-output regression tests.

### 25.1 Test matrix

| # | PDB ID | Description | Expected `status` | Expected `ribosome_classification` | Notable expectations |
|---|--------|-------------|--------------------|-------------------------------------|----------------------|
| 1 | `5J7L` | *E. coli* 70S with mRNA + A/P/E tRNAs (classic) | `annotated` | `bacterial_ribosome` | All four functional chains assigned; states `A/A`, `P/P`, `E/E` |
| 2 | `5TBW` | Yeast 80S with tRNAs | `annotated` | `eukaryotic_ribosome` | mRNA + tRNAs assigned; states from `A/A`/`P/P`/`E/E` family |
| 3 | `6ZMI` | Human 80S | `annotated` | `eukaryotic_ribosome` | Same shape as case 2 with human chain IDs (e.g. `S2`, `L5`) |
| 4 | `6ZM5` or `6ZM6` | Human mitoribosome | `annotated` | `eukaryotic_organellar_ribosome` | bacterial-like Rfam (`RF00177`/`RF02541` analogues), Eukaryota protein superkingdom |
| 5 | An isolated 30S entry (e.g. `4V42` 30S-only assembly, or any entry where assembly contains only `RF00177` and no `RF02541`) | Isolated SSU | `skipped` | n/a | `skip_reason="partial_ribosome_missing_ssu_or_lsu"` |
| 6 | An isolated 50S entry (analogous: LSU-only) | Isolated LSU | `skipped` | n/a | `skip_reason="partial_ribosome_missing_ssu_or_lsu"` |
| 7 | An archaeal ribosome (e.g. *Haloarcula* or *Pyrococcus* 70S) | Archaeal | `skipped` | n/a | `skip_reason="archaeal_ribosome_not_supported"` |
| 8 | An entry with two complete biological assemblies (e.g. `5FDV`, which deposits two ribosome copies as assembly1 and assembly2) | Multi-assembly | `annotated` × 2 | bacterial or eukaryotic per the entry | Both assemblies emitted as separate result objects |

For cases 5, 6, and 7 the implementation should not hard-code PDB IDs in
tests — pick representative entries when implementing and document the
chosen IDs in `tests/fixtures/acceptance_set.md`.

### 25.2 Acceptance smoke test

A single CLI invocation exercises the whole pipeline end-to-end:

```bash
ribostate annotate-batch tests/fixtures/acceptance_pdbs.txt \
    -o /tmp/acceptance.json \
    --continue-on-error \
    --verbose
```

The test harness loads `/tmp/acceptance.json` and verifies the expected
`status` / `ribosome_classification` / `skip_reason` for each entry. Exact
chain IDs may shift slightly between PDB revisions, so the assertions
should match on **roles and shapes**, not specific letter codes.

## 26. Glossary

Brief definitions for terms used throughout this spec. They are not
intended as definitive biological references — just as anchors so that
readers (and Claude Code) interpret the spec consistently.

- **SSU** — small ribosomal subunit. 16S rRNA in bacteria (Rfam `RF00177`);
  18S in eukaryotes (`RF01960`). Contains the decoding centre.
- **LSU** — large ribosomal subunit. 23S in bacteria (`RF02541`); 25S
  (yeast) / 28S (most other eukaryotes) (`RF02543`). Contains the
  peptidyl transferase centre.
- **A-site / P-site / E-site** — aminoacyl, peptidyl, and exit sites of
  the ribosome. tRNAs occupy these in order during the elongation cycle:
  A → P → E.
- **CCA end** — the universally conserved 3′-CCA sequence of mature tRNAs.
  The aminoacyl/peptidyl moiety is attached here. Used in §12.4 as a
  geometric anchor for finding bound protein factors (EF-Tu, release
  factors, etc.).
- **ASL** — anticodon stem-loop. A truncated tRNA construct often used
  in crystallography; <30 nt. Triggers the `**` state suffix.
- **Rfam** — a database of RNA family covariance models. Used here to
  identify rRNA roles in a structure-agnostic way.
- **BGSU unit ID** — the BGSU RNA Bioinformatics Group's identifier for
  a single nucleotide in a 3D structure. Format
  `PDB|model|chain|residue_name|residue_number`. Example
  `5J7L|1|AA|G|926`. Chain and residue identifiers are **author-style**.
- **Biological assembly** — the functional macromolecular complex as
  determined or inferred by the depositor. A single PDB entry can have
  multiple biological assemblies (e.g. `5FDV` has two ribosome copies).
- **Auth asym ID** — the chain identifier assigned by the structure's
  author (e.g. `AA`, `L1`, `S2`). This is what users see in molecular
  viewers and in BGSU unit IDs.
- **Label asym ID** — the mmCIF-internal chain identifier. Always
  present, sometimes differs from auth asym ID.
- **Superkingdom** — top-level NCBI taxonomic rank for cellular
  organisms: `Bacteria`, `Eukaryota`, or `Archaea`. Surfaced by RCSB
  as `rcsb_entity_source_organism.ncbi_parent_scientific_name`.
- **Contact-transfer annotation** — the central scientific abstraction
  of this package (see §3.3). Reference nucleotides known to contact
  mRNA / tRNAs in a reference ribosome are mapped to the query ribosome
  via BGSU correspondence, and the query ribosome's mRNA / tRNA chains
  are then identified by which chains physically contact those mapped
  nucleotides.

## 27. What to delete from the original prototype

The original prototype (`process_annotation.py` + `ribosome_annotation_new.py`)
relies on a local BGSU-style relational database (`UnitPairsInteractions2024`,
`UnitInfo`, `ChainInfo`, `ChainPropertyValue`) and on a species-name
matching layer that does not scale beyond four hard-coded organisms.
v1 is a **rewrite**, not a refactor, and must not preserve compatibility
shims for any of the following:

1. **The entire `database.py` and `models.py` SQLAlchemy layer.** No
   `db_session()`, no `UnitPairsInteractions2024`, no `UnitInfo`, no
   `ChainInfo`, no `ChainPropertyValue`. All contact information is
   derived from Gemmi at run time.
2. **`get_neighboring_chains.test_run()`** and any FR3D-derived pre-computed
   contact files. Replaced by Gemmi neighbour search.
3. **`get_common_species_name` + `get_formatted_species_name`** and the
   `CURRENT_SPECIES_NAMES` list. Replaced by NCBI superkingdom
   classification (§8.3). The reference ribosome (E. coli vs yeast) is
   chosen from the classification, not from the query's species name.
4. **The `REF_UNITS_NEW["Thermus thermophilus"]` and
   `REF_UNITS_NEW["Homo sapiens"]` reference sets.** v1 uses only the
   *E. coli* (`5J7L`) and yeast (`5TBW`) reference sets — BGSU
   correspondence covers transfer to all other organisms. Including more
   species references is a v2 concern.
5. **All `sys.exit(...)` calls outside the CLI entrypoint**, including
   the multiple-subunit and missing-rRNA exits in
   `ribosome_annotation_new.py`. Replaced by structured `status="skipped"`
   results.
6. **Python 2 syntax.** The original has `print e` (line 24 of
   `process_annotation.py`), `print "Usage: ..."` (line 461 of
   `ribosome_annotation_new.py`). v1 is Python 3.10+ only.
7. **Hard-coded `LSU_large_chain` / `LSU_medium_chain` /
   `LSU_small_chain` keys** as the **canonical output model**. The
   canonical model is role-based (§6.2, §9.1). However, the equivalent
   columns are **still emitted in the CSV** (§15.3) — they're derived
   from the role-based fields. What is deleted is using them as the
   primary in-memory shape, not the CSV column.
8. **The global `HAS_MULTIPLE_ASSEMBLIES` / `HAS_rRNA_COMPONENT` /
   `HAS_MULTIPLE_SUBUNITS` module-level flags.** No global mutable state
   (§23.2).
9. **CSV file mutation by append (`open(..., 'a')`)** in
   `generate_ribosome_chain_annotation_csv` and
   `generate_ribosome_assembly_annotation_csv`. v1 CSV output is
   write-mode and idempotent per invocation.
10. **The `_new`-suffix shadow functions** (`infer_etrna_ife_new`,
    `infer_atrna_ife_new`, `infer_interacting_chain`, etc.). Pick the
    cleaner of each pair and drop the other; do not carry both forward.

Useful logic to **preserve** (reimplemented, not copied):

- The biological inference rules in `infer_tRNA_state` (`A/A`, `A/P`,
  `ap/AP`, `**`, `*` semantics).
- The reference unit IDs for E. coli (5J7L) and yeast (5TBW) in
  `REF_UNITS_NEW`. These are the curated functional-site anchors of
  §3.3 and §6.3 and are the most valuable scientific contribution the
  prototype carries forward.
- The exclusion logic for already-assigned chains (`remove_ife` →
  pool-based assignment in §11.1).

## 28. v3.1 Live-API addendum

This section consolidates every load-bearing way the v1 implementation
diverges from the spec body above. Each item exists because the live
external API (RCSB / BGSU / PDBe) or a real-world PDB entry surfaced a
constraint that the spec didn't anticipate. The body sections above
have been amended in place to point here; this section is the
authoritative summary.

### 28.1 RCSB GraphQL: field rename + missing rRNA Rfam

- The polymer-entity annotation field is named
  `rcsb_polymer_entity_annotation` (singular) in the live schema, not
  `rcsb_polymer_entity_annotations` (plural). The implementation
  accepts both forms defensively (see §5.1.1).
- The live RCSB schema returns Rfam annotations for tRNA / SRP RNA /
  other small ncRNAs, but **not for rRNA polymer entities**. Because
  every rRNA-typing rule in §6.1 / §6.2 / §8.2 depends on those Rfam
  accessions, the implementation augments the per-chain Rfam set with
  PDBe's REST endpoint `/api/v2/nucleic_mappings/rfam/<pdb_id>` (see
  §5.3). Cache namespace: `pdbe/`.

### 28.2 BGSU correspondence: four live-API quirks

The spec §5.2 describes a clean JSON shape and a single batch query;
the live endpoint has four behaviours that the implementation must
respect.

1. **`format=json` is required** — without it the response is
   tab-delimited `text/plain` regardless of the `Accept` header.
2. **`depth=700` is required** — the default depth caps the response
   at ~120 NR equivalence-class representatives and excludes common
   ribosome PDBs (5J7L, 4V5Q). At depth 700 the response includes the
   major bacterial / eukaryotic / archaeal ribosome deposits we care
   about.
3. **NR-representative chain only** — even with high depth, BGSU
   returns one chain per PDB. Multi-chain / multi-assembly PDBs (e.g.
   4V5Q's 16S on `AA`/`CA`, 5TBW's 18S on `A`/`sR`) only have one
   chain in any response. The `correspondence._try_chain_substitution`
   fallback rewrites the chain segment of each mapped unit ID from
   the NR representative to the assembly's equivalent chain (looked up
   via `polymer_entity_instances`) before reapplying the §5.2.2
   filter, so multi-assembly entries still get complete contact
   transfer (see §5.2.2 v3.1 chain-substitution fallback paragraph).
4. **Batch semantic is all-or-nothing + intersection** — when ANY
   query unit references a residue not resolved in the deposited
   structure, the entire batch returns zero mappings. When the batch
   succeeds, the response is the **intersection** of PDBs across all
   query positions: a target organism whose Rfam alignment lacks one
   of the query positions silently disappears from the response. This
   constraint shapes the yeast reference choice (§6.3 / §28.4).

### 28.3 §12.4 polymer filter for CCA-end selector

mmCIF chains routinely include HETATM cofactors (Mg ions, waters,
modified ligands) deposited under the same chain ID as the polymer.
For a tRNA chain that nominally contains residues 1..76, the same
chain often also contains `MG/101..` and `HOH/201..`. Naïvely picking
the three highest-seqid residues gives Mg / waters whose distance to
nearby factors is meaningless, and the §12.4 protein-factor LSU label
collapses to `*`. The fix (and the requirement) is to filter to
`gemmi.EntityType.Polymer` before sorting by seqid — see §12.4 step 1.
Modified nucleotides (4SU, PSU, MIA, 7MG, 5MU, H2U) keep
`EntityType.Polymer` and are correctly retained.

### 28.4 §6.3 yeast reference: switch from 5TBW to 7ZW0

The spec originally specified *S. cerevisiae* 5TBW as the eukaryotic
cytoplasmic reference. The implementation switched to 7ZW0 because of
three independent live-API constraints stacking on top of each other:

1. **5TBW chain `A` is not a BGSU NR representative.** BGSU picked
   assembly-2 chain `sR` as the SSU representative, so the spec's
   natural assembly-1 chain `A` returns 0 mappings (per §28.2 quirk
   3). 7ZW0 chain `2` IS the BGSU representative, so queries work
   directly.
2. **5TBW is missing residue 2454 on its LSU chain** — one of the
   original `lsu_etrna` anchors. Per §28.2 quirk 4, this means the
   batch returns 0 mappings. 7ZW0 has 2454 resolved, but see #3.
3. **Even on 7ZW0, position 2454 is in a less-conserved region of
   the 25S rRNA Rfam alignment**, and including it in the
   `lsu_etrna` batch causes the response to lose cross-organism
   targets such as 6Y57 (human 80S) entirely — collapsing the P-tRNA
   state from `P/E` to `P/*`. The `lsu_etrna` anchor tuple was
   therefore trimmed from 3 to 2 well-conserved positions
   (`G|2793`, `G|2794`). §5.2.3 explicitly permits one or more
   anchors per site as sufficient.

### 28.5 §6.1 archaeal Rfam additions

The role table in §6.1 must include the archaeal SSU/LSU Rfam
accessions (RF01959 for the archaeal 16S-like SSU rRNA, RF02540 for
the archaeal 23S-like LSU rRNA) in `SSU_MAIN_RRNA` and `LSU_MAIN_RRNA`
respectively. Without these the §8.2 rRNA-core determination returns
`ambiguous_rrna_core` for any archaeal ribosome and the §8.5 dispatch
cannot find a reference. The implementation buckets the
archaeal-pair (RF01959 + RF02540) as `bacterial_like` for the purpose
of selecting the *E. coli* anchor set in §8.5 — archaeal ribosomes are
phylogenetically closer to bacteria than to eukaryotes at the
secondary-structure level the contact-transfer workflow operates on.

### 28.6 §11.2 mRNA assignment: tRNA-Rfam exclusion

The mRNA-assignment rule in §11.2 must, **before** picking the closest
SSU-contacting chain, exclude any candidate chain whose Rfam mapping
classifies it as a tRNA. In 6Y57 (human 80S) the P-site tRNA's CCA
end is close to the SSU mRNA decoding centre and was being
mis-assigned as the mRNA chain. The implementation rule is:

```python
trna_rfam_ifes = {c.ife for c in by_role.get("trna", [])}
mrna_candidates = [c for c in candidate_pool if c.ife not in trna_rfam_ifes]
```

This is load-bearing for the 6Y57 / 6ZMI human-80S acceptance cases.

### 28.7 Multi-assembly PDB entries

When a PDB contains multiple biological assemblies (e.g. 4V5Q, 5TBW),
two acceptance-test conventions must be respected:

- Use a per-assembly lookup (the test helper `_by_assembly_id`) rather
  than the `_only` accessor, because `_only` raises when there are
  multiple `AnnotatedAssembly` results in one `AnnotationResult`.
- Apply the §5.2.2 step-2 chain filter per-assembly. Without this,
  rows whose chain belongs to a *different* assembly of the same PDB
  leak through and corrupt the contact-transfer result. The
  per-assembly filter is the implementation of the spec's "biological
  assembly currently being processed" wording.

### 28.8 Cache namespaces

The v1 cache layer uses four namespaces (one more than the spec
implied):

- `rcsb/` — GraphQL responses, keyed by query hash + variables.
- `bgsu/` — Correspondence-mapping responses, keyed by sorted unit ID
  list + scope + resolution + depth + format.
- `pdbe/` — PDBe Rfam-mapping responses (per §28.1), keyed by PDB ID.
- `coords/` — mmCIF coordinate files, keyed by PDB ID + assembly ID.

The `CacheInfo` summary surface (§17) reports the per-namespace entry
count for all four.
