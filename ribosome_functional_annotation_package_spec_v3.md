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
   13. Attach RADdb large-scale movement metrics
       (`intersubunit_rotation`, `ssu_head_rotation`) for the
       assembly's `(pdb_id, lsu_chain_id, ssu_chain_id)` triple —
       see §29.
   14. Extract per-site codon/anticodon evidence
       (`trna_mrna_interactions`) from FR3D base pairs — see §30.
   15. Emit structured annotation output.

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

### 5.3 EBI Rfam ``pdb_full_region`` flat file — rRNA Rfam source (v3.2)

> **Replaces the v3.1 PDBe REST endpoint.** The v3.1 spec used PDBe's
> per-entry `/nucleic_mappings/rfam/<pdb_id>` REST endpoint as the
> Rfam-by-chain source. v3.2 replaces it with a single weekly-refreshed
> flat file from EBI (see §32) that delivers better data with one
> background download instead of an HTTP call per assembly. The
> `pdbe_client.py` module and the `pdbe/` cache namespace are
> deprecated — kept on disk only for backward compatibility with old
> caches.

The live RCSB GraphQL schema returns no `Rfam` annotation entries for
rRNA polymer entities (it still returns Rfam for small ncRNAs such as
tRNA / SRP RNA). The package therefore augments its Rfam-by-chain map
with values pulled from EBI's

```text
https://ftp.ebi.ac.uk/pub/databases/Rfam/.preview/pdb_full_region.txt.gz
```

This is a single tab-separated file (~200 KB compressed, ~28 k rows,
~4400 distinct PDBs) listing every PDB chain's Rfam hits with cmsearch
alignment bit-scores. The package caches it under `~/.cache/
ribosome-state-annotator/rfam/` and refreshes it weekly via the
upstream `Last-Modified` header.

For each `(pdb_id, auth_chain_id)` the package picks the **single
highest-bit-score** Rfam accession. This collapses the cross-family
HMM hits that PDBe's REST endpoint surfaces (e.g. a eukaryotic 18S
chain being tagged simultaneously with RF00177 bacterial 16S +
RF01959 archaeal SSU + RF01960 eukaryotic 18S) down to the single
biologically-correct family — which eliminates the MIXED rRNA-core
edge case at its source.

Treatment:

- Parsed into `{auth_chain_id: [single_best_Rfam_accession]}`.
- Replaces (not merges with) the chain-level Rfam set produced from
  the RCSB `rcsb_polymer_entity_annotation` field — the file is
  authoritative when it has a hit. Chains not in the file keep their
  RCSB-supplied Rfam tags (typically none for rRNA).
- The resulting map is what §6.1 (Rfam role table) and §8.2 (rRNA-core
  determination) consult.
- A network failure or missing local file is non-fatal — the
  RCSB-only map is used unchanged.
- See §32 for the full data-source, refresh-policy, and CLI surface
  specification.

### 5.4 PDB/mmCIF coordinate files

Use biological assembly mmCIF files where possible. The coordinate source should be configurable:

- `auto`: download from RCSB/PDB if not already cached.
- `local`: read from a local file or local directory.
- `url`: read from a user-provided coordinate URL.

Gemmi should be used for parsing coordinates and performing neighbour searches.

### 5.5 RADdb LSU↔SSU CSV (large-scale movement metrics)

A weekly snapshot of inter-subunit + SSU-head rotation /
translation parameters for every public ribosome assembly,
published at
`https://radtool.rc.northeastern.edu/public_database/RADdb.<YYYYMMDD>.LSUSSU.csv`.
Used to attach two motion metrics
(`intersubunit_rotation`, `ssu_head_rotation`) to each annotated
assembly. See §29 for the full local-cache + refresh policy, the
match-key rules, and the JSON output schema.

### 5.6 FR3D base-pair CSV (codon-anticodon evidence)

Curated per-PDB base-pair annotations from the FR3D pipeline,
published at
`http://rna.bgsu.edu/rna3dhub/pdb/<pdb_id>/interactions/fr3d/basepairs/csv`.
Three columns: `unit_id_1, interaction_type, unit_id_2`. Used to
populate per-site `trna_mrna_interactions` evidence (codon residues,
anticodon residues, FR3D base-pair labels). See §30 for the full
extraction pipeline.

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
    reference_units = ECOLI_REFERENCE_UNITS_ORGANELLAR   # see §31.4
elif ribosome_classification == eukaryotic_ribosome:
    reference_units = YEAST_REFERENCE_UNITS
else:
    skip as unsupported
```

Rationale:

- bacterial ribosomes use bacterial-like rRNA site references;
- eukaryotic cytoplasmic ribosomes use yeast-like rRNA site references;
- eukaryotic organellar ribosomes use a **filtered** variant of the
  E. coli reference set (§31.4) — four anchors are removed because the
  corresponding E. coli residues have no mt-rRNA equivalent and trigger
  BGSU's intersection-semantic to drop every mt deposit from a batched
  response. BGSU's structural alignment cross-walks the remaining 16
  anchors onto mt-12S / mt-16S residues in any deposit's native
  numbering.

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

If both SSU and LSU rRNA roles are present but with **asymmetric chain
counts** (e.g. 1 SSU chain + 3 LSU chains, as in Tetrahymena /
Chlamydomonas chloroplast / fragmented human mitoribosome deposits
where the 28S or 23S rRNA is biologically cleaved into multiple
chains), skip with:

```json
{
  "status": "skipped",
  "skip_reason": "fragmented_ribosome_not_supported"
}
```

Matched-count multi-chain assemblies (`n_ssu == n_lsu >= 2`) are
**not** fragmented — they are di-ribosomes / in-situ polysomes and
are split into per-ribosome sub-contexts (§31.1) rather than skipped.
The fragmented-skip is `(n_ssu >= 2 OR n_lsu >= 2) AND n_ssu != n_lsu`.

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

`"mixed"` → defensive demotion to `"bacterial_like"` with warning
`mixed_rrna_core_treated_as_bacterial_like`.

This branch is **unreachable from real data** in v3.2 because the
production Rfam source (§32, EBI `pdb_full_region.txt.gz`) selects
the single highest-bit-score Rfam accession per chain, so each
chain's `rfam_accessions` list contains at most one entry and the
`bacterial_like AND eukaryotic_like` predicate cannot both be true.
The demotion is preserved as a belt-and-braces fallback for
synthetic / hand-crafted `ChainRef`s in tests that supply multi-tag
annotations directly. See §31.5 for the v3.2-rc1 protein-vote
resolution that was used briefly before §32 made the workaround
unnecessary.

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
    # RADdb-derived large-scale movement metrics (§29). Always emitted
    # for annotated assemblies — null when RADdb is unavailable or the
    # (pdb, lsu, ssu) triple doesn't match a RADdb row.
    large_scale_movements: LargeScaleMovements | None = None
    # Per-site (A/P/E) codon/anticodon base-pair evidence from FR3D (§30).
    # Always emitted (possibly empty) for stable consumer schema.
    trna_mrna_interactions: list[TRNAmRNAInteraction] = []
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
>
> **v3.2 extension — description-based fallback (§31.6).** Some
> deposits omit the `RF00005` Rfam annotation on tRNA chains (notably
> mitoribosome A/A-tRNAs like 7QI5 chain Aw, and the synthetic
> CCA-end tRNA analogs in 7RQA / 8T8C). To catch these, any chain
> whose `description` contains the substring `"tRNA"`
> (case-insensitive) is *also* excluded from the mRNA candidate pool.
> The risk of a false-positive (an mRNA labelled "tRNA-binding-site
> mimic") is acceptable because mRNA assignment then just skips and
> the chain falls into `other_rna_chains`.

### 11.3 P-site tRNA assignment

Find RNA chains neighbouring the mapped `ssu_ptrna` reference site. Exclude core rRNAs and assigned mRNA.

### 11.4 A-site tRNA assignment

Find RNA chains neighbouring the mapped `ssu_atrna` reference site. Exclude core rRNAs and assigned mRNA. If both mRNA and tRNA neighbour A1492/A1493-like positions, prefer the non-mRNA chain as A-site tRNA.

### 11.5 E-site tRNA assignment

Find RNA chains neighbouring the mapped `lsu_etrna` reference site. Exclude core rRNAs and assigned P-site tRNA. If the only candidate is the P-site tRNA, do not assign an independent E-site tRNA.

### 11.6 LSU-based fallback for A and P (v3.2)

After the canonical SSU-anchor pass above (§11.2–§11.5) finishes, any
A or P slot that is still empty is filled by a closer-LSU-site
fallback using `lsu_atrna` / `lsu_ptrna` proximity. See §31.2 for the
full algorithm. The fallback recovers tRNA analogs that engage only
the PTC (CCA-end fragments in 7RQA, 8T8C), pre-accommodation full
tRNAs whose acceptor end sits at the PTC before anticodon-mRNA
pairing (3JAG, 7O7Z, 7OSM, 7UG7), and mitoribosome P-tRNAs whose SSU
contact is displaced (8OIR P-Met-tRNA).

## 12. tRNA state inference

Keep the existing state logic but implement it with Gemmi-derived contacts rather than database unit-pair interactions.

### 12.1 A-site / aminoacyl tRNA

For assigned A-tRNA:

- Check whether it contacts `lsu_atrna` site.
- Check whether it contacts neighbouring P-site reference regions:
  - `ssu_ptrna`
  - `lsu_ptrna`

Rules (v3.2):

```text
SSU state:
- if contacts SSU A-site AND contacts SSU P-site: ap
- if contacts SSU A-site only:                    A
- if contacts neither (LSU-fallback assignment, §11.6):
    - if chain length < 30: **
    - else:                  *

LSU state:
- if contacts LSU A-site AND LSU P-site:          AP
- if contacts LSU A-site only:                    A
- if contacts LSU P-site only:                    P
- if contacts neither:
    - if chain length < 30: **
    - else compute the protein-factor LSU label (see §12.4);
      if none found, use *
```

The `*` / `**` distinction is symmetric: `**` always means "chain
shorter than 30 nt, physically cannot reach this subunit", `*` always
means "full-length polymer that doesn't make a canonical contact at
this subunit" (positionally displaced, or no protein factor at CCA).
See §31.7 for the rationale and the `**/**` safeguard (§31.8) that
demotes a chain with no anchor contact on either side to
`unmapped_rna_chains`.

Final state:

```text
<SSU state>/<LSU state>
```

Examples:

```text
A/A         canonical A site
A/P         classical A on SSU, displaced toward P on LSU
ap/AP       fully chimeric (translocation intermediate)
A/*         no LSU contact, no factor
A/**        ASL fragment (no LSU contact, chain < 30 nt)
*/AP        full-length pre-accommodation A-tRNA at PTC (3JAG, 7O7Z)
**/A        CCA-end-only tRNA analog at PTC A-site (7RQA)
```

### 12.2 P-site / peptidyl tRNA

For assigned P-tRNA:

- Check whether it contacts `ssu_ptrna` site.
- Check whether it contacts neighbouring E-site reference regions:
  - `ssu_etrna`
  - `lsu_etrna`

Rules (v3.2):

```text
SSU state:
- if contacts SSU P-site AND contacts SSU E-site: pe
- if contacts SSU P-site only:                    P
- if contacts neither (LSU-fallback assignment, §11.6):
    - if chain length < 30: **
    - else:                  *

LSU state:
- if contacts LSU P-site AND LSU E-site:          PE
- if contacts LSU P-site only:                    P
- if contacts LSU E-site only:                    E
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
*/P         CCA-end tRNA analog at PTC P-site (7RQA chain 1x)
**/P        same, when chain < 30 nt
```

### 12.3 E-site / exit tRNA

For assigned E-tRNA, the LSU state is always `E` because assignment
itself was made by `lsu_etrna` contact. The SSU side uses the same
fragment-vs-displaced convention as A- and P-tRNA:

```text
SSU state:
- if contacts SSU E-site:                E
- if contacts neither SSU E nor P:
    - if chain length < 30: **
    - else:                  *

LSU state: always E (assignment is by lsu_etrna contact)
```

Final state:

```text
<SSU state>/E
```

Examples:

```text
E/E         canonical E site
*/E         full polymer at LSU E only (mitoribosome E-tRNA where the
            mt-12S exit-site anchor is filtered out, §31.4)
**/E        CCA-tripeptide fragment at lsu_etrna (rare; safeguard
            §31.8 then demotes the chain to unmapped_rna_chains)
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
- PDBe Rfam-mapping responses (§28.1).
- Downloaded coordinate (assembly mmCIF) files.
- The RADdb LSU↔SSU CSV plus its metadata sidecar (§29 — separate
  refresh policy: 7-day weekly window, not content-addressed).
- FR3D per-PDB base-pair CSVs (§30 — content-addressed).
- Per-component PDB Chemical Component Dictionary CIFs (§30 — used to
  resolve authoritative parent-base info for modified nucleotides that
  Gemmi's built-in tabulated table doesn't fully describe).

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
pdbe/<pdb_id_lower>.json
coords/<pdb_id_lower>-assembly<n>.cif.gz
fr3d/<pdb_id_lower>.csv
ccd/<comp_id_upper>.cif
raddb/RADdb.LSUSSU.csv
raddb/RADdb.LSUSSU.metadata.json
```

The `raddb/` namespace uses a fixed filename (not content-addressed)
because the package keeps a single rolling copy of the upstream
weekly release. See §29.3 for the metadata sidecar schema and §29.4
for the refresh policy.

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
│       ├── cache.py
│       ├── raddb.py
│       ├── trna_mrna.py
│       └── ccd_client.py
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
- `raddb/` — Cached RADdb LSU↔SSU CSV plus its metadata sidecar
  (per §29). Refreshed weekly on a separate policy from the other four
  namespaces.

The `CacheInfo` summary surface (§17) reports the per-namespace entry
count for all five.

## 29. RADdb integration (large-scale ribosome movement metrics)

This section is the authoritative spec for the RADdb integration that
annotates each ribosome assembly with the canonical inter-subunit and
SSU-head rotation angles. RADdb is the only data source in the
package that is **periodic** (not content-addressed): the upstream
file is regenerated weekly and the package transparently keeps a
local copy fresh.

### 29.1 Goal

For every annotated assembly, surface the two most widely-used
RADdb-derived motion metrics:

- `intersubunit_rotation` — the canonical ratchet-like inter-subunit
  rotation (RADdb column `body rot.`).
- `ssu_head_rotation` — the SSU-head swivel that gates translocation
  (RADdb column `head rot.`).

These two angles are sufficient to compare ribosome conformations
across structures and are biologically interpretable on their own.
RADdb's remaining columns (tilt, translation, directionality) are
loaded into memory but **not** exposed in v1 — they may surface in a
future version for advanced analysis / clustering / ML use cases.

### 29.2 Source

```text
https://radtool.rc.northeastern.edu/public_database/RADdb.<YYYYMMDD>.LSUSSU.csv
```

The `<YYYYMMDD>` segment is the release date and changes weekly. The
implementation must NOT hard-code any single date. See §29.5 for the
latest-URL discovery procedure.

### 29.3 Local cache layout

The RADdb cache lives under the existing user-cache root introduced
in §17:

```text
~/.cache/ribosome-state-annotator/raddb/
├── RADdb.LSUSSU.csv
└── RADdb.LSUSSU.metadata.json
```

The sidecar metadata file uses the following schema:

```json
{
  "source_url": "https://radtool.rc.northeastern.edu/public_database/RADdb.20260508.LSUSSU.csv",
  "downloaded_at": "2026-05-17T16:42:16+00:00",
  "rad_date": "20260508",
  "local_file": "RADdb.LSUSSU.csv"
}
```

`downloaded_at` is ISO-8601 UTC with second resolution; `rad_date` is
the `YYYYMMDD` substring of the source URL.

### 29.4 Refresh policy

At the start of every annotation run:

1. If the local CSV is missing → attempt to download the latest
   available release (see §29.5).
2. If the local CSV is present and its `downloaded_at` is **less than
   7 days** old → use the cached file; no network call.
3. If the local CSV is present and **≥ 7 days** old → probe the
   upstream server for a newer release. Download if newer; otherwise
   keep using the cached file and log a warning.
4. If the user passes `--refresh-raddb` (CLI) or `refresh_raddb=True`
   (library) → unconditionally probe for a newer release regardless
   of cache age.
5. If any network step fails:
   - With a cached file present → fall back to the cached file and
     log a warning; the annotation pipeline continues with the
     potentially stale metrics.
   - Without a cached file → continue with `large_scale_movements`
     emitted in the **unavailable** shape (see §29.8).

**Hard requirement:** RADdb integration MUST NEVER crash the
annotation pipeline. Every failure path degrades to null metrics.

### 29.5 Latest-URL discovery

To find the newest available RADdb release without an authoritative
index:

1. Generate candidate dates by walking backwards from today (UTC)
   for at most **60 days**.
2. Format each candidate as `YYYYMMDD` and substitute into the URL
   template above.
3. Issue an HTTP `HEAD` request for each candidate URL (cheap, no
   body transfer).
4. The first URL that responds `200 OK` is the latest release.
5. If no candidate within the window responds `200` → return `None`
   and degrade per §29.4 step 5.

`HEAD` is preferred over `GET` because RADdb files are ~500 KB and
the probe is run at every annotation startup beyond the 7-day
window.

### 29.6 CSV parsing

The CSV has a single header row with the following columns relevant
to v1:

| Column | Type | Maps to |
|--------|------|---------|
| `RCSB` | string (uppercase 4-character PDB ID) | match key |
| `LSU chain ID` | string (case-sensitive) | match key |
| `SSU chain ID` | string (case-sensitive) | match key |
| `body rot.` | float (degrees) | `intersubunit_rotation` |
| `head rot.` | float (degrees) | `ssu_head_rotation` |

Note the trailing periods on `body rot.` / `head rot.` — preserve
them verbatim when reading the column.

Required columns: all five above. If any is missing from the header,
the implementation logs a warning and treats the dataset as
unavailable (null metrics for every assembly).

Numeric coercion: each value is stripped of whitespace, then
`float()`-converted. Blank or non-numeric values yield `None` for
that metric (the lookup still returns a row, but the per-metric
field is null).

### 29.7 Lookup keying and matching

#### 29.7.1 Match key

Each RADdb row is keyed by the triple:

```text
(pdb_id.upper(), lsu_chain_id, ssu_chain_id)
```

Rules:

- **PDB ID** is uppercased before comparison (RADdb stores upper).
- **Chain IDs** preserve exact case — `A` and `a` are distinct.
- The implementation must NOT match by `assembly_id`. A PDB entry
  with multiple biological assemblies has one RADdb row per
  assembly, distinguished by the assembly's specific LSU + SSU
  chain IDs. Example: 5J7L produces two rows, `(5J7L, DA, AA)` for
  assembly 1 and `(5J7L, CA, BA)` for assembly 2, each with
  different rotation values.

#### 29.7.2 Lookup table

When loading the CSV, build a dictionary keyed by the triple. This
avoids repeated pandas-style filtering during batch runs:

```python
{
    ("5J7L", "DA", "AA"): {
        "intersubunit_rotation": 5.8,
        "ssu_head_rotation": 9.4,
        ...
    },
    ("5J7L", "CA", "BA"): {
        "intersubunit_rotation": -0.6,
        "ssu_head_rotation": 5.9,
        ...
    },
}
```

Rows whose key is missing the PDB ID, LSU chain, or SSU chain are
skipped silently.

#### 29.7.3 Duplicate keys

If two or more rows share the same `(pdb, lsu, ssu)` key:

- A warning is logged at load time (`raddb_duplicate_keys`).
- Subsequent lookups for that key return `None` (null metrics)
  rather than silently picking one row.

#### 29.7.4 Per-assembly lookup

For each annotated assembly:

1. Take the assembly's SSU main rRNA chain `auth_asym_id` and LSU
   main rRNA chain `auth_asym_id`.
2. Look up the triple `(annotation.pdb_id.upper(),
   lsu_chain.auth_asym_id, ssu_chain.auth_asym_id)`.
3. On a unique match → populate `large_scale_movements` with the
   row's metrics + the cached file's `rad_date`.
4. On a miss → emit `large_scale_movements` with non-null
   `rad_date` (since the dataset is loaded) and null metrics. Log
   an `info` message.
5. **Multi-chain assemblies** (more than one SSU or LSU main rRNA
   chain — rare; organellar fragments or double-ribosome packing):
   try every (lsu, ssu) pair the assembly contains and return
   metrics only when **exactly one** pair matches a RADdb row.
   Otherwise emit null metrics and add a
   `raddb_ambiguous_chain_pair` warning.

### 29.8 JSON output schema

Add a `large_scale_movements` block to each annotated assembly's
JSON object. The block is **always emitted** for annotated assemblies
— never absent — so the consumer schema is stable across runs.

Three shapes:

```jsonc
// Match found
"large_scale_movements": {
  "source": "RADdb",
  "rad_date": "20260508",
  "intersubunit_rotation": 5.8,
  "ssu_head_rotation": 9.4
}

// RADdb dataset loaded, but no row matches this (pdb, lsu, ssu) triple
"large_scale_movements": {
  "source": "RADdb",
  "rad_date": "20260508",
  "intersubunit_rotation": null,
  "ssu_head_rotation": null
}

// RADdb dataset entirely unavailable (no cache + download failed)
"large_scale_movements": {
  "source": "RADdb",
  "rad_date": null,
  "intersubunit_rotation": null,
  "ssu_head_rotation": null
}
```

Skipped / failed annotations do NOT carry the
`large_scale_movements` block (it is irrelevant when no rRNA chains
were assigned).

`large_scale_movements` is **not** mirrored into the chain-level or
assembly-level CSV outputs in v1. It is JSON-only.

### 29.9 API surface

#### 29.9.1 Library

```python
annotate_pdb(
    pdb_id: str,
    *,
    # ... existing params ...
    raddb_dataset: RADdbDataset | None = None,
    refresh_raddb: bool = False,
    no_raddb: bool = False,
) -> list[RibosomeAnnotation]
```

- `raddb_dataset` — pre-loaded dataset; takes precedence over the
  refresh path. Useful for batch runs that want to share one load
  across many calls.
- `refresh_raddb` — force an online check at the start of the call
  (overrides the 7-day window). No-op when `raddb_dataset` is
  supplied.
- `no_raddb` — skip RADdb integration entirely. The JSON still
  emits `large_scale_movements` in the unavailable shape so the
  schema stays stable.

`annotate_many` accepts the same three parameters; it loads RADdb
once at the start of the batch and passes the dataset through to
every per-entry `annotate_pdb` call (cost: one load per batch, not
one per entry).

#### 29.9.2 CLI

Add to both `annotate` and `annotate-batch`:

```text
--refresh-raddb     Force an online check for a newer RADdb release.
```

Add a new `raddb` subcommand group with two commands:

```bash
ribostate raddb info       # show cached rad_date, downloaded_at, file size
ribostate raddb refresh    # force a check + download if newer; --force re-downloads
```

The `ribostate cache info` output gains a `raddb` row showing the
file count in the RADdb cache namespace.

### 29.10 Module layout

The integration lives in one new module:

```text
src/ribosome_state_annotator/raddb.py
```

Recommended public surface:

```python
get_raddb_cache_dir(cache_root: Path | None = None) -> Path
get_local_raddb_csv_path(cache_root: Path | None = None) -> Path
get_local_raddb_metadata_path(cache_root: Path | None = None) -> Path

load_raddb_metadata(cache_root=None) -> RADdbMetadata | None
save_raddb_metadata(metadata, cache_root=None) -> Path

find_latest_raddb_url(*, client=None, today=None, lookback_days=60)
    -> tuple[str, str] | None    # (url, rad_date_YYYYMMDD)

download_raddb_csv(url, rad_date, *, cache_root=None, client=None)
    -> tuple[Path, RADdbMetadata]

ensure_raddb_available(*, cache_root=None, force_refresh=False,
                       client=None, now=None) -> RADdbMetadata | None

load_raddb_table(*, cache_root=None, csv_path=None)
    -> list[dict[str, str]] | None

build_raddb_lookup(rows) -> tuple[dict, frozenset]
load_raddb_dataset(*, cache_root=None, metadata=None) -> RADdbDataset | None
get_motion_metrics(dataset, pdb_id, lsu_chain_id, ssu_chain_id)
    -> dict | None
```

The `RADdbDataset` dataclass carries the metadata, the lookup dict,
and the frozenset of duplicate keys.

The corresponding Pydantic model lives in `models.py`:

```python
class LargeScaleMovements(BaseModel):
    source: Literal["RADdb"] = "RADdb"
    rad_date: str | None = None
    intersubunit_rotation: float | None = None
    ssu_head_rotation: float | None = None
```

and `RibosomeAnnotation` gains:

```python
large_scale_movements: LargeScaleMovements | None = None
```

### 29.11 Logging

Required `INFO`-level events:

- RADdb file exists / missing
- checking latest RADdb version
- downloading RADdb (with size + destination)
- using cached RADdb (with rad_date + downloaded date)
- RADdb match found / not found per assembly

Required `WARNING`-level events:

- RADdb download / probe failure
- stale RADdb cache but no newer release in the lookback window
- duplicate `(pdb, lsu, ssu)` keys in the CSV
- missing required columns in the CSV
- ambiguous chain-pair match (multi-chain assemblies, §29.7.4 step 5)

### 29.12 Implementation requirements

- Use the package's existing `httpx` client for all network calls
  (do not introduce `requests`).
- Use stdlib `csv` for parsing (the file is ~500 KB / ~2k rows;
  pandas would add a heavy dependency for no benefit). The optional
  `pandas` extra remains available for downstream consumers.
- All public helpers must have docstrings and full typing.
- All failure paths must return `None` / null metrics rather than
  raising. RADdb integration MUST remain optional and robust.

### 29.13 Tests

Coverage requirements (all unit tests, no live network calls):

1. No local file → triggers download.
2. Local file < 7 days old → uses cache; no HTTP calls.
3. Local file ≥ 7 days old → probes for newer; downloads if newer.
4. Probe / download failure with cached file present → returns
   cached metadata; pipeline continues.
5. Probe / download failure without a cached file → returns `None`;
   downstream emits the "unavailable" JSON shape.
6. Column mapping verifies `body rot.` → `intersubunit_rotation`
   and `head rot.` → `ssu_head_rotation`.
7. Match key uses `(pdb_id.upper(), lsu_chain_id, ssu_chain_id)`,
   chain case preserved.
8. Duplicate keys → warning at load + null metrics at lookup.
9. Missing required columns → warning + dataset treated as
   unavailable.
10. PDB ID uppercase normalization works for lowercased input.

Additionally, `tests/conftest.py` must redirect the RADdb cache
root to an isolated tmp directory via an autouse fixture so a
developer's pre-existing
`~/.cache/ribosome-state-annotator/raddb/` does not contaminate
the test suite.

### 29.14 References

See `REFERENCES.md` for the RADtool / RADdb methodology citation
(Mears *et al.* 2002 and the RADtool project page). The metric
naming follows the established ribosome-conformation literature:
`body rot.` is the canonical inter-subunit ratchet rotation;
`head rot.` is the SSU-head swivel.

## 30. tRNA-mRNA codon/anticodon extraction (FR3D)

This section is the authoritative spec for the per-site codon ↔
anticodon evidence extraction. The feature surfaces, for every
A/P/E-site tRNA assigned by the upstream pipeline, the three anticodon
residues, the three codon residues, and the FR3D base-pair
interactions that connect them — as raw evidence only.

### 30.1 Goal

Add a `trna_mrna_interactions` field on `RibosomeAnnotation` (one
entry per A/P/E site present) that captures:

- The three anticodon residues at biological positions 34, 35, 36.
- The three codon residues at codon positions 1, 2, 3.
- Each FR3D base-pair interaction (`cWW`, `tHS`, etc.) observed
  between a codon residue and an anticodon residue.
- Per-residue modification metadata (the observed CCD code,
  parent base, and `is_modified` flag).
- Per-codon assignment status (`complete` / `partial` / `missing`)
  and per-pair assignment status (`assigned` / `ambiguous`).

This is **evidence-only**. The module MUST NOT:

- classify cognate / near-cognate / non-cognate.
- infer tRNA type or amino acid identity.
- infer canonical Watson-Crick status from `cWW`.

Consumers receive raw FR3D labels and decide.

### 30.2 Run condition

Run only when the annotation already has:

- a non-null `mrna_chain`, AND
- at least one of `aminoacyl_trna_chain`, `peptidyl_trna_chain`,
  `exit_trna_chain` populated.

When the condition fails, the field is emitted as an empty list
(stable consumer schema).

Each site (A, P, E) is processed independently. Missing sites are
skipped silently.

### 30.3 Source

```text
http://rna.bgsu.edu/rna3dhub/pdb/<pdb_id_lower>/interactions/fr3d/basepairs/csv
```

Three-column CSV: `unit_id_1, interaction_type, unit_id_2`. Each
physical base pair appears twice (once per direction); the
implementation must deduplicate.

The endpoint is per-PDB and content-stable, so the raw CSV bytes are
cached under the existing user-cache root:

```text
~/.cache/ribosome-state-annotator/fr3d/<pdb_id_lower>.csv
```

Cache invalidation follows the same content-addressed never-expires
policy as the other §17 namespaces.

### 30.4 Anticodon residue identification

The three anticodon residues are the polymer residues at
`auth_seq_id` **34, 35, 36** of the tRNA chain. The pick is anchored
on the canonical first residue (auth_seq_id 1) rather than on the
34th element of the polymer list. This matters when a deposited
chain has a pre-residue numbered `0` (or negative) at the 5′ end —
for example 5UYM chain W, whose polymer runs 0..75: the "34th
polymer residue" would land on auth_seq_id 33 and silently miss the
true Sprinzl-34 wobble residue.

Procedure:

1. Filter the Gemmi chain to `EntityType.Polymer` (skips Mg / HOH /
   modified-nucleotide-as-ligand HETATMs that share the chain ID).
2. From that polymer set, take the residues whose `seqid.num` is
   exactly 34, 35, and 36.
3. If any of those is missing (chain truncated, organellar
   mt-tRNA with renumbering, anticodon-stem-loop fragment), skip the
   site entirely and add a warning to the assembly.

Record the source as `anticodon_position_source: "polymer_sequence_index"`
(or `"auth_seq_id_fallback"` for a future Sprinzl-aware version that
falls back when its lookup fails). v1 always uses the auth_seq_id
34/35/36 anchor described above.

### 30.5 Unit-ID construction

For both anticodon residues and codon residues, construct
BGSU-style unit IDs as:

```text
<pdb_id_upper>|1|<chain>|<observed_chem_comp_id>|<auth_seq_id>
```

Critical: the unit ID's residue number is the residue's **author
sequence number** (`seqid.num`), NOT the polymer-sequence-index. FR3D
unit IDs embed deposited residue numbers; the biological position
(34/35/36) lives separately in the `trna_position` field. The
chem_comp_id is the *observed* CCD code, so modified residues keep
their modified name (`PSU`, `5MU`, `MIA`, etc.) — required for matching
FR3D rows verbatim.

### 30.6 Residue chemistry metadata

For each anticodon residue, record:

- `parent_base`: the canonical base (`"A" / "C" / "G" / "U"`).
- `trna_chem_comp_id`: the observed CCD code, verbatim.
- `is_modified`: True iff the residue is not a canonical base.

Resolution order:

1. **Gemmi tabulated dictionary**
   (`gemmi.find_tabulated_residue(comp).one_letter_code`). Canonical
   bases surface as uppercase (`"A"` → `"A"`); modified nucleotides
   Gemmi knows about surface as lowercase (`"PSU"` → `"u"`).
2. **Per-component CCD fetch** when Gemmi returns blank /
   whitespace-only / missing — the package downloads
   `https://files.rcsb.org/ligands/view/<comp>.cif`, parses the
   ``_chem_comp.mon_nstd_parent_comp_id`` (preferred) or
   ``_chem_comp.one_letter_code`` field, and uses that as the
   authoritative parent base. Cached on disk under the `ccd/`
   namespace; one network call per unrecognized comp_id per cache
   lifetime. The CCD is also the source of truth for unusual
   modifications that Gemmi ships with a blank one-letter code (e.g.
   `U8U` = 5-methylaminomethyl-2-thiouridine-5′-monophosphate, a
   *Thermus* tRNA wobble modification).
3. **First-character heuristic** as a last-resort fallback (network
   failure, parse failure): use `comp_id[0]` and flag `is_modified =
   len(comp_id) > 1`. Works for most CCD codes because the PDB
   convention names modified nucleotides starting with the parent
   base letter (`PSU` → U, `7MG` → G).

For each codon residue:

- `base`: the residue name (canonical bases are typically what mRNA
  carries; modified bases are rare on coding mRNA).
- `source`: `"fr3d_observed"`, `"mmcif_reconstructed"`, or
  `"mrna_frame_inference"` per §30.9.

### 30.7 FR3D codon-pairing fallback for missed sites

Contact-transfer (`infer.assign_functional_chains`, §11) compares
candidate tRNA chains against the canonical SSU decoding-centre
monitor bases (E. coli A1492/A1493 → yeast 18S A1824/A1825). Those
monitor bases only flip out and contact the codon-anticodon duplex
once the aa-tRNA has been **accommodated** — in pre-accommodation
states (eEF1A·GTP·aa-tRNA decoding intermediates, e.g. PDB **3JAG**),
the anticodon is base-paired to the mRNA codon but the rRNA
fingerprint hasn't formed yet. Contact-transfer correctly refuses to
claim the chain, but the result is an unassigned A-site even though
codon-anticodon pairing is unambiguous in FR3D.

The fallback closes this gap. Before per-site extraction, scan every
unassigned tRNA-Rfam chain (those in `other_rna_chains` with Rfam
RF00005) for cWW FR3D pairs to the mRNA chain at the anticodon
residues (auth_seq_id 34/35/36 per §30.4). Disambiguate which empty
site each candidate fills using the **mRNA codon position**:

- The mRNA is translated 5' → 3': A-codon is the most-3'
  (largest auth_seq_id), P is one codon (3 residues) upstream, E is
  two codons upstream.
- Score each candidate by the **mean `auth_seq_id`** of the mRNA
  residues it pairs with.
- Sort candidates descending by that score.
- Assign in that order to the remaining empty slots in canonical
  A → P → E order — so the most-downstream candidate goes to the
  most-A-side empty slot, the next-downstream to the next-A-side, etc.

Each fallback assignment:

- Mutates `annotation.aminoacyl_trna_chain` /
  `peptidyl_trna_chain` / `exit_trna_chain` in place.
- Removes the chain from `other_rna_chains`.
- Appends a warning of the form
  `<site>trna_assigned_from_fr3d_codon_pairing_<chain_id>`.
- Records the assignment under
  `classification_evidence["fr3d_codon_pairing_fallback"][site]` so
  consumers can distinguish FR3D-derived assignments from canonical
  contact-transfer assignments.

The fallback runs only when:

- An mRNA chain is present.
- At least one of A / P / E is empty.
- An unassigned tRNA-Rfam chain (`RF00005` in `rfam_accessions`)
  exists in `other_rna_chains`.
- That candidate makes ≥ 1 cWW pair with the mRNA at one of the
  anticodon residues (non-cWW pairs alone don't trigger assignment).

The fallback never overrides an existing assignment. It does NOT
recompute tRNA states for the newly-assigned chain in v1 — the state
inference is left to the existing §12 contact-transfer rules, which
will typically report `*/...` (no canonical SSU contact) for
fallback-assigned chains. A future version may add state inference
informed by the FR3D codon-pairing evidence.

### 30.8 FR3D pair filtering

Keep only FR3D rows where:

- One side's chain == the mRNA chain ID.
- The other side's chain == the current tRNA chain ID.
- The tRNA side's `auth_seq_id` is one of the three anticodon
  residues' `auth_seq_id`s (computed in §30.4).

Each physical pair appears twice in FR3D (`A → B` and `B → A`).
Canonicalise as `(mrna_unit_id, interaction, trna_unit_id)` and
deduplicate.

### 30.9 Codon assignment

#### 30.8.1 Anticodon ↔ codon position mapping

Antiparallel base pairing maps tRNA positions to codon positions:

| tRNA position | codon position | Notes |
|---------------|----------------|-------|
| 34            | 3              | wobble position — `is_wobble_position = true` on this pair |
| 35            | 2              | |
| 36            | 1              | |

#### 30.8.2 Primary source: FR3D-observed pairs

For each codon position, the mRNA residue comes from the FR3D pair
with the matching anticodon residue.

If multiple mRNA residues pair with the same anticodon position:

1. Prefer the `cWW` pair (canonical Watson-Crick geometry).
2. If still ambiguous (no `cWW`, or multiple `cWW`), surface one
   candidate and mark `assignment_status = "ambiguous"` on the
   resulting `BasePair` entry.

#### 30.8.3 Secondary source: mmCIF polymer-order reconstruction

If at least one codon position is FR3D-observed but the others
aren't, fill the missing positions from the immediately-adjacent
residues in mmCIF polymer order. The triplet is contiguous on the
mRNA polymer.

Examples (highest seeded codon position drives the others):

- Position 3 observed → positions 1, 2 = the two residues immediately
  preceding position 3 in polymer order.
- Position 2 observed → position 1 = immediately preceding, position
  3 = immediately following.
- Position 1 observed → positions 2, 3 = the two residues immediately
  following position 1.

Mark reconstructed residues as `source: "mmcif_reconstructed"`.

#### 30.8.4 Tertiary source: mRNA frame inference

If **any** site (A, P, or E) is fully resolved (either FR3D-observed
or FR3D + locally reconstructed), the other two sites can be
inferred by stepping the mRNA polymer-order frame. Along the mRNA
read 5' → 3', the three sites are ordered **E < P < A**, with one
codon (3 residues) between each adjacent pair. The polymer-index
offset from the anchor site's codon position 1 to the target site's
codon position 1 is:

```text
offset = 3 * (target_index - anchor_index)
```

where `_SITE_FRAME_INDEX = {"E": 0, "P": 1, "A": 2}`. Examples:

- Anchor A → P: offset = `3 * (1 - 2) = -3` (3 residues upstream).
- Anchor A → E: offset = `3 * (0 - 2) = -6` (6 residues upstream).
- Anchor P → A: offset = `3 * (2 - 1) = +3` (3 residues downstream).
- Anchor P → E: offset = `3 * (0 - 1) = -3` (3 residues upstream).

The implementation selects the anchor in canonical A → P → E priority
order — the most-A-side resolved site is preferred because in routine
elongation states the A-site codon is most likely to carry direct
FR3D evidence. The P-anchored case rescues pre-accommodation states
(e.g. PDB 5UYN, where the A-tRNA is held by EF-Tu and lacks
Watson-Crick contacts to the codon while the P-site is firmly
resolved): anchoring on P, the inferred A and E codons are
3 residues downstream and upstream of the P-codon respectively.

Frame inference is permitted only if:

- The mRNA polymer run spanning the anchor and target codons is
  **continuous** by `auth_seq_id` (no gaps).
- The mRNA chain has at least **9 consecutive polymer residues**
  somewhere in the chain (a conservative check that an E/P/A frame
  can plausibly be spanned).

Otherwise mark the unresolved sites' codons as incomplete and add the
warning `insufficient_mrna_for_frame_inference`. Frame-inferred
residues are tagged with `source: "mrna_frame_inference"`.

#### 30.8.5 Assignment status

Per-codon `assignment_status`:

- `"complete"` — all three positions assigned (by any combination of
  sources). The `sequence` field is populated as the 3-letter parent-
  base string.
- `"partial"` — one or two positions assigned. `sequence` is `null`.
- `"missing"` — no positions assigned. `sequence` is `null`.

Per-pair `assignment_status`:

- `"assigned"` — single unique pair for this codon position.
- `"ambiguous"` — multiple FR3D pairs at this position and no `cWW`
  disambiguator.

### 30.10 JSON output schema

Add `trna_mrna_interactions` (an ordered list, A → P → E) to each
annotated assembly:

```jsonc
"trna_mrna_interactions": [
  {
    "site": "A",
    "mrna_chain_id": "V",
    "trna_chain_id": "Y",
    "anticodon_position_source": "polymer_sequence_index",

    "codon": {
      "sequence": "UUC",
      "assignment_status": "complete",
      "residues": [
        { "codon_position": 1, "unit_id": "5UYM|1|V|U|19", "base": "U",
          "source": "mrna_frame_inference" },
        { "codon_position": 2, "unit_id": "5UYM|1|V|U|20", "base": "U",
          "source": "mrna_frame_inference" },
        { "codon_position": 3, "unit_id": "5UYM|1|V|C|21", "base": "C",
          "source": "fr3d_observed" }
      ]
    },

    "anticodon": {
      "sequence_parent": "GAA",
      "residues": [
        { "trna_position": 34, "unit_id": "5UYM|1|Y|G|34", "parent_base": "G",
          "trna_chem_comp_id": "G", "is_modified": false },
        { "trna_position": 35, "unit_id": "5UYM|1|Y|A|35", "parent_base": "A",
          "trna_chem_comp_id": "A", "is_modified": false },
        { "trna_position": 36, "unit_id": "5UYM|1|Y|A|36", "parent_base": "A",
          "trna_chem_comp_id": "A", "is_modified": false }
      ]
    },

    "pairs": [
      {
        "codon_position": 3, "trna_position": 34,
        "codon_unit_id": "5UYM|1|V|C|21", "trna_unit_id": "5UYM|1|Y|G|34",
        "codon_base": "C", "trna_parent_base": "G",
        "trna_chem_comp_id": "G", "trna_is_modified": false,
        "fr3d_interaction": "cWW", "basepair": "C-G",
        "is_wobble_position": true, "assignment_status": "assigned"
      }
    ],

    "warnings": []
  }
]
```

When the run condition fails, `trna_mrna_interactions: []`. The block
is JSON-only — it is **not** mirrored into the chain-level or
assembly-level CSV outputs in v1.

### 30.11 API surface

#### 30.10.1 Library

```python
annotate_pdb(
    pdb_id: str,
    *,
    # ... existing params ...
    no_fr3d: bool = False,
) -> list[RibosomeAnnotation]
```

`no_fr3d` skips the extraction entirely; the field is still emitted
as an empty list. `annotate_many` accepts the same parameter.

There is no `refresh_fr3d` flag in v1 — the cache is content-
addressed and never expires (same policy as `rcsb/`, `bgsu/`, etc.).
Users wanting fresh data delete the cache or pass `--no-cache`.

#### 30.10.2 CLI

No new flags are required. `ribostate cache info` gains an `fr3d`
row showing the cache entry count.

### 30.12 Module layout

The integration lives in one new module:

```text
src/ribosome_state_annotator/trna_mrna.py
```

Recommended public surface:

```python
fetch_fr3d_basepairs(pdb_id, *, cache=None, client=None)
    -> list[_ParsedFr3dRow] | None

extract_trna_mrna_interactions(annotation, structure, *, cache=None, client=None)
    -> list[TRNAmRNAInteraction]
```

The corresponding Pydantic models live in `models.py`:

```python
class AnticodonResidue(BaseModel): ...
class Anticodon(BaseModel):  ...
class CodonResidue(BaseModel):  ...
class Codon(BaseModel):  ...
class BasePair(BaseModel):  ...
class TRNAmRNAInteraction(BaseModel):  ...
```

and `RibosomeAnnotation` gains:

```python
trna_mrna_interactions: list[TRNAmRNAInteraction] = []
```

### 30.13 Logging

Required `INFO`-level events:

- detected mRNA chain
- detected A/P/E-site tRNA chains
- current site being processed
- anticodon residues identified
- FR3D fetch URL (with cache hit/miss)
- retained mRNA-anticodon basepair count per site
- codon reconstruction / inference status per site

Required `WARNING`-level events:

- tRNA chain has fewer than 36 polymer residues
- FR3D download / parse failure
- ambiguous codon assignment (per codon position)
- insufficient mRNA residues for frame inference
- no FR3D pairs found for a site that has a tRNA chain
- FR3D codon-pairing fallback assigned a previously empty site
  (recorded both in `warnings` and in `classification_evidence`)
- per-component CCD fetch / parse failure (the package falls back to
  the first-character heuristic in that case)

### 30.14 Implementation requirements

- Use the package's existing `httpx` client for the FR3D fetch.
- Use stdlib `csv` for parsing.
- Every failure path must degrade to an empty list / partial result
  rather than raising — the codon/anticodon extraction must never
  block the surrounding annotation pipeline.

### 30.15 Tests

Coverage requirements:

1. Polymer-sequence-index 34/35/36 picker filters HETATM ions.
2. Picker uses polymer index, not `auth_seq_id`, when the tRNA is
   renumbered (chains starting at residue 101 still give a hit at
   biological position 34).
3. Picker returns `None` for chains with fewer than 36 polymer
   residues.
4. Unit ID embeds `auth_seq_id` and observed chem_comp (`PSU`, not
   parent `U`).
5. Parent-base / `is_modified` lookup correctly handles canonical
   (`A` → `("A", False)`) and modified (`PSU` → `("U", True)`)
   residues.
6. FR3D CSV parser tolerates malformed rows.
7. FR3D fetch caches on success and reuses the cache on the second
   call (one HTTP call).
8. Cross-chain pair filter keeps both orientations and deduplicates.
9. Codon assignment maps tRNA positions to codon positions correctly
   and prefers `cWW` when multiple pairs exist for the same anticodon
   position.
10. Local reconstruction fills missing codon positions from polymer
    order in both directions.
11. Frame inference fills P/E codons from a complete A-site codon.
12. Frame inference returns `None` when the mRNA polymer run is too
    short or has gaps.
13. End-to-end extraction returns `[]` when there is no mRNA chain,
    no tRNA chain, or FR3D is unreachable.
14. End-to-end extraction populates all three sites correctly when
    given a synthetic FR3D CSV.
15. Per-site warning is added when frame inference is impossible.
16. Fallback assigns A-site from FR3D evidence when contact-transfer
    left A empty but a tRNA-Rfam chain has cWW codon-anticodon pairs.
17. Fallback distinguishes A vs P vs E using the mRNA codon position
    (descending `auth_seq_id` → A first).
18. Fallback only triggers when ≥ 1 cWW pair exists (non-cWW pairs
    alone do not trigger assignment).
19. Fallback ignores short tRNA chains (< 36 polymer residues) and
    non-tRNA-Rfam chains.
20. Per-component CCD fetch returns parent base for residues Gemmi's
    tabulated dictionary doesn't fully describe (e.g. `U8U`).
21. Per-component CCD fetch caches on disk; second call for the same
    comp_id is served from the cache (no second HTTP request).
22. CCD parse failure / network failure falls back to the
    first-character heuristic without raising.
23. Frame inference anchors on the P-site codon and fills A and E
    correctly when only the P-site has FR3D pairs (e.g. 5UYN
    pre-accommodation state). The most-A-side complete site is used
    as the anchor in A → P → E priority order.

### 30.16 References

See `REFERENCES.md` for the FR3D citation (Sarver et al. 2008) and
the Leontis-Westhof base-pair geometric nomenclature (2001).


## 31. v3.2 robustness and edge-case extensions

This section collects the v3.2 design extensions that emerged from
running the package on a random 200-PDB sample. Each subsection
describes one extension; the in-place updates in earlier sections
(§6.4, §7.4, §8.2, §11.2, §11.6, §12) reference back here.

### 31.1 Multi-ribosome bundle splitting

**Problem.** Some deposits pack two complete ribosomes into a single
biological assembly — Mycobacterium 70S di-ribosomes (8R3V, 8RCL),
plant bacterial dimers (9GFT, 9O3L), CCA-end-fragment dimers (8T8C),
in-situ human di-ribosomes (9B0S). The classification truth table
(§8) and the per-assembly contact-transfer (§11) treat the whole
bundle as one ribosome, filling only one A/P/E set and leaving the
second ribosome's tRNAs in `other_rna_chains`.

**Detection.** A multi-ribosome bundle is defined by:

```text
n_ssu := len(by_role["ssu_main_rrna"])
n_lsu := len(by_role["lsu_main_rrna"])

is_multi_ribosome := (n_ssu >= 2) AND (n_lsu >= 2) AND (n_ssu == n_lsu)
```

The equal-counts requirement distinguishes a multi-ribosome bundle
from a fragmented-LSU deposit (§7.4 / §31.3) where `n_ssu != n_lsu`.

**Pairing rule (greedy nearest-centroid).** For each SSU and LSU
chain, compute the centroid of all atom positions. Greedily pair the
SSU↔LSU pair with the smallest centroid distance, remove both from
contention, repeat. The result is a deterministic list of
`(ssu_chain, lsu_chain)` pairs whose length equals
`min(n_ssu, n_lsu)`.

The greedy rule is robust for the cases observed because the two
ribosomes' centroids are typically 100–200 Å apart while
within-ribosome SSU↔LSU centroid distance is ~50 Å — the "right"
pairings beat the "wrong" cross-pairings by an order of magnitude.

**Partitioning other chains.** Every non-main-rRNA chain (mRNA, tRNAs,
5S/5.8S, non-ribosomal proteins) is assigned to the ribosome whose
combined SSU+LSU centroid it is closest to. Each ribosome receives:

- its paired SSU + LSU main rRNA chains;
- the non-main-rRNA chains in its centroid partition;
- the non-ribosomal protein chains in its centroid partition.

**Per-ribosome correspondence.** BGSU's NR correspondence response
typically returns mappings to **one** chain per PDB per equivalence
class — the canonical NR representative (see the live-API behaviour
documented in our examples to BGSU). For the bundle this means BGSU
returns anchors only for ribosome 1's chains; ribosome 2's chains
have no native mapping in the response. The package rebuilds a
per-ribosome `CorrespondenceResult` by re-applying the §5.2.2
chain-substitution fallback with the ribosome's own SSU/LSU chains as
the substitution targets. For same-organism / same-numbering cases
(all bacterial dimers in the test set) this produces valid anchors;
for cross-organism cases (9B0S queried with yeast anchors) the
substitution preserves the reference residue number which may not
exist in the target chain — a limitation noted as a follow-up that
resolves once BGSU returns all chain instances per PDB.

**Output shape.** Each ribosome emits its own `RibosomeAnnotation`
with a suffixed `assembly_id`:

```text
8R3V assembly 1 → annotations with assembly_id "1-1" and "1-2"
```

Single-ribosome assemblies (the 99% case) keep their plain
`assembly_id`. A warning `multi_ribosome_bundle_split` is recorded on
each sub-annotation. `annotate_assembly(pdb_id, "1-1")` retrieves a
specific sub-ribosome; `annotate_assembly(pdb_id, "1")` returns the
first sub-ribosome of the bundle for backward compatibility.

**Module.** `src/ribosome_state_annotator/multiribo.py` owns the
detection (`detect_multi_ribosome`), pairing
(`pair_ssu_lsu_by_centroid`), and partitioning
(`partition_chains_by_ribosome`).

### 31.2 LSU-based A/P fallback

**Problem.** The canonical SSU-anchor pass (§11.2–§11.4) misses two
classes of tRNA:

1. **CCA-end-only analogs.** Synthetic peptidyl-tRNA / aminoacyl-tRNA
   mimics (~10–30 nt) that engage only the PTC and lack the anticodon
   end. Examples: 7RQA chains 1w / 1x / 2w / 2x; 8T8C chains 1w / 1x.
2. **Pre-accommodation full tRNAs.** Full-length tRNAs delivered into
   the PTC before anticodon-codon pairing in the SSU decoding centre.
   Examples: 3JAG chain 2 (tRNA-Val); 7O7Z chain AT; 7OSM chain ASIT;
   7UG7 chain Pt.

Both have strong PTC contacts (`lsu_atrna` 2–3 Å, `lsu_ptrna` 2–3 Å)
but `ssu_atrna` / `ssu_ptrna` distances of 10–60 Å — outside the 5 Å
cutoff.

**Algorithm.** After the canonical SSU pass fills mRNA → P → A → E,
collect each remaining candidate's distances to `lsu_atrna` and
`lsu_ptrna`. For each candidate that contacts either site within
`cutoff`:

```text
preferred_site := A if d(lsu_atrna) <= d(lsu_ptrna) else P
```

Process candidates in best-first order (smallest LSU contact first).
For each:

1. If the preferred site is still unfilled, assign the candidate there.
2. Else if the alternative LSU site is still unfilled AND the
   candidate also contacts that site within cutoff, assign there.
3. Else skip.

The closer-site-first ordering is essential: with strict A-before-P
ordering, 8OIR's P-Met-tRNA (`lsu_atrna` 3.2 Å, `lsu_ptrna` 1.8 Å)
would be mis-assigned to A because the A-fallback runs first.

E-tRNA already uses LSU anchors (`lsu_etrna`) as its primary anchor
in §11.5, so no separate E-tRNA fallback is needed.

**Helper.** `infer._compute_chain_distances` returns the full
`{chain_id: min_distance}` map for an anchor site; the fallback uses
it twice (once each for `lsu_atrna` and `lsu_ptrna`) and combines.

### 31.3 Fragmented-ribosome skip

**Problem.** Some organisms biologically fragment their LSU rRNA into
multiple chains: Tetrahymena (8OVA, 8OVJ, 8QHU, 8QIE, 8RXH),
Chlamydomonas chloroplast (28JW, 28LU, 9TVU, 9HL9), some fragmented
human mitoribosomes (7A5K). The deposit's `lsu_main_rrna` role
collects all fragments (e.g. 1 SSU + 3 LSU). Canonical anchor
residue numbers don't transfer onto fragments — each fragment uses
its own per-fragment numbering, and BGSU's NR correspondence rarely
covers these unusual rRNAs.

**Detection and skip.** At the start of `_annotate_one_assembly`,
after `partition_rna_chains_by_role` produces `by_role` and before
the multi-ribo splitter:

```text
n_ssu := len(by_role["ssu_main_rrna"])
n_lsu := len(by_role["lsu_main_rrna"])

is_fragmented := (n_ssu >= 2 OR n_lsu >= 2) AND n_ssu != n_lsu
```

Skip with:

```json
{
  "status": "skipped",
  "skip_reason": "fragmented_ribosome_not_supported"
}
```

An `INFO` log line `"fragmented ribosome detected for {pdb} assembly
{aid} (SSU chains=N, LSU chains=M); skipping"` is emitted so batch
runs surface the reason without forcing the caller to inspect each
annotation.

Multi-ribosome bundles (`n_ssu == n_lsu >= 2`) bypass this check
because they're handled by §31.1.

### 31.4 Filtered E. coli reference set for organellar

**Problem.** Routing organellar classifications through the full
`ECOLI_REFERENCE_UNITS` triggers BGSU's intersection semantic on
batched correspondence queries — four E. coli anchors (`G|693`,
`A|694`, `G|2112`, `G|2421`) have no mt-rRNA equivalent (residues in
helices that mt-rRNA lost during reduction). Their inclusion in a
batched query drops every mt deposit from the response.

**Solution.** `ECOLI_REFERENCE_UNITS_ORGANELLAR` is the same as
`ECOLI_REFERENCE_UNITS` minus those four anchors:

| Site         | Anchors retained               | Anchors removed             |
|--------------|--------------------------------|-----------------------------|
| `ssu_mrna`   | G\|926, 4OC\|1402, C\|1403     | —                           |
| `ssu_ptrna`  | G\|1338, A\|1339, C\|1400      | —                           |
| `ssu_atrna`  | G\|530, A\|1492, A\|1493       | —                           |
| `ssu_etrna`  | *(empty)*                      | G\|693, A\|694              |
| `lsu_atrna`  | G\|2553, G\|2583, U\|2585      | —                           |
| `lsu_ptrna`  | OMG\|2251, G\|2252, G\|2253    | —                           |
| `lsu_etrna`  | C\|2422                        | G\|2112, G\|2421            |

With the filtered set, BGSU's batched query returns mt mappings for
all 16 remaining anchors across all known human / Toxoplasma /
*spinach* / *Chlamydomonas* mitoribosome and chloroplast ribosome
deposits in the sample, in each deposit's native residue numbering
(8ANY-style 1557 / canonical 909 / 7A5G-style 1561, etc.) — BGSU's
R3D structural alignment cross-walks correctly across both Rfam
families and numbering schemes.

The empty `ssu_etrna` site means the SSU E-site half-state for
organellar P-tRNAs collapses (the `pe` SSU hybrid label cannot be
detected); E-tRNA *assignment* (§11.5) is unaffected because it uses
`lsu_etrna` only.

### 31.5 MIXED rRNA-core handling — superseded by §32

**Historical context.** PDBe's HMM-based Rfam mapping over-annotates
some rRNA chains with multiple cross-family hits. 9B0S' 18S rRNA
chain was tagged with RF00177 (bacterial 16S) + RF01959 (archaeal
SSU) + RF01960 (eukaryotic 18S); its 28S had RF02540 (archaeal LSU)
+ RF02541 (bacterial 23S) + RF02543 (eukaryotic 28S). This triggered
`determine_rrna_core` to return `"mixed"`. The v3.0 spec demoted
MIXED to `"bacterial_like"` unconditionally, which combined with the
Eukaryota protein vote routed 9B0S into
`eukaryotic_organellar_ribosome` — the wrong classification (it's a
human cytoplasmic 80S di-ribosome).

An interim v3.2-rc1 fix delayed the MIXED demotion until after the
§8.3 protein vote so 9B0S would resolve to eukaryotic based on the
protein evidence. That fix was reverted in v3.2 — see §32. The
production Rfam source (EBI `pdb_full_region.txt.gz`, single best
bit-score per chain) yields at most one Rfam accession per chain,
which makes MIXED unreachable from real data. The defensive demotion
in `classify_assembly` is preserved as a fallback (returns
`bacterial_like`) for synthetic / hand-crafted `ChainRef`s in tests
that supply multi-tag annotations directly.

### 31.6 Description-based mRNA-pool exclusion

**Problem.** Some deposits omit `RF00005` on tRNA chains despite the
chain being a tRNA. Two patterns observed:

- Mitoribosome A/A-tRNAs (e.g. 7QI5 chain Aw, labelled "A/A-tRNA"
  with empty `rfam_accessions`).
- Synthetic CCA-end tRNA analogs (e.g. 7RQA chains 1w/1x, labelled
  "A-site Aminoacyl-tRNA Analog" with no Rfam tag).

Without the `RF00005` tag, these chains enter the mRNA candidate
pool and can win the mRNA slot via geometric proximity to the SSU
mRNA anchors (which sit at the decoding-centre).

**Fix.** Extend the §11.2 mRNA-pool exclusion to also drop any chain
whose `description` contains the substring `"tRNA"`
(case-insensitive). The check runs alongside the `RF00005` filter,
before the per-pool ranking.

```python
def _looks_like_trna(chain: ChainRef) -> bool:
    return "trna" in (chain.description or "").lower()

mrna_candidates = [
    c for c in candidate_pool
    if c.ife not in trna_rfam_ifes and not _looks_like_trna(c)
]
```

False-positive risk (an mRNA chain literally described as
"tRNA-binding-site mimic") is acceptable: mRNA assignment then
silently skips and the chain falls into `other_rna_chains`, where
the consumer can re-inspect.

### 31.7 Fragment-vs-displaced state vocabulary

The half-state placeholders in §12 are unified across SSU and LSU
sides with the following convention:

| Symbol | Meaning |
|--------|---------|
| `A` / `P` / `E` | canonical contact at that site |
| `AP` / `PE` | LSU hybrid (uppercase) |
| `ap` / `pe` | SSU hybrid (lowercase) |
| `*` | **full polymer (≥ 30 nt) but no canonical contact at this subunit** — positionally displaced (pre-accommodation A-tRNA, mt P-tRNA with displaced SSU contact, …) |
| `**` | **chain shorter than 30 nt** — physically cannot reach this subunit (anticodon-stem-loop fragments at SSU, CCA-end analogs at PTC) |
| `<factor name>` (LSU only) | full polymer, no LSU rRNA contact, but a non-ribosomal protein factor at the CCA end (§12.4) |

Concrete real-world state strings the package now produces:

| State | Example deposits |
|-------|------------------|
| `A/A`, `P/P`, `E/E` | classical states (5UYM, 7K00) |
| `A/<factor>` | A-tRNA in EF-Tu / RF complex (5UYM A-tRNA Tu, 5KPX RelA) |
| `ap/AP` | EF-G-trapped translocation intermediate (5UYN) |
| `*/AP` | pre-accommodation full A-tRNA at PTC (3JAG, 7O7Z, 7OSM, 7UG7) |
| `**/A`, `**/P` | CCA-end-only tRNA analog at PTC (7RQA, 8T8C) |
| `*/E` | mt-rRNA E-site tRNA (no SSU exit-site anchor available in §31.4's filtered set) |
| `**/**` | (rejected — §31.8) |

### 31.8 No-contact-fragment safeguard

A state of `**/**` would mean a chain shorter than 30 nt that makes
no anchor contact on either subunit. The §11.2–§11.5 / §11.6
assignment passes require a ≤cutoff anchor contact for any role, so
this state shouldn't appear from a normal run. A defensive
post-processing step in `api._demote_no_contact_fragments` checks
each assigned tRNA's state after `compute_trna_states` and demotes
to `other_rna_chains` any chain whose state begins with `**/` and
ends with `**` (i.e. `**/**`), and any E-tRNA whose state begins
with `**/`.

The check guards against hand-crafted `ChainAssignments` (used in
unit tests) and any future code path that bypasses the contact
requirement.

### 31.9 Open follow-up — BGSU all-chains-per-PDB

The §31.1 chain-substitution path for multi-ribosome bundles is a
workaround for a BGSU API behaviour: `map_across_chains` currently
returns **one chain instance per PDB per Rfam equivalence class**
(the NR representative), even when the deposit contains multiple
instances of that chain. For bacterial same-organism dimers the
workaround succeeds because reference and target share residue
numbering; for cross-organism queries (yeast anchor → human 9B0S)
the workaround fails because the preserved reference seqid doesn't
exist in the target chain's numbering. When BGSU is updated to
return all chain instances per PDB, the chain-substitution shim and
the per-ribosome correspondence rebuild can both be deleted; the
existing per-assembly chain filter (§5.2.2) will then pick out the
correct chain for each ribosome automatically.

Examples of the one-chain-per-PDB pattern (BGSU input → returned mt
mapping):

| PDB  | Query                  | Returned                  | Missing chain         |
|------|------------------------|---------------------------|-----------------------|
| 8R3V | `5J7L\|1\|AA\|A\|1492` | `8R3V\|1\|A1\|A\|1492`    | A2 (second ribosome)  |
| 8R3V | `5J7L\|1\|DA\|G\|2553` | `8R3V\|1\|72\|G\|2553`    | 71                    |
| 9O3L | `5J7L\|1\|AA\|A\|1492` | `9O3L\|1\|1a\|A\|1492`    | 2a                    |
| 9O3L | `5J7L\|1\|DA\|G\|2553` | `9O3L\|1\|1A\|G\|2553`    | 2A                    |
| 9B0S | `7ZW0\|1\|2\|G\|1150`  | `9B0S\|1\|s2\|G\|1207`    | S2                    |
| 9B0S | `7ZW0\|1\|LA\|G\|2922` | `9B0S\|1\|l5\|G\|4499`    | L5                    |
| 7AZO | `5J7L\|1\|AA\|A\|1492` | `7AZO\|1\|16SB\|A\|2115`  | 16SA (assembly 1 has 16SA, BGSU returned a chain from assembly 2) |

## 32. Rfam ``pdb_full_region`` integration (rRNA Rfam source)

### 32.1 Goal

Provide a single best-bit-score Rfam accession per `(pdb_id, chain)`
for every PDB rRNA / tRNA chain, replacing the per-entry PDBe REST
endpoint specified in §5.3 v3.1. The new source resolves the MIXED
rRNA-core edge case (§31.5 / §8.2) at the data layer: single-best
selection means each chain carries at most one Rfam family, so the
`determine_rrna_core` function cannot produce a MIXED verdict from
real data.

### 32.2 Source

EBI publishes the file at:

```
https://ftp.ebi.ac.uk/pub/databases/Rfam/.preview/pdb_full_region.txt.gz
```

Format: tab-separated, gzipped, ~200 KB compressed, ~28 k rows
spanning ~4400 distinct PDBs. Columns:

| # | Field | Notes |
|--:|-------|-------|
| 1 | `rfam_acc`         | e.g. `RF00177` (bacterial 16S) |
| 2 | `pdb_id`           | lowercase 4-character PDB ID |
| 3 | `auth_chain_id`    | case-sensitive author chain identifier |
| 4 | `seq_start`        | start residue (sequence-based) |
| 5 | `seq_end`          | end residue |
| 6 | `bit_score`        | cmsearch alignment bit-score |
| 7 | `e_value`          | cmsearch e-value |
| 8 | `cm_start`         | covariance-model alignment start |
| 9 | `cm_end`           | covariance-model alignment end |
| 10 | `color`           | display colour (ignored) |
| 11 | `strand`          | orientation flag (ignored) |

The package reads columns 1, 2, 3, 6 and ignores the rest.

### 32.3 Local cache layout

```
~/.cache/ribosome-state-annotator/rfam/
├── pdb_full_region.txt.gz                    # the latest file
└── pdb_full_region.metadata.json             # sidecar
```

The metadata sidecar records:

```json
{
  "source_url": "https://ftp.ebi.ac.uk/pub/databases/Rfam/.preview/pdb_full_region.txt.gz",
  "downloaded_at": "2026-05-19T12:00:00+00:00",
  "last_modified": "Mon, 18 May 2026 09:37:40 GMT",
  "local_file": "pdb_full_region.txt.gz"
}
```

### 32.4 Refresh policy

Mirrors the §29.4 RADdb policy:

1. **Cache missing** → download immediately.
2. **Cache present + force-refresh flag** → HEAD-probe the upstream
   URL; compare `Last-Modified` against the cached value; download if
   they differ.
3. **Cache present + age ≥ 7 days** (`STALENESS_DAYS`) → same HEAD
   probe + conditional re-download.
4. **Cache present + age < 7 days** → use the cached file, no network.
5. **Any network failure with a local file** → keep the local file and
   log a warning.
6. **Any network failure without a local file** → return `None`; the
   pipeline continues with RCSB-supplied Rfam tags only.

The HEAD probe + `Last-Modified` comparison is cheaper than a full
download. EBI's `Last-Modified` is a stable identifier of the file
version.

### 32.5 Parsing and best-score selection

```python
def build_rfam_pdb_region_lookup(
    rows: list[tuple[str, str, str, float]],   # (rfam, pdb, chain, bit_score)
) -> dict[tuple[str, str], str]:
    """Pick the single best-score Rfam per (pdb_id_lower, chain_id).

    Ties are broken by input row order (rare in practice).
    """
    best_score: dict[tuple[str, str], float] = {}
    best_acc:   dict[tuple[str, str], str]   = {}
    for rfam_acc, pdb_id, chain_id, bit_score in rows:
        key = (pdb_id, chain_id)
        if key not in best_score or bit_score > best_score[key]:
            best_score[key] = bit_score
            best_acc[key]   = rfam_acc
    return best_acc
```

### 32.6 Why single-best disambiguates MIXED at the source

PDBe's HMM-based scanner hits multiple cross-family Rfam models on the
same rRNA chain because the models share covariance ancestry. Example:
9B0S chain `s2` (human cytoplasmic 18S, 1740 nt) reports four hits:

| Rfam   | family            | bit_score |
|--------|-------------------|----------:|
| RF00177 | bacterial 16S    |    456.3 |
| RF01959 | archaeal SSU     |    581.4 |
| **RF01960** | **eukaryotic 18S** | **1669.0** |
| RF02542 | (broad rRNA CM)  |    822.0 |

The bit-score is dominated by the biologically-correct family.
Best-score selection picks `RF01960` cleanly. With the v3.1 PDBe REST
source all four Rfam accessions were merged onto the chain, which
combined with the matching pattern on the 28S triggered `MIXED` in
`determine_rrna_core` (§8.2) and routed 9B0S into
`eukaryotic_organellar_ribosome`. The MIXED-protein-vote fix in §31.5
was a workaround for that downstream symptom; the v3.2 single-best
selection eliminates the symptom at the source.

### 32.7 API surface

```python
from ribosome_state_annotator.rfam_pdb_region import (
    RfamPdbRegionDataset,
    RfamPdbRegionMetadata,
    ensure_rfam_pdb_region_available,
    load_rfam_pdb_region_dataset,
    get_rfam_for_chain,
    get_rfam_mapping_for_pdb,
)
```

`get_rfam_mapping_for_pdb(dataset, pdb_id) -> dict[str, list[str]]`
returns `{auth_chain_id: [single_rfam_accession]}` — wire-compatible
with the legacy PDBe REST parsed shape, so `api.py`'s
`_apply_rfam_pdb_region` is a one-line replacement of the old
`_augment_chain_rfam`.

`annotate_pdb` and `annotate_many` accept:

- `rfam_dataset: RfamPdbRegionDataset | None` — pre-loaded dataset to
  thread through a batch.
- `refresh_rfam: bool` — force an upstream check on this run.
- `no_rfam: bool` — skip the file entirely (RCSB tags only).

### 32.8 CLI

Mirrors the RADdb CLI surface:

```text
ribostate rfam info       # show cache location, last_modified, downloaded_at, bytes
ribostate rfam refresh    # force HEAD-probe + conditional download
```

The `annotate` and `annotate-batch` subcommands accept
`--refresh-rfam` (same semantics as `--refresh-raddb`).

### 32.9 Tests

`tests/test_rfam_pdb_region.py` covers:

- Path helpers and metadata round-trip.
- Parser tolerance for malformed rows (skip silently).
- Best-bit-score selection across overlapping Rfam hits (9B0S vignette).
- `get_rfam_for_chain` lookup behaviour (PDB ID case-insensitive,
  chain case-sensitive, unknown PDB returns `None`).
- `ensure_rfam_pdb_region_available` flow: download-when-missing,
  reuse-when-unchanged, refresh-when-stale, return-None-when-offline.

### 32.10 Migration note

The v3.1 PDBe REST machinery (module `pdbe_client.py`, cache namespace
`pdbe/`, helper `api._augment_chain_rfam`) is removed from the
annotation pipeline. The `pdbe/` cache namespace remains in `cache.py`
so existing user caches don't crash on the new version; the entries
inside it are no longer consulted and can be deleted at the user's
convenience via `ribostate cache clear`.
