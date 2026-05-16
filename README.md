# ribosome-state-annotator

Python package that annotates the functional state of complete
bacterial and eukaryotic ribosome assemblies from PDB entries. Usable
as a CLI (`ribostate`) or as a library
(`import ribosome_state_annotator`).

## Status

v1, pre-release. See
[`ribosome_functional_annotation_package_spec_v3.md`](./ribosome_functional_annotation_package_spec_v3.md)
for the full design.

## Short summary

### Biological context

The ribosome is the macromolecular machine that translates messenger
RNA into protein. During translation it cycles through structurally
distinct functional states — initiation, decoding, peptide-bond
formation, translocation, termination, and recycling — each
characterised by a specific arrangement of mRNA, transfer RNAs, and
auxiliary protein factors (e.g. EF-Tu, EF-G, RF1/2, RRF) bound to the
small (SSU) and large (LSU) subunits. The three tRNA binding sites —
**A** (aminoacyl), **P** (peptidyl), and **E** (exit) — act as a
fingerprint for where in the elongation cycle a given structure was
captured.

A tRNA's functional state is reported as `<SSU>/<LSU>` — its contact
with the small-subunit site followed by its contact with the
large-subunit site:

- **Classical states** (`A/A`, `P/P`, `E/E`) — the tRNA occupies the
  same site on both subunits. These are the stable resting positions
  before and after each round of elongation.
- **Hybrid states** (`A/P`, `P/E`) — after peptide-bond formation but
  before translocation, the acceptor end of each tRNA spontaneously
  swings into the next LSU site while the anticodon end stays put on
  the SSU. Hybrid configurations are direct evidence of a
  post-peptidyl-transfer, pre-translocation intermediate.
- **Chimeric / intermediate states** (lowercase letters, e.g.
  `ap/AP`) — the tRNA simultaneously contacts two adjacent sites,
  catching it mid-translocation. Lowercase signals partial contact;
  uppercase signals full contact.
- **Factor-bound A-site** (`A/Elongation factor Tu`,
  `A/Release factor 1`, …) — when the A-tRNA's CCA end is held by a
  protein factor rather than docked into the LSU peptidyl-transferase
  centre, the LSU half of the label is replaced by the factor's name.
  Typical of decoding-state ternary complexes (EF-Tu·GTP·aa-tRNA) and
  termination complexes.
- **`*` and `**`** — `*` means no contact was found at that subunit;
  `**` means contact was found but couldn't be labelled.

The PDB now contains thousands of ribosome structures spanning
bacteria, archaea, and eukaryotes (cytoplasmic and organellar). Chain
naming conventions vary widely across deposits, and identifying which
chain is the A-site tRNA vs the P-site tRNA vs the mRNA from a raw
mmCIF is non-trivial. This package automates that annotation.

### How the package does it

Functional-site nucleotides on the rRNA — the residues that directly
contact mRNA, A-tRNA, P-tRNA, and E-tRNA — are among the most
phylogenetically conserved positions in all of biology. The package
exploits that conservation:

1. **Fetch + classify.** Pull entry metadata from RCSB GraphQL
   (augmented with PDBe REST Rfam mappings for rRNA chains) and
   classify the assembly as bacterial, eukaryotic cytoplasmic, or
   eukaryotic organellar.
2. **Project reference anchors.** Curated functional-site nucleotides
   from a reference ribosome (*E. coli* 5J7L for bacterial /
   organellar, yeast 7ZW0 for eukaryotic cytoplasmic) are mapped onto
   the query's rRNA via BGSU's RNA correspondence API, which uses
   Rfam-based covariance alignments to find evolutionarily equivalent
   positions even across distant organisms.
3. **Detect contacts.** Download the biological-assembly mmCIF and use
   Gemmi neighbour search to find which chains in the query
   physically contact each set of projected anchors.
4. **Assign + label.** The contacting chains become the mRNA / A-tRNA
   / P-tRNA / E-tRNA assignments. The SSU/LSU contact pattern for
   each tRNA determines its functional state (classical, hybrid, or
   chimeric — see vocabulary above), and a neighbour search around
   the A-tRNA's CCA end identifies the bound elongation/release
   factor when one is present.

This is the **contact-transfer annotation** workflow: instead of
relying on chain names or sequence-only heuristics, the package
transfers functional identity through three-dimensional contacts to
conserved positions.

## Workflow

```text
                       ┌──────────────────────────┐
                       │  PDB ID (e.g. 5UYM)      │
                       └────────────┬─────────────┘
                                    ▼
        ┌───────────────────────────────────────────────────┐
        │  RCSB GraphQL  +  PDBe REST (Rfam fallback)       │
        │  → assemblies, chains, Rfam, taxonomy, ligands    │
        └────────────┬──────────────────────────────────────┘
                     ▼
        ┌───────────────────────────────────────────────────┐
        │  Classify                                         │
        │  rRNA core + ribosomal-protein superkingdom vote  │
        │  → bacterial / eukaryotic / organellar / skip     │
        └────────────┬──────────────────────────────────────┘
                     ▼
        ┌───────────────────────────────────────────────────┐
        │  Pick reference                                   │
        │  5J7L (bacterial/organellar) or 7ZW0 (eukaryotic) │
        └────────────┬──────────────────────────────────────┘
                     ▼
        ┌───────────────────────────────────────────────────┐
        │  BGSU correspondence                              │
        │  reference anchors → query-PDB equivalent units   │
        │  per-assembly chain-substitution fallback applied │
        └────────────┬──────────────────────────────────────┘
                     ▼
        ┌───────────────────────────────────────────────────┐
        │  Download assembly mmCIF + Gemmi neighbour search │
        │  → chains contacting each anchor set              │
        └────────────┬──────────────────────────────────────┘
                     ▼
        ┌───────────────────────────────────────────────────┐
        │  Assign mRNA / A-tRNA / P-tRNA / E-tRNA            │
        │  Infer tRNA states                                │
        │  Label A-site factor by CCA-end neighbour         │
        └────────────┬──────────────────────────────────────┘
                     ▼
        ┌───────────────────────────────────────────────────┐
        │  JSON  +  ribosome_chain_annotation.csv           │
        │        +  ribosome_assembly_annotation.csv        │
        └───────────────────────────────────────────────────┘
```

Multi-assembly PDBs (e.g. 4V5Q) run the classify-through-assign block
**once per biological assembly**, independently, so each assembly gets
its own functional-chain assignment and tRNA states.

## Installation

Requires Python ≥3.10.

```bash
git clone https://github.com/BGSU-RNA/ribosome-functional-annotation.git
cd ribosome-functional-annotation
python3 -m venv .venv
source .venv/bin/activate
pip install -e .          # runtime only
pip install -e ".[dev]"   # add ruff / mypy / pytest
```

## Running the package

The package ships a CLI (`ribostate`) and a Python API. Pick whichever
fits your workflow.

### From the command line

**One PDB → JSON + companion CSVs:**

```bash
ribostate annotate 5UYM -o 5uym.json
# writes 5uym.json, ribosome_chain_annotation.csv, ribosome_assembly_annotation.csv
```

**One PDB → JSON to stdout (no files):**

```bash
ribostate annotate 5UYM --stdout
```

**One biological assembly only:**

```bash
ribostate annotate 4V5Q --assembly-id 1 -o 4v5q_a1.json
```

**Batch over many PDBs:**

```bash
echo -e "5J7L\n7ZW0\n6ZMI" > pdb_ids.txt
ribostate annotate-batch pdb_ids.txt -o batch.json --continue-on-error
```

**Local mmCIF (skip RCSB download):**

```bash
ribostate annotate 5UYM --input-file ./5uym-assembly1.cif -o 5uym.json
```

**Cache maintenance:**

```bash
ribostate cache info        # entry counts per namespace
ribostate cache clear --yes # wipe everything
```

Other useful flags:

| Flag | Effect |
|------|--------|
| `--assembly-id N` | Restrict to one biological assembly. |
| `--stdout` | JSON to stdout, no files written. |
| `--no-csv` | JSON only — suppress the two companion CSVs. |
| `--cutoff 5.0` | Gemmi neighbour-search cutoff (Å). |
| `--cache-dir PATH` | Override `~/.cache/ribosome-state-annotator`. |
| `--no-cache` | Disable caching for this invocation. |
| `--strict` | Skip (don't just warn) assemblies with low ribosomal-protein counts. |
| `--input-file PATH` | Parse a local mmCIF instead of downloading from RCSB. |
| `--verbose` / `--debug` | Increase log verbosity. |

Full help: `ribostate --help`, `ribostate annotate --help`,
`ribostate annotate-batch --help`, `ribostate cache --help`.

## Library usage

```python
from ribosome_state_annotator import annotate_pdb, annotate_assembly, annotate_many

# All assemblies of one entry.
annotations = annotate_pdb("5J7L")

# A specific assembly.
annotation = annotate_assembly("5J7L", "1")
print(annotation.ribosome_classification)        # "bacterial_ribosome"
print(annotation.aminoacyl_trna_chain.ife)       # "5J7L|1|V"
print(annotation.aminoacyl_trna_state)           # "A/A"

# Batch over many PDB IDs (don't abort the batch on per-entry errors).
results = annotate_many(["5J7L", "5TBW", "6ZMI"], continue_on_error=True)
```

Each result is a `RibosomeAnnotation` Pydantic model — see
`src/ribosome_state_annotator/models.py` for the full field list. The role-based rRNA outputs (`ssu_main_rrna_chains` etc.)
are canonical; `ssu_chain` and `lsu_chain` are convenience aliases that
resolve to `None` whenever the underlying list isn't exactly one entry.

JSON / CSV output is also exposed via the `output` module:

```python
from ribosome_state_annotator.output import write_json, write_chain_csv

write_json(annotations, Path("out.json"))
write_chain_csv(annotations, Path("chain.csv"))
```

## Output formats

| Format | When emitted | Notes |
|--------|--------------|-------|
| JSON | Default, or `--output foo.json` | Full `RibosomeAnnotation` list. Field layout in `models.py`. |
| JSONL | `--output foo.jsonl` | One annotation per line. For streaming consumers. |
| `ribosome_chain_annotation.csv` | Default companion (suppress with `--no-csv`) | One row per annotated assembly, 13 columns. Matches the legacy prototype byte-for-byte. |
| `ribosome_assembly_annotation.csv` | Default companion | One row per `(property, value)` tuple: species, non-ribosomal proteins, bound ligands, unmapped RNA chains, plus v1 extensions (classification, superkingdom, warnings). |

Skipped and failed annotations appear in JSON but are omitted from CSV.

## Known limitations

Behaviours v1 accepts but you should know about:

- **`ribosome_state_annotator.classify.matches_ribosomal_protein_narrow`** is
  a plain substring check for `"ribosomal protein"`. It produces:
  - A false positive: `"ribosomal recycling factor"` contains the
    substring and gets flagged. Acceptable in v1 because RRF's
    geometry near the CCA end is similar enough that excluding it from
    the factor search doesn't lose anything important.
  - A false negative: chains named just `"S1"` / `"L11"` without the
    `"ribosomal"` substring are missed. Rare in modern PDB depositions;
    common in older entries. The broader superkingdom-vote rule catches
    these via an anchored regex, but the `is_ribosomal_protein` flag
    used by the factor search stays on the narrow rule.
- The allow-list / deny-list overrides for ribosomal-protein detection
  (`extra_ribosomal_descriptions`, `non_ribosomal_overrides`) are NOT
  exposed in the v1 API. Workaround: rebuild the `ChainRef`s with
  edited `is_ribosomal_protein` flags before passing them to the
  inference layer.

v1 out-of-scope:

- Archaeal ribosomes (skipped with `archaeal_ribosome_not_supported`).
- Ribosome fragments and partial assemblies missing either SSU or LSU.
- Entries solved only by NMR.
- Multiple complete ribosomes packed into one biological assembly.
- Complex symmetry / operator expansion beyond what the deposited
  biological assembly already contains.
- Advanced tRNA state labels beyond the existing contact-based scheme.

## Caching

Every external call is cached on disk at
`~/.cache/ribosome-state-annotator/` (override with `--cache-dir`,
disable with `--no-cache`). Four namespaces: `rcsb/`, `bgsu/`, `pdbe/`,
`coords/`. The cache is content-addressed and never expires — to
refresh, use `ribostate cache clear` or delete the cache root.

## Package layout

The source tree is split by responsibility — each module owns one stage
of the contact-transfer workflow above.

| File | Responsibility |
|------|----------------|
| `api.py` | Top-level orchestration. `annotate_pdb` / `annotate_assembly` / `annotate_many` live here; this is the file most library callers reach for. |
| `cli.py` | Typer CLI (`ribostate annotate`, `…-batch`, `cache`). Thin wrapper over `api.py`. |
| `models.py` | Pydantic v2 models: `ChainRef`, `LigandRef`, `AssemblyContext`, `CorrespondenceResult`, `RibosomeAnnotation`. All JSON / CSV output round-trips through these. |
| `constants.py` | Curated reference unit IDs (`BACTERIAL_REFERENCE_UNITS` for 5J7L, `YEAST_REFERENCE_UNITS` for 7ZW0) — the functional-site anchors. |
| `rcsb_client.py` | RCSB GraphQL fetcher + assembly parser. Reads `assemblies → polymer_entity_instances` and produces per-assembly `AssemblyContext` records. |
| `pdbe_client.py` | PDBe REST Rfam-mappings fetcher. Augments rRNA chains with Rfam accessions RCSB no longer ships. |
| `bgsu_client.py` | BGSU correspondence HTTP client + tolerant JSON parser. Handles the live `mappings` and idealised `alignment` response shapes. |
| `correspondence.py` | PDB-prefix + assembly-chain filter, plus the multi-assembly chain-substitution fallback. |
| `coordinates.py` | mmCIF download + Gemmi parsing. Caches under `coords/`. |
| `gemmi_contacts.py` | Gemmi `NeighborSearch` wrapper restricted to the active assembly. |
| `classify.py` | rRNA-core determination, ribosomal-protein detection, dominant-superkingdom vote, final classification rule. |
| `infer.py` | Functional chain assignment + tRNA-state inference. Owns the polymer-filtered CCA-end selector for the A-site factor label. |
| `cache.py` | Content-addressed on-disk cache with four namespaces: `rcsb/`, `bgsu/`, `pdbe/`, `coords/`. |
| `config.py` | Tunables: contact cutoff, completeness thresholds, network timeouts. |
| `exceptions.py` | Typed exception hierarchy (`ApiRequestError`, `CorrespondenceMappingError`, `CoordinateDownloadError`, …). |
| `output.py` | JSON / JSONL / CSV writers. |

For the full design rationale, classification rules, output schemas,
and the v3.1 live-API addendum, see
[`ribosome_functional_annotation_package_spec_v3.md`](./ribosome_functional_annotation_package_spec_v3.md).

## Development

```bash
ruff format src tests
ruff check src tests
mypy --strict src/
pytest                  # default — excludes network-marked tests
pytest --cov            # with coverage
pytest -m network       # opt in to live RCSB + BGSU + PDBe tests
```

The integration tests under `tests/integration/` hit the live
RCSB / BGSU / PDBe APIs and download mmCIF files; they're excluded from
the default `pytest` run. See `tests/fixtures/acceptance_set.md` for
the PDB IDs used by the acceptance tests.

## Contributors

- Sri Devan Appasamy
- Claude (Anthropic) — pair-programmed the v1 build under direction

## License

MIT.
