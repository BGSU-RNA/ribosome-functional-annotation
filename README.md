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

Given a PDB ID, the package pulls entry metadata from RCSB GraphQL
(augmented with PDBe REST Rfam mappings for rRNA chains), classifies
the assembly as bacterial / eukaryotic / eukaryotic-organellar, and
projects curated **functional-site anchor nucleotides** from a
reference ribosome (*E. coli* 5J7L for bacterial / organellar, yeast
7ZW0 for eukaryotic cytoplasmic) onto the query via BGSU correspondence.
It then downloads the biological-assembly mmCIF and uses Gemmi to find
which chains in the query physically contact those mapped anchors —
giving the mRNA / A-tRNA / P-tRNA / E-tRNA assignments and the SSU/LSU
contact pattern that determines the tRNA functional state (`A/A`,
`ap/AP`, `A/Elongation factor Tu`, etc.). This is the
"contact-transfer annotation" workflow — spec §3.3 has the long
version.

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
        │  Classify (§7, §8)                                │
        │  rRNA core + ribosomal-protein superkingdom vote  │
        │  → bacterial / eukaryotic / organellar / skip     │
        └────────────┬──────────────────────────────────────┘
                     ▼
        ┌───────────────────────────────────────────────────┐
        │  Pick reference (§6.4)                            │
        │  5J7L (bacterial/organellar) or 7ZW0 (eukaryotic) │
        └────────────┬──────────────────────────────────────┘
                     ▼
        ┌───────────────────────────────────────────────────┐
        │  BGSU correspondence (§5.2)                       │
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
        │  Assign mRNA / A-tRNA / P-tRNA / E-tRNA (§11)     │
        │  Infer tRNA states (§12)                          │
        │  Label A-site factor by CCA-end neighbour (§12.4) │
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
`src/ribosome_state_annotator/models.py` (or spec §9.1) for the full
field list. The role-based rRNA outputs (`ssu_main_rrna_chains` etc.)
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
| JSON | Default, or `--output foo.json` | Full `RibosomeAnnotation` list. Field layout in spec §15.1 and `models.py`. |
| JSONL | `--output foo.jsonl` | One annotation per line. For streaming consumers. |
| `ribosome_chain_annotation.csv` | Default companion (suppress with `--no-csv`) | One row per annotated assembly, 13 columns. Matches the legacy prototype byte-for-byte (spec §15.3). |
| `ribosome_assembly_annotation.csv` | Default companion | One row per `(property, value)` tuple: species, non-ribosomal proteins, bound ligands, unmapped RNA chains, plus v1 extensions (classification, superkingdom, warnings). |

Skipped and failed annotations appear in JSON but are omitted from CSV.

## Known limitations (spec §13.1 + §3.2)

Behaviours v1 accepts but you should know about:

- **`ribosome_state_annotator.classify.matches_ribosomal_protein_narrow`** is
  a plain substring check for `"ribosomal protein"`. It produces:
  - A false positive: `"ribosomal recycling factor"` contains the
    substring and gets flagged. Acceptable in v1 because RRF's
    geometry near the CCA end is similar enough that excluding it from
    the §12.4 factor search doesn't lose anything important.
  - A false negative: chains named just `"S1"` / `"L11"` without the
    `"ribosomal"` substring are missed. Rare in modern PDB depositions;
    common in older entries. The §8.4 broader rule (used only for the
    superkingdom vote) catches these via an anchored regex, but the
    `is_ribosomal_protein` flag used by §12.4 stays on the narrow rule.
- The allow-list / deny-list overrides documented in spec §13.1
  (`extra_ribosomal_descriptions`, `non_ribosomal_overrides`) are NOT
  exposed in the v1 API. Workaround: rebuild the `ChainRef`s with
  edited `is_ribosomal_protein` flags before passing them to the
  inference layer.

v1 out-of-scope (spec §3.2):

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
refresh, use `ribostate cache clear` or delete the cache root. See spec
§17 for the key format.

## Package layout

The source tree is split by responsibility — each module owns one stage
of the contact-transfer workflow above.

| File | Responsibility |
|------|----------------|
| `api.py` | Top-level orchestration. `annotate_pdb` / `annotate_assembly` / `annotate_many` live here; this is the file most library callers reach for. |
| `cli.py` | Typer CLI (`ribostate annotate`, `…-batch`, `cache`). Thin wrapper over `api.py`. |
| `models.py` | Pydantic v2 models: `ChainRef`, `LigandRef`, `AssemblyContext`, `CorrespondenceResult`, `RibosomeAnnotation`. All JSON / CSV output round-trips through these. |
| `constants.py` | Curated reference unit IDs (`BACTERIAL_REFERENCE_UNITS` for 5J7L, `YEAST_REFERENCE_UNITS` for 7ZW0). The curated functional-site anchors of §6.3. |
| `rcsb_client.py` | RCSB GraphQL fetcher + assembly parser. Reads `assemblies → polymer_entity_instances` and produces per-assembly `AssemblyContext` records. |
| `pdbe_client.py` | PDBe REST Rfam-mappings fetcher. Augments rRNA chains with Rfam accessions RCSB no longer ships (spec §5.3). |
| `bgsu_client.py` | BGSU correspondence HTTP client + tolerant JSON parser. Handles `format=json`, `depth=700`, and the live `mappings` / spec `alignment` response shapes (spec §5.2, §28.2). |
| `correspondence.py` | §5.2.2 PDB-prefix + assembly-chain filter, and the multi-assembly chain-substitution fallback (§28.2). |
| `coordinates.py` | mmCIF download + Gemmi parsing. Caches under `coords/`. |
| `gemmi_contacts.py` | Gemmi `NeighborSearch` wrapper restricted to the active assembly (§10). |
| `classify.py` | rRNA-core determination, ribosomal-protein detection, dominant-superkingdom vote, final classification rule (§8). |
| `infer.py` | Functional chain assignment + tRNA-state inference (§11, §12). Owns the polymer-filtered CCA-end selector for §12.4. |
| `cache.py` | Content-addressed on-disk cache with four namespaces: `rcsb/`, `bgsu/`, `pdbe/`, `coords/`. |
| `config.py` | Tunables: contact cutoff, completeness thresholds, network timeouts. |
| `exceptions.py` | Typed exception hierarchy (`ApiRequestError`, `CorrespondenceMappingError`, `CoordinateDownloadError`, …). |
| `output.py` | JSON / JSONL / CSV writers (spec §15). |

Spec sections referenced above all live in
[`ribosome_functional_annotation_package_spec_v3.md`](./ribosome_functional_annotation_package_spec_v3.md);
§28 in that document is the v3.1 addendum consolidating every place the
implementation diverges from the original spec body.

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
the PDB IDs the spec §25.1 acceptance tests use.

## Contributors

- Sri Devan Appasamy
- Claude (Anthropic) — pair-programmed the v1 build under direction

## License

MIT.
