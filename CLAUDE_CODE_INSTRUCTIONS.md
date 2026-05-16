# Initial instructions for Claude Code: ribosome-state-annotator

## What you're building

A production-quality Python package called `ribosome-state-annotator` (import name `ribosome_state_annotator`) that annotates functional states of complete bacterial and eukaryotic ribosome assemblies from PDB entries. Usable both as a CLI tool (`ribostate`) and as a Python library.

This is a rewrite of an existing prototype, **not a refactor**. The prototype relied on a local SQLAlchemy database (`UnitPairsInteractions2024` etc.) that is being thrown away. The new package retrieves data from public APIs (RCSB GraphQL, BGSU correspondence) and detects atomic contacts locally with Gemmi.

## Files in this directory

Before you do anything else, read these in order:

1. **`ribosome_functional_annotation_package_spec_v3.md`** — The full specification. This is your source of truth. Every design decision, threshold, output field, API endpoint, and edge case is documented there. Read all of it. It is ~2,200 lines but most sections are short and self-contained.

2. **`process_annotation.py`** and **`ribosome_annotation_new.py`** — The original prototype. **Do not copy code from these files**, but read them to understand the biological inference logic (especially `infer_tRNA_state`). Section 27 of the spec lists exactly what to delete from the prototype and what scientific logic to preserve.

3. **`ribosome_chain_annotation.csv`** (1,104 rows) and **`ribosome_assembly_annotation.csv`** (6,043 rows) — Real output files from the prototype. These define the de-facto CSV contract that downstream consumers depend on. The new package must produce CSVs with the same column layout (spec §15.3). Move these into `tests/fixtures/legacy_csv/` once you've created the package skeleton — they are your scale regression corpus (spec §22.3).

## Before you start coding

Confirm to me that you have read the spec by answering these five questions in your first reply. Don't write any code until I've confirmed your answers. If any answer is "I'd need to check" or you're guessing, stop and re-read the relevant section.

1. What is the central scientific abstraction of this package, and why is it different from chain-name-based or sequence-alignment-based annotation? (Spec §3.3.)

2. The output model has both `ssu_main_rrna_chains: list[ChainRef]` (canonical) and `ssu_chain: ChainRef | None` (derived). Why is the list-based form canonical, and what triggers the derived alias to be `None`? (Spec §6.2, §9.1.)

3. An assembly has a bacterial-like rRNA Rfam (`RF00177` + `RF02541`) but the dominant ribosomal protein superkingdom vote returns `Eukaryota` with 78 votes vs `Bacteria` with 2. What classification does the package emit, and what string goes in the `rule` field of `classification_evidence`? (Spec §8.3, §8.5.)

4. A full-length A-site tRNA does not contact either the LSU A-site or the LSU P-site reference nucleotides. How does the package decide what LSU label to put in `aminoacyl_trna_state`? Specifically: what is the rule for which protein chain is chosen, and is the resulting label a canonical short token or the raw `pdbx_description`? (Spec §12.1, §12.4.)

5. A user runs `ribostate annotate 5FDV` and 5FDV has two complete biological assemblies. How many JSON result objects are emitted? What's in the `assembly_id` field of each? (Spec §4, §14.1, §25.1 row 8.)

## How to work

### Engineering priorities (build in this order — spec §23.1)

1. Package skeleton (`pyproject.toml`, `src/ribosome_state_annotator/`, `tests/`, CLI scaffolding via Typer). Get `ribostate --help` working with a stub.
2. Models and constants (`models.py`, `constants.py`, `exceptions.py`). Make `ChainRef`, `LigandRef`, `RibosomeAnnotation`, `AssemblyContext`, `CorrespondenceResult` pass `mypy --strict` before moving on.
3. RCSB client (`rcsb_client.py`). Implement the GraphQL query from spec §5.1.1 verbatim. Write unit tests that parse the example response from §5.1.2.
4. Complete-ribosome filter + classification (`classify.py`). This is the hardest pure-logic component. Cover every branch of §8.5 with unit tests (each rule has a test case; ties have test cases; insufficient-evidence has a test case).
5. BGSU correspondence client + parser (`bgsu_client.py`, `correspondence.py`). Mock the HTTP layer with `respx`. Test the two-step filtering rule from §5.2.2.
6. Coordinate cache + Gemmi loader (`coordinates.py`, `cache.py`). Test against a small mmCIF fixture.
7. Unit ID parser (`correspondence.py`). Handle modified nucleotides (`4OC`, `OMG`).
8. Gemmi neighbour-search wrapper (`gemmi_contacts.py`).
9. Chain assignment + tRNA-state inference (`infer.py`). This carries the most prototype logic; cross-reference `infer_tRNA_state` in `process_annotation.py` for the SSU/LSU state rules but **rewrite, do not copy**.
10. JSON / JSONL / CSV output (`output.py`). Match the prototype CSV column layout exactly.
11. Tests, docs, and the legacy-CSV regression suite (§22.3).

Do not jump ahead. Each step has tests that gate the next step.

### Working style

- **One PR-sized chunk at a time.** After completing each numbered step above, stop and tell me what you did, what tests pass, and any open questions. Don't proceed to the next step without explicit confirmation.
- **Tests before integration.** Every module gets unit tests using `pytest`. Mock all network calls — no test should hit the real RCSB or BGSU API during `pytest`. Integration tests that do hit the network live in `tests/integration/` and are marked `@pytest.mark.network` so they can be excluded by default.
- **No `print()` outside the CLI rendering layer.** Use `logging.getLogger(__name__)` per module (spec §18).
- **No `sys.exit()` outside the CLI entrypoint** (spec §16).
- **No global mutable state.** The prototype has module-level flags like `HAS_MULTIPLE_ASSEMBLIES` — these are explicitly forbidden (spec §27 item 8).
- **Type hints everywhere.** Run `mypy --strict src/` and `ruff check src/` clean before declaring a step done.
- **Pydantic v2** for all data models. Use `model_config = ConfigDict(frozen=True)` where appropriate.

### When the spec is silent or ambiguous

Ask me. Do not invent. The spec has been iterated specifically to remove ambiguity in the high-risk areas (RCSB API contract, classification rule, tRNA state labelling, output model). If you find yourself reaching for "reasonable default" reasoning on something biologically significant, stop and ask.

Three specific things that look like spec gaps but aren't:

- **The protein-factor LSU label uses the raw `pdbx_description`, not a canonical token.** An earlier draft proposed canonical tokens (`T`, `R`, `G`, etc.); §12.4 explicitly rejects them. Do not reintroduce a token table.
- **The CSV format reproduces the prototype's column layout exactly.** Including `ssu_chain`, `lsu_large_chain`, `lsu_medium_chain`, `lsu_small_chain` as columns. These are derived from the role-based canonical fields; the derivation is in §15.3 and §11.1. The role-based form is canonical *in JSON and in Python*; the legacy form is canonical *in CSV*.
- **Chain identifiers in CSV output are IFE strings** (`PDB|1|chain`, e.g. `5J7L|1|AA`), not bare auth_asym_id. `ChainRef` has a `.ife` computed field for this.

### What "done" looks like for v1

- All eight acceptance-test cases in spec §25.1 pass.
- The legacy-CSV regression test (spec §22.3) passes on the sample of ~20 entries.
- `mypy --strict src/` clean.
- `ruff check src/` and `ruff format src/` clean.
- `pytest` (without `network` marker) achieves ≥85% coverage.
- README documents installation, CLI usage, library usage, and lists the §13.1 known limitations.
- `ribostate annotate 5J7L -o /tmp/out.json` works end-to-end against the live APIs and produces a sensible result.

## First task

Once you've answered the five questions above and I've confirmed, start with **step 1 only**: create the package skeleton. `pyproject.toml` with the dependencies from spec §21, the `src/ribosome_state_annotator/` layout from §20, a stubbed Typer CLI that prints "not yet implemented" for `annotate` and `annotate-batch`, and a working `pytest` invocation that finds zero tests but exits 0. That's it. Then stop and tell me what you did.
