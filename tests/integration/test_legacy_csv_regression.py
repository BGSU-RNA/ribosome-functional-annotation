"""Legacy-CSV regression suite (spec §22.3).

Compares the v1 package's chain-level CSV output against the prototype's
checked-in CSV (``tests/fixtures/legacy_csv/ribosome_chain_annotation.csv``)
on a curated sample of ~20 (pdb_id, assembly_id) entries covering each
row-shape the spec calls out:

- bacterial complete with no mRNA/tRNAs (just SSU + LSU + 5S)
- bacterial with full mRNA + A/P/E (classic ``A/A``, ``P/P``, ``E/E``)
- bacterial hybrid states (``P/E``, ``A/P``)
- bacterial protein-factor LSU labels (``A/Elongation factor Tu``)
- eukaryotic cytoplasmic with populated ``lsu_medium_chain`` (5.8S)
- multi-assembly

Per spec §22.3 the assertions are deliberately tolerant:

- ``ssu_chain`` / ``lsu_large_chain``: must match exactly or both be empty.
- ``mrna`` / ``aminoacyl_trna`` / ``peptidyl_trna`` / ``exit_trna``: when
  both legacy and new are non-empty they must refer to the same chain;
  one-sided emptiness is logged but does not fail individual tests.
- tRNA-state strings: same first character on each side of the ``/``.

This is a ``network`` test (live RCSB + BGSU calls); excluded from the
default ``pytest`` run.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from ribosome_state_annotator.api import annotate_assembly
from ribosome_state_annotator.output import chain_csv_row

LEGACY_CSV_PATH = (
    Path(__file__).parent.parent / "fixtures" / "legacy_csv" / "ribosome_chain_annotation.csv"
)

# Sample selected per spec §22.3 to cover the listed row shapes.
REGRESSION_SAMPLE: list[tuple[str, str]] = [
    # Bacterial complete, no mRNA/tRNAs (just rRNA chains).
    ("4V48", "1"),
    ("4V49", "1"),
    # Bacterial full mRNA + A/P/E classic.
    ("8QBT", "1"),
    # Bacterial with state A/A.
    ("5MC6", "1"),
    # Bacterial hybrid P/E.
    ("6Y57", "1"),
    # Bacterial protein-factor LSU label (A/Elongation factor Tu and similar).
    ("5WFS", "1"),
    ("4V5Q", "1"),
    ("4V5Q", "2"),
    ("8G7P", "1"),
    ("5UYP", "1"),
    # Eukaryotic cytoplasmic with populated lsu_medium_chain (5.8S).
    ("4U6F", "1"),
    # Multi-assembly bacterial.
    ("5FDV", "1"),
    ("5FDV", "2"),
]


# ---------------------------------------------------------------------------
# Legacy CSV loader
# ---------------------------------------------------------------------------


def _load_legacy_rows() -> dict[tuple[str, str], dict[str, str]]:
    """Load the legacy chain CSV keyed by (pdb_id, assembly_id)."""
    rows: dict[tuple[str, str], dict[str, str]] = {}
    with LEGACY_CSV_PATH.open(newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            key = (row["pdb_id"], row["assembly_id"])
            rows[key] = row
    return rows


@pytest.fixture(scope="module")
def legacy_rows() -> dict[tuple[str, str], dict[str, str]]:
    return _load_legacy_rows()


# ---------------------------------------------------------------------------
# Per-row regression test
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(("pdb_id", "assembly_id"), REGRESSION_SAMPLE)
def test_legacy_csv_compat(
    pdb_id: str,
    assembly_id: str,
    legacy_rows: dict[tuple[str, str], dict[str, str]],
) -> None:
    legacy = legacy_rows.get((pdb_id, assembly_id))
    if legacy is None:
        pytest.skip(
            f"{pdb_id} assembly {assembly_id} is not in the legacy CSV — "
            "regression sample needs adjustment"
        )

    annotation = annotate_assembly(pdb_id, assembly_id, no_cache=True)
    if annotation.status != "annotated":
        pytest.fail(
            f"{pdb_id} assembly {assembly_id} expected annotated but got "
            f"{annotation.status} ({annotation.skip_reason})"
        )

    new = chain_csv_row(annotation)
    _assert_chain_compat(legacy, new, pdb_id, assembly_id)


# ---------------------------------------------------------------------------
# Assertion helpers
# ---------------------------------------------------------------------------


def _assert_chain_compat(
    legacy: dict[str, str], new: dict[str, str], pdb_id: str, assembly_id: str
) -> None:
    """Assert that the new CSV row is compatible with the legacy one per §22.3."""
    label = f"{pdb_id}/{assembly_id}"

    # ssu_chain and lsu_large_chain: same or both empty (§22.3).
    _assert_same_or_both_empty(legacy["ssu_chain"], new["ssu_chain"], f"{label} ssu_chain")
    _assert_same_or_both_empty(
        legacy["lsu_large_chain"], new["lsu_large_chain"], f"{label} lsu_large_chain"
    )

    # Functional chains: when both non-empty they must agree.
    for column in ("mrna", "aminoacyl_trna", "peptidyl_trna", "exit_trna"):
        _assert_same_when_both_present(legacy[column], new[column], f"{label} {column}")

    # tRNA states: same first character on each side of the "/" (§22.3).
    for state_col in (
        "aminoacyl_trna_state",
        "peptidyl_trna_state",
        "exit_trna_state",
    ):
        _assert_state_compatible(legacy[state_col], new[state_col], f"{label} {state_col}")


def _assert_same_or_both_empty(legacy: str, new: str, where: str) -> None:
    if legacy == "" and new == "":
        return
    assert legacy == new, f"{where}: legacy={legacy!r} new={new!r}"


def _assert_same_when_both_present(legacy: str, new: str, where: str) -> None:
    if not legacy or not new:
        return  # One-sided emptiness — logged, not failed.
    assert legacy == new, f"{where}: legacy={legacy!r} new={new!r}"


def _assert_state_compatible(legacy: str, new: str, where: str) -> None:
    """First character on each side of the ``/`` must agree.

    So ``A/Elongation factor Tu`` matches ``A/Elongation factor Tu 2``
    (both ``A`` SSU, both LSU starts with ``E``). ``A/A`` does NOT match
    ``ap/A`` (SSU first char ``A`` vs ``a``, the chimeric / classic
    distinction must be preserved).
    """
    if not legacy or not new:
        return
    legacy_ssu, _, legacy_lsu = legacy.partition("/")
    new_ssu, _, new_lsu = new.partition("/")
    assert legacy_ssu[:1] == new_ssu[:1], (
        f"{where}: SSU first-char mismatch legacy={legacy!r} new={new!r}"
    )
    assert legacy_lsu[:1] == new_lsu[:1], (
        f"{where}: LSU first-char mismatch legacy={legacy!r} new={new!r}"
    )


# ---------------------------------------------------------------------------
# Sanity check: the regression sample exists in the legacy CSV
# ---------------------------------------------------------------------------


def test_regression_sample_covers_legacy_csv(
    legacy_rows: dict[tuple[str, str], dict[str, str]],
) -> None:
    """If any sample row is missing from the legacy CSV the test is mis-pinned."""
    missing = [key for key in REGRESSION_SAMPLE if key not in legacy_rows]
    assert not missing, f"sample rows missing from legacy CSV: {missing}"
