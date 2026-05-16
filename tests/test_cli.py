"""Unit tests for the Typer CLI (spec §14.3, §17.3).

The CLI body delegates to :mod:`api`; tests mock the api functions at
the module boundary so we don't re-test the api orchestrator here. CSV
companion-file emission and stdout/file routing ARE tested end-to-end.
"""

from __future__ import annotations

import json
from collections.abc import Iterable
from pathlib import Path

import pytest
from typer.testing import CliRunner

from ribosome_state_annotator import cli
from ribosome_state_annotator.cache import Cache
from ribosome_state_annotator.models import (
    ChainRef,
    RibosomeAnnotation,
)

runner = CliRunner()


# ---------------------------------------------------------------------------
# Test fixtures
# ---------------------------------------------------------------------------


def _annotated() -> RibosomeAnnotation:
    ssu = ChainRef(
        pdb_id="5J7L",
        assembly_id="1",
        auth_asym_id="AA",
        rfam_accessions=["RF00177"],
        scientific_name="Escherichia coli",
    )
    lsu = ChainRef(
        pdb_id="5J7L",
        assembly_id="1",
        auth_asym_id="DA",
        rfam_accessions=["RF02541"],
    )
    return RibosomeAnnotation(
        pdb_id="5J7L",
        assembly_id="1",
        status="annotated",
        ribosome_classification="bacterial_ribosome",
        ssu_main_rrna_chains=[ssu],
        lsu_main_rrna_chains=[lsu],
        aminoacyl_trna_state="A/A",
    )


def _patch_annotate_pdb(monkeypatch: pytest.MonkeyPatch, results: list[RibosomeAnnotation]) -> None:
    """Replace cli.annotate_pdb with a stub returning ``results``."""

    def _stub(pdb_id: str, **kwargs: object) -> list[RibosomeAnnotation]:
        return list(results)

    monkeypatch.setattr(cli, "annotate_pdb", _stub)


def _patch_annotate_many(
    monkeypatch: pytest.MonkeyPatch, results: list[RibosomeAnnotation]
) -> None:
    def _stub(pdb_ids: Iterable[str], **kwargs: object) -> list[RibosomeAnnotation]:
        return list(results)

    monkeypatch.setattr(cli, "annotate_many", _stub)


# ---------------------------------------------------------------------------
# Top-level
# ---------------------------------------------------------------------------


def test_top_level_help() -> None:
    result = runner.invoke(cli.app, ["--help"])
    assert result.exit_code == 0
    assert "annotate" in result.stdout
    assert "annotate-batch" in result.stdout
    assert "cache" in result.stdout


def test_version_flag() -> None:
    result = runner.invoke(cli.app, ["--version"])
    assert result.exit_code == 0
    assert "ribostate" in result.stdout


# ---------------------------------------------------------------------------
# annotate
# ---------------------------------------------------------------------------


def test_annotate_writes_json_to_output(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    _patch_annotate_pdb(monkeypatch, [_annotated()])
    out = tmp_path / "out.json"
    result = runner.invoke(cli.app, ["annotate", "5J7L", "--output", str(out), "--no-cache"])
    assert result.exit_code == 0
    assert out.exists()
    parsed = json.loads(out.read_text())
    assert isinstance(parsed, list)
    assert parsed[0]["pdb_id"] == "5J7L"


def test_annotate_writes_companion_csvs_alongside_json(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    _patch_annotate_pdb(monkeypatch, [_annotated()])
    out = tmp_path / "out.json"
    result = runner.invoke(cli.app, ["annotate", "5J7L", "--output", str(out), "--no-cache"])
    assert result.exit_code == 0
    assert (tmp_path / cli.CHAIN_CSV_FILENAME).exists()
    assert (tmp_path / cli.ASSEMBLY_CSV_FILENAME).exists()


def test_annotate_no_csv_suppresses_companion_files(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    _patch_annotate_pdb(monkeypatch, [_annotated()])
    out = tmp_path / "out.json"
    result = runner.invoke(
        cli.app,
        ["annotate", "5J7L", "--output", str(out), "--no-cache", "--no-csv"],
    )
    assert result.exit_code == 0
    assert out.exists()
    assert not (tmp_path / cli.CHAIN_CSV_FILENAME).exists()
    assert not (tmp_path / cli.ASSEMBLY_CSV_FILENAME).exists()


def test_annotate_stdout_writes_json_to_stdout(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _patch_annotate_pdb(monkeypatch, [_annotated()])
    result = runner.invoke(cli.app, ["annotate", "5J7L", "--stdout", "--no-cache"])
    assert result.exit_code == 0
    # Strip rich/ANSI noise — the JSON is the only printable content.
    assert "5J7L" in result.stdout
    assert "bacterial_ribosome" in result.stdout


def test_annotate_without_output_or_stdout_auto_names_in_cwd(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Omitting both --output and --stdout is now valid: the CLI auto-names
    ``<PDB>.json`` (uppercase) in the current working directory and writes
    the two companion CSVs alongside."""
    _patch_annotate_pdb(monkeypatch, [_annotated()])
    monkeypatch.chdir(tmp_path)
    result = runner.invoke(cli.app, ["annotate", "5j7l", "--no-cache"])
    assert result.exit_code == 0, result.stderr
    assert (tmp_path / "5J7L.json").exists()
    assert (tmp_path / "ribosome_chain_annotation.csv").exists()
    assert (tmp_path / "ribosome_assembly_annotation.csv").exists()


def test_annotate_output_and_stdout_mutually_exclusive(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    _patch_annotate_pdb(monkeypatch, [_annotated()])
    out = tmp_path / "out.json"
    result = runner.invoke(
        cli.app,
        ["annotate", "5J7L", "--output", str(out), "--stdout", "--no-cache"],
    )
    assert result.exit_code == 2


def test_annotate_jsonl_extension_writes_jsonl(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    _patch_annotate_pdb(monkeypatch, [_annotated(), _annotated()])
    out = tmp_path / "out.jsonl"
    result = runner.invoke(cli.app, ["annotate", "5J7L", "--output", str(out), "--no-cache"])
    assert result.exit_code == 0
    text = out.read_text()
    # Two annotations → two lines.
    non_empty_lines = [line for line in text.splitlines() if line.strip()]
    assert len(non_empty_lines) == 2
    for line in non_empty_lines:
        json.loads(line)


def test_annotate_csv_extension_writes_chain_csv_only(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    _patch_annotate_pdb(monkeypatch, [_annotated()])
    out = tmp_path / "out.csv"
    result = runner.invoke(cli.app, ["annotate", "5J7L", "--output", str(out), "--no-cache"])
    assert result.exit_code == 0
    assert out.exists()
    text = out.read_text()
    assert "ssu_chain" in text
    # No companion files for --output=*.csv (companions only attach to JSON output).
    assert (
        not (tmp_path / cli.CHAIN_CSV_FILENAME).exists()
        or (tmp_path / cli.CHAIN_CSV_FILENAME).resolve() == out.resolve()
    )


def test_annotate_no_create_dirs_fails_when_parent_missing(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    _patch_annotate_pdb(monkeypatch, [_annotated()])
    out = tmp_path / "nope" / "out.json"
    result = runner.invoke(
        cli.app,
        [
            "annotate",
            "5J7L",
            "--output",
            str(out),
            "--no-cache",
            "--no-create-dirs",
        ],
    )
    assert result.exit_code == 2


# ---------------------------------------------------------------------------
# annotate-batch
# ---------------------------------------------------------------------------


def test_annotate_batch_reads_pdb_ids_and_writes_combined_output(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    _patch_annotate_many(monkeypatch, [_annotated(), _annotated()])
    pdb_ids_file = tmp_path / "ids.txt"
    pdb_ids_file.write_text("5J7L\n5TBW\n# a comment line\n\n")
    out = tmp_path / "batch.json"
    result = runner.invoke(
        cli.app,
        ["annotate-batch", str(pdb_ids_file), "--output", str(out), "--no-cache"],
    )
    assert result.exit_code == 0
    parsed = json.loads(out.read_text())
    assert isinstance(parsed, list)
    assert len(parsed) == 2


def test_annotate_batch_strips_inline_comments(tmp_path: Path) -> None:
    ids_file = tmp_path / "ids.txt"
    ids_file.write_text("5J7L  # E. coli\n# full comment\n5TBW\n")
    ids = cli._read_pdb_ids(ids_file)
    assert ids == ["5J7L", "5TBW"]


def test_annotate_batch_empty_file_exits_with_error(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    _patch_annotate_many(monkeypatch, [])
    pdb_ids_file = tmp_path / "empty.txt"
    pdb_ids_file.write_text("# only comments\n\n")
    out = tmp_path / "out.json"
    result = runner.invoke(
        cli.app,
        ["annotate-batch", str(pdb_ids_file), "--output", str(out)],
    )
    assert result.exit_code == 2


# ---------------------------------------------------------------------------
# cache info / cache clear
# ---------------------------------------------------------------------------


def test_cache_info_on_missing_cache(tmp_path: Path) -> None:
    result = runner.invoke(cli.app, ["cache", "info", "--cache-dir", str(tmp_path / "absent")])
    assert result.exit_code == 0
    assert "missing" in result.stdout


def test_cache_info_on_populated_cache(tmp_path: Path) -> None:
    cache = Cache(tmp_path / "filled")
    cache.put_rcsb_payload("5J7L", {"x": 1})
    result = runner.invoke(cli.app, ["cache", "info", "--cache-dir", str(cache.root)])
    assert result.exit_code == 0
    assert "rcsb" in result.stdout
    assert "1" in result.stdout  # one entry


def test_cache_clear_with_yes_removes_root(tmp_path: Path) -> None:
    cache = Cache(tmp_path / "to-clear")
    cache.put_rcsb_payload("5J7L", {"x": 1})
    assert cache.root.exists()
    result = runner.invoke(cli.app, ["cache", "clear", "--cache-dir", str(cache.root), "--yes"])
    assert result.exit_code == 0
    assert not cache.root.exists()


def test_cache_clear_on_missing_cache_is_noop(tmp_path: Path) -> None:
    target = tmp_path / "absent"
    result = runner.invoke(cli.app, ["cache", "clear", "--cache-dir", str(target), "--yes"])
    assert result.exit_code == 0
    # Click 8.3 separates stdout and stderr; status messages go to stderr.
    assert "nothing to clear" in result.stderr


def test_cache_clear_prompts_without_yes(tmp_path: Path) -> None:
    """Without --yes, the user is prompted; declining cancels."""
    cache = Cache(tmp_path / "filled")
    cache.put_rcsb_payload("5J7L", {"x": 1})
    result = runner.invoke(cli.app, ["cache", "clear", "--cache-dir", str(cache.root)], input="n\n")
    # Non-zero exit because the user cancelled.
    assert result.exit_code == 1
    # Cache still there.
    assert cache.root.exists()
