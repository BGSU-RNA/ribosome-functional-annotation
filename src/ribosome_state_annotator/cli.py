"""Typer CLI for ``ribostate`` (spec §14.3, §17.3).

The CLI is the only module allowed to call :func:`sys.exit` (via
:exc:`typer.Exit`) per spec §16. All work is delegated to
:mod:`.api` and rendered through :mod:`.output`.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console
from rich.table import Table

from ribosome_state_annotator import __version__
from ribosome_state_annotator.api import annotate_many, annotate_pdb
from ribosome_state_annotator.cache import Cache
from ribosome_state_annotator.coordinates import CoordinateSource
from ribosome_state_annotator.models import RibosomeAnnotation
from ribosome_state_annotator.output import (
    render_json,
    render_jsonl,
    write_assembly_csv,
    write_chain_csv,
    write_json,
    write_jsonl,
)

logger = logging.getLogger(__name__)

app = typer.Typer(
    name="ribostate",
    help="Annotate functional states of complete ribosome assemblies from PDB entries.",
    no_args_is_help=True,
    add_completion=False,
    context_settings={"help_option_names": ["-h", "--help"]},
)

cache_app = typer.Typer(
    name="cache",
    help="Inspect or clear the on-disk cache (spec §17.3).",
    no_args_is_help=True,
    add_completion=False,
)
app.add_typer(cache_app, name="cache")

stdout_console = Console()


def _err(message: str) -> None:
    """Emit a Rich-markup message to stderr.

    A fresh ``Console`` is constructed per call so that the binding picks
    up whatever ``sys.stderr`` currently is — important for Click's
    ``CliRunner``, which redirects stderr at test time. Rich auto-strips
    ANSI sequences when stderr is not a TTY, so markup like
    ``[red]...[/red]`` renders as colour on terminals and as plain text
    under capture.
    """
    Console(stderr=True).print(message)


CHAIN_CSV_FILENAME = "ribosome_chain_annotation.csv"
ASSEMBLY_CSV_FILENAME = "ribosome_assembly_annotation.csv"


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------


def _configure_logging(verbose: bool, debug: bool) -> None:
    """Configure root logging level for the CLI. Library code does not call this."""
    if debug:
        level = logging.DEBUG
    elif verbose:
        level = logging.INFO
    else:
        level = logging.WARNING
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )


def _ensure_parent(output_path: Path, no_create_dirs: bool) -> None:
    """Make sure ``output_path.parent`` exists, unless ``--no-create-dirs``."""
    parent = output_path.parent
    if not parent.exists():
        if no_create_dirs:
            _err(f"[red]parent directory does not exist: {parent}[/red]")
            raise typer.Exit(code=2)
        parent.mkdir(parents=True, exist_ok=True)


def _validate_output_flags(output: Path | None, stdout: bool) -> None:
    """Enforce --output / --stdout mutex BEFORE any expensive work runs.

    Called at the top of each command body so flag misuse fails fast
    without first hitting the live RCSB / BGSU / PDBe / RCSB-file APIs.
    """
    if stdout and output is not None:
        _err("[red]--stdout and --output are mutually exclusive[/red]")
        raise typer.Exit(code=2)
    if not stdout and output is None:
        _err("[red]one of --output PATH or --stdout is required[/red]")
        raise typer.Exit(code=2)


def _emit_annotations(
    annotations: list[RibosomeAnnotation],
    *,
    output: Path | None,
    stdout: bool,
    no_create_dirs: bool,
    no_csv: bool,
) -> None:
    """Route annotations to stdout or to ``output`` per CLI flags.

    Format inference from ``output`` extension:

    - ``.jsonl`` → JSON-Lines
    - ``.csv`` → chain-level CSV (no companion assembly CSV)
    - anything else → JSON

    When writing JSON to a file and ``--no-csv`` is *not* set, two
    companion CSVs are written alongside (``ribosome_chain_annotation.csv``
    and ``ribosome_assembly_annotation.csv``) per spec §15.3.

    Callers are expected to have already run :func:`_validate_output_flags`;
    re-checking here for defence in depth.
    """
    _validate_output_flags(output, stdout)

    if stdout:
        # Stream JSON to stdout. JSONL alternative is available via
        # explicit shell-side redirection of a .jsonl --output instead.
        stdout_console.print(render_json(annotations), highlight=False)
        return

    assert output is not None
    _ensure_parent(output, no_create_dirs)
    suffix = output.suffix.lower()
    if suffix == ".jsonl":
        write_jsonl(annotations, output)
    elif suffix == ".csv":
        write_chain_csv(annotations, output)
    else:
        write_json(annotations, output)
        if not no_csv:
            chain_csv_path = output.parent / CHAIN_CSV_FILENAME
            assembly_csv_path = output.parent / ASSEMBLY_CSV_FILENAME
            write_chain_csv(annotations, chain_csv_path)
            write_assembly_csv(annotations, assembly_csv_path)
            _err(
                f"[green]wrote {chain_csv_path.name} and {assembly_csv_path.name}"
                f" alongside {output.name}[/green]"
            )

    n_annotated = sum(1 for a in annotations if a.status == "annotated")
    n_skipped = sum(1 for a in annotations if a.status == "skipped")
    n_failed = sum(1 for a in annotations if a.status == "failed")
    _err(
        f"[green]wrote {len(annotations)} annotation(s) to {output} "
        f"(annotated={n_annotated}, skipped={n_skipped}, failed={n_failed})[/green]"
    )


def _read_pdb_ids(pdb_ids_file: Path) -> list[str]:
    """Read PDB IDs one per line, stripping blanks and ``#`` comments."""
    ids: list[str] = []
    for raw_line in pdb_ids_file.read_text(encoding="utf-8").splitlines():
        stripped = raw_line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        # Allow inline trailing comments: "5J7L  # E. coli 70S"
        if "#" in stripped:
            stripped = stripped.split("#", 1)[0].strip()
        if stripped:
            ids.append(stripped)
    return ids


def _resolve_coordinate_source(
    local_file: Path | None,
) -> tuple[CoordinateSource, Path | None]:
    if local_file is not None:
        return ("local", local_file)
    return ("auto", None)


# ---------------------------------------------------------------------------
# Top-level callback
# ---------------------------------------------------------------------------


def _version_callback(value: bool) -> None:
    if value:
        typer.echo(f"ribostate {__version__}")
        raise typer.Exit(code=0)


@app.callback()
def main(
    version: Annotated[
        bool,
        typer.Option(
            "--version",
            help="Show the ribostate version and exit.",
            callback=_version_callback,
            is_eager=True,
        ),
    ] = False,
) -> None:
    """Top-level callback. Sub-commands carry the real flag surface."""
    return None


# ---------------------------------------------------------------------------
# annotate
# ---------------------------------------------------------------------------


@app.command("annotate")
def annotate(
    pdb_id: Annotated[str, typer.Argument(help="Four-character PDB accession, e.g. 5J7L.")],
    output: Annotated[
        Path | None,
        typer.Option("--output", "-o", help="Path to write the annotation result."),
    ] = None,
    stdout: Annotated[
        bool,
        typer.Option("--stdout", help="Write JSON to stdout instead of a file."),
    ] = False,
    assembly_id: Annotated[
        str | None,
        typer.Option("--assembly-id", help="Restrict to a single biological assembly ID."),
    ] = None,
    cutoff: Annotated[
        float,
        typer.Option("--cutoff", help="Gemmi neighbour-search cutoff in angstroms."),
    ] = 5.0,
    cache_dir: Annotated[
        Path | None,
        typer.Option(
            "--cache-dir",
            help="Override the default ~/.cache/ribosome-state-annotator.",
        ),
    ] = None,
    no_cache: Annotated[
        bool,
        typer.Option("--no-cache", help="Disable caching for this invocation."),
    ] = False,
    strict: Annotated[
        bool,
        typer.Option(
            "--strict",
            help="Skip (don't just warn) assemblies with low ribosomal protein counts.",
        ),
    ] = False,
    no_csv: Annotated[
        bool,
        typer.Option(
            "--no-csv",
            help="Suppress the companion chain/assembly CSV files (JSON only).",
        ),
    ] = False,
    no_create_dirs: Annotated[
        bool,
        typer.Option(
            "--no-create-dirs",
            help="Fail rather than auto-creating missing parent directories for --output.",
        ),
    ] = False,
    input_file: Annotated[
        Path | None,
        typer.Option(
            "--input-file",
            help="Local mmCIF (.cif or .cif.gz) instead of downloading from RCSB.",
        ),
    ] = None,
    verbose: Annotated[bool, typer.Option("--verbose", help="INFO-level logging.")] = False,
    debug: Annotated[bool, typer.Option("--debug", help="DEBUG-level logging.")] = False,
) -> None:
    """Annotate one PDB entry (every biological assembly, or one --assembly-id)."""
    _configure_logging(verbose=verbose, debug=debug)
    _validate_output_flags(output, stdout)
    coordinate_source, local_path = _resolve_coordinate_source(input_file)
    annotations = annotate_pdb(
        pdb_id,
        assembly_id=assembly_id,
        contact_cutoff_angstrom=cutoff,
        cache_dir=cache_dir,
        no_cache=no_cache,
        strict_complete_check=strict,
        coordinate_source=coordinate_source,
        local_coordinate_path=local_path,
    )
    _emit_annotations(
        annotations,
        output=output,
        stdout=stdout,
        no_create_dirs=no_create_dirs,
        no_csv=no_csv,
    )


# ---------------------------------------------------------------------------
# annotate-batch
# ---------------------------------------------------------------------------


@app.command("annotate-batch")
def annotate_batch(
    pdb_ids_file: Annotated[
        Path,
        typer.Argument(
            help="Plain-text file with one PDB ID per line. Blank/# comment lines ignored.",
            exists=True,
            dir_okay=False,
            readable=True,
        ),
    ],
    output: Annotated[
        Path | None,
        typer.Option("--output", "-o", help="Path to write the combined annotation file."),
    ] = None,
    stdout: Annotated[
        bool,
        typer.Option("--stdout", help="Write JSON to stdout instead of a file."),
    ] = False,
    continue_on_error: Annotated[
        bool,
        typer.Option(
            "--continue-on-error",
            help="Don't abort the batch on per-entry errors; emit a failed record and continue.",
        ),
    ] = False,
    cutoff: Annotated[
        float,
        typer.Option("--cutoff", help="Gemmi neighbour-search cutoff in angstroms."),
    ] = 5.0,
    cache_dir: Annotated[
        Path | None,
        typer.Option("--cache-dir", help="Override the default cache directory."),
    ] = None,
    no_cache: Annotated[
        bool,
        typer.Option("--no-cache", help="Disable caching for this invocation."),
    ] = False,
    strict: Annotated[
        bool,
        typer.Option("--strict", help="Strict completeness filtering."),
    ] = False,
    no_csv: Annotated[
        bool,
        typer.Option(
            "--no-csv",
            help="Suppress the companion chain/assembly CSV files (JSON only).",
        ),
    ] = False,
    no_create_dirs: Annotated[
        bool,
        typer.Option(
            "--no-create-dirs",
            help="Fail rather than auto-creating missing parent directories for --output.",
        ),
    ] = False,
    verbose: Annotated[bool, typer.Option("--verbose", help="INFO-level logging.")] = False,
    debug: Annotated[bool, typer.Option("--debug", help="DEBUG-level logging.")] = False,
) -> None:
    """Annotate every PDB ID listed in PDB_IDS_FILE."""
    _configure_logging(verbose=verbose, debug=debug)
    _validate_output_flags(output, stdout)
    pdb_ids = _read_pdb_ids(pdb_ids_file)
    if not pdb_ids:
        _err(f"[yellow]no PDB IDs found in {pdb_ids_file}[/yellow]")
        raise typer.Exit(code=2)
    annotations = annotate_many(
        pdb_ids,
        continue_on_error=continue_on_error,
        contact_cutoff_angstrom=cutoff,
        cache_dir=cache_dir,
        no_cache=no_cache,
        strict_complete_check=strict,
    )
    _emit_annotations(
        annotations,
        output=output,
        stdout=stdout,
        no_create_dirs=no_create_dirs,
        no_csv=no_csv,
    )


# ---------------------------------------------------------------------------
# cache info / cache clear
# ---------------------------------------------------------------------------


@cache_app.command("info")
def cache_info(
    cache_dir: Annotated[
        Path | None,
        typer.Option("--cache-dir", help="Cache directory to inspect."),
    ] = None,
) -> None:
    """Print on-disk cache size and entry count (spec §17.3)."""
    cache = Cache(cache_dir) if cache_dir is not None else Cache()
    info = cache.info()
    table = Table(title=f"ribostate cache @ {info.root}")
    table.add_column("Namespace")
    table.add_column("Entries", justify="right")
    if not info.exists:
        table.add_row("status", "[yellow]missing[/yellow]")
    else:
        table.add_row("rcsb", str(info.rcsb_entries))
        table.add_row("bgsu", str(info.bgsu_entries))
        table.add_row("pdbe", str(info.pdbe_entries))
        table.add_row("coords", str(info.coords_entries))
        table.add_row("[bold]total entries", f"[bold]{info.total_entries}[/bold]")
        table.add_row("total bytes", f"{info.total_bytes:,}")
    stdout_console.print(table)


@cache_app.command("clear")
def cache_clear(
    cache_dir: Annotated[
        Path | None,
        typer.Option("--cache-dir", help="Cache directory to clear."),
    ] = None,
    yes: Annotated[
        bool,
        typer.Option("--yes", "-y", help="Skip the interactive confirmation prompt."),
    ] = False,
) -> None:
    """Delete the cache directory (spec §17.3)."""
    cache = Cache(cache_dir) if cache_dir is not None else Cache()
    if not cache.root.exists():
        _err(f"[yellow]nothing to clear at {cache.root}[/yellow]")
        return
    if not yes:
        confirmed = typer.confirm(f"Clear cache at {cache.root}?")
        if not confirmed:
            _err("[yellow]cache clear cancelled[/yellow]")
            raise typer.Exit(code=1)
    cache.clear()
    _err(f"[green]cleared cache at {cache.root}[/green]")


# Silence "unused JSONL helper" warnings — render_jsonl is currently
# unused by the CLI but is part of the public output API surface and
# is exercised by tests/test_output.py.
_ = render_jsonl


if __name__ == "__main__":  # pragma: no cover
    app()
