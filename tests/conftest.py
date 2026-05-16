"""Shared pytest fixtures.

The minimal Gemmi-built fixtures are reused by tests for
:mod:`coordinates` (step 6) and :mod:`gemmi_contacts` (step 8). They
deliberately stay tiny — a few residues per chain with one atom each —
so the structure builder is fast and the assertions stay focused on the
loader/neighbour-search code rather than on biological realism.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import gemmi
import pytest


def _add_residue(chain: gemmi.Chain, name: str, seqid: int, x: float, y: float, z: float) -> None:
    """Append a single-atom residue to ``chain`` at the given coordinates."""
    residue = gemmi.Residue()
    residue.name = name
    residue.seqid = gemmi.SeqId(seqid, " ")
    # The §12.4 CCA-end selector in :mod:`infer` filters to
    # ``EntityType.Polymer`` so chain-included waters / ions don't get
    # mistaken for tRNA acceptor-end residues. Synthetic test residues
    # default to ``Unknown``; mark them as polymer so they survive the
    # filter.
    residue.entity_type = gemmi.EntityType.Polymer
    atom = gemmi.Atom()
    # C1' for nucleotides, CA for amino acids; both are valid Gemmi elements.
    atom.name = "C1'" if len(name) == 1 else "CA"
    atom.element = gemmi.Element("C")
    atom.pos = gemmi.Position(x, y, z)
    residue.add_atom(atom)
    chain.add_residue(residue)


def build_minimal_structure() -> gemmi.Structure:
    """Build a small in-memory :class:`gemmi.Structure` for tests.

    Layout:

    - Model ``"1"`` (the only model — biological-assembly mmCIFs since the
      May 2022 wwPDB change have a single model with operator-expanded
      chains, per spec §10.4).
    - RNA chain ``"AA"`` with 5 nucleotides (G/C/A/U/G) at residue numbers
      1..5, single C1' atom each, x-coordinates 1.0..5.0.
    - RNA chain ``"DA"`` with 5 nucleotides at residue numbers 100..104,
      x-coordinates 10.0..14.0.
    - Protein chain ``"P1"`` with 3 amino acids (ALA/GLY/VAL) at residue
      numbers 1..3, single CA atom each, x-coordinates 20.0..22.0.
    """
    structure = gemmi.Structure()
    structure.name = "TEST"

    model = gemmi.Model("1")

    chain_aa = gemmi.Chain("AA")
    for i, name in enumerate(["G", "C", "A", "U", "G"], start=1):
        _add_residue(chain_aa, name, seqid=i, x=float(i), y=0.0, z=0.0)
    model.add_chain(chain_aa)

    chain_da = gemmi.Chain("DA")
    for i, name in enumerate(["G", "C", "U", "U", "G"], start=1):
        _add_residue(chain_da, name, seqid=100 + i - 1, x=float(10 + i - 1), y=0.0, z=0.0)
    model.add_chain(chain_da)

    chain_p1 = gemmi.Chain("P1")
    for i, name in enumerate(["ALA", "GLY", "VAL"], start=1):
        _add_residue(chain_p1, name, seqid=i, x=float(20 + i - 1), y=0.0, z=0.0)
    model.add_chain(chain_p1)

    structure.add_model(model)
    # NeighborSearch (step 8) requires a unit cell. Use a large P 1 cell
    # so the spatial grid is sized correctly for the structure's extent
    # without periodic-boundary surprises.
    structure.cell = gemmi.UnitCell(100.0, 100.0, 100.0, 90.0, 90.0, 90.0)
    return structure


def build_contact_structure() -> gemmi.Structure:
    """Build a structure with deliberate inter-chain contact pairs (step 8 tests).

    Layout (all atoms on the x-axis at y=z=0):

    - Chain ``AA`` (RNA, 3 nts): residue 1 G at x=0, residue 2 C at x=1,
      residue 3 A at x=2.
    - Chain ``BB`` (RNA, 2 nts): residue 1 U at x=2.5, residue 2 G at x=10.0.
    - Chain ``CC`` (protein, 1 aa): residue 1 ALA at x=3.5.
    - Chain ``FAR`` (protein, 1 aa): residue 1 GLY at x=50.0 (no contacts).

    At a 5.0 Å cutoff every nearest-pair distance is computable analytically:

    - ``AA``↔``BB``: min 0.5 Å (A3↔U1) — contact.
    - ``AA``↔``CC``: min 1.5 Å (A3↔ALA1) — contact.
    - ``BB``↔``CC``: min 1.0 Å (U1↔ALA1) — contact.
    - ``FAR`` is >5 Å from every other chain — no contacts.
    """
    structure = gemmi.Structure()
    structure.name = "CONTACTS"

    model = gemmi.Model("1")

    chain_aa = gemmi.Chain("AA")
    for i, (name, x) in enumerate([("G", 0.0), ("C", 1.0), ("A", 2.0)], start=1):
        _add_residue(chain_aa, name, seqid=i, x=x, y=0.0, z=0.0)
    model.add_chain(chain_aa)

    chain_bb = gemmi.Chain("BB")
    for i, (name, x) in enumerate([("U", 2.5), ("G", 10.0)], start=1):
        _add_residue(chain_bb, name, seqid=i, x=x, y=0.0, z=0.0)
    model.add_chain(chain_bb)

    chain_cc = gemmi.Chain("CC")
    _add_residue(chain_cc, "ALA", seqid=1, x=3.5, y=0.0, z=0.0)
    model.add_chain(chain_cc)

    chain_far = gemmi.Chain("FAR")
    _add_residue(chain_far, "GLY", seqid=1, x=50.0, y=0.0, z=0.0)
    model.add_chain(chain_far)

    structure.add_model(model)
    structure.cell = gemmi.UnitCell(200.0, 200.0, 200.0, 90.0, 90.0, 90.0)
    return structure


@pytest.fixture
def minimal_structure() -> gemmi.Structure:
    return build_minimal_structure()


@pytest.fixture
def contact_structure() -> gemmi.Structure:
    """Structure with deliberate inter-chain close contacts; see
    :func:`build_contact_structure` for the layout and pre-computed distances."""
    return build_contact_structure()


# ---------------------------------------------------------------------------
# Ribosome-shaped fixture for step 9 (chain assignment + tRNA state inference)
# ---------------------------------------------------------------------------


def _add_chain_of_residues(
    model: gemmi.Chain,
    chain_name: str,
    residues: list[tuple[str, int, tuple[float, float, float]]],
) -> None:
    """(Unused; see :func:`_make_chain`.) Retained as documentation of the
    shape ``residues`` takes — (residue_name, seqid_num, (x, y, z)) tuples."""


def _make_chain(
    name: str,
    residues: list[tuple[str, int, tuple[float, float, float]]],
) -> gemmi.Chain:
    chain = gemmi.Chain(name)
    for resname, seqid, (x, y, z) in residues:
        _add_residue(chain, resname, seqid=seqid, x=x, y=y, z=z)
    return chain


def build_ribosome_fixture() -> gemmi.Structure:
    """A synthetic structure shaped like a ribosome for step-9 inference tests.

    Chains (all on a deterministic Cartesian grid; one atom per residue):

    - ``S`` (SSU rRNA, 100 nts at residue numbers 1..100, x=0..99 along x-axis).
      Anchor positions inside this chain that the tests use:
      * residue 10 = "ssu_mrna" anchor (mRNA chain ``M`` sits 0.5 Å away)
      * residue 30 = "ssu_ptrna" anchor (P-tRNA ``TP`` sits 0.5 Å away)
      * residue 50 = "ssu_atrna" anchor (A-tRNA ``TA`` sits 0.5 Å away)
      * residue 70 = "ssu_etrna" anchor (E-tRNA ``TE`` sits 0.5 Å away)
    - ``L`` (LSU rRNA, 100 nts on a parallel x-axis offset y=20.0).
      Anchor positions:
      * residue 20 = "lsu_atrna" anchor (A-tRNA acceptor end neighbours)
      * residue 40 = "lsu_ptrna" anchor (P-tRNA acceptor end neighbours)
      * residue 60 = "lsu_etrna" anchor (E-tRNA acceptor end neighbours)
    - ``M`` (mRNA candidate, 6 nts): one residue at the SSU mrna anchor,
      five more strung out along x just outside the SSU envelope.
    - ``TA`` (A-tRNA candidate, 76 nts): residue 1 next to the SSU
      atrna anchor; residue 76 (CCA end) next to the LSU atrna anchor
      AND next to the factor protein ``EFTU``.
    - ``TP`` (P-tRNA candidate, 76 nts): residue 1 next to SSU ptrna
      anchor; residue 76 next to LSU ptrna anchor.
    - ``TE`` (E-tRNA candidate, 76 nts): residue 1 next to SSU etrna
      anchor; residue 76 next to LSU etrna anchor.
    - ``ASL`` (anticodon-stem-loop fragment, 20 nts only — under
      :data:`constants.ASL_FRAGMENT_MAX_LENGTH`). Sits far from
      everything; used to exercise the ``"**"`` state branch.
    - ``EFTU`` (elongation factor Tu, single protein residue, near TA's
      CCA-end residue 76). Used to exercise the §12.4 protein-factor
      label branch on the A-tRNA side.
    - ``IF2`` (initiation factor 2, single protein residue, near TP's
      CCA-end residue 76). Used to exercise the §12.4 protein-factor
      label branch on the P-tRNA side. Sits 20 Å from EFTU so it
      can't be picked up by the A-tRNA factor search.
    - ``L1`` (ribosomal protein, single protein residue, also near TA's
      CCA-end). The §12.4 algorithm must exclude this chain because the
      caller passes a ChainRef with ``is_ribosomal_protein=True``.

    All chain positions are pre-computed in Cartesian coordinates so the
    inference truth table for the fixture is deterministic.
    """
    structure = gemmi.Structure()
    structure.name = "RIBOFIXTURE"

    model = gemmi.Model("1")

    # SSU rRNA spans x = 0..99 at y = z = 0.
    ssu_residues: list[tuple[str, int, tuple[float, float, float]]] = [
        ("U", i, (float(i), 0.0, 0.0)) for i in range(1, 101)
    ]
    model.add_chain(_make_chain("S", ssu_residues))

    # LSU rRNA spans x = 0..99 at y = 20, z = 0.
    lsu_residues = [("U", i, (float(i), 20.0, 0.0)) for i in range(1, 101)]
    model.add_chain(_make_chain("L", lsu_residues))

    # mRNA: residue 1 sits 0.5 Å from SSU residue 10 (anchor for ssu_mrna).
    # Five more residues sit out at x = 110..114 (well outside the SSU).
    mrna_residues = [("A", 1, (10.5, 0.0, 0.0))] + [
        ("A", i, (110.0 + (i - 2), 0.0, 50.0)) for i in range(2, 7)
    ]
    model.add_chain(_make_chain("M", mrna_residues))

    # A-tRNA: 76 nts. Residue 1 near SSU atrna anchor (S/50, x=50, y=0).
    # Residue 76 (CCA end) near LSU atrna anchor (L/20, x=20, y=20).
    # Residues 2..75 strung along y=10 at x=50 → keeps them away from the
    # other functional sites so they don't make spurious contacts.
    atrna_residues = (
        [("A", 1, (50.5, 0.0, 0.0))]
        + [("A", i, (60.0 + i, 10.0, 30.0)) for i in range(2, 76)]
        + [("A", 76, (20.5, 20.0, 0.0))]
    )
    model.add_chain(_make_chain("TA", atrna_residues))

    # P-tRNA: 76 nts. Residue 1 near SSU ptrna anchor (S/30). Residue 76
    # near LSU ptrna anchor (L/40).
    ptrna_residues = (
        [("A", 1, (30.5, 0.0, 0.0))]
        + [("A", i, (60.0 + i, 10.0, 60.0)) for i in range(2, 76)]
        + [("A", 76, (40.5, 20.0, 0.0))]
    )
    model.add_chain(_make_chain("TP", ptrna_residues))

    # E-tRNA: 76 nts. Residue 1 near SSU etrna anchor (S/70). Residue 76
    # near LSU etrna anchor (L/60).
    etrna_residues = (
        [("A", 1, (70.5, 0.0, 0.0))]
        + [("A", i, (60.0 + i, 10.0, 90.0)) for i in range(2, 76)]
        + [("A", 76, (60.5, 20.0, 0.0))]
    )
    model.add_chain(_make_chain("TE", etrna_residues))

    # ASL fragment: 20 nts placed far from any anchor. Used for the
    # length-based ``"**"`` LSU-state branch.
    asl_residues = [("A", i, (200.0 + i, 0.0, 0.0)) for i in range(1, 21)]
    model.add_chain(_make_chain("ASL", asl_residues))

    # Elongation factor Tu: one residue 0.5 Å from A-tRNA's CCA-end
    # residue 76 at (20.5, 20.0, 0.0). Tests can swap A-tRNA into a
    # "no LSU contact" geometry by removing the L/20 anchor; the §12.4
    # rule should then pick EFTU and produce the "A/Elongation factor Tu"
    # state.
    eftu_residues = [("ALA", 1, (21.0, 20.0, 0.0))]
    model.add_chain(_make_chain("EFTU", eftu_residues))

    # Initiation factor 2: one residue 0.5 Å from P-tRNA's CCA-end
    # residue 76 at (40.5, 20.0, 0.0). Mirrors EFTU for the P-tRNA
    # factor-label path — biologically the P-site case is the
    # IF2 / eIF5B initiator-tRNA configuration (initiation factors
    # latching the acceptor stem of the initiator tRNA at the P site;
    # see REFERENCES.md for the canonical structures).
    if2_residues = [("ALA", 1, (41.0, 20.0, 0.0))]
    model.add_chain(_make_chain("IF2", if2_residues))

    # Ribosomal protein L1: also near A-tRNA's CCA-end (so it would win the
    # "minimum distance" race if not excluded). Tests pass
    # is_ribosomal_protein=True on the ChainRef so §12.4 skips it.
    l1_residues = [("ALA", 1, (20.7, 20.0, 0.0))]
    model.add_chain(_make_chain("L1", l1_residues))

    structure.add_model(model)
    structure.cell = gemmi.UnitCell(500.0, 500.0, 500.0, 90.0, 90.0, 90.0)
    return structure


@pytest.fixture
def ribosome_fixture() -> gemmi.Structure:
    """Synthetic ribosome-shaped structure for chain assignment + tRNA-state
    inference tests; see :func:`build_ribosome_fixture` for the layout."""
    return build_ribosome_fixture()


@pytest.fixture
def minimal_cif_path(tmp_path: Path) -> Path:
    """Write the minimal structure to ``tmp_path/test.cif`` and return the path."""
    structure = build_minimal_structure()
    path = tmp_path / "test.cif"
    structure.make_mmcif_document().write_file(str(path))
    return path


@pytest.fixture
def minimal_cif_gz_path(tmp_path: Path, minimal_cif_path: Path) -> Path:
    """Gzip the minimal mmCIF fixture and return the ``.cif.gz`` path."""
    gz_path = tmp_path / "test.cif.gz"
    gz_path.write_bytes(gzip.compress(minimal_cif_path.read_bytes()))
    return gz_path


@pytest.fixture
def minimal_cif_bytes(minimal_cif_path: Path) -> bytes:
    """Gzipped mmCIF bytes — the shape :func:`download_assembly_cif` returns."""
    return gzip.compress(minimal_cif_path.read_bytes())
