"""Curated biological reference data.

This module is the single source of truth for the values the contact-transfer
workflow (spec §3.3) depends on:

- Rfam accessions and their ribosomal roles (§6.1).
- Functional-site anchor nucleotide sets for the bacterial (E. coli / 5J7L)
  and eukaryotic cytoplasmic (yeast / 5TBW) reference ribosomes (§6.3).
- Organellar-classification supporting keywords (§8.6).
- Default neighbour-search cutoff, ASL fragment threshold, and ligand
  exclusions (§10.2, §12.1, §13.3).

None of these values should be modified at runtime; the module-level
constants are typed with ``Final`` and use immutable containers
(``frozenset``, ``tuple``) so any accidental mutation surfaces as a type
error or attribute error.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Final

# ---------------------------------------------------------------------------
# Rfam accessions and role map (spec §6.1)
# ---------------------------------------------------------------------------

SSU_MAIN_RRNA: Final[frozenset[str]] = frozenset({"RF00177", "RF01960", "RF01959"})
"""SSU main rRNA Rfam accessions: bacterial-like 16S (RF00177), eukaryotic-like
18S (RF01960), and archaeal SSU (RF01959). Archaeal SSU is included so that
archaeal deposits successfully reach the §8.3 superkingdom vote and trigger
the proper `archaeal_ribosome_not_supported` skip rather than being dropped
upstream as `partial_ribosome_missing_ssu_or_lsu`. Spec §3.2 still lists
archaeal ribosomes as out-of-scope for v1 annotation."""

LSU_MAIN_RRNA: Final[frozenset[str]] = frozenset({"RF02541", "RF02543", "RF02540"})
"""LSU main rRNA Rfam accessions: bacterial-like 23S (RF02541), eukaryotic-like
25S/28S (RF02543), and archaeal LSU (RF02540). See SSU_MAIN_RRNA for the
archaeal-inclusion rationale."""

LSU_ASSOCIATED_RRNA: Final[frozenset[str]] = frozenset({"RF00001", "RF00002"})
"""Associated LSU rRNAs: 5S (RF00001) and 5.8S-like (RF00002)."""

TRNA: Final[frozenset[str]] = frozenset({"RF00005"})
"""tRNA Rfam accession."""

RFAM_ROLE_MAP: Final[Mapping[str, str]] = {
    "RF00177": "ssu_main_rrna",  # bacterial-like 16S SSU rRNA
    "RF01960": "ssu_main_rrna",  # eukaryotic-like 18S SSU rRNA
    "RF01959": "ssu_main_rrna",  # archaeal SSU rRNA (kept for archaeal-skip detection)
    "RF02541": "lsu_main_rrna",  # bacterial-like 23S LSU rRNA
    "RF02543": "lsu_main_rrna",  # eukaryotic-like 25S/28S LSU rRNA
    "RF02540": "lsu_main_rrna",  # archaeal LSU rRNA (kept for archaeal-skip detection)
    "RF00002": "lsu_associated_rrna",  # 5.8S rRNA-like component
    "RF00001": "lsu_associated_rrna",  # 5S rRNA-like component
    "RF00005": "trna",
}
"""Mapping of Rfam accession to ribosomal role string used by the inference layer."""

# Rfam accessions that signal the bacterial vs eukaryotic rRNA core (spec §8.2).
BACTERIAL_RRNA_CORE_SSU: Final[str] = "RF00177"
BACTERIAL_RRNA_CORE_LSU: Final[str] = "RF02541"
EUKARYOTIC_RRNA_CORE_SSU: Final[str] = "RF01960"
EUKARYOTIC_RRNA_CORE_LSU: Final[str] = "RF02543"

# Rfam accessions corresponding to the legacy CSV columns (spec §11.1, §15.3).
LSU_MEDIUM_RFAM: Final[str] = "RF00002"  # 5.8S-like
LSU_SMALL_RFAM: Final[str] = "RF00001"  # 5S-like

# ---------------------------------------------------------------------------
# Functional-site anchor nucleotide sets (spec §6.3)
# ---------------------------------------------------------------------------
#
# These are BGSU unit IDs (PDB|model|chain|residue_name|residue_number) for
# conserved rRNA nucleotides known to contact mRNA / A-tRNA / P-tRNA /
# E-tRNA in well-characterised reference ribosomes. They are the input to
# the BGSU correspondence API (§5.2), which transfers them to equivalent
# positions in any homologous query ribosome.

# The dictionary keys are the site labels used throughout the codebase.
REFERENCE_SITE_KEYS: Final[tuple[str, ...]] = (
    "ssu_mrna",
    "ssu_atrna",
    "ssu_ptrna",
    "ssu_etrna",
    "lsu_atrna",
    "lsu_ptrna",
    "lsu_etrna",
)

ECOLI_REFERENCE_UNITS: Final[Mapping[str, tuple[str, ...]]] = {
    "ssu_mrna": (
        "5J7L|1|AA|G|926",
        "5J7L|1|AA|4OC|1402",
        "5J7L|1|AA|C|1403",
    ),
    "ssu_ptrna": (
        "5J7L|1|AA|G|1338",
        "5J7L|1|AA|A|1339",
        "5J7L|1|AA|C|1400",
    ),
    "ssu_atrna": (
        "5J7L|1|AA|G|530",
        "5J7L|1|AA|A|1492",
        "5J7L|1|AA|A|1493",
    ),
    "ssu_etrna": (
        "5J7L|1|AA|G|693",
        "5J7L|1|AA|A|694",
    ),
    "lsu_atrna": (
        "5J7L|1|DA|G|2553",
        "5J7L|1|DA|G|2583",
        "5J7L|1|DA|U|2585",
    ),
    "lsu_ptrna": (
        "5J7L|1|DA|OMG|2251",
        "5J7L|1|DA|G|2252",
        "5J7L|1|DA|G|2253",
    ),
    "lsu_etrna": (
        "5J7L|1|DA|G|2112",
        "5J7L|1|DA|G|2421",
        "5J7L|1|DA|C|2422",
    ),
}
"""E. coli 70S (PDB ``5J7L``) functional-site anchors. Used for bacterial
ribosomes and, in v1, for eukaryotic organellar ribosomes (spec §6.4)."""

YEAST_REFERENCE_UNITS: Final[Mapping[str, tuple[str, ...]]] = {
    "ssu_mrna": (
        "7ZW0|1|2|G|1150",
        "7ZW0|1|2|C|1639",
        "7ZW0|1|2|C|1640",
    ),
    "ssu_atrna": (
        "7ZW0|1|2|G|577",
        "7ZW0|1|2|A|1755",
        "7ZW0|1|2|A|1756",
    ),
    "ssu_ptrna": (
        "7ZW0|1|2|G|1575",
        "7ZW0|1|2|A|1576",
        "7ZW0|1|2|C|1637",
    ),
    "ssu_etrna": (
        "7ZW0|1|2|G|904",
        "7ZW0|1|2|A|905",
    ),
    "lsu_atrna": (
        "7ZW0|1|LA|G|2922",
        "7ZW0|1|LA|G|2952",
        "7ZW0|1|LA|U|2954",
    ),
    "lsu_ptrna": (
        "7ZW0|1|LA|G|2619",
        "7ZW0|1|LA|G|2620",
        "7ZW0|1|LA|G|2621",
    ),
    "lsu_etrna": (
        "7ZW0|1|LA|G|2793",
        "7ZW0|1|LA|G|2794",
    ),
}
"""Yeast 80S functional-site anchors. Used for eukaryotic cytoplasmic
ribosomes (spec §6.4).

**Diverges from spec §6.3 — uses 7ZW0 instead of 5TBW.** Three reasons:

1. **5TBW assembly-1 chain ``A`` isn't BGSU's NR representative.**
   Querying BGSU with ``5TBW|1|A|...`` returns zero mappings; the
   representative is the assembly-2 chain ``sR``. 7ZW0 sidesteps the
   issue — its SSU chain ``2`` IS its own BGSU representative and
   queries work directly.
2. **5TBW is missing residue 2454 on chain ``1``** (LSU ``G|2454``,
   one of the lsu_etrna anchors). With BGSU's all-or-nothing batch
   semantic this poisons the entire lsu_etrna query. 7ZW0 has 2454
   resolved as ``G``.
3. **7ZW0 is an active 80S** (FAP-80S complex, rotated state) with
   bound mRNA + two tRNAs — a more representative reference for the
   contact-transfer workflow than 5TBW's empty 80S.

**`lsu_etrna` drops position 2454** (kept only ``G|2793`` + ``G|2794``).
Position 2454 sits in a less-conserved region of the LSU that isn't
mapped across all eukaryotic Rfam representatives; including it in a
batch query causes BGSU's intersection-semantic to silently drop
target PDBs (e.g. 6Y57) that lack a 2454 equivalent. The two
remaining well-conserved anchors are sufficient — spec §5.2.3
explicitly says "a site with at least one mapped nucleotide in the
assembly is sufficient for that site's contact search".

All anchor positions and residue identities verified against the live
7ZW0 mmCIF. Same S. cerevisiae residue numbering as 5TBW so the
positions match the spec values one-to-one."""

# ---------------------------------------------------------------------------
# Classification (spec §8)
# ---------------------------------------------------------------------------

ORGANELLAR_KEYWORDS: Final[tuple[str, ...]] = (
    "mitochondrial",
    "chloroplast",
    "plastid",
    "organellar",
)
"""Case-insensitive supporting-evidence keywords for organellar classification
(§8.6). Their presence/absence does not change the classification — they only
suppress the ``organellar_classification_without_keyword_support`` warning."""

SUPERKINGDOM_TIEBREAK_ORDER: Final[tuple[str, ...]] = (
    "Bacteria",
    "Eukaryota",
    "Archaea",
)
"""Tiebreak order for the dominant-superkingdom vote (§8.3). Bacteria beats
Eukaryota on a tie because organellar ribosomes typically have a clean
Eukaryota plurality; a tie usually means thin/chimeric data and
bacterial-like is the safer default."""

MIN_RIBOSOMAL_PROTEIN_VOTES: Final[int] = 3
"""Minimum number of ribosomal protein chains with non-null taxonomy required
for the §8.3 superkingdom vote to count. Below this, return ``unknown`` and
emit ``insufficient_taxonomy_evidence``."""

# ---------------------------------------------------------------------------
# Inference (spec §10, §12, §13)
# ---------------------------------------------------------------------------

DEFAULT_CONTACT_CUTOFF_ANGSTROM: Final[float] = 5.0
"""Default Gemmi neighbour-search cutoff (§10.2). Configurable via CLI."""

ASL_FRAGMENT_MAX_LENGTH: Final[int] = 30
"""tRNA chains shorter than this many nucleotides are treated as anticodon
stem-loop fragments and trigger the ``**`` LSU-state suffix (§12.1, §12.2)."""

DEFAULT_LIGAND_EXCLUSIONS: Final[frozenset[str]] = frozenset({"HOH"})
"""Default comp_ids excluded from the ``bound_ligands`` output (§13.3).
Water only — magnesium and other biologically meaningful cofactors are kept."""

# ---------------------------------------------------------------------------
# Status and skip-reason vocabularies
# ---------------------------------------------------------------------------
#
# Defining these as constants (rather than scattering string literals
# throughout the codebase) prevents typos and gives callers a stable enum
# to switch on.

SKIP_NMR: Final[str] = "nmr_structure_not_supported"
SKIP_ARCHAEAL: Final[str] = "archaeal_ribosome_not_supported"
SKIP_PARTIAL_MISSING_SSU_OR_LSU: Final[str] = "partial_ribosome_missing_ssu_or_lsu"
SKIP_LOW_PROTEIN_COUNT: Final[str] = "likely_partial_ribosome_low_ribosomal_protein_count"
SKIP_AMBIGUOUS_CLASSIFICATION: Final[str] = "ambiguous_ribosome_classification"
SKIP_AMBIGUOUS_RRNA_CORE: Final[str] = "ambiguous_rrna_core"
SKIP_INSUFFICIENT_TAXONOMY: Final[str] = "insufficient_taxonomy_evidence"
SKIP_FRAGMENTED_RIBOSOME: Final[str] = "fragmented_ribosome_not_supported"
"""Emitted when the assembly carries asymmetric SSU/LSU rRNA chain
counts (e.g. one SSU chain with a 28S rRNA split across three LSU
chains, as in Tetrahymena / Trypanosoma / human mitoribosome deposits).
The canonical BGSU anchors are numbered against a single reference
rRNA molecule and don't transfer onto fragmented chains, so the
contact-transfer step cannot resolve A/P/E sites. Multi-ribosome
bundles (matched SSU/LSU counts) are still supported."""

CLASSIFICATION_RULE_BACTERIAL: Final[str] = "bacterial_like_rfam_plus_bacterial_proteins"
CLASSIFICATION_RULE_EUKARYOTIC: Final[str] = "eukaryotic_like_rfam_plus_eukaryotic_proteins"
CLASSIFICATION_RULE_ORGANELLAR: Final[str] = "bacterial_like_rfam_plus_eukaryotic_proteins"
