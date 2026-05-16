"""Complete-ribosome filter and three-way classification (spec §7, §8).

This module is the package's hardest pure-logic component. It owns:

- the two ribosomal-protein predicates (§13.1 narrow, §8.4 broader) — kept
  here rather than in :mod:`rcsb_client` so the rule lives next to the
  classification logic that consumes it;
- :func:`partition_rna_chains_by_role`, which fans RNA chains out into
  the canonical Rfam role buckets (§6.1, §11.1);
- :func:`determine_rrna_core` (§8.2), :func:`compute_dominant_superkingdom`
  (§8.3), :func:`extract_supporting_terms` (§8.6);
- :func:`classify_assembly`, the orchestrator that runs the §7 complete-
  ribosome filter, the §8.5 truth table, and the §7.4 completeness check
  in the correct order.

The module is pure: it has no network or filesystem dependencies. All
inputs are :class:`~.models.AssemblyContext` / :class:`~.models.ChainRef`
objects produced upstream by :mod:`rcsb_client`.
"""

from __future__ import annotations

import logging
import re
from collections import Counter
from typing import Any

from pydantic import BaseModel, Field

from ribosome_state_annotator import constants as C
from ribosome_state_annotator.config import CompletenessThresholds
from ribosome_state_annotator.models import (
    AssemblyContext,
    ChainRef,
    RibosomeClassification,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Result model
# ---------------------------------------------------------------------------


class ClassificationResult(BaseModel):
    """Outcome of :func:`classify_assembly` for one assembly.

    ``classification`` is ``None`` whenever the assembly was rejected by
    the §7 filter or the §8.5 truth-table fall-through cases; in those
    situations ``skip_reason`` carries one of the ``SKIP_*`` constants
    from :mod:`constants`. ``evidence`` is the §8.7 dict (always populated
    with as much as is known by the time of the decision). ``warnings``
    are surfaced verbatim on the final :class:`RibosomeAnnotation`.
    """

    classification: RibosomeClassification | None = None
    skip_reason: str | None = None
    evidence: dict[str, Any] = Field(default_factory=dict)
    warnings: list[str] = Field(default_factory=list)

    @property
    def is_skip(self) -> bool:
        return self.classification is None


# ---------------------------------------------------------------------------
# Ribosomal-protein predicates
# ---------------------------------------------------------------------------

_RIBOSOMAL_PROTEIN_SUBSTRING = "ribosomal protein"
_RIBOSOMAL_SUBUNIT_PROTEIN_SUBSTRING = "ribosomal subunit protein"

# Spec §8.4 — anchored regex matching common ribosomal-protein chain names.
# Verbatim from the spec; see the spec note for the deliberate trade-offs
# (e.g. it matches bare "S1" / "L11" / "uS3" but does not match "RPS3" with
# trailing digits because the anchor + word boundary disallow it).
_RIBOSOMAL_PROTEIN_NAME_REGEX = re.compile(
    r"^(MRP|MRPS|MRPL|RPS|RPL|"
    r"S\d{1,2}[A-Za-z]?|L\d{1,2}[A-Za-z]?|"
    r"uS\d{1,2}|uL\d{1,2}|bS\d{1,2}|bL\d{1,2}|"
    r"eS\d{1,2}|eL\d{1,2}|mS\d{1,2}|mL\d{1,2})\b"
)

_SUBUNIT_NUMBERS = ("30s", "40s", "50s", "60s")
_SUBUNIT_WORDS = ("subunit", "protein", "component")


def matches_ribosomal_protein_narrow(
    description: str | None,
    uniprot_name: str | None,
) -> bool:
    """§13.1 narrow rule: ``"ribosomal protein"`` substring (case-insensitive).

    Also matches the Ban-nomenclature pattern ``"ribosomal subunit protein"``
    (e.g. UniProt ``"Large ribosomal subunit protein uL2m"`` for human
    mitoribosomal uL2m), which is the canonical form for all modern
    mitoribosomal and cytoplasmic ribosomal proteins in UniProt.

    Used by :mod:`rcsb_client` to set ``ChainRef.is_ribosomal_protein`` at
    parse time, and by §12.4's factor search to exclude ribosomal proteins
    from the "nearest non-ribosomal protein at the CCA end" lookup.
    """
    for candidate in (description, uniprot_name):
        if not candidate:
            continue
        lowered = candidate.lower()
        if (
            _RIBOSOMAL_PROTEIN_SUBSTRING in lowered
            or _RIBOSOMAL_SUBUNIT_PROTEIN_SUBSTRING in lowered
        ):
            return True
    return False


def matches_ribosomal_protein_broad(chain: ChainRef) -> bool:
    """§8.4 broader rule used **only** for the §8.3 superkingdom vote.

    A chain qualifies if any of the following holds:

    1. Description contains ``"ribosomal protein"`` or
       ``"ribosomal subunit protein"`` (case-insensitive).
    2. Description matches the §8.4 anchored regex
       (``^(MRP|MRPS|MRPL|RPS|RPL|S\\d{1,2}…)\\b``).
    3. UniProt name contains ``"ribosomal protein"``.
    4. Description contains one of ``"30S"|"40S"|"50S"|"60S"`` together
       with one of ``"subunit"|"protein"|"component"``.

    This is intentionally broader than :func:`matches_ribosomal_protein_narrow`;
    do not unify them (spec §13.1 explicitly keeps them separate to avoid
    voting-side false positives leaking into the §12.4 factor exclusion).
    """
    description = chain.description
    uniprot = chain.uniprot_name

    # Sub-rule 3 (UniProt) is identical to the narrow rule's UniProt branch.
    if uniprot and _RIBOSOMAL_PROTEIN_SUBSTRING in uniprot.lower():
        return True

    if not description:
        return False

    desc_lower = description.lower()
    if (
        _RIBOSOMAL_PROTEIN_SUBSTRING in desc_lower
        or _RIBOSOMAL_SUBUNIT_PROTEIN_SUBSTRING in desc_lower
    ):
        return True

    if _RIBOSOMAL_PROTEIN_NAME_REGEX.match(description):
        return True

    has_subunit_number = any(n in desc_lower for n in _SUBUNIT_NUMBERS)
    has_subunit_word = any(w in desc_lower for w in _SUBUNIT_WORDS)
    return has_subunit_number and has_subunit_word


# ---------------------------------------------------------------------------
# Rfam role partitioning
# ---------------------------------------------------------------------------

# Role bucket keys produced by :func:`partition_rna_chains_by_role`. The
# four ribosomal-role keys mirror :data:`constants.RFAM_ROLE_MAP` values;
# ``"unmapped"`` holds RNA chains whose Rfam annotations don't match any
# role (typically mRNA, tRNA-without-RF00005, or unannotated chains).
RNA_ROLE_KEYS: tuple[str, ...] = (
    "ssu_main_rrna",
    "lsu_main_rrna",
    "lsu_associated_rrna",
    "trna",
    "unmapped",
)


def role_for_chain(chain: ChainRef) -> str | None:
    """Return the single Rfam role (§6.1) for a chain, or ``None`` if unmapped.

    Chains with multiple matched roles are extremely rare in practice — at
    most one ribosomal Rfam family should ever be annotated per chain. When
    it does happen, log a warning and pick the alphabetically-first role
    for determinism; downstream callers can re-classify if needed.
    """
    matched = {C.RFAM_ROLE_MAP[acc] for acc in chain.rfam_accessions if acc in C.RFAM_ROLE_MAP}
    if not matched:
        return None
    if len(matched) > 1:
        sorted_roles = sorted(matched)
        logger.warning(
            "chain %s has multiple ribosomal Rfam roles %s; using %r",
            chain.ife,
            sorted_roles,
            sorted_roles[0],
        )
        return sorted_roles[0]
    return next(iter(matched))


def partition_rna_chains_by_role(
    rna_chains: list[ChainRef],
) -> dict[str, list[ChainRef]]:
    """Bucket RNA chains by Rfam role (§6.1, §11.1).

    Returns a dict with keys :data:`RNA_ROLE_KEYS`. Every chain lands in
    exactly one bucket; ``"unmapped"`` holds anything not in the role map.
    """
    buckets: dict[str, list[ChainRef]] = {key: [] for key in RNA_ROLE_KEYS}
    for chain in rna_chains:
        role = role_for_chain(chain) or "unmapped"
        buckets[role].append(chain)
    return buckets


# ---------------------------------------------------------------------------
# rRNA core determination (§8.2)
# ---------------------------------------------------------------------------

RRNA_CORE_BACTERIAL = "bacterial_like"
RRNA_CORE_EUKARYOTIC = "eukaryotic_like"
RRNA_CORE_MIXED = "mixed"
RRNA_CORE_AMBIGUOUS = "ambiguous"

WARNING_MIXED_RRNA_CORE = "mixed_rrna_core_treated_as_bacterial_like"


# Archaeal Rfam accessions (added to RFAM_ROLE_MAP so SSU/LSU detection
# succeeds for archaeal entries). The rRNA-core classifier treats an
# archaeal pair as ``bacterial_like`` so the §8.5 truth table reaches
# the archaeal-protein-superkingdom branch and correctly emits
# ``archaeal_ribosome_not_supported`` rather than ``ambiguous_rrna_core``.
# Biological justification: archaeal ribosomes are structurally and
# evolutionarily closer to bacterial than to eukaryotic, so the
# bacterial_like label is the right v1 bucket for the §8.5 dispatch.
_ARCHAEAL_SSU_RFAM = "RF01959"
_ARCHAEAL_LSU_RFAM = "RF02540"


def determine_rrna_core(
    ssu_chains: list[ChainRef],
    lsu_chains: list[ChainRef],
) -> str:
    """Determine the rRNA core flavour (spec §8.2).

    Returns one of :data:`RRNA_CORE_BACTERIAL`, :data:`RRNA_CORE_EUKARYOTIC`,
    :data:`RRNA_CORE_MIXED`, or :data:`RRNA_CORE_AMBIGUOUS`. Callers handle
    the "mixed" → bacterial-like demotion and the "ambiguous" → skip
    themselves (see :func:`classify_assembly`).

    An archaeal SSU+LSU pair (Rfam ``RF01959`` + ``RF02540``) is treated
    as ``bacterial_like`` so the §8.5 truth table can route it through
    the dominant-superkingdom check and emit ``archaeal_ribosome_not_supported``
    rather than ``ambiguous_rrna_core``. Spec §3.2 keeps archaea
    out-of-scope for v1; this routing only ensures the right skip
    reason is reported.
    """
    ssu_accessions = {acc for chain in ssu_chains for acc in chain.rfam_accessions}
    lsu_accessions = {acc for chain in lsu_chains for acc in chain.rfam_accessions}

    bacterial_like = (
        C.BACTERIAL_RRNA_CORE_SSU in ssu_accessions and C.BACTERIAL_RRNA_CORE_LSU in lsu_accessions
    )
    eukaryotic_like = (
        C.EUKARYOTIC_RRNA_CORE_SSU in ssu_accessions
        and C.EUKARYOTIC_RRNA_CORE_LSU in lsu_accessions
    )
    archaeal_like = (
        _ARCHAEAL_SSU_RFAM in ssu_accessions and _ARCHAEAL_LSU_RFAM in lsu_accessions
    )
    if bacterial_like and eukaryotic_like:
        return RRNA_CORE_MIXED
    if bacterial_like:
        return RRNA_CORE_BACTERIAL
    if eukaryotic_like:
        return RRNA_CORE_EUKARYOTIC
    if archaeal_like:
        # Bucket with bacterial_like so the §8.5 archaeal-protein check fires.
        return RRNA_CORE_BACTERIAL
    return RRNA_CORE_AMBIGUOUS


# ---------------------------------------------------------------------------
# Dominant protein superkingdom (§8.3)
# ---------------------------------------------------------------------------

SUPERKINGDOM_UNKNOWN = "unknown"


def compute_dominant_superkingdom(
    voting_chains: list[ChainRef],
    *,
    min_voters: int = C.MIN_RIBOSOMAL_PROTEIN_VOTES,
) -> tuple[str, dict[str, int]]:
    """Compute the dominant superkingdom vote (spec §8.3).

    Returns ``(dominant_superkingdom, vote_counts)`` where
    ``dominant_superkingdom`` is one of ``"Bacteria"``, ``"Eukaryota"``,
    ``"Archaea"``, or :data:`SUPERKINGDOM_UNKNOWN`. Ties are broken in the
    :data:`constants.SUPERKINGDOM_TIEBREAK_ORDER` order
    (Bacteria > Eukaryota > Archaea).

    Returns ``("unknown", counts)`` when fewer than ``min_voters`` chains
    carry a non-null superkingdom — the spec mandates ≥3 voters before the
    classification is trusted.
    """
    counts: Counter[str] = Counter()
    for chain in voting_chains:
        sk = chain.superkingdom
        if sk:
            counts[sk] += 1

    total_votes = sum(counts.values())
    if total_votes < min_voters:
        return (SUPERKINGDOM_UNKNOWN, dict(counts))

    if not counts:
        # All voting chains had null superkingdom.
        return (SUPERKINGDOM_UNKNOWN, {})

    max_votes = max(counts.values())
    leaders = {sk for sk, v in counts.items() if v == max_votes}
    for sk in C.SUPERKINGDOM_TIEBREAK_ORDER:
        if sk in leaders:
            return (sk, dict(counts))
    # Leader is outside the canonical superkingdom set (e.g. a clade name
    # for virus-source proteins, which are out of scope for v1). Fall back
    # to the leader as-is; the §8.5 truth table will route it to
    # ``insufficient_taxonomy_evidence``.
    return (sorted(leaders)[0], dict(counts))


# ---------------------------------------------------------------------------
# Organellar supporting evidence (§8.6)
# ---------------------------------------------------------------------------

WARNING_ORGANELLAR_NO_KEYWORDS = "organellar_classification_without_keyword_support"

# Compile word-boundary patterns once. "exact word match" per §8.6 means
# whole-word, not substring (so "mitochondrially" does NOT count for
# "mitochondrial").
_ORGANELLAR_KEYWORD_PATTERNS: dict[str, re.Pattern[str]] = {
    kw: re.compile(rf"\b{re.escape(kw)}\b", re.IGNORECASE) for kw in C.ORGANELLAR_KEYWORDS
}


def extract_supporting_terms(chains: list[ChainRef]) -> list[str]:
    """Return the deduplicated, sorted organellar keywords (§8.6) found in
    the descriptions of ``chains``.

    Whole-word match, case-insensitive. Each keyword appears at most once
    in the returned list regardless of how many chains mention it.
    """
    found: set[str] = set()
    for chain in chains:
        description = chain.description
        if not description:
            continue
        for keyword, pattern in _ORGANELLAR_KEYWORD_PATTERNS.items():
            if pattern.search(description):
                found.add(keyword)
    return sorted(found)


# ---------------------------------------------------------------------------
# Orchestrator (§7 + §8)
# ---------------------------------------------------------------------------

WARNING_LOW_PROTEIN_COUNT = "low_ribosomal_protein_count"


def _classification_rule(rrna_core: str, dominant_sk: str) -> RibosomeClassification | None:
    """Return the §8.5 classification or ``None`` if the (rrna_core, sk)
    combination falls into the ambiguous skip branch.

    Callers handle the Archaea / unknown / insufficient short-circuits
    BEFORE calling this; this function is only the bacterial/eukaryotic
    truth table.
    """
    if rrna_core == RRNA_CORE_BACTERIAL and dominant_sk == "Bacteria":
        return "bacterial_ribosome"
    if rrna_core == RRNA_CORE_EUKARYOTIC and dominant_sk == "Eukaryota":
        return "eukaryotic_ribosome"
    if rrna_core == RRNA_CORE_BACTERIAL and dominant_sk == "Eukaryota":
        return "eukaryotic_organellar_ribosome"
    return None


def _rule_string_for(classification: RibosomeClassification) -> str:
    if classification == "bacterial_ribosome":
        return C.CLASSIFICATION_RULE_BACTERIAL
    if classification == "eukaryotic_ribosome":
        return C.CLASSIFICATION_RULE_EUKARYOTIC
    return C.CLASSIFICATION_RULE_ORGANELLAR


def _is_nmr(experimental_methods: list[str]) -> bool:
    return any("NMR" in method.upper() for method in experimental_methods)


def _collect_rfam_in_role(chains: list[ChainRef], role_set: frozenset[str]) -> list[str]:
    """Sorted unique list of Rfam accessions on ``chains`` that fall in ``role_set``."""
    return sorted({acc for chain in chains for acc in chain.rfam_accessions if acc in role_set})


def classify_assembly(
    assembly: AssemblyContext,
    *,
    thresholds: CompletenessThresholds | None = None,
    strict_complete_check: bool = False,
) -> ClassificationResult:
    """Apply the §7 complete-ribosome filter and the §8 three-way classifier.

    Order of decisions (spec §7 then §8 then §7.4):

    1. NMR check (§7.3) — skip with ``nmr_structure_not_supported``.
    2. SSU/LSU presence (§7.4) — skip with ``partial_ribosome_missing_ssu_or_lsu``.
    3. rRNA core (§8.2) — "ambiguous" skips with ``ambiguous_rrna_core``;
       "mixed" emits a warning and is treated as bacterial-like.
    4. Voting eligible ribosomal protein chains (§8.4).
    5. Dominant superkingdom vote (§8.3) — Archaea skips with
       ``archaeal_ribosome_not_supported``; "unknown" skips with
       ``insufficient_taxonomy_evidence``.
    6. Primary classification truth table (§8.5) — eukaryotic-like rRNA
       with Bacteria-dominant proteins skips with
       ``ambiguous_ribosome_classification``.
    7. Organellar supporting evidence (§8.6) — warn-only.
    8. Completeness threshold (§7.4) — strict mode skips with
       ``likely_partial_ribosome_low_ribosomal_protein_count``; non-strict
       emits ``low_ribosomal_protein_count`` and continues.

    Returns a :class:`ClassificationResult` for the assembly. The function
    never raises for "this isn't supported" cases — every such case is
    represented as a structured skip, per spec §16.
    """
    thresholds = thresholds or CompletenessThresholds()
    warnings: list[str] = []
    evidence: dict[str, Any] = {}

    # 1. NMR
    if _is_nmr(assembly.experimental_methods):
        return ClassificationResult(
            skip_reason=C.SKIP_NMR,
            evidence={"experimental_methods": list(assembly.experimental_methods)},
            warnings=warnings,
        )

    # 2. SSU/LSU presence
    by_role = partition_rna_chains_by_role(assembly.rna_chains)
    ssu_chains = by_role["ssu_main_rrna"]
    lsu_chains = by_role["lsu_main_rrna"]
    ssu_rfam = _collect_rfam_in_role(ssu_chains, C.SSU_MAIN_RRNA)
    lsu_rfam = _collect_rfam_in_role(lsu_chains, C.LSU_MAIN_RRNA)
    evidence["ssu_rfam"] = ssu_rfam
    evidence["lsu_rfam"] = lsu_rfam

    if not ssu_chains or not lsu_chains:
        return ClassificationResult(
            skip_reason=C.SKIP_PARTIAL_MISSING_SSU_OR_LSU,
            evidence=evidence,
            warnings=warnings,
        )

    # 3. rRNA core
    rrna_core = determine_rrna_core(ssu_chains, lsu_chains)
    if rrna_core == RRNA_CORE_AMBIGUOUS:
        evidence["rrna_core"] = rrna_core
        return ClassificationResult(
            skip_reason=C.SKIP_AMBIGUOUS_RRNA_CORE,
            evidence=evidence,
            warnings=warnings,
        )
    if rrna_core == RRNA_CORE_MIXED:
        warnings.append(WARNING_MIXED_RRNA_CORE)
        rrna_core = RRNA_CORE_BACTERIAL
    evidence["rrna_core"] = rrna_core

    # 4. Voting eligible chains (§8.4)
    voting_chains = [c for c in assembly.protein_chains if matches_ribosomal_protein_broad(c)]
    voters_with_taxonomy = [c for c in voting_chains if c.superkingdom]

    # 5. Dominant superkingdom
    dominant_sk, vote_counts = compute_dominant_superkingdom(voting_chains)
    evidence["dominant_ribosomal_protein_superkingdom"] = dominant_sk
    evidence["ribosomal_protein_superkingdom_votes"] = vote_counts
    evidence["ribosomal_protein_chains_voting"] = len(voters_with_taxonomy)

    if dominant_sk == "Archaea":
        return ClassificationResult(
            skip_reason=C.SKIP_ARCHAEAL,
            evidence=evidence,
            warnings=warnings,
        )
    if dominant_sk == SUPERKINGDOM_UNKNOWN:
        return ClassificationResult(
            skip_reason=C.SKIP_INSUFFICIENT_TAXONOMY,
            evidence=evidence,
            warnings=warnings,
        )

    # 6. Primary classification truth table
    classification = _classification_rule(rrna_core, dominant_sk)
    if classification is None:
        return ClassificationResult(
            skip_reason=C.SKIP_AMBIGUOUS_CLASSIFICATION,
            evidence=evidence,
            warnings=warnings,
        )
    rule = _rule_string_for(classification)
    evidence["rule"] = rule

    # 7. Organellar supporting evidence (warn-only — does NOT change classification)
    supporting_terms = extract_supporting_terms(voting_chains)
    evidence["supporting_terms"] = supporting_terms
    if classification == "eukaryotic_organellar_ribosome" and not supporting_terms:
        warnings.append(WARNING_ORGANELLAR_NO_KEYWORDS)

    # 8. Completeness threshold
    threshold = thresholds.for_classification(classification)
    if len(voting_chains) < threshold:
        if strict_complete_check:
            return ClassificationResult(
                skip_reason=C.SKIP_LOW_PROTEIN_COUNT,
                evidence=evidence,
                warnings=warnings,
            )
        warnings.append(WARNING_LOW_PROTEIN_COUNT)

    return ClassificationResult(
        classification=classification,
        evidence=evidence,
        warnings=warnings,
    )
