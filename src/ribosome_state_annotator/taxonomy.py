"""Per-assembly NCBI taxonomy aggregation (spec §34).

The :func:`aggregate_assembly_lineage` function takes the rRNA chains of
one biological assembly and returns the assembly's consensus
:class:`AssemblyTaxonomy`. Voting is restricted to rRNA chains because
heterologous in-vitro reconstitutions routinely mix ribosomes from one
organism with tRNAs from another — tRNAs, mRNAs, and bound factors
appear in :attr:`AssemblyTaxonomy.source_organisms` for visibility but
do not contribute to the consensus lineage.

Algorithm — hybrid strict-LCA then majority-mode constrained to branch:

1. Drop voting chains with an empty :attr:`ChainRef.taxonomy_lineage`.
2. Walk lineage depth-by-depth (depth 1 = "cellular organisms").
3. **Strict phase** — while every voting chain has the *same*
   :attr:`TaxonNode.tax_id` at the current depth, append that node to
   the consensus and continue. ``strict_lca_depth`` records the last
   depth at which this held.
4. **Majority phase** — once chains disagree, drop chains that have
   already terminated *or* that diverged in the strict phase. Among the
   remaining chains, take the mode at each subsequent depth. This
   guarantees the consensus lineage is a valid root-to-leaf path
   through the NCBI tree, because each majority pick is constrained to
   descendants of the last strict node.
5. Stop when no voting chain extends to the current depth.

Edge cases:

- All voting chains share a single full-length lineage → ``is_mixed=False``,
  ``strict_lca_depth == len(lineage)``.
- All voting chains lack lineage (synthetic constructs only) → returns
  ``None`` and the caller emits ``no_source_organism_taxonomy``.
- One voting chain → its lineage is returned verbatim.

The function is pure: no network or filesystem access.
"""

from __future__ import annotations

from collections import Counter

from ribosome_state_annotator.models import (
    AssemblyTaxonomy,
    ChainRef,
    OrganismRef,
    TaxonNode,
)

# NCBI domain (a.k.a. superkingdom) ``tax_id``s. Used to populate the
# convenience ``AssemblyTaxonomy.domain`` field. Viruses are included
# even though we don't classify viral particles as ribosomes — better
# to surface the label than to silently mis-tag a viral entry.
_DOMAIN_TAX_IDS: dict[int, str] = {
    2: "Bacteria",
    2157: "Archaea",
    2759: "Eukaryota",
    10239: "Viruses",
}


def aggregate_assembly_lineage(
    rrna_chains: list[ChainRef],
    *,
    all_chains: list[ChainRef] | None = None,
) -> AssemblyTaxonomy | None:
    """Aggregate the assembly's NCBI taxonomy from its rRNA chains.

    Args:
        rrna_chains: All rRNA chains (``ssu_main_rrna`` +
            ``lsu_main_rrna`` + ``lsu_associated_rrna``) — the voting
            set. Chains without a ``taxonomy_lineage`` are silently
            skipped.
        all_chains: Every chain in the assembly (rRNA + protein + tRNA
            + mRNA + factor). Used only to populate
            :attr:`AssemblyTaxonomy.source_organisms`. ``None`` falls
            back to ``rrna_chains``.

    Returns:
        :class:`AssemblyTaxonomy` when at least one rRNA chain has a
        non-empty lineage, otherwise ``None`` so the caller can emit
        the ``no_source_organism_taxonomy`` warning.
    """
    voting_chains = [c for c in rrna_chains if c.taxonomy_lineage]
    if not voting_chains:
        return None

    consensus, strict_depth, is_mixed = _hybrid_consensus(voting_chains)
    organisms = _collect_source_organisms(all_chains if all_chains is not None else rrna_chains)

    domain = _resolve_domain(consensus)
    species = _resolve_species(consensus, domain)

    return AssemblyTaxonomy(
        lineage=tuple(consensus),
        domain=domain,
        species=species,
        voting_chains=tuple(c.ife for c in voting_chains),
        n_voting_chains=len(voting_chains),
        strict_lca_depth=strict_depth,
        is_mixed=is_mixed,
        source_organisms=organisms,
    )


def _hybrid_consensus(
    voting_chains: list[ChainRef],
) -> tuple[list[TaxonNode], int, bool]:
    """Run the strict-LCA → majority-mode walk.

    Returns ``(consensus_lineage, strict_lca_depth, is_mixed)``.
    """
    # Index each chain's lineage by depth for O(1) per-depth lookups.
    lineages: list[dict[int, TaxonNode]] = [
        {node.depth: node for node in chain.taxonomy_lineage} for chain in voting_chains
    ]
    if not lineages:
        return [], 0, False

    # Strict phase: walk from depth 1 upwards while every chain agrees.
    all_depths = sorted({d for lineage in lineages for d in lineage})
    consensus: list[TaxonNode] = []
    strict_depth = 0
    diverge_depth: int | None = None
    for depth in all_depths:
        nodes_at_depth = [lineage.get(depth) for lineage in lineages]
        if any(n is None for n in nodes_at_depth):
            # At least one chain ran out of lineage; strict can't
            # continue beyond every-chain-present. Stop strict here.
            diverge_depth = depth
            break
        unique_ids = {n.tax_id for n in nodes_at_depth if n is not None}
        if len(unique_ids) == 1:
            consensus.append(nodes_at_depth[0])  # type: ignore[arg-type]
            strict_depth = depth
            continue
        diverge_depth = depth
        break

    if diverge_depth is None:
        # Every chain agreed all the way to its leaf — fully resolved.
        return consensus, strict_depth, False

    # Majority phase: filter to chains that agreed up to strict_depth
    # (i.e. their depth-strict_depth node matches our consensus). These
    # are the chains rooted in the strict-LCA subtree; they get to vote
    # at depth diverge_depth and onward.
    if strict_depth > 0:
        anchor_id = consensus[-1].tax_id
        eligible = [
            lineage for lineage in lineages
            if lineage.get(strict_depth) is not None
            and lineage[strict_depth].tax_id == anchor_id
        ]
    else:
        # No strict-LCA reached — everything diverges at the root. This
        # is essentially impossible for biological deposits (every
        # cellular chain has "cellular organisms" at depth 1).
        eligible = lineages

    current_depth = diverge_depth
    while True:
        nodes_at_depth = [lineage.get(current_depth) for lineage in eligible]
        present = [n for n in nodes_at_depth if n is not None]
        if not present:
            break
        # Mode by tax_id.
        counts = Counter(n.tax_id for n in present)
        winning_id, _ = counts.most_common(1)[0]
        # Pick the canonical TaxonNode for the winning id (first chain
        # whose depth-current_depth node has that id).
        winning_node = next(n for n in present if n.tax_id == winning_id)
        consensus.append(winning_node)
        # Constrain future picks to chains that voted for the winner at
        # this depth — keeps the lineage on a single tree branch.
        eligible = [
            lineage for lineage in eligible
            if lineage.get(current_depth) is not None
            and lineage[current_depth].tax_id == winning_id
        ]
        if not eligible:
            break
        current_depth += 1

    return consensus, strict_depth, True


def _resolve_domain(consensus: list[TaxonNode]) -> str | None:
    """Find the domain (Bacteria/Archaea/Eukaryota/Viruses) in the
    consensus lineage."""
    for node in consensus:
        if node.tax_id in _DOMAIN_TAX_IDS:
            return _DOMAIN_TAX_IDS[node.tax_id]
    return None


def _resolve_species(consensus: list[TaxonNode], domain: str | None) -> str | None:
    """Best-effort species name.

    NCBI doesn't tag species with a stable rank label in the lineage;
    we pick the deepest node whose name looks like a binomial
    ("Genus species"), falling back to the leaf node when no node
    matches. Returns ``None`` for empty lineages or pure-domain
    results.
    """
    if not consensus:
        return None
    for node in reversed(consensus):
        if " " in node.name and not node.name.startswith(domain or ""):
            return node.name
    # Final fallback — the leaf, but only if it's not the domain itself.
    leaf = consensus[-1]
    if domain and leaf.name == domain:
        return None
    return leaf.name


def _collect_source_organisms(chains: list[ChainRef]) -> tuple[OrganismRef, ...]:
    """Every unique ``(tax_id, scientific_name)`` across all chains.

    Synthetic chains (no source organism) collapse to a single
    ``(None, None)`` entry.
    """
    seen: list[OrganismRef] = []
    seen_keys: set[tuple[int | None, str | None]] = set()
    for chain in chains:
        key = (chain.tax_id, chain.scientific_name)
        if key in seen_keys:
            continue
        seen_keys.add(key)
        seen.append(OrganismRef(tax_id=chain.tax_id, scientific_name=chain.scientific_name))
    return tuple(seen)
