# References

Citations for the biological data sources, classification rules, and
naming conventions the package relies on.

## tRNA functional-state notation

The `<SSU>/<LSU>` state-string notation used throughout this package
(e.g. `A/A`, `A/P`, `ap/AP`, `pe/E`) is the standard convention from
two key papers:

- **Moazed, D. & Noller, H. F. (1989). Intermediate states in the
  movement of transfer RNA in the ribosome. *Nature*, 342, 142–148.**
  DOI: [10.1038/342142a0](https://doi.org/10.1038/342142a0). The
  foundational hybrid-states paper, introducing the
  `<SSU-site>/<LSU-site>` notation. Showed via chemical-protection
  footprinting that after peptide-bond formation but before
  translocation, deacylated and peptidyl tRNAs occupy hybrid `P/E`
  and `A/P` configurations — the acceptor end has swung into the
  next LSU site while the anticodon end remains on the SSU. All
  single-letter both-uppercase states (`A/A`, `P/P`, `E/E`, `A/P`,
  `P/E`) originate here.

- **Zhou, J., Lancaster, L., Donohue, J. P. & Noller, H. F. (2014).
  How the ribosome hands the A-site tRNA to the P site during
  EF-G-catalyzed translocation. *Science*, 345, 1188–1191.** DOI:
  [10.1126/science.1255030](https://doi.org/10.1126/science.1255030).
  Defined the **chimeric** hybrid states and the lowercase/uppercase
  case convention used by this package: **lowercase letters denote
  SSU (small-subunit) contacts; uppercase letters denote LSU
  (large-subunit) contacts**. A doubled letter on either side
  (`ap`, `pe`, `AP`, `PE`) means the tRNA simultaneously contacts
  two adjacent sites on that subunit. The fully chimeric `ap/AP`
  and `pe/PE` states correspond to the EF-G-trapped translocation
  intermediate (PDB **4W29** is the canonical structure of this
  state).

## Factor-bound tRNA states

The §12.4 protein-factor LSU label (a tRNA's CCA end held by a
non-ribosomal protein rather than docked into a canonical LSU site)
applies symmetrically to both A- and P-site tRNAs. Canonical
biological examples:

### A-site factor binding (elongation, termination)

- **EF-Tu·GTP·aa-tRNA decoding complex** — *Schmeing et al.*,
  *The crystal structure of the ribosome bound to EF-Tu and aminoacyl-tRNA*,
  *Science* **326**, 688–694 (2009). DOI:
  [10.1126/science.1179700](https://doi.org/10.1126/science.1179700).
  EF-Tu's domain II contacts the acceptor stem of the A-site
  aminoacyl-tRNA during codon decoding.
- **eEF1A·GTP·aa-tRNA** (eukaryotic homolog) — *Shao et al.*,
  *Decoding mammalian ribosome–mRNA states by translational GTPase
  complexes*, *Cell* **167**, 1229–1240 (2016). DOI:
  [10.1016/j.cell.2016.10.046](https://doi.org/10.1016/j.cell.2016.10.046).
- **Termination complexes** with release factors (RF1/RF2 in bacteria,
  eRF1·eRF3 in eukaryotes) contact the A-site tRNA region.

### P-site factor binding (initiation)

- **Initiation factor IF2 + initiator fMet-tRNAfMet** — *Sprink, Ramrath
  et al.*, *Structures of ribosome-bound initiation factor 2 reveal the
  mechanism of subunit association*, *Science Advances* **2**, e1501502
  (2016). DOI:
  [10.1126/sciadv.1501502](https://doi.org/10.1126/sciadv.1501502).
  Bacterial IF2's C-terminal domain (C2) latches onto the acceptor stem
  of the initiator tRNA at the P site. Canonical structure: PDB **1ZO1**
  (Allen et al. *Cell* 2005, E. coli 70S initiation complex).
- **Compact IF2 gating P-site accommodation** — *Basu, Yusupova et al.*,
  *Compact IF2 allows initiator tRNA accommodation into the P site and
  gates the ribosome to elongation*, *Nature Communications* **13**,
  3388 (2022). DOI:
  [10.1038/s41467-022-31129-2](https://doi.org/10.1038/s41467-022-31129-2).
- **eIF5B + Met-tRNAi (eukaryotic homolog)** — *Wang et al.*, *eIF5B
  and eIF1A reorient initiator tRNA to allow ribosomal subunit
  joining*, *Nature* **607**, 165–170 (2022). DOI:
  [10.1038/s41586-022-04858-z](https://doi.org/10.1038/s41586-022-04858-z).
  Domain IV of eIF5B forms extensive contacts with the methionylated
  3′ CCA end of the initiator tRNA, stabilising a non-canonical P/I
  orientation prior to GTP hydrolysis.
- **80S–eIF5B initiation→elongation transition** — *Lapointe et al.*,
  *Structural basis for the transition from translation initiation to
  elongation by an 80S-eIF5B complex*, *Nature Communications* **11**,
  4839 (2020). DOI:
  [10.1038/s41467-020-18829-3](https://doi.org/10.1038/s41467-020-18829-3).

## Ribosomal-protein nomenclature

- **Ban, N., Beckmann, R., Cate, J. H. D., Dinman, J. D., Dragon, F.,
  Ellis, S. R., Lafontaine, D. L. J., Lindahl, L., Liljas, A., Lipton,
  J. M., McAlear, M. A., Moore, P. B., Noller, H. F., Ortega, J.,
  Panse, V. G., Ramakrishnan, V., Spahn, C. M. T., Steitz, T. A.,
  Tchorzewski, M., Tollervey, D., Warren, A. J., Williamson, J. R.,
  Wilson, D., Yonath, A., & Yusupov, M. (2014). A new system for
  naming ribosomal proteins. *Current Opinion in Structural Biology*,
  24, 165–169.** DOI:
  [10.1016/j.sbi.2014.01.002](https://doi.org/10.1016/j.sbi.2014.01.002).
  Defines the unified **uS / uL / bS / bL / eS / eL / mS / mL**
  systematic naming used throughout this package's classification
  rules. `u` = universal (present in bacteria, archaea, eukaryotes),
  `b` = bacteria-specific, `e` = eukaryote-specific, `m` =
  mitochondria-specific. The short-form regex in
  `classify._BAN_NOMENCLATURE_REGEX` matches names following this
  convention and surfaces them as ribosomal proteins even when the
  RCSB description is just the bare token (e.g. `"bL28m"`) and no
  UniProt cross-reference is supplied.

## Curated reference ribosomes for contact-transfer

- **5J7L** — *Escherichia coli* 70S ribosome with mRNA and tRNAs;
  used as the bacterial / organellar reference. Noeske *et al.*
  (2015). *Nat. Struct. Mol. Biol.* 22, 336–341.
- **7ZW0** — *Saccharomyces cerevisiae* 80S ribosome; used as the
  eukaryotic cytoplasmic reference (chosen over the original spec's
  5TBW because 5TBW chain `A` is not a BGSU NR representative and
  5TBW is missing residue 2454 needed for the `lsu_etrna` anchors —
  see the v3.1 spec addendum for the full history).

## Mitoribosome structural biology

- **Amunts, A., Brown, A., Toots, J., Scheres, S. H. W., &
  Ramakrishnan, V. (2015). The structure of the human mitochondrial
  ribosome. *Science*, 348, 95–98.** DOI:
  [10.1126/science.aaa1193](https://doi.org/10.1126/science.aaa1193).
  Original report of PDB **3J9M** — the human 55S mitoribosome,
  showing that mt-tRNA-Val occupies the central protuberance as a
  structural surrogate for the 5S rRNA that mammalian mitoribosomes
  lack. Used as a regression case in this package because it
  exercises the organellar-classification path, the missing-5S case,
  and the "tRNA without a canonical contact site" case.

## Large-scale ribosome movements (RADdb)

- **Mears, J. A., Cannone, J. J., Stagg, S. M., Gutell, R. R.,
  Agrawal, R. K., & Harvey, S. C. (2002). Modeling a minimal ribosome
  based on comparative sequence analysis.** *Journal of Molecular
  Biology*, 321, 215–234. DOI:
  [10.1016/S0022-2836(02)00568-5](https://doi.org/10.1016/S0022-2836(02)00568-5).
  Foundational work for the rotation/tilt/translation decomposition
  later refined by RADtool.
- **Bonin, J. P., Aramini, J. M., Dunkle, J. A., & Harvey, S. C.
  RADtool / RADdb (https://radtool.rc.northeastern.edu/).** The
  package consumes RADdb's weekly LSUSSU CSV release, which decomposes
  every PDB ribosome assembly's inter-subunit geometry into body and
  head rotation/tilt/translation parameters. Two metrics are surfaced
  in v1: **`body rot.` → `intersubunit_rotation`** (the canonical
  ratchet-like inter-subunit rotation) and **`head rot.` →
  `ssu_head_rotation`** (the SSU-head swivel that gates translocation).
  Rows are joined to annotations by the `(pdb_id, lsu_chain_id,
  ssu_chain_id)` triple, which is unique per ribosome assembly in
  RADdb. See `raddb.py` for the local-cache / weekly-refresh policy.

## External APIs

- **RCSB Data API (GraphQL).** <https://data.rcsb.org/index.html#data-api>.
  Used for entry / assembly / polymer-entity-instance metadata.
- **PDBe REST API — `/nucleic_mappings/rfam/<pdb_id>`.**
  <https://www.ebi.ac.uk/pdbe/api/doc/>. Used to recover rRNA Rfam
  accessions, which the live RCSB GraphQL schema no longer ships for
  rRNA polymer entities.
- **BGSU RNA 3D Hub correspondence API —
  `/correspondence/map_across_chains`.**
  <https://rna.bgsu.edu/main/data-and-services/>. Used to project
  curated functional-site anchor nucleotides from a reference
  ribosome onto the equivalent positions in any query ribosome via
  Rfam covariance-model alignment.
- **RADdb LSU↔SSU CSV.**
  <https://radtool.rc.northeastern.edu/public_database/>. Weekly
  release of inter-subunit + SSU-head rotation/tilt/translation
  parameters for every public ribosome assembly; see the citation
  above for the underlying methodology.
- **Rfam.** <https://rfam.org/>. Source of the rRNA / tRNA family
  accessions consulted by `RFAM_ROLE_MAP` in `constants.py`.
