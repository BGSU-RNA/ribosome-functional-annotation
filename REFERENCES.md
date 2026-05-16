# References

Citations for the biological data sources, classification rules, and
naming conventions the package relies on.

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
- **Rfam.** <https://rfam.org/>. Source of the rRNA / tRNA family
  accessions consulted by `RFAM_ROLE_MAP` in `constants.py`.
