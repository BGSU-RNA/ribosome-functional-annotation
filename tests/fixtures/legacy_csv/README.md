# Legacy CSV regression corpus

Real outputs from the original prototype (`process_annotation.py` +
`ribosome_annotation_new.py`), checked in as the de-facto scale regression
contract for the rewrite (spec §22.3).

- `ribosome_chain_annotation.csv` — 1,104 rows, one per (pdb_id, assembly_id).
  Columns:
  `pdb_id,assembly_id,ssu_chain,lsu_large_chain,lsu_medium_chain,lsu_small_chain,mrna,aminoacyl_trna,peptidyl_trna,exit_trna,aminoacyl_trna_state,peptidyl_trna_state,exit_trna_state`.
- `ribosome_assembly_annotation.csv` — 6,043 rows, one per (property, value)
  pair. Columns: `pdb_id,assembly_id,chain,property,value`.

Chain identifiers in both files are BGSU-style IFE strings
(`<pdb_id>|1|<auth_asym_id>`, e.g. `5J7L|1|AA`). Do not edit these files;
they are the fixed comparison target for `tests/test_legacy_csv_regression.py`
(implemented in step 11 of the build plan).
