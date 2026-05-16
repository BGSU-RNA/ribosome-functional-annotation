# Acceptance test PDB selection (spec §25.1)

The spec calls out three slots that require picking representative PDB
entries: isolated SSU (case 5), isolated LSU (case 6), and archaeal
(case 7). The picks below are tentative — verify against the live API
on the first network run and adjust if any entry has been redeposited or
no longer matches the expected shape.

| Case | PDB ID | Reason for pick |
|------|--------|-----------------|
| 1 | `5J7L` | Canonical *E. coli* 70S with mRNA + A/P/E tRNAs. The reference structure for the BGSU bacterial anchor set (§6.3). |
| 2 | `5TBW` | Reference yeast 80S. Also the eukaryotic anchor set. |
| 3 | `6ZMI` | Human cytoplasmic 80S — exercises multi-character chain IDs (`S2`, `L5`). |
| 4 | `6ZM6` | Human mitochondrial ribosome — bacterial-like rRNA core + Eukaryota protein superkingdom → `eukaryotic_organellar_ribosome`. (Alternative: `6ZM5`.) |
| 5 | `1J5E` | *Thermus thermophilus* 30S subunit, bacterial, deposited as 30S-only — exercises the `partial_ribosome_missing_ssu_or_lsu` skip path with a single SSU. |
| 6 | `1NKW` | *Deinococcus radiodurans* 50S subunit, bacterial, 50S-only deposit. Alternative: `4WSM`. |
| 7 | `4V6U` | *Pyrococcus furiosus* 70S, archaeal. Carries archaeal SSU (Rfam `RF01959`) + LSU (`RF02540`); both are in `RFAM_ROLE_MAP` so SSU/LSU detection succeeds and the §8.3 vote correctly identifies the Archaea taxonomy, triggering `archaeal_ribosome_not_supported`. (Note: `8HKV` looks archaeal but is actually 50S-only — it would trip the partial-skip first.) |
| 8 | `5FDV` | Multi-assembly bacterial — assembly 1 and assembly 2 each carry a separate ribosome. Should produce two annotated `RibosomeAnnotation` objects. |

If any of cases 5/6/7 turn out not to match (e.g. RCSB redeposits with a
combined assembly), pick a close substitute, update the table, and
update the test's PDB ID — the assertion shape stays the same.
