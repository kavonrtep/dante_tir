# short dataset

`short.fasta` is `Chr1[1..20 000 000]` of `test-data/tiny_pea.fasta`
(20 000 000 bp, single record under the original `Chr1` name).

`short.gff3` contains every DANTE record fully inside that interval
(6 702 records, 220 `Subclass_1`).

`tests/short.sh` on 4 CPUs:

- `dante_tir.py` ~ 67 s, produces `DANTE_TIR_final.gff3` with **14 TIR
  records** at the time the dataset was committed.
- `dante_tir_summary.R` ~ 20 s, produces `report.html`,
  `all_representative_elements_min3.fasta`, per-superfamily fasta/csv,
  and PNG plots under `img/`.

The `≥ 5` GFF-record assertion in `tests/short.sh` is set well below 14
to tolerate cross-version drift in detection sensitivity.

Regeneration: see `tests/data/smoke/README.md`; substitute `20000000`
for the slice length.
