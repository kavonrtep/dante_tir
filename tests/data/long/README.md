# long dataset

`long.fasta` is the full Chr1 of `test-data/tiny_pea.fasta`
(30 014 000 bp, single record under the original `Chr1` name).

`long.gff3` contains every DANTE record on Chr1 (10 162 records, 332
`Subclass_1`).

`tests/long.sh` on 4 CPUs:

- `dante_tir.py` ~ 85 s, produces `DANTE_TIR_final.gff3` with **29 TIR
  records** at the time the dataset was committed.
- `dante_tir_summary.R` ~ 28 s, produces the full summary tree.

The `≥ 10` GFF-record assertion in `tests/long.sh` is set well below 29
to tolerate cross-version drift.

Local deep mode: export `LONG_FASTA` and `LONG_DANTE` to point at
external (e.g. HPC `/mnt/raid/...`) paths to override the committed
inputs without touching the script.

Regeneration: see `tests/data/smoke/README.md`; substitute `30000000`
for the slice length.
