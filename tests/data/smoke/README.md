# smoke dataset

`smoke.fasta` is `Chr1[1..2 000 000]` of `test-data/tiny_pea.fasta`
(2 000 000 bp, single record kept under its original `Chr1` name).

`smoke.gff3` contains every DANTE record from
`test-data/DANTE_tiny_pea.gff3` whose `Chr1` coordinates fall fully
inside that interval (400 records, 5 of which are `Subclass_1` TIR
domains).

`tests/smoke.sh` runs in ~ 15 s on 4 CPUs. Five Subclass_1 records is
below dante_tir's per-superfamily assembly threshold, so smoke
deliberately produces **zero detected TIRs** — its job is structural:
verify the python + R + CAP3 + mmseqs2 + BLAST integration runs cleanly
and that both entry-point CLIs (`dante_tir.py`, `dante_tir_summary.R`)
work.

To regenerate from a fresh `test-data/`, run
`tmp/smoke_dataset_eval/make_slice.py 2000000 tmp/.../slice_2000000`
and copy `slice.{fasta,gff3}` over `smoke.{fasta,gff3}`.
