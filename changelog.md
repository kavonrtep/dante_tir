# Changelog

## Unreleased

Round-3 scaling: the self-BLAST table is never parsed in R at all.

- The `awk` reduction added in 0.2.7 kept the Round-3 table off the heap but
  could not shrink it enough: on an 89 Gb genome with 168k EnSpm/CACTA domains
  the *filtered* table was still 41.5 GB / 2.5e9 rows and `read.table` aborted
  with `long vectors not supported yet`. The row count is intrinsic to an
  all-vs-all self-BLAST, so no per-row filter can fix it.
- `run_blast_tir_analysis` now streams `blastn` straight into `blast_cp.py`,
  which applies `filter_blast3`'s predicates, folds each surviving hit into a
  per-subject coverage profile and computes the switch points, returning only
  the `id`/`cp` table Round 3 actually consumes. Nothing proportional to the
  hit count is ever stored: memory is a dense int32 profile array
  (`n_subjects x max_length`, 4.2 GB for the class above) and the BLAST stream
  itself is discarded as it is read. Round-3 outputs are now
  `*_upstream3.cp.tsv` / `*_downstream3.cp.tsv`.
- Ask `blastn` for only the six columns the predicates need (it was formatting
  twelve and `awk` discarded nine) and push the identity cut-off down into
  `-perc_identity`. `filter_blast3`'s evalue predicate was already unreachable
  behind `-evalue 1e-10` and is no longer evaluated.
- Results are unchanged: `tests/test_round3_cp.R` (replacing
  `tests/test_blast_reduce.R`) asserts the new path reproduces the R
  implementation's switch points exactly, across predicate edge cases, random
  coverage profiles for both switch-point methods, and a real self-BLAST.
- New runtime dependency: `numpy`.
- New `--max_round3_queries N` (default 0 = off, i.e. the exact path). Above
  that many copies, a class runs its Round-3 self-BLAST on a random N-query
  sample with the coverage scaled back up, and every subject whose sampled
  support falls below `min_support` (default 30) is then re-resolved exactly
  against the full query set in a cheap second pass over just those subjects.
  BLAST cost for the CACTA class is linear in query count (measured), so at
  the recommended `--max_round3_queries 20000` that class costs ~23% of the
  exact run (~2.3 h instead of ~10 h). Sampling is not free: measured on
  MuDR/Mutator at a comparable sampling fraction, switch points move by <=10 bp
  for 96% of subjects but ~4% change between called and not-called. Numbers,
  calibration and the known limitation of the support gate are in
  `docs/round3_scaling.md` 11.
- Sampled runs are deterministic: the draw is a function of `--seed` and the
  input alone. The generator is pinned explicitly rather than inherited, which
  closes a real hole -- under `RNGkind("L'Ecuyer-CMRG")`, which `parallel` code
  commonly sets, `set.seed()` draws a different sample -- and the caller's RNG
  stream is left untouched. Each sampled class logs its seed and a fingerprint
  of the draw to `log/stderr.txt` (which survives without `--debug`, unlike the
  sampled query FASTA). Thread count does not affect results either: the
  coverage profile is an order-independent integer accumulation, verified
  byte-identical at `-num_threads 1` vs `8` on real data.
- Stop loading the whole genome into R. `detect_tirs.R` and
  `cluster_tir_sequences` both did `readDNAStringSet(genome)` to pull out a few
  thousand short TIR ranges; a `DNAStringSet` costs ~1 byte per base, so for
  run-000129's 94.3 Gbp assembly that is ~94 GB of RAM, twice. Neither call had
  ever been reached on that genome because Round 3 failed first. Both now use
  `Rsamtools::FaFile` random access through the `.fai`, which reads only the
  requested ranges (`genome_fa_handle`, which indexes the genome if needed and
  fails loudly if it cannot). Output is unchanged: `tests/test_genome_access.R`
  asserts the two paths return identical sequences including
  reverse-complementing of minus-strand ranges, and `tests.sh short` produces a
  byte-identical `DANTE_TIR_final.fasta`. `bioconductor-rsamtools` is now a
  declared dependency (it arrived via BSgenome before).
- Known issue documented, not changed: `find_switch_point_from_blast_coverage3`
  (used for EnSpm/CACTA) reads its QC thresholds 200 bp before the switch point
  it detected. `blast_cp.py` reproduces the behaviour deliberately; see
  `docs/round3_scaling.md` 8.1.

- Fix a hang on large inputs: `cluster_aa_sequences_mmseqs2` and
  `make_blast_db` ran their child processes with `subprocess.check_call(...,
  stdout=PIPE, stderr=PIPE)`. `check_call` waits without ever reading those
  pipes, so the child blocks forever once it writes past the 64 KB pipe buffer
  -- `mmseqs easy-cluster` reached that on real data and deadlocked the whole
  run in `pipe_write`. mmseqs now goes through `subprocess.run` (which drains
  while waiting) and reports its stderr when clustering fails; `makeblastdb`
  discards its output via `DEVNULL`. Added `tests/test_subprocess_pipes.py`
  (in `tests/unit.sh`) to keep the pattern from coming back.

## 0.2.8 — 2026-07-21

Round-3 memory fix plus a release guard.

- Round 3 no longer reads the full self-BLAST table into R. On large,
  high-copy genomes the all-vs-all self-BLAST can exceed hundreds of GB, and
  `read.table` aborted with `long vectors not supported yet` (R's 2^31-element
  limit) before any element was called. `run_blast_tir_analysis` now streams
  `blastn` through an `awk` filter (`blast3_reduce_awk`) that reproduces
  `filter_blast3`'s predicates and keeps only the `saccver/sstart/send` columns
  coverage needs, so the full table is never written to disk nor loaded into R.
  Results are unchanged — verified identical at the filter, coverage, switch-
  point, and final-GFF3 level. Round-3 BLAST outputs are now named
  `*_upstream3.filtered.tsv` / `*_downstream3.filtered.tsv`.
- Add `tests/test_blast_reduce.R` (run from `tests/smoke.sh`) proving the
  streamed reduction is identical to the previous read.table/filter path across
  predicate edge cases and a real self-BLAST.
- Release plumbing: add `dev_scripts/check_release_version.sh`, a guard that
  refuses to (re)release a version already published on the conda channel (and,
  with `--require-untagged`, one that already has a local git tag). Wired into
  `conda-release.yml` as a fail-fast preflight so a duplicate version bump can
  no longer waste a full build and then fail at upload.

## 0.2.7 — 2026-07-15

Memory-usage fix plus CI plumbing.

- `extract_flanking_regions` no longer loads the whole genome into RAM
  (previously ~genome_bp, causing OOM on large assemblies such as a 90 Gbp
  genome). It now takes the FASTA path and streams one sequence at a time
  via a new `fasta_record_generator`, so peak memory is the largest single
  sequence rather than the whole assembly. Output is byte-identical to the
  previous dict-based logic.
- Add `tests/test_extract_flanking_regions.py` (plus `tests/unit.sh` and a
  `unit` level in `tests.sh`) verifying the streaming output matches the
  original logic across both strands, window clamping, and line-wrapped
  multi-sequence FASTA.
- CI: switch the conda-release workflow from the removed `conda mambabuild`
  to the standalone `conda-build` executable.

## 0.2.6 — 2026-06-30

- `dante_tir_summary.R`: reorder the TIR consensus logo figure in the HTML
  report so the 5' TIR is shown first (top panel) and the 3' TIR second
  (bottom panel).

## 0.2.5 — 2026-04-28

CI / release-plumbing migration with a small batch of robustness fixes.

CI / release plumbing:

- Add in-repo conda recipe under `conda/dante_tir/`.
- Add tag-driven GitHub Actions release pipeline under
  `.github/workflows/conda-release.yml`; tests on every push/PR via
  `.github/workflows/tests.yml`.
- Refactor `tests.sh` into a `{smoke|short|long|all}` dispatcher; add
  tiered scripts under `tests/`.
- Commit tiered datasets under `tests/data/{smoke,short,long}/` derived
  from the first 2 / 20 / 30 Mb of `tiny_pea` Chr1.
- Move old developer harnesses (`tests2.sh`, `tests3.sh`, `tests2.py`,
  `run_parameter_tests.sh`) under `dev_scripts/` — they reference HPC
  paths and are not run in CI.
- Add `requirements.txt` for use by both the recipe and the CI workflows.

Robustness:

- When no TIRs are detected, write valid empty `DANTE_TIR_final.gff3`
  (with `##gff-version 3` + `##DANTE_TIR version` banner) and
  `DANTE_TIR_final.fasta` instead of leaving the output dir bare. Print
  a final summary line so users can tell apart "ran cleanly, found 0"
  from "wrote no output / crashed".
- `dante_tir_summary.R`: capture per-class processing errors and emit
  a "Per-class summary unavailable: <reason>" placeholder in the HTML
  report instead of blowing up the whole report.
- `dt_utils.R`: keep the consensus matrix in matrix shape via
  `drop = FALSE`, return an empty `DNAStringSet` when every column has
  zero counts.

## 0.2.4

(Pre-migration release. See git log for details.)
