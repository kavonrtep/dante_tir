# Changelog

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
