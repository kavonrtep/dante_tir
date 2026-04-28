# Changelog

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
