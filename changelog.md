# Changelog

## 0.2.5 — 2026-04-28

CI / release-plumbing only. No user-visible behaviour change.

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

## 0.2.4

(Pre-migration release. See git log for details.)
