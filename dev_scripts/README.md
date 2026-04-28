# dev_scripts/

Developer-only harnesses kept for HPC use. **Not run by CI.**

| script | what it does | requires |
|---|---|---|
| `tests2.sh` | runs `dante_tir.py` on a JunEffu_240903 DANTE GFF + genome | `/mnt/raid/454_data/DToL/Henderson_paper/...` |
| `tests3.sh` | another large-genome run | external paths under `/mnt/raid/...` |
| `run_parameter_tests.sh` | sweeps `n_beast_iter` × `max_class_size` across 4 genomes | external `test-data/` tree |
| `tests2.py` | small ad-hoc benchmark | none |

These were previously at the repo root. They were moved here when
`tests.sh` became the `{smoke|short|long|all}` dispatcher in 0.2.5.
