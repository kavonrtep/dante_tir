# DANTE_TIR conda-release migration plan

Implementation plan for migrating `kavonrtep/dante_tir` from the centralised
`kavonrtep/recipes/dante_tir` workflow to in-repo, tag-driven GitHub-Actions
conda releases, following `conda_release_migration_playbook.md`.

This document is the working plan referenced in playbook §5. Keep it updated
as steps are completed.

---

## 1. Repo facts (playbook §0 answers)

| Question | Answer |
|---|---|
| Repo | `kavonrtep/dante_tir` |
| Conda package name | `dante_tir` |
| Anaconda channel | `petrnovak` |
| CLIs | `dante_tir.py` (supports `--version` / `--help`), `dante_tir_summary.R` (optparse `--help`, no `--version`) |
| Existing recipe location | `github.com/kavonrtep/recipes/tree/main/dante_tir` (top level — NOT under `recipes/`) |
| Current released version | `0.2.4` |
| Existing run-deps | `r-optparse`, `blast`, `cap3`, `r-rbeast`, `bioconductor-biostrings`, `bioconductor-bsgenome`, `bioconductor-rtracklayer`, `python>=3.6.0`, `mmseqs2>=14` |
| Existing `build.sh` symlinks | `dante_tir.py`, `dante_tir_summary.R` |
| Version source | `version.py` at root: `__version__ = '0.2.4'` (matches playbook default — regex template works as-is) |
| Tag style | 3-component, no `v` prefix (`0.2.4`, `0.2.3`, …) — playbook tag filter `['[0-9]+.[0-9]+.[0-9]+*']` works as-is |
| License | GPL-3.0 (file `LICENSE` at root) |
| Test data status | `test-data/` is huge (`pearl_v2.fa.split/scaffold_1.fa` ≈ 461 MB, `tiny_pea.fasta` ≈ 78 MB). **Cannot be pushed to GitHub or run in CI as-is.** Tiered subsets curated from `tiny_pea` (see §4). |

---

## 2. Target layout (after migration)

```
dante_tir/
├── .github/
│   └── workflows/
│       ├── tests.yml
│       └── conda-release.yml
├── conda/
│   └── dante_tir/
│       ├── meta.yaml
│       └── build.sh
├── tests/
│   ├── data/
│   │   ├── smoke/                       (~ 0.4 MB committed, gz)
│   │   │   ├── smoke.fasta              # Chr1[1..2 000 000] of tiny_pea
│   │   │   ├── smoke.gff3               # matching DANTE annotation
│   │   │   └── README.md                # provenance
│   │   ├── short/                       (~ 6.6 MB committed, gz)
│   │   │   ├── short.fasta              # Chr1[1..20 000 000] of tiny_pea
│   │   │   ├── short.gff3
│   │   │   └── README.md
│   │   └── long/                        (~ 10 MB committed, gz)
│   │       ├── long.fasta               # Chr1 (full, 30 014 000 bp) of tiny_pea
│   │       ├── long.gff3
│   │       └── README.md
│   ├── smoke.sh                          # ~ 15 s, structural only
│   ├── short.sh                          # ~ 90 s, dante_tir + summary
│   └── long.sh                           # ~ 115 s, release gate
├── tests.sh                              # dispatcher (rewritten)
├── version.py                            # already exists, keep
├── requirements.txt                      # NEW — conda-format run deps
├── changelog.md                          # NEW
├── LICENSE                               # already exists
├── README.md                             # already exists, install snippet bumped
├── dante_tir.py / detect_tirs.R / dante_tir_summary.R / dt_utils.{py,R}   # unchanged
├── dev_scripts/                          # NEW — old shell harnesses kept for HPC use
│   ├── tests2.sh
│   ├── tests3.sh
│   ├── tests2.py
│   └── run_parameter_tests.sh
└── docs/
    └── conda_release_migration_plan.md   # this file
```

Existing `tests.sh`, `tests2.sh`, `tests3.sh`, `tests2.py`, and
`run_parameter_tests.sh` are absorbed/replaced as follows:

- `tests.sh` (root) → rewritten as the `{smoke|short|long|all}` dispatcher.
- The current `tests.sh` body (scaffold_1, 461 MB) becomes the local-only
  long path under `tests/long.sh`, gated on the file actually being present.
- `tests2.sh` / `tests3.sh` reference `/mnt/raid/...`. They stay in the repo
  for the developer's HPC use **only if** they are moved to `dev_scripts/`
  or kept as-is — they don't run in CI. Recommendation: move into
  `dev_scripts/` and add a top-of-file note that they need external data.
- `tests2.py` — small ad-hoc benchmark; either deleted or moved into
  `dev_scripts/`. Decision before commit (default: move).

---

## 3. New / changed files — concrete drafts

### 3.1 `requirements.txt`

Conda-format dependency list that the workflow's
`mamba install --file requirements.txt` consumes. Must cover everything
the recipe `run:` block lists.

```
python >=3.6.0
mmseqs2 >=14
blast
cap3
r-optparse
r-rbeast
bioconductor-biostrings
bioconductor-bsgenome
bioconductor-rtracklayer
```

Note: the recipe does not currently list `r-base`, but
`bioconductor-*` will pull it in. We keep the same minimal set as the
existing recipe to avoid resolver surprises.

### 3.2 `conda/dante_tir/meta.yaml`

Playbook §3.1 verbatim, filled in. Run-deps mirror the existing recipe
with one additional sanity item (`r-base` left implicit). Test commands
cover both CLIs.

```yaml
{% set m = load_file_regex(load_file='version.py',
                           regex_pattern="__version__\\s*=\\s*['\"]([^'\"]+)['\"]",
                           from_recipe_dir=True) %}
{% set version = m.group(1) if m else '0.0.0.dev' %}

package:
  name: dante_tir
  version: {{ version }}

source:
  path: ../..

build:
  number: 0
  noarch: generic

requirements:
  run:
    - python >=3.6.0
    - mmseqs2 >=14
    - blast
    - cap3
    - r-optparse
    - r-rbeast
    - bioconductor-biostrings
    - bioconductor-bsgenome
    - bioconductor-rtracklayer

test:
  commands:
    - dante_tir.py --version
    - dante_tir.py --help | grep -qi 'usage'
    - dante_tir_summary.R --help | grep -qi 'usage'

about:
  home: https://github.com/kavonrtep/dante_tir
  dev_url: https://github.com/kavonrtep/dante_tir
  license: GPL-3.0-only
  license_family: GPL3
  license_file: LICENSE
  summary: >
    DANTE_TIR identifies full-length DNA transposons with Terminal Inverted
    Repeats (TIRs), starting from DANTE annotations of conserved transposase
    domains.

extra:
  copy_test_source_files: true
  final: true
  recipe-maintainers:
    - petrnovak
```

### 3.3 `conda/dante_tir/build.sh`

```sh
#!/bin/sh
set -x -e

PKG_DIR="${PREFIX}/share/dante_tir"
mkdir -p "${PREFIX}/bin" "${PKG_DIR}"
cp -r . "${PKG_DIR}"

ln -s "${PKG_DIR}/dante_tir.py"          "${PREFIX}/bin/dante_tir.py"
ln -s "${PKG_DIR}/dante_tir_summary.R"   "${PREFIX}/bin/dante_tir_summary.R"
```

Note: `detect_tirs.R` is *not* on PATH in the existing recipe — it is
called via `Rscript` from `dante_tir.py` using a path computed relative
to the python script. That continues to work because `cp -r .` puts it
next to `dante_tir.py` inside `$PKG_DIR`. **Verify** at local-build time
that `dante_tir.py` resolves `detect_tirs.R` correctly.

### 3.4 `.github/workflows/tests.yml`

Playbook §3.3 with `<PKG>=dante_tir`, `<USER>=petrnovak`. Runs smoke +
short on push and PR to `main`.

### 3.5 `.github/workflows/conda-release.yml`

Playbook §3.4 verbatim with the same substitutions. **Do not** modify
the `setuptools<81` pin, the `Sync version.py into recipe dir` step, or
the `load_file_regex` (playbook §7 pitfalls 1, 3, 4).

### 3.6 `tests.sh` (new dispatcher)

Playbook §3.5 verbatim with `<PKG>=dante_tir`. Backwards-compat: an
integer first arg routes to `long` with that NCPU.

### 3.7 `tests/smoke.sh` — structural, ~ 15 s

The smoke slice (Chr1[1..2 Mb] of tiny_pea, 0.42 MB committed gz)
deliberately produces **zero detected TIRs** — there aren't enough
Subclass_1 copies in 2 Mb to cross the assembly cutoff. Smoke verifies
that the python + R + CAP3 + mmseqs2 + BLAST integration runs cleanly
and the entry-point CLIs work. Functional assertions move to short/long.

```bash
#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DATA="$ROOT/tests/data/smoke"
OUT="$ROOT/tmp/tests/smoke"
NCPU="${NCPU:-2}"

rm -rf "$OUT"; mkdir -p "$OUT"; cd "$ROOT"

echo "=== CLI checks ==="
./dante_tir.py --version
./dante_tir.py --help > /dev/null
./dante_tir_summary.R --help > /dev/null

echo "=== dante_tir on smoke subset (~ 15 s, expect 0 TIRs detected) ==="
./dante_tir.py -g "$DATA/smoke.gff3" -f "$DATA/smoke.fasta" \
               -o "$OUT/out" -c "$NCPU"

# Structural: pipeline must exit 0 and create the final GFF/FASTA
# files, even if they are empty (no TIRs detected on a 2 Mb slice).
[ -f "$OUT/out/DANTE_TIR_final.gff3" ]  || { echo "FAIL: DANTE_TIR_final.gff3 missing"; exit 1; }
[ -f "$OUT/out/DANTE_TIR_final.fasta" ] || { echo "FAIL: DANTE_TIR_final.fasta missing"; exit 1; }
[ -d "$OUT/out/log" ]                   || { echo "FAIL: log/ dir missing"; exit 1; }

echo "smoke PASSED"
```

### 3.8 `tests/short.sh` — functional, ~ 90 s

Short slice: Chr1[1..20 Mb] of tiny_pea (~ 6.6 MB committed gz).
Benchmark on 4 CPUs: dante_tir 67 s + summary 20 s = ~ 87 s; expect
≥ 14 TIR records and a non-empty `report.html`.

```bash
#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DATA="$ROOT/tests/data/short"
OUT="$ROOT/tmp/tests/short"
NCPU="${NCPU:-2}"

rm -rf "$OUT"; mkdir -p "$OUT"; cd "$ROOT"

echo "=== dante_tir on short subset (~ 70 s) ==="
./dante_tir.py -g "$DATA/short.gff3" -f "$DATA/short.fasta" \
               -o "$OUT/out" -c "$NCPU"

GFF="$OUT/out/DANTE_TIR_final.gff3"
N=$(grep -c -v '^#' "$GFF" || true)
[ "$N" -ge 5 ] || { echo "FAIL: too few TIR records in $GFF: $N"; exit 1; }

echo "=== dante_tir_summary.R on short output (~ 20 s) ==="
./dante_tir_summary.R -g "$GFF" -f "$DATA/short.fasta" \
                      -o "$OUT/summary" -t "$NCPU"

[ -s "$OUT/summary/report.html" ] || { echo "FAIL: report.html missing"; exit 1; }
[ -s "$OUT/summary/all_representative_elements_min3.fasta" ] || \
  { echo "FAIL: all_representative_elements_min3.fasta missing/empty"; exit 1; }

echo "short PASSED ($N TIR records)"
```

### 3.9 `tests/long.sh` — release gate, ~ 115 s

Long slice: full Chr1 of tiny_pea, 30 Mb (~ 10 MB committed gz).
Benchmark on 4 CPUs: dante_tir 85 s + summary 28 s = ~ 113 s; expect
≥ 29 TIR records.

Two modes:

1. **CI / default.** Use the committed `tests/data/long/long.{fasta,gff3}`.
2. **Local deep mode.** If `LONG_FASTA` and `LONG_DANTE` are exported and
   point to existing files, use them instead. Keeps the spirit of the
   old `tests2.sh` / `tests3.sh` (HPC `/mnt/raid/...` runs) without
   breaking CI (playbook §7 pitfall 5).

```bash
#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DEFAULT_FA="$ROOT/tests/data/long/long.fasta"
DEFAULT_GFF="$ROOT/tests/data/long/long.gff3"
LONG_FASTA="${LONG_FASTA:-$DEFAULT_FA}"
LONG_DANTE="${LONG_DANTE:-$DEFAULT_GFF}"
OUT="$ROOT/tmp/tests/long"
NCPU="${NCPU:-2}"

[ -s "$LONG_FASTA" ] || { echo "FAIL: $LONG_FASTA missing"; exit 1; }
[ -s "$LONG_DANTE" ] || { echo "FAIL: $LONG_DANTE missing"; exit 1; }

rm -rf "$OUT"; mkdir -p "$OUT"; cd "$ROOT"

echo "=== dante_tir on long input ==="
./dante_tir.py -g "$LONG_DANTE" -f "$LONG_FASTA" \
               -o "$OUT/out" -c "$NCPU"

GFF="$OUT/out/DANTE_TIR_final.gff3"
N=$(grep -c -v '^#' "$GFF" || true)
# committed long input gives 29 TIR records on tiny_pea Chr1; allow drift
[ "$N" -ge 10 ] || { echo "FAIL: too few TIR records in $GFF: $N"; exit 1; }

echo "=== dante_tir_summary.R on long output ==="
./dante_tir_summary.R -g "$GFF" -f "$LONG_FASTA" \
                      -o "$OUT/summary" -t "$NCPU"

[ -s "$OUT/summary/report.html" ] || { echo "FAIL: report.html missing"; exit 1; }
[ -s "$OUT/summary/all_representative_elements_min3.fasta" ] || \
  { echo "FAIL: all_representative_elements_min3.fasta missing/empty"; exit 1; }

echo "long PASSED ($N TIR records)"
```

### 3.10 `.gitignore`

Add:

```
/conda/dante_tir/version.py
/tmp/
/__pycache__/
*.pyc
/test2/
/old_stuff/
/.idea/
/CLAUDE.md
/MAX_CLASS_SIZE_PLAN.md
/conda_release_migration_playbook.md
/hermit/
```

(Audit list against `git status` before committing — keep only the
entries that actually correspond to local-only files. Do **not** add
`tests/data/` here.)

### 3.11 `changelog.md`

Seed file:

```
## 0.2.5 — <date>
- Add in-repo conda recipe under `conda/dante_tir/`.
- Add tag-driven GitHub Actions release pipeline.
- Refactor `tests.sh` into `{smoke|short|long|all}` dispatcher.
- Commit a tiny smoke dataset under `tests/data/smoke/`.

## 0.2.4
- (existing release — pre-migration)
```

### 3.12 `README.md` install snippet

Update the existing line:

```
conda install -c conda-forge -c r -c bioconda -c petrnovak  dante_tir
```

— no change strictly needed; the channel hasn't moved. Add one paragraph
after it noting that the package is now built from this repository's
tags via GitHub Actions, and that anaconda.org/petrnovak/dante_tir is
canonical going forward.

---

## 4. Test datasets — sources, sizes, and benchmark

All three slices are derived from `test-data/tiny_pea.fasta` +
`test-data/DANTE_tiny_pea.gff3` by taking the first N bp of Chr1 and
filtering the DANTE GFF3 to records fully within that range. No
coordinate remapping is needed because the slice keeps the original
chromosome name and starts at position 1.

Generator: `tmp/smoke_dataset_eval/make_slice.py` (gitignored). Full
benchmark + reasoning: `tmp/smoke_dataset_eval/REPORT.md`.

### Detection threshold (the reason smoke is structural-only)

dante_tir filters DANTE for `Subclass_1` TIR domains and assembles
≥ 4 copies per superfamily before R/BEAST detection. Slices below
~ 15 Mb of contiguous Chr1 don't reach that threshold for any
superfamily. Sharp transition between 15 Mb (0 TIRs) and 20 Mb
(14 TIRs).

Carving a denser non-contiguous region works but reorganises
coordinates and is fragile across DANTE re-annotations. Decision:
keep smoke structural, push functional assertions into short/long.

### Final selection (4 CPU on host with SSD)

| tier  | slice                          | dante_tir | summary | total | size committed (gz) | TIRs detected | assertion                                                           |
|-------|--------------------------------|----------:|--------:|------:|--------------------:|--------------:|---------------------------------------------------------------------|
| smoke | Chr1[1..2 000 000]             |     14 s  |   skip  | ~15 s | ~ 0.42 MB           |             0 | structural — exit 0; final files exist; log/ dir exists             |
| short | Chr1[1..20 000 000]            |     67 s  |   20 s  | ~87 s | ~ 6.6 MB            |            14 | `≥ 5` GFF records; summary `report.html` non-empty; combined `min3.fasta` non-empty |
| long  | Chr1 (full, 30 014 000 bp)     |     85 s  |   28 s  | ~113 s | ~ 10 MB            |            29 | `≥ 10` GFF records; summary as above                                |

50 % of tiny_pea (40 Mb, Chr1 + 10 Mb of Chr2) was tested too — works
fine but adds ~ 3 MB of committed data for only +6 TIRs over full Chr1.
Not worth it.

---

## 5. Step-by-step execution checklist

### Pre-work
- [ ] Read `conda_release_migration_playbook.md` end-to-end (already done).
- [ ] Confirm anaconda channel still `petrnovak` (yes).
- [ ] Confirm `version.py` regex works (yes — already in playbook default form).
- [ ] Inventory CLIs against existing `build.sh`: `dante_tir.py`,
      `dante_tir_summary.R`. Match.

### Scaffolding (one commit per logical unit)
- [ ] Commit 1 — *recipe*: add `conda/dante_tir/meta.yaml`,
      `conda/dante_tir/build.sh` (chmod +x), and `requirements.txt`.
- [ ] Commit 2 — *smoke data*: add `tests/data/smoke/{smoke.fasta,smoke.gff3}`
      (and a small `README.md` explaining provenance).
- [ ] Commit 3 — *test scaffolding*: add `tests/smoke.sh`,
      `tests/short.sh`, `tests/long.sh` (chmod +x); replace root
      `tests.sh` with the dispatcher; move `tests2.sh`, `tests3.sh`,
      `tests2.py`, `run_parameter_tests.sh` to `dev_scripts/` (or leave
      in place with documentation — decide before commit).
- [ ] Commit 4 — *workflows*: add `.github/workflows/tests.yml` and
      `.github/workflows/conda-release.yml`.
- [ ] Commit 5 — *housekeeping*: update `.gitignore`, add `changelog.md`,
      bump README install paragraph.

Per playbook §5, do NOT bundle these into a single blob commit.

### Local verification (mandatory — playbook §4 + §5)
- [ ] `mamba create -n relcheck -c conda-forge python=3.11 conda-build boa
      anaconda-client "setuptools<81" -y && conda activate relcheck`.
- [ ] `anaconda --version` and `anaconda upload --help | head -1` both
      print without `pkg_resources` traceback.
- [ ] Local meta.yaml render check (small Python script in `tmp/` that
      replicates `load_file_regex` and confirms the rendered version is
      `0.2.5`, NOT `None`, NOT `0.0.0.dev`).
- [ ] `./tests.sh smoke` exits 0 in < 30 s.
- [ ] `./tests.sh short` exits 0 in < 3 min.
- [ ] `cp version.py conda/dante_tir/ && conda mambabuild -c conda-forge
      -c bioconda -c petrnovak --output-folder /tmp/bld conda/dante_tir`
      produces `dante_tir-0.2.5-0.tar.bz2` (NOT `dante_tir-None-0.tar.bz2`,
      NOT `dante_tir-0.0.0.dev-0.tar.bz2`).
- [ ] Optional — install the freshly built `.tar.bz2` into a scratch
      env and run `dante_tir.py --version` plus `dante_tir_summary.R --help`
      to verify both CLIs end up on PATH and `detect_tirs.R` resolves.

### Secrets + first push
- [ ] User creates anaconda token in anaconda.org → Settings → Access:
      scopes "Allow write access to the API site" + "Allow all operations
      on Conda repositories", name `dante_tir_github_release`, expiration
      1 year.
- [ ] User adds it as repo secret `ANACONDA_API_TOKEN` at
      github.com/kavonrtep/dante_tir/settings/secrets/actions.
- [ ] Push commits 1–5 to `main`. Wait for `tests.yml` to go green.
- [ ] Bump `version.py` to `0.2.5`, append changelog entry, commit.
- [ ] Tag: `git tag -a 0.2.5 -m "release 0.2.5"`.
- [ ] Push tag: `git push origin 0.2.5` — fires `conda-release.yml`.
- [ ] Watch the Actions log:
      - tag↔version.py assertion passes
      - long test passes
      - `conda mambabuild` produces `dante_tir-0.2.5-0.tar.bz2`
      - `anaconda upload` succeeds (no `pkg_resources` traceback)
- [ ] Verify https://anaconda.org/petrnovak/dante_tir shows `0.2.5`.

### Cleanup in `kavonrtep/recipes`
- [ ] After successful first release, replace
      `kavonrtep/recipes/dante_tir/meta.yaml` and `…/build.sh` with a
      short `README.md` redirecting to
      `github.com/kavonrtep/dante_tir/tree/main/conda/dante_tir`, **or**
      delete the `dante_tir/` folder entirely. Commit in the recipes
      repo.

### Stability watch
- [ ] After three consecutive successful tag-driven releases (e.g.
      0.2.5, 0.2.6, 0.2.7 or whatever sequence), declare stable and
      remove the external recipe directory if it was only redirected.

---

## 6. Risks and gotchas specific to this repo

1. **R-package resolution time.** This recipe has the heaviest deps of
   the kavonrtep recipes (Bioconductor + Rbeast). Channel order
   `conda-forge,bioconda,petrnovak` with `channel-priority: strict` is
   essential — playbook §7 pitfall 8.
2. **`detect_tirs.R` is invoked as a sibling file, not a CLI.** The
   build.sh `cp -r .` ensures it lands in `$PKG_DIR` next to
   `dante_tir.py`. Verify the python invokes it via a relative-to-script
   path (it does — `dt_utils.py` uses `os.path.dirname(__file__)`).
   Fastest verification: install the locally-built tarball and run
   the smoke test inside that env.
3. **Test data scale.** Existing `test-data/` is dev-only (461 MB
   scaffold, etc.). The whole `test-data/` must NOT be added to git in
   this migration if it isn't already tracked. Verify with
   `git ls-files test-data/ | head` before committing — only the curated
   `tests/data/smoke/` (and possibly `tests/data/short/`) goes into the
   repo.
4. **Tag style.** All historical tags are 3-component (`0.2.4`). The
   playbook tag filter `['[0-9]+.[0-9]+.[0-9]+*']` already matches; do
   not switch to a `v`-prefix style.
5. **`dante_tir_summary.R --help`.** optparse prints `Usage: …` on
   `--help`, which satisfies the `grep -qi 'usage'` test command.
   Confirm at local-build time before pushing.
6. **First post-migration version.** Bump from `0.2.4` to `0.2.5` (not
   `0.3.0`) — this is a CI/release-plumbing change, not a feature
   release. Any user-visible changes can ride the next normal version.

---

## 7. Decisions (resolved)

1. **Test datasets:** structural smoke at 2 Mb, functional short at
   20 Mb, long (release gate) at full Chr1 = 30 Mb. All three are
   contiguous Chr1[1..N] slices of `tiny_pea`. See §4.
2. **First release version:** `0.2.5` — CI/release-plumbing only, no
   user-visible behaviour change.
3. **Old shell harnesses:** move `tests2.sh`, `tests3.sh`, `tests2.py`,
   `run_parameter_tests.sh` to `dev_scripts/` (kept for HPC use; they
   reference `/mnt/raid/...`, not run in CI).
4. **Smoke-data licensing:** `tiny_pea` slices are derived from the
   committed `test-data/` already in the dev tree; assumed OK to
   redistribute as part of GPL-3 source. **Re-confirm with user before
   first push if `test-data/` has any redistribution restriction.**
