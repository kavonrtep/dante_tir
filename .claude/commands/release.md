---
description: Cut a dante_tir release (guard, bump version, changelog, commit, tag) — stops before push
argument-hint: [X.Y.Z]
---

Prepare a new dante_tir release. The target version is `$ARGUMENTS` (if empty,
propose the next patch version after the latest git tag and confirm it).

Follow this sequence exactly. Tools (Rscript, blastn, etc.) are NOT on the base
PATH — prepend the conda env when you need them:
`export PATH="$(pwd)/hermit/envs/conda/envs/dante_tir_test/bin:$PATH"`.

1. **Pick the version.** Use `$ARGUMENTS` if given; otherwise take the highest
   `git tag`, bump the patch component, and state your choice.

2. **Guard FIRST — before writing anything.** A version already published on
   the conda channel cannot be re-released (anaconda.org rejects duplicate
   uploads — this bit us with 0.2.7). Confirm the target is free:
   - `RELEASE_VERSION=<X.Y.Z> ./dev_scripts/check_release_version.sh --require-untagged`
   - It checks the `petrnovak/dante_tir` conda channel and local git tags.
   - If it fails, STOP and report — do not proceed. Pick a higher version.

3. **Run the release gate locally — `./tests.sh long`.** This is the step that
   `release.yml` runs before it builds anything, and `smoke`/`short` do
   NOT cover it: `long` is the only test that exercises Round 4, the mmseqs
   clustering of TIR sequences and the final genome extraction. 0.3.0 burned
   three tag cycles because this was skipped. It takes ~90 s:
   `NCPU=2 ./tests.sh long`
   Also run `unit`, `smoke` and `short`. If any fail, STOP.

4. **Bump** `version.py` to `<X.Y.Z>`.

5. **Changelog.** Add a `## <X.Y.Z> — <today's date>` section at the top of
   `changelog.md`, summarizing the commits since the last release tag
   (`git log <last-tag>..HEAD`). Match the prose style of existing entries.

6. **Commit** exactly as `release <X.Y.Z>` (subject line), with a short body and
   the trailer:
   `Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>`
   If `git commit` complains about identity, set it for this repo:
   `git config user.name "Petr Novak" && git config user.email "petr@umbr.cas.cz"`.

7. **Tag** `git tag <X.Y.Z>`.

8. **STOP. Do NOT push.** The user pushes themselves — the tag-driven
   `release.yml` CI publishes on tag push (conda package, SIF on GHCR, GitHub
   release), and they control timing. Report
   the state and give the exact push command, e.g.
   `git push origin main && git push origin <X.Y.Z>`.
   (Only mention `--force`/`--force-with-lease` if the branch/tag was rewritten.)

Related: see the `release-workflow` and `env-tools` memories.

## If the release gate fails after the tag is pushed

A failed gate publishes **nothing** — `release.yml` runs the gate before
`conda-build` and the SIF build, so there is no artifact and no anaconda upload. Confirm with
`RELEASE_VERSION=<X.Y.Z> ./dev_scripts/check_release_version.sh`; if it still
reports the version as free, the tag may be moved rather than burning a version:

```bash
git commit ...                       # the fix
git tag -f <X.Y.Z> HEAD
# then the user pushes:
git push origin main
git push --delete origin <X.Y.Z>
git push origin <X.Y.Z>
```

Point the tag at the *fix* commit, not the `release <X.Y.Z>` commit — the
workflow only asserts that `version.py` matches the tag name, which still holds.

## Hard-won notes

- **CI environment ≠ a fresh conda env.** The 0.3.0 gate failed with
  `there is no package called 'GenomeInfoDbData'` while the same
  `requirements.txt` resolved perfectly in a locally created env. Installing the
  runtime stack *into* an env that already held `conda-build` produced a broken
  Bioconductor installation. The gate now builds its own env; do not move the
  runtime deps back into the build env.
- **Do not debug CI blind.** `dante_tir.py` prints R's stderr on failure and
  both workflows load the R stack right after install. If a failure ever again
  reports only an exit status, fix the visibility first — two cycles were spent
  guessing before that was added.
- **The GitHub release is automatic** (`github-release` job, notes taken from
  the changelog section for the tag). It runs only after both the conda upload
  and the GHCR push succeed. Before 0.3.1 this was manual.
- **SIF failed but conda succeeded?** Run the workflow manually:
  `gh workflow run release.yml -f tag=<X.Y.Z>`. It rebuilds, tests and pushes
  the image, skips conda, and creates the GitHub release if it is missing.
  It builds from the tag, so a transient failure needs nothing else; a fix to
  `Singularity.def` needs the tag moved to the fix commit first. The push run
  the moved tag triggers fails at the conda guard, which is expected.
- **Runtime pins** live in `requirements.txt` and are mirrored in
  `conda/dante_tir/meta.yaml`. Change both together.
