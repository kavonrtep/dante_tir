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

3. **Bump** `version.py` to `<X.Y.Z>`.

4. **Changelog.** Add a `## <X.Y.Z> — <today's date>` section at the top of
   `changelog.md`, summarizing the commits since the last release tag
   (`git log <last-tag>..HEAD`). Match the prose style of existing entries.

5. **Commit** exactly as `release <X.Y.Z>` (subject line), with a short body and
   the trailer:
   `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`
   If `git commit` complains about identity, set it for this repo:
   `git config user.name "Petr Novak" && git config user.email "petr@umbr.cas.cz"`.

6. **Tag** `git tag <X.Y.Z>`.

7. **STOP. Do NOT push.** The user pushes themselves — the tag-driven
   `conda-release.yml` CI publishes on tag push, and they control timing. Report
   the state and give the exact push command, e.g.
   `git push origin main && git push origin <X.Y.Z>`.
   (Only mention `--force`/`--force-with-lease` if the branch/tag was rewritten.)

Related: see the `release-workflow` and `env-tools` memories.
