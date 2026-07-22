#!/bin/bash
# dev_scripts/check_release_version.sh
#
# Release guard: refuse to (re)release a version that is already published.
# Anaconda.org rejects re-uploading an existing version/build, so a duplicate
# version bump wastes a full CI build and then fails at the upload step. This
# script catches the collision up front.
#
# Checks the version in version.py against:
#   1. the conda channel (anaconda.org/<CHANNEL>/<PKG>)  -- always
#   2. existing local git tags                            -- only with
#      --require-untagged (use before creating a new release tag; skip in the
#      release CI, where the tag necessarily already exists)
#
# Usage:
#   dev_scripts/check_release_version.sh                 # conda-channel guard (CI)
#   dev_scripts/check_release_version.sh --require-untagged  # local pre-tag guard
#
# Exit status 0 = version is free to release; non-zero = collision / error.
set -euo pipefail

CHANNEL="petrnovak"
PKG="dante_tir"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"

require_untagged=0
[ "${1:-}" = "--require-untagged" ] && require_untagged=1

# RELEASE_VERSION overrides version.py (used by tests to probe the guard
# without touching the real version.py).
if [ -n "${RELEASE_VERSION:-}" ]; then
  V="$RELEASE_VERSION"
  echo "version (RELEASE_VERSION override) = $V"
else
  V=$(python3 -c "exec(open('$ROOT/version.py').read()); print(__version__)")
  echo "version.py = $V"
fi

fail=0

# 1) conda channel -----------------------------------------------------------
api="https://api.anaconda.org/package/$CHANNEL/$PKG"
json=$(curl -fsSL "$api" 2>/dev/null || true)
if [ -z "$json" ]; then
  echo "WARNING: could not reach $api (skipping conda-channel check)" >&2
else
  published=$(printf '%s' "$json" | python3 -c \
    "import sys,json; print('\n'.join(json.load(sys.stdin).get('versions',[])))")
  if printf '%s\n' "$published" | grep -qx "$V"; then
    echo "::error::version $V is already published on anaconda.org/$CHANNEL/$PKG" >&2
    echo "  -> bump version.py to a new version before releasing." >&2
    fail=1
  else
    echo "ok: $V not present on anaconda.org/$CHANNEL/$PKG"
  fi
fi

# 2) local git tag -----------------------------------------------------------
if [ "$require_untagged" -eq 1 ]; then
  if git -C "$ROOT" rev-parse -q --verify "refs/tags/$V" >/dev/null 2>&1; then
    echo "::error::git tag $V already exists locally" >&2
    fail=1
  else
    echo "ok: no local git tag $V yet"
  fi
fi

if [ "$fail" -ne 0 ]; then
  echo "RELEASE GUARD FAILED for version $V" >&2
  exit 1
fi
echo "release guard OK: $V is free to release"
