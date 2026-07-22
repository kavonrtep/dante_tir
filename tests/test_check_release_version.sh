#!/bin/bash
# tests/test_check_release_version.sh
#
# Regression test for the release guard (dev_scripts/check_release_version.sh):
# a version already published on the conda channel must be rejected, and a
# fresh version must be accepted. Uses the RELEASE_VERSION override so the
# repo's real version.py is never touched.
#
# Requires network access to anaconda.org. If the channel is unreachable the
# guard prints a warning and skips its conda check; this test then SKIPs rather
# than producing a false failure.
set -uo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
GUARD="$ROOT/dev_scripts/check_release_version.sh"

echo "=== guard rejects an already-published version (0.2.7) ==="
out=$(RELEASE_VERSION=0.2.7 bash "$GUARD" 2>&1); rc=$?
echo "$out"
if printf '%s\n' "$out" | grep -q "could not reach"; then
  echo "SKIP: anaconda.org unreachable"; exit 0
fi
[ "$rc" -ne 0 ] || { echo "FAIL: guard accepted already-published 0.2.7"; exit 1; }
printf '%s\n' "$out" | grep -q "already published" || \
  { echo "FAIL: guard failed for the wrong reason"; exit 1; }
echo "  ok: 0.2.7 rejected"

echo
echo "=== guard accepts a fresh, unpublished version (99.99.99) ==="
out=$(RELEASE_VERSION=99.99.99 bash "$GUARD" --require-untagged 2>&1); rc=$?
echo "$out"
[ "$rc" -eq 0 ] || { echo "FAIL: guard rejected fresh version 99.99.99"; exit 1; }
echo "  ok: 99.99.99 accepted"

echo
echo "check_release_version guard test OK"
