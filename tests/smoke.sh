#!/bin/bash
# tests/smoke.sh — structural sanity check (~ 15 s on 4 CPU).
#
# The smoke slice (Chr1[1..2 Mb] of tiny_pea) deliberately produces zero
# detected TIRs — there aren't enough Subclass_1 copies in 2 Mb to clear
# dante_tir's per-superfamily assembly threshold. Smoke verifies that
# the python + R + CAP3 + mmseqs2 + BLAST integration runs cleanly and
# that both entry-point CLIs work.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DATA="$ROOT/tests/data/smoke"
OUT="$ROOT/tmp/tests/smoke"
NCPU="${NCPU:-2}"

rm -rf "$OUT"
mkdir -p "$OUT"
cd "$ROOT"

echo "=== CLI checks ==="
./dante_tir.py --version
./dante_tir.py --help > /dev/null
./dante_tir_summary.R --help > /dev/null

echo
echo "=== dante_tir on smoke subset (~ 15 s, expect 0 TIRs detected) ==="
./dante_tir.py -g "$DATA/smoke.gff3" -f "$DATA/smoke.fasta" \
               -o "$OUT/out" -c "$NCPU"

# When zero TIRs are detected, dante_tir does not write
# DANTE_TIR_final.{gff3,fasta} — only the log/ directory. The smoke
# slice is sized to hit that path on purpose.
[ -d "$OUT/out" ]      || { echo "FAIL: output dir missing"; exit 1; }
[ -d "$OUT/out/log" ]  || { echo "FAIL: log/ dir missing"; exit 1; }
[ -s "$OUT/out/log/stderr.txt" ] || \
  { echo "FAIL: log/stderr.txt missing or empty"; exit 1; }
grep -q "No TIRs found\|Exiting" "$OUT/out/log/stderr.txt" || \
  { echo "FAIL: R log doesn't show expected 'No TIRs found' marker"; \
    tail -20 "$OUT/out/log/stderr.txt"; exit 1; }

echo
echo "smoke PASSED"
