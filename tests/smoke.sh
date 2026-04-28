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

# Even when zero TIRs are detected, the pipeline now emits valid empty
# DANTE_TIR_final.{gff3,fasta} with proper headers (gff3 has the
# ##gff-version 3 + ##DANTE_TIR version banner) so downstream tools
# don't have to special-case "ran cleanly, found nothing".
[ -d "$OUT/out" ]                         || { echo "FAIL: output dir missing"; exit 1; }
[ -d "$OUT/out/log" ]                     || { echo "FAIL: log/ dir missing"; exit 1; }
[ -f "$OUT/out/DANTE_TIR_final.gff3" ]    || { echo "FAIL: DANTE_TIR_final.gff3 missing"; exit 1; }
[ -f "$OUT/out/DANTE_TIR_final.fasta" ]   || { echo "FAIL: DANTE_TIR_final.fasta missing"; exit 1; }
[ -f "$OUT/out/TIR_classification_summary.txt" ] || \
  { echo "FAIL: TIR_classification_summary.txt missing"; exit 1; }
grep -q '^##gff-version 3' "$OUT/out/DANTE_TIR_final.gff3" || \
  { echo "FAIL: DANTE_TIR_final.gff3 missing gff-version header"; exit 1; }
N=$(grep -c -v '^#' "$OUT/out/DANTE_TIR_final.gff3" || true)
[ "$N" -eq 0 ] || { echo "FAIL: smoke unexpectedly detected $N TIRs"; exit 1; }

echo
echo "smoke PASSED"
