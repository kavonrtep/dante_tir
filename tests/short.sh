#!/bin/bash
# tests/short.sh — functional check on Chr1[1..20 Mb] of tiny_pea
# (~ 90 s on 4 CPU: dante_tir 67 s + dante_tir_summary 20 s).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DATA="$ROOT/tests/data/short"
OUT="$ROOT/tmp/tests/short"
NCPU="${NCPU:-2}"

rm -rf "$OUT"
mkdir -p "$OUT"
cd "$ROOT"

echo "=== dante_tir on short subset (~ 70 s) ==="
./dante_tir.py -g "$DATA/short.gff3" -f "$DATA/short.fasta" \
               -o "$OUT/out" -c "$NCPU"

GFF="$OUT/out/DANTE_TIR_final.gff3"
[ -s "$GFF" ] || { echo "FAIL: DANTE_TIR_final.gff3 missing or empty"; exit 1; }
N=$(grep -c -v '^#' "$GFF" || true)
[ "$N" -ge 5 ] || { echo "FAIL: too few TIR records in $GFF: $N (expected >= 5)"; exit 1; }

echo
echo "=== dante_tir_summary.R on short output (~ 20 s) ==="
./dante_tir_summary.R -g "$GFF" -f "$DATA/short.fasta" \
                      -o "$OUT/summary" -t "$NCPU"

[ -s "$OUT/summary/report.html" ] || \
  { echo "FAIL: report.html missing or empty"; exit 1; }
[ -s "$OUT/summary/all_representative_elements_min3.fasta" ] || \
  { echo "FAIL: all_representative_elements_min3.fasta missing or empty"; exit 1; }

echo
echo "short PASSED ($N TIR records)"
