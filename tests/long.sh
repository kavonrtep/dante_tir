#!/bin/bash
# tests/long.sh — release gate. Default input is full Chr1 of tiny_pea
# (~ 115 s on 4 CPU: dante_tir 85 s + summary 28 s).
#
# Local deep mode: export LONG_FASTA and LONG_DANTE to override the
# committed inputs. Used on HPC hosts that have larger genomes available
# under /mnt/raid/... etc.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
LONG_FASTA="${LONG_FASTA:-$ROOT/tests/data/long/long.fasta}"
LONG_DANTE="${LONG_DANTE:-$ROOT/tests/data/long/long.gff3}"
OUT="$ROOT/tmp/tests/long"
NCPU="${NCPU:-2}"

[ -s "$LONG_FASTA" ] || { echo "FAIL: $LONG_FASTA missing or empty"; exit 1; }
[ -s "$LONG_DANTE" ] || { echo "FAIL: $LONG_DANTE missing or empty"; exit 1; }

rm -rf "$OUT"
mkdir -p "$OUT"
cd "$ROOT"

echo "=== dante_tir on long input ==="
echo "    fasta: $LONG_FASTA"
echo "    gff:   $LONG_DANTE"
./dante_tir.py -g "$LONG_DANTE" -f "$LONG_FASTA" \
               -o "$OUT/out" -c "$NCPU"

GFF="$OUT/out/DANTE_TIR_final.gff3"
[ -s "$GFF" ] || { echo "FAIL: DANTE_TIR_final.gff3 missing or empty"; exit 1; }
N=$(grep -c -v '^#' "$GFF" || true)
[ "$N" -ge 10 ] || { echo "FAIL: too few TIR records in $GFF: $N (expected >= 10)"; exit 1; }

echo
echo "=== dante_tir_summary.R on long output ==="
./dante_tir_summary.R -g "$GFF" -f "$LONG_FASTA" \
                      -o "$OUT/summary" -t "$NCPU"

[ -s "$OUT/summary/report.html" ] || \
  { echo "FAIL: report.html missing or empty"; exit 1; }
[ -s "$OUT/summary/all_representative_elements_min3.fasta" ] || \
  { echo "FAIL: all_representative_elements_min3.fasta missing or empty"; exit 1; }

echo
echo "long PASSED ($N TIR records)"
