#!/bin/bash
# Re-run TIR detection (rounds 1-4) on an existing DANTE_TIR working_dir.
#
# dante_tir.py's Python stages -- domain extraction, flanking regions, CAP3 --
# are deterministic and already succeeded in the run being re-checked, so this
# reuses their output and re-runs only detect_tirs.R. That is where every
# recent change lives (streaming Round 3, query sampling, indexed genome
# access), and it turns a multi-day full re-run into a couple of hours.
#
# The source working_dir is only ever read; everything is written to --work.
#
# Usage:
#   rerun_from_working_dir.sh --src <existing working_dir> --work <fresh dir> \
#                             --genome <genome.fasta> [--threads N] \
#                             [--max-round3-queries N] [--seed N] [--code DIR]
set -euo pipefail

SRC=""; WORK=""; GENOME=""; THREADS=8; MAXQ=0; SEED=42
CODE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MIN_FREE_GB=${MIN_FREE_GB:-20}

die() { echo "ERROR: $*" >&2; exit 1; }
note() { echo "[$(date +%H:%M:%S)] $*"; }

while [ $# -gt 0 ]; do
  case "$1" in
    --src) SRC="$2"; shift 2 ;;
    --work) WORK="$2"; shift 2 ;;
    --genome) GENOME="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --max-round3-queries) MAXQ="$2"; shift 2 ;;
    --seed) SEED="$2"; shift 2 ;;
    --code) CODE="$(cd "$2" && pwd)"; shift 2 ;;
    -h|--help) sed -n '2,20p' "$0"; exit 0 ;;
    *) die "unknown argument: $1" ;;
  esac
done

[ -n "$SRC" ]    || die "--src is required (an existing DANTE_TIR working_dir)"
[ -n "$WORK" ]   || die "--work is required (a fresh, writable directory)"
[ -n "$GENOME" ] || die "--genome is required"
[ -d "$SRC" ]    || die "--src is not a directory: $SRC"
[ -f "$GENOME" ] || die "genome not found: $GENOME"

# ---------------------------------------------------------------- preflight
note "preflight"

[ -f "$SRC/tir_flank_coords.txt" ] || die "$SRC has no tir_flank_coords.txt -- is it a DANTE_TIR working_dir?"
ls "$SRC"/*_regions.fasta >/dev/null 2>&1 || die "$SRC has no *_regions.fasta"

for tool in Rscript blastn makeblastdb python3; do
  command -v "$tool" >/dev/null || die "$tool not on PATH"
done

# blast_cp.py is new in this version and needs numpy; failing here beats
# failing hours later in the middle of Round 3.
[ -f "$CODE/detect_tirs.R" ] || die "detect_tirs.R not found in --code dir: $CODE"
[ -f "$CODE/blast_cp.py" ]   || die "blast_cp.py not found in $CODE (pre-0.2.9 checkout?)"
python3 -c 'import numpy' 2>/dev/null || \
  die "python3 on PATH has no numpy -- blast_cp.py requires it (conda install numpy)"

# Round 4 and the final extraction read the genome through its .fai. If it is
# missing, R will try to create one next to the genome, which needs write
# access there.
if [ ! -f "$GENOME.fai" ]; then
  if [ -w "$(dirname "$GENOME")" ]; then
    note "note: $GENOME.fai missing; it will be created on first use"
  else
    die "$GENOME.fai is missing and $(dirname "$GENOME") is not writable -- run 'samtools faidx $GENOME' somewhere writable first"
  fi
fi

mkdir -p "$WORK"
[ -w "$WORK" ] || die "--work is not writable: $WORK"
FREE_GB=$(df -BG --output=avail "$WORK" | tail -1 | tr -dc '0-9')
[ "$FREE_GB" -ge "$MIN_FREE_GB" ] || \
  die "only ${FREE_GB}G free in $WORK; this needs about ${MIN_FREE_GB}G"

note "code:    $CODE ($(cd "$CODE" && git rev-parse --short HEAD 2>/dev/null || echo 'not a git checkout'))"
note "source:  $SRC"
note "work:    $WORK (${FREE_GB}G free)"
note "genome:  $GENOME"
note "threads: $THREADS   max_round3_queries: $MAXQ   seed: $SEED"

# ---------------------------------------------------------------------- copy
# Only what detect_tirs.R reads: round-1 contigs, the region FASTAs and their
# BLAST databases, and the flank coordinates. Not the CAP3 inputs, not the
# per-class fragment FASTAs -- on run-000129 that is 3.9 GB instead of 52 GB.
# Idempotent, so a re-invocation after a failure skips what is already there.
note "copying inputs (skipping anything already present)"
copy_glob() {   # copy_glob <pattern> <label>
  local n
  n=$(find "$SRC" -maxdepth 1 -name "$1" | wc -l)
  [ "$n" -gt 0 ] || { note "  $2: none found"; return 0; }
  find "$SRC" -maxdepth 1 -name "$1" -print0 | xargs -0 -r cp -n -t "$WORK"
  note "  $2: $n files"
}
copy_glob 'tir_flank_coords.txt'  'flank coordinates'
copy_glob '*_regions.fasta'       'region FASTAs'
copy_glob '*_regions.fasta.n*'    'region BLAST databases'
copy_glob '*upstream_Contig*'     'upstream contigs'
copy_glob '*downstream_Contig*'   'downstream contigs'
note "copied $(du -sh "$WORK" | cut -f1) into $WORK"

# ----------------------------------------------------------------------- run
mkdir -p "$WORK/log"
note "running detect_tirs.R (logs: $WORK/log/{stdout,stderr}.txt)"
START=$(date +%s)
set +e
"$CODE/detect_tirs.R" \
    --contig_dir "$WORK" \
    --output "$WORK" \
    --threads "$THREADS" \
    --genome "$GENOME" \
    --seed "$SEED" \
    --max_round3_queries "$MAXQ" \
    > "$WORK/log/stdout.txt" 2> >(tee "$WORK/log/stderr.txt" >&2)
RC=$?
set -e
ELAPSED=$(( $(date +%s) - START ))
note "detect_tirs.R exited $RC after $((ELAPSED / 3600))h $(((ELAPSED % 3600) / 60))m"

# ------------------------------------------------------------------- summary
echo
echo "================= summary ================="
if [ "$RC" -ne 0 ]; then
  echo "FAILED (exit $RC). Last 20 lines of stderr:"
  tail -20 "$WORK/log/stderr.txt"
  echo
  echo "An RData workspace may have been saved to $WORK/DANTE_TIR.RData"
  exit "$RC"
fi

GFF="$WORK/DANTE_TIR_final.gff3"
if [ -s "$GFF" ]; then
  echo "TIR records: $(grep -c -v '^#' "$GFF" || true)"
  echo "per classification:"
  grep -v '^#' "$GFF" | sed 's/.*Classification=\([^;]*\).*/  \1/' | sort | uniq -c | sort -rn
else
  echo "DANTE_TIR_final.gff3 is missing or empty"
fi
for f in DANTE_TIR_final.fasta TIR_classification_summary.txt; do
  [ -f "$WORK/$f" ] && echo "$f: $(du -h "$WORK/$f" | cut -f1)"
done
echo
echo "per-round element counts:"
grep -E "Number of elements found" "$WORK/log/stderr.txt" || true
echo
echo "Round-3 sampling (empty if --max-round3-queries 0):"
grep -E "round3 (sampling|sample fingerprint|pass 2)" "$WORK/log/stderr.txt" || echo "  (none -- exact path)"
echo
echo "blast_cp.py:"
grep -E "^blast_cp:" "$WORK/log/stderr.txt" | grep -v "^blast_cp: note" || true
echo "==========================================="
