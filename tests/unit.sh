#!/bin/bash
# tests/unit.sh — fast pure-python unit tests (no external tools / data).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"

echo "=== unit: extract_flanking_regions streaming == dict (OOM fix) ==="
python3 "$ROOT/tests/test_extract_flanking_regions.py"

echo "unit tests OK"
