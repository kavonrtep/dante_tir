#!/bin/bash
# tests/unit.sh — fast pure-python unit tests (no external tools / data).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"

echo "=== unit: extract_flanking_regions streaming == dict (OOM fix) ==="
python3 "$ROOT/tests/test_extract_flanking_regions.py"

echo
echo "=== unit: no blocking subprocess call pipes its child's output ==="
python3 "$ROOT/tests/test_subprocess_pipes.py"

echo
echo "=== unit: CAP3 failures are reported, not silently swallowed ==="
python3 "$ROOT/tests/test_cap3_guard.py"

echo "unit tests OK"
