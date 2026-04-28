#!/bin/sh
# Install dante_tir into $PREFIX. noarch package; copies the tree into
# $PREFIX/share/dante_tir and symlinks the user-facing CLIs into
# $PREFIX/bin. detect_tirs.R is a sibling helper invoked by dante_tir.py
# via os.path.dirname(__file__); leaving it inside share/ keeps that
# resolution working.
set -x -e

PKG_DIR="${PREFIX}/share/dante_tir"
mkdir -p "${PREFIX}/bin" "${PKG_DIR}"
cp -r . "${PKG_DIR}"

# tests/data/ is a multi-megabyte CI fixture; not useful to end users
# and would inflate the installed package. Keep tests/*.sh on disk in
# case someone wants to re-run the test scripts after installing.
rm -rf "${PKG_DIR}/tests/data"

ln -s "${PKG_DIR}/dante_tir.py"        "${PREFIX}/bin/dante_tir.py"
ln -s "${PKG_DIR}/dante_tir_summary.R" "${PREFIX}/bin/dante_tir_summary.R"
