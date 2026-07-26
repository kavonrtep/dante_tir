#!/usr/bin/env python3
"""Stream a BLAST tabular table into per-subject coverage profiles and switch points.

Round 3 of DANTE_TIR self-BLASTs each superfamily's flanking regions against
themselves and, for every subject, looks for the position where coverage jumps
from "flank" to "element" (the switch point / cp). On large, high-copy genomes
that BLAST table is enormous -- for a 89 Gb genome with 168k EnSpm/CACTA
domains the *already filtered* table was 41.5 GB / 2.5e9 rows, which no R
data.frame can hold (`long vectors not supported yet`).

Nothing downstream needs the table: `round3()` consumes only the named vector
`saccver -> cp`. This tool therefore replaces the read-then-filter-then-split
path entirely. It reads outfmt-6 from stdin, applies the row predicates, folds
each surviving hit straight into a per-subject coverage profile, and emits one
small `id<TAB>cp` table.

Memory is a dense int32 difference array of `n_subjects x (max_len + 2)`
(4.2 GB for the CACTA class above); the BLAST stream itself is never stored.

The filter predicates mirror `filter_blast3()` and the two switch-point
functions mirror `find_switch_point_from_blast_coverage2/3()` in dt_utils.R --
including their quirks, see cp_cumsum(). tests/test_round3_cp.R asserts the
equivalence against the R implementations.
"""

import argparse
import gzip
import os
import resource
import sys

import numpy as np

# Column names understood in --columns. Only saccver/sstart/send are required;
# the rest enable the corresponding predicate.
KNOWN_COLUMNS = ("qaccver", "saccver", "pident", "length", "sstart", "send",
                 "evalue", "bitscore", "mismatch", "gapopen", "qstart", "qend",
                 "slen", "qlen")

# The default layout Round 3 asks blastn for: everything the predicates need
# and nothing else (blastn's own 12-column default is ~2x the text volume).
DEFAULT_COLUMNS = "qaccver,saccver,pident,length,sstart,send"

# The full `-outfmt 6` default, for reprocessing tables produced elsewhere.
STD_COLUMNS = ("qaccver,saccver,pident,length,mismatch,gapopen,"
               "qstart,qend,sstart,send,evalue,bitscore")


def log(msg):
    sys.stderr.write("blast_cp: %s\n" % msg)
    sys.stderr.flush()


def scan_fasta(path):
    """Return (names, max_length) for a FASTA file.

    Names are the first whitespace-delimited token of each header, in file
    order; that order defines the subject index used internally.
    """
    names = []
    max_len = 0
    cur = 0
    have = False
    with open(path, "rb") as fh:
        for line in fh:
            if line[:1] == b">":
                if have and cur > max_len:
                    max_len = cur
                names.append(line[1:].split(None, 1)[0].decode("utf-8", "replace"))
                cur = 0
                have = True
            else:
                cur += len(line.rstrip(b"\r\n"))
    if have and cur > max_len:
        max_len = cur
    return names, max_len


def read_blocks(fh, chunk_bytes):
    """Yield byte blocks from fh, each ending on a line boundary."""
    rem = b""
    while True:
        block = fh.read(chunk_bytes)
        if not block:
            if rem:
                yield rem
            return
        if rem:
            block = rem + block
        cut = block.rfind(b"\n")
        if cut < 0:
            rem = block
            continue
        yield block[:cut + 1]
        rem = block[cut + 1:]


class SubjectIndex:
    """Maps subject accessions (as bytes) to dense row indices.

    Subject names in Round 3 are bare integers, which allows a vectorised
    lookup table. Anything else falls back to a dict over the distinct names
    seen in each chunk (fewer lookups than rows, so still cheap).
    """

    def __init__(self, names):
        self.names = names
        self.by_name = {n.encode("utf-8", "replace"): i for i, n in enumerate(names)}
        self.numeric = all(n.isdigit() and (n == "0" or n[0] != "0") for n in names)
        self.lut = None
        if self.numeric and names:
            ids = np.array([int(n) for n in names], dtype=np.int64)
            self.lut = np.full(int(ids.max()) + 1, -1, dtype=np.int32)
            self.lut[ids] = np.arange(len(names), dtype=np.int32)

    def lookup(self, col):
        """Vectorised name -> index; -1 for names absent from the subject set."""
        if self.lut is not None:
            try:
                vals = col.astype(np.int64)
            except ValueError:
                pass
            else:
                out = np.full(vals.shape, -1, dtype=np.int32)
                ok = (vals >= 0) & (vals < self.lut.size)
                out[ok] = self.lut[vals[ok]]
                return out
        uniq, inv = np.unique(col, return_inverse=True)
        mapped = np.array([self.by_name.get(u, -1) for u in uniq], dtype=np.int32)
        return mapped[inv]


def cp_win200(cov):
    """find_switch_point_from_blast_coverage2() -- sliding 200 bp windows.

    R (dt_utils.R:478):
        swp <- seq(W, L-200, by = 1)
        m1  <- window mean ending at swp; m2 <- window mean starting at swp+1
        cp  <- which.max((m2 + 1)/(m1 + 2)) + W
        keep if mcov1 < 3 & mcov2 > 8 | mcov1 < 2 & mcov2 > 6
    """
    W = 200
    L = cov.size
    if L < 500:
        return None
    # csum[i] = sum of positions 1..i, so a window (a..b] is csum[b] - csum[a].
    csum = np.empty(L + 1, dtype=np.float64)
    csum[0] = 0.0
    np.cumsum(cov, dtype=np.float64, out=csum[1:])
    swp = np.arange(W, L - 200 + 1, dtype=np.int64)   # 1-based positions
    if swp.size == 0:
        return None
    m1 = (csum[swp] - csum[swp - W]) / W
    m2 = (csum[swp + W] - csum[swp]) / W
    k = int(np.argmax((m2 + 1.0) / (m1 + 2.0)))       # R which.max: first max
    cp = k + 1 + W
    mcov1 = m1[k]
    mcov2 = m2[k]
    if (mcov1 < 3 and mcov2 > 8) or (mcov1 < 2 and mcov2 > 6):
        return cp
    return None


def cp_cumsum(cov):
    """find_switch_point_from_blast_coverage3() -- cumulative means (CACTA).

    Replicates dt_utils.R:499 exactly, including two quirks that change which
    subjects pass. Do not "fix" them here; they are what the R path does and
    tests/test_round3_cp.R asserts the match. See docs/round3_scaling.md 8.1.

      * m1/m2 are indexed by position, but the QC thresholds are read at
        `cp - W`, i.e. 200 bp before the detected switch point (in cp_win200
        the same expression is correct, because there the vectors are indexed
        by window offset).
      * m2 divides by (L - i) rather than (L - i + 1), so it is the mean of
        positions i..L over one element too few, and is Inf at i = L. The edge
        zeroing below hides that from which.max.
    """
    W = 200
    L = cov.size
    if L < 500:
        return None
    pos = np.arange(1, L + 1, dtype=np.float64)
    sl = np.cumsum(cov, dtype=np.float64)
    sr = np.cumsum(cov[::-1], dtype=np.float64)[::-1]
    m1 = sl / pos
    with np.errstate(divide="ignore", invalid="ignore"):
        m2 = sr / (L - pos)
        m12 = (m2 + cov.mean() * 0.04) / (m1 + cov.mean() * 0.04)
    m12[:W] = 0.0            # R: m12[1:W] <- 0
    m12[L - W - 1:] = 0.0    # R: m12[(L-W):L] <- 0
    if not np.any(np.isfinite(m12)):
        return None
    cp = int(np.nanargmax(m12)) + 1   # R which.max skips NA/NaN
    j = cp - W                        # R: m1[cp - W], 1-based
    if j < 1:
        # R would index with a negative subscript here and error out; the edge
        # zeroing above makes this unreachable for L >= 500.
        return None
    mcov1 = m1[j - 1]
    mcov2 = m2[j - 1]
    if (mcov1 < 3 and mcov2 > 20) or (mcov1 < 2 and mcov2 > 10):
        return cp
    return None


METHODS = {"win200": cp_win200, "cumsum": cp_cumsum}


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Stream a BLAST tabular table into per-subject coverage "
                    "profiles and switch points.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--subjects", required=True,
                   help="FASTA of the BLAST subject database (defines the "
                        "subject set, their order and the profile width)")
    p.add_argument("--method", required=True, choices=sorted(METHODS),
                   help="switch-point method: win200 = "
                        "find_switch_point_from_blast_coverage2, cumsum = ...3")
    p.add_argument("--columns", default=DEFAULT_COLUMNS,
                   help="comma-separated outfmt-6 column names, or 'std' for "
                        "the blastn default 12-column layout")
    p.add_argument("--input", default="-",
                   help="input table ('-' = stdin, '.gz' supported)")
    p.add_argument("--out", default="-", help="output cp table ('-' = stdout)")
    p.add_argument("--keep-hits", default=None,
                   help="also write the surviving saccver/sstart/send rows "
                        "here ('.gz' compresses)")
    p.add_argument("--min-length", type=float, default=150,
                   help="keep hits with length > this")
    p.add_argument("--min-identity", type=float, default=80,
                   help="keep hits with pident > this (needs the pident column)")
    p.add_argument("--max-evalue", type=float, default=None,
                   help="keep hits with evalue < this (needs the evalue column)")
    p.add_argument("--min-sstart", type=float, default=12,
                   help="keep hits with sstart > this")
    p.add_argument("--no-self-filter", action="store_true",
                   help="keep hits whose query and subject accessions match "
                        "(filter_blast3 drops them)")
    p.add_argument("--chunk-bytes", type=int, default=64 * 1024 * 1024,
                   help="input block size")
    p.add_argument("--progress", type=int, default=0,
                   help="report progress every N million input rows (0 = off)")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    cols = STD_COLUMNS if args.columns == "std" else args.columns
    cols = [c.strip() for c in cols.split(",") if c.strip()]
    for c in cols:
        if c not in KNOWN_COLUMNS:
            sys.exit("unknown column name: %s" % c)
    idx_of = {c: i for i, c in enumerate(cols)}
    for required in ("saccver", "sstart", "send"):
        if required not in idx_of:
            sys.exit("--columns must include %s" % required)
    if args.min_identity is not None and "pident" not in idx_of:
        # Identity is then expected to have been pushed into blastn itself
        # (-perc_identity); saying so beats silently dropping the predicate.
        log("note: no pident column - identity filtering is left to BLAST")
    if args.max_evalue is not None and "evalue" not in idx_of:
        sys.exit("--max-evalue needs the evalue column in --columns")

    names, max_len = scan_fasta(args.subjects)
    if not names:
        sys.exit("no sequences in %s" % args.subjects)
    if max_len < 1:
        sys.exit("no sequence data in %s" % args.subjects)
    subjects = SubjectIndex(names)

    # Difference array: +1 at sstart, -1 at send+1, one row per subject.
    # Row 0 of each subject is unused so positions stay 1-based.
    stride = max_len + 2
    nbytes = len(names) * stride * 4
    log("%d subjects, max length %d, profile array %.2f GB"
        % (len(names), max_len, nbytes / 1e9))
    flat = np.zeros(len(names) * stride, dtype=np.int32)

    from_stdin = args.input == "-"
    fh = (sys.stdin.buffer if from_stdin
          else (gzip.open(args.input, "rb") if args.input.endswith(".gz")
                else open(args.input, "rb")))
    keep_fh = None
    if args.keep_hits:
        keep_fh = (gzip.open(args.keep_hits, "wb") if args.keep_hits.endswith(".gz")
                   else open(args.keep_hits, "wb"))

    ncol = len(cols)
    i_sacc, i_sstart, i_send = idx_of["saccver"], idx_of["sstart"], idx_of["send"]
    i_q = idx_of.get("qaccver")
    n_in = n_kept = n_unknown = n_oob = 0
    next_report = args.progress * 1000000

    try:
        for block in read_blocks(fh, args.chunk_bytes):
            fields = block.split()
            if len(fields) % ncol:
                sys.exit("input has %d fields, not a multiple of the %d "
                         "columns in --columns" % (len(fields), ncol))
            table = np.array(fields).reshape(-1, ncol)
            n_in += table.shape[0]

            sstart = table[:, i_sstart].astype(np.int64)
            send = table[:, i_send].astype(np.int64)
            # filter_blast3(): sstart > 12 and sstart < send (drops the
            # minus-strand orientation as well as zero-length hits)
            keep = (sstart > args.min_sstart) & (sstart < send)
            if "length" in idx_of:
                keep &= table[:, idx_of["length"]].astype(np.int64) > args.min_length
            if "pident" in idx_of and args.min_identity is not None:
                keep &= table[:, idx_of["pident"]].astype(np.float64) > args.min_identity
            if args.max_evalue is not None:
                keep &= table[:, idx_of["evalue"]].astype(np.float64) < args.max_evalue
            if i_q is not None and not args.no_self_filter:
                # filter_blast3(): gsub("_+", "", qaccver) != saccver
                q = table[:, i_q]
                if b"_" in block:
                    q = np.char.replace(q, b"_", b"")
                keep &= q != table[:, i_sacc]

            # A coordinate past the subject length can only mean --subjects is
            # not the database that was searched; folding it in would corrupt
            # the neighbouring subject's row, so drop it loudly instead.
            too_long = send > max_len
            if too_long.any():
                n_oob += int((keep & too_long).sum())
                keep &= ~too_long

            if not keep.all():
                sacc = table[:, i_sacc][keep]
                sstart = sstart[keep]
                send = send[keep]
            else:
                sacc = table[:, i_sacc]

            sidx = subjects.lookup(sacc) if sacc.size else np.empty(0, np.int32)
            if sidx.size and not (sidx >= 0).all():
                known = sidx >= 0
                n_unknown += int((~known).sum())
                sidx = sidx[known]
                sacc = sacc[known]
                sstart = sstart[known]
                send = send[known]
            n_kept += sidx.size

            if sidx.size:
                base = sidx.astype(np.int64) * stride
                starts = base + sstart
                ends = base + np.minimum(send + 1, stride - 1)
                # Fold the chunk down to unique cells first: a scattered
                # read-modify-write over a multi-GB array is dominated by
                # cache misses, and there are far fewer distinct cells than
                # hits.
                u, c = np.unique(starts, return_counts=True)
                flat[u] += c.astype(np.int32)
                u, c = np.unique(ends, return_counts=True)
                flat[u] -= c.astype(np.int32)
                if keep_fh is not None:
                    keep_fh.write(b"".join(
                        b"%s\t%d\t%d\n" % (a, s, e)
                        for a, s, e in zip(sacc, sstart, send)))

            if next_report and n_in >= next_report:
                log("%d M rows in, %d M kept" % (n_in // 1000000, n_kept // 1000000))
                next_report += args.progress * 1000000
    finally:
        if not from_stdin:
            fh.close()
        if keep_fh is not None:
            keep_fh.close()

    if n_unknown:
        log("warning: %d hits had a subject accession absent from %s"
            % (n_unknown, os.path.basename(args.subjects)))
    if n_oob:
        log("warning: %d hits ended past the subject length in %s and were "
            "dropped -- is it the database that was searched?"
            % (n_oob, os.path.basename(args.subjects)))

    method = METHODS[args.method]
    out = sys.stdout if args.out == "-" else open(args.out, "w")
    n_cp = n_profiles = 0
    try:
        out.write("id\tcp\n")
        for i, name in enumerate(names):
            row = flat[i * stride:(i + 1) * stride]
            if not row.any():
                continue
            cov = np.cumsum(row[1:])
            nz = np.flatnonzero(cov)
            if nz.size == 0:
                continue
            # coverage(GRanges) in R yields a vector of length max(send), so
            # the profile stops at the last covered base, not the subject length.
            cov = cov[:nz[-1] + 1]
            n_profiles += 1
            cp = method(cov)
            if cp is None:
                out.write("%s\tNA\n" % name)
            else:
                out.write("%s\t%d\n" % (name, cp))
                n_cp += 1
    finally:
        if out is not sys.stdout:
            out.close()

    log("%d rows read, %d kept, %d subjects with coverage, %d switch points, "
        "peak RSS %.2f GB"
        % (n_in, n_kept, n_profiles, n_cp,
           resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6))
    return 0


if __name__ == "__main__":
    sys.exit(main())
