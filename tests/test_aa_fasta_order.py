#!/usr/bin/env python3
"""The amino-acid FASTA handed to mmseqs2 must have a stable record order.

--max_class_size splits classes along mmseqs2 clusters, and mmseqs2 clustering
is *order-sensitive*: measured on run-000129's MuDR domains (11,260 sequences,
pipeline parameters), reordering the input alone moved 258 of 5,983 clusters
and changed the cluster count to 6,016. The clustering *result* is otherwise
stable -- rerunning gives the same partition, and 1 thread matches 4.

Note what that stability does not cover: mmseqs does not write that result in a
stable *order*. Two maize runs over a bit-identical FASTA produced an identical
partition but a different line order in clusters_cluster.tsv. Nothing downstream
may inherit that order -- see test_cluster_grouping_order.py.

So the determinism of a split run rests on this FASTA coming out the same way
every time. Order here is the GFF3's, carried through dicts and lists; the risk
would be a `set` creeping into that chain, whose iteration order varies between
processes because Python randomises string hashing. These tests therefore build
the FASTA in separate interpreters with different PYTHONHASHSEED values and
compare bytes.
"""

import hashlib
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FAILS = 0


def check(cond, msg):
    global FAILS
    if cond:
        print("  ok: " + msg)
    else:
        print("FAIL: " + msg)
        FAILS += 1


GFF_HEADER = "##gff-version 3\n"
# Deliberately interleaved classifications and non-monotonic IDs, so any
# regrouping or sorting shows up as a different order.
ROWS = [
    ("Chr1", 100, 400, "Class_II|Subclass_1|TIR|hAT", "MKVLAAGIVR"),
    ("Chr1", 900, 1200, "Class_II|Subclass_1|TIR|EnSpm/CACTA", "MTTQWWPLAS"),
    ("Chr2", 50, 350, "Class_II|Subclass_1|TIR|hAT", "MGGKLPQRTY"),
    ("Chr1", 500, 800, "Class_II|Subclass_1|TIR|MuDR/Mutator", "MPPLKDDSCV"),
    ("Chr2", 700, 1000, "Class_II|Subclass_1|TIR|EnSpm/CACTA", "MRRAYYNETL"),
    ("Chr3", 10, 300, "Class_II|Subclass_1|TIR|hAT", "MFFSDDQPLK"),
    ("Chr2", 1500, 1800, "Class_II|Subclass_1|TIR|MuDR/Mutator", "MNNCVVTREE"),
]


def write_gff(path):
    with open(path, "w") as f:
        f.write(GFF_HEADER)
        for seqid, start, end, cls, aa in ROWS:
            # the quality attributes gff3_quality_ok() requires, all passing
            attrs = ("Final_Classification=%s;Region_Seq=%s;Identity=0.85;"
                     "Similarity=0.92;Relat_Interruptions=0.0;Relat_Length=1.0"
                     % (cls, aa))
            f.write("\t".join([
                seqid, "dante", "protein_match", str(start), str(end), ".",
                "+", ".", attrs]) + "\n")
    return path


BUILD = r"""
import sys, os, hashlib
sys.path.insert(0, %r)
import dt_utils as dt
gff, outdir = sys.argv[1], sys.argv[2]
domains = dt.get_tir_records_from_dante(gff)
aa = dt.extract_aa_sequences_from_dante(domains)
files = dt.save_aa_sequences_by_superfamily(aa, outdir)
for cls in sorted(files):
    with open(files[cls], 'rb') as fh:
        print("%%s %%s" %% (os.path.basename(files[cls]),
                          hashlib.md5(fh.read()).hexdigest()))
"""


def build(gff, outdir, hashseed):
    env = dict(os.environ, PYTHONHASHSEED=str(hashseed))
    os.makedirs(outdir, exist_ok=True)
    res = subprocess.run([sys.executable, "-c", BUILD % ROOT, gff, outdir],
                         capture_output=True, text=True, env=env)
    if res.returncode != 0:
        print(res.stderr)
        raise SystemExit("building the aa FASTA failed")
    return res.stdout.strip()


def main():
    wd = tempfile.mkdtemp()
    gff = write_gff(os.path.join(wd, "dante.gff3"))

    # Same input, different string-hash seeds, separate interpreters.
    a = build(gff, os.path.join(wd, "a"), 0)
    b = build(gff, os.path.join(wd, "b"), 12345)
    c = build(gff, os.path.join(wd, "c"), 999)
    check(a == b == c,
          "the aa FASTA is byte-identical across PYTHONHASHSEED values")
    if a != b:
        print("    seed 0    : " + a.replace("\n", " | "))
        print("    seed 12345: " + b.replace("\n", " | "))

    # And the order is the GFF3's, which is what makes it reproducible.
    hat = os.path.join(wd, "a", "Class_II_Subclass_1_TIR_hAT_aa_sequences.fasta")
    if not os.path.exists(hat):
        check(False, "expected per-class FASTA was not written")
    else:
        ids = [l[1:].strip() for l in open(hat) if l.startswith(">")]
        expected = [str(i + 1) for i, r in enumerate(ROWS)
                    if r[3] == "Class_II|Subclass_1|TIR|hAT"]
        check(ids == expected,
              "records keep GFF3 order (%s)" % ",".join(ids))

    if FAILS:
        print("test_aa_fasta_order: %d FAILURES" % FAILS)
        return 1
    print("test_aa_fasta_order: PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
