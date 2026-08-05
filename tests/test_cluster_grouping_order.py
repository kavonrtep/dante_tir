#!/usr/bin/env python3
"""Grouping must not inherit mmseqs2's output order.

mmseqs2 clusters *deterministically* but does not write its result in a stable
order: measured on maize, two runs over a bit-identical amino-acid FASTA gave an
identical partition, an identical representative set and an identical cluster
count, but a different line order in `clusters_cluster.tsv` on every rerun (the
`--threads N` completion order).

That order used to reach the output. `group_sequences_by_clusters()` sorted
clusters by size with a *stable* sort, so equal-sized clusters kept mmseqs'
order -- and with many singleton clusters, ties are the common case. The group
member order then became the fragment order handed to CAP3, and CAP3's assembly
is order-sensitive, so the whole pipeline inherited run-to-run churn from it.

These tests feed the same clustering in many different line orders and require
the grouping to come out identical.
"""

import os
import random
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import dt_utils as dt  # noqa: E402

FAILS = 0


def check(cond, msg):
    global FAILS
    if cond:
        print(f"  ok   {msg}")
    else:
        FAILS += 1
        print(f"  FAIL {msg}")


def write_tsv(path, pairs):
    with open(path, 'w') as f:
        for rep, member in pairs:
            f.write(f"{rep}\t{member}\n")


def build_clustering():
    """A clustering with the shapes that matter: ties, singletons, big clusters."""
    clusters = {}
    seq = 1
    # several clusters of identical size -> ties in the size sort
    for _ in range(6):
        rep = seq
        clusters[str(rep)] = [str(seq + i) for i in range(4)]
        seq += 4
    # a large cluster that triggers the splitting path
    rep = seq
    clusters[str(rep)] = [str(seq + i) for i in range(50)]
    seq += 50
    # many singletons -> ties everywhere
    for _ in range(30):
        clusters[str(seq)] = [str(seq)]
        seq += 1
    return clusters


def as_pairs(clusters):
    return [(rep, m) for rep, members in clusters.items() for m in members]


def group_from(pairs, tmpdir, max_class_size):
    d = os.path.join(tmpdir, 'mm')
    os.makedirs(d, exist_ok=True)
    write_tsv(os.path.join(d, 'clusters_cluster.tsv'), pairs)
    return dt.group_sequences_by_clusters({'CLS': d}, max_class_size)['CLS']


def main():
    clusters = build_clustering()
    pairs = as_pairs(clusters)

    for max_class_size in (10000, 20):
        print(f"\n--- max_class_size={max_class_size} ---")
        with tempfile.TemporaryDirectory() as tmp:
            ref = group_from(pairs, tmp, max_class_size)
            check(bool(ref), "reference grouping is non-empty")

            # every member appears exactly once, whatever the order
            flat = [m for members in ref.values() for m in members]
            check(len(flat) == len(set(flat)) == len(pairs),
                  f"grouping covers each of the {len(pairs)} sequences exactly once")

            for seed in range(12):
                shuffled = list(pairs)
                random.Random(seed).shuffle(shuffled)
                got = group_from(shuffled, tmp, max_class_size)
                if got != ref:
                    check(False, f"line order seed={seed} changed the grouping")
                    break
            else:
                check(True, "grouping is identical across 12 shuffled line orders")

            # cluster blocks reversed (mmseqs emits whole clusters, not stray lines)
            rev = as_pairs(dict(reversed(list(clusters.items()))))
            check(group_from(rev, tmp, max_class_size) == ref,
                  "grouping is identical with cluster blocks in reverse order")

            # members reversed within every cluster
            memrev = as_pairs({r: list(reversed(m)) for r, m in clusters.items()})
            check(group_from(memrev, tmp, max_class_size) == ref,
                  "grouping is identical with members reversed within clusters")

    print()
    if FAILS:
        print(f"test_cluster_grouping_order: {FAILS} FAILED")
        sys.exit(1)
    print("test_cluster_grouping_order OK")


if __name__ == '__main__':
    main()
