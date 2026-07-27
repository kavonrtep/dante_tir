#!/usr/bin/env python3
"""A region's fragments must depend on the region alone.

Until 0.3.0 the fragmentation jitter came from one global RNG, so a sequence's
fragments depended on how many draws had happened before it -- i.e. on how many
other sequences were processed first, and how long they were. Regrouping the
input therefore changed the fragments even when it split nothing, which is why
enabling --max_class_size moved tests/data/short from 14 records to 9 with a
single part per class, and why the flag could never be turned on by default.

Each sequence now gets its own generator seeded from (seed, salt, seq_id).
These tests pin the properties that buys, with the subset invariance (I2) being
the one that makes splitting a genuine no-op.
"""

import os
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import dt_utils as dt

FAILS = 0


def check(cond, msg):
    global FAILS
    if cond:
        print("  ok: " + msg)
    else:
        print("FAIL: " + msg)
        FAILS += 1


def seqs(n=6, length=900):
    """Deterministic pseudo-sequences of differing lengths."""
    out = {}
    for i in range(n):
        L = length + 137 * i
        out["seq%d" % i] = "".join("ACGT"[(i * 7 + j * 3) % 4] for j in range(L))
    return out


def main():
    data = seqs()

    # --- I1: order invariance --------------------------------------------
    forward = dt.dict_fasta_to_dict_fragments(data)
    reversed_order = dt.dict_fasta_to_dict_fragments(
        {k: data[k] for k in reversed(list(data))})
    check(forward == reversed_order,
          "I1 fragmenting in a different order gives identical fragments")

    # --- I2: subset invariance -------------------------------------------
    # The property that makes --max_class_size a no-op when it does not split:
    # fragmenting a part must give what that part had inside the whole.
    part = {k: data[k] for k in ("seq1", "seq3")}
    part_frags = dt.dict_fasta_to_dict_fragments(part)
    expected = {k: v for k, v in forward.items()
                if k.rsplit("_", 2)[0] in part}
    check(part_frags == expected,
          "I2 fragments of a subset match that subset inside the whole")

    # --- I3: independent of Python's per-process string hashing -----------
    prog = (
        "import sys; sys.path.insert(0, %r);"
        "import dt_utils as dt, hashlib;"
        "d={'a'*1:'ACGT'*400, 'b':'TTTTGGGGCCCCAAAA'*70, 'seq_x':'ACGTAC'*200};"
        "f=dt.dict_fasta_to_dict_fragments(d);"
        "print(hashlib.md5(repr(sorted(f.items())).encode()).hexdigest())"
        % os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    )
    digests = set()
    for hashseed in ("0", "12345", "999"):
        env = dict(os.environ, PYTHONHASHSEED=hashseed)
        res = subprocess.run([sys.executable, "-c", prog],
                             capture_output=True, text=True, env=env)
        if res.returncode != 0:
            print(res.stderr)
            check(False, "I3 subprocess failed")
            break
        digests.add(res.stdout.strip())
    check(len(digests) == 1,
          "I3 fragments are identical across PYTHONHASHSEED values")

    # --- I4: the seed still does something --------------------------------
    other_seed = dt.dict_fasta_to_dict_fragments(data, seed=7)
    check(other_seed != forward, "I4 a different seed gives different fragments")
    check(dt.dict_fasta_to_dict_fragments(data, seed=7) == other_seed,
          "I4 the same seed reproduces them")

    # --- I5: salt separates upstream from downstream ----------------------
    up = dt.dict_fasta_to_dict_fragments(data, salt="upstream")
    down = dt.dict_fasta_to_dict_fragments(data, salt="downstream")
    check(up != down,
          "I5 the same ids under a different salt are jittered independently")

    # --- I6: the fragmentation scheme itself is unchanged -----------------
    step, length, jitter = 70, 100, 10
    bad_len = bad_bounds = bad_offset = bad_content = 0
    for frag_id, frag in forward.items():
        seq_id, start, end = frag_id.rsplit("_", 2)
        start, end = int(start), int(end)
        seq = data[seq_id]
        if len(frag) != length:
            bad_len += 1
        if start < 0 or end > len(seq):
            bad_bounds += 1
        if frag != seq[start:end]:
            bad_content += 1
        # start is a multiple of step nudged by at most `jitter`, unless it was
        # clamped to the last valid offset
        nearest = round(start / step) * step
        if abs(start - nearest) > jitter and start != len(seq) - length:
            bad_offset += 1
    check(bad_len == 0, "I6 every fragment is exactly %d bp" % length)
    check(bad_bounds == 0, "I6 every fragment lies inside its sequence")
    check(bad_content == 0, "I6 every fragment matches the source sequence")
    check(bad_offset == 0,
          "I6 every offset is within +/-%d of a %d bp step" % (jitter, step))

    # rough count check: one fragment per step, minus those clamped together
    for seq_id, seq in data.items():
        n = sum(1 for k in forward if k.rsplit("_", 2)[0] == seq_id)
        upper = len(range(0, len(seq), step))
        check(0 < n <= upper,
              "I6 %s yields %d fragments (<= %d step positions)"
              % (seq_id, n, upper))
        break   # one is enough; the rest follow the same code path

    if FAILS:
        print("test_fragmentation: %d FAILURES" % FAILS)
        return 1
    print("test_fragmentation: PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
