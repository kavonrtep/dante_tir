#!/usr/bin/env python3
"""cap3assembly must not turn a CAP3 crash into a silent empty assembly.

CAP3 indexes its concatenated forward + reverse-complement sequence with a
signed 32-bit int and segfaults past ~1.07 Gbp of input -- measured on
cap3 10.2011: 1.060 Gbp assembles, 1.073 Gbp dies in ~26 s, driven by base
count rather than read count or content. On run-000129 that killed the
168,012-copy EnSpm/CACTA assembly in about a minute, and because the exit
status was never checked, the pipeline wrote a zero-byte .cap.aln, carried on,
and lost that superfamily's entire Round 1.

These tests use a stub `cap3` on PATH, so nothing here needs a real assembler
or a gigabase of input.
"""

import os
import shutil
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


def write_fasta(path, nseq=3, seqlen=50):
    with open(path, "w") as f:
        for i in range(nseq):
            f.write(">r%d\n%s\n" % (i, "ACGT" * (seqlen // 4)))
    return path


def stub_cap3(directory, script):
    """Put a fake `cap3` at the front of PATH and return the new PATH."""
    binp = os.path.join(directory, "bin")
    os.makedirs(binp, exist_ok=True)
    p = os.path.join(binp, "cap3")
    with open(p, "w") as f:
        f.write(script)
    os.chmod(p, 0o755)
    return binp + os.pathsep + os.environ["PATH"]


def main():
    wd = tempfile.mkdtemp()
    try:
        # --- fasta_total_bases counts sequence only ---------------------------
        fa = write_fasta(os.path.join(wd, "a.fasta"), nseq=3, seqlen=48)
        check(dt.fasta_total_bases(fa) == 144,
              "fasta_total_bases counts bases, not headers or newlines (%d)"
              % dt.fasta_total_bases(fa))

        # --- oversized input is refused without invoking CAP3 -----------------
        # A stub that would report success if it ever ran, so reaching CAP3 at
        # all is detectable.
        old_path = os.environ["PATH"]
        os.environ["PATH"] = stub_cap3(wd, "#!/bin/sh\necho ran > %s/CAP3_RAN\nexit 0\n" % wd)
        try:
            res = dt.cap3assembly(fa, max_bases=100)   # 144 bases > 100
            check(res is None, "oversized input returns None instead of a bogus path")
            check(not os.path.exists(os.path.join(wd, "CAP3_RAN")),
                  "oversized input is refused without running CAP3")
            check(not os.path.exists(fa + ".cap.aln"),
                  "oversized input leaves no .cap.aln behind")
            check(os.path.exists(fa + ".cap.err")
                  and "max_class_size" in open(fa + ".cap.err").read(),
                  ".cap.err explains the refusal and names the knob to turn")

            # --- a crashing CAP3 is reported, not swallowed -------------------
            fb = write_fasta(os.path.join(wd, "b.fasta"))
            os.environ["PATH"] = stub_cap3(
                wd, "#!/bin/sh\necho 'Segmentation fault' >&2\nexit 139\n")
            res = dt.cap3assembly(fb)
            check(res is None, "a crashing CAP3 returns None")
            check(not os.path.exists(fb + ".cap.aln"),
                  "a crashing CAP3 leaves no zero-byte .cap.aln (the run-000129 bug)")
            check(os.path.exists(fb + ".cap.err")
                  and "Segmentation fault" in open(fb + ".cap.err").read(),
                  "CAP3's stderr is preserved for diagnosis")

            # --- a zero-byte .cap.aln must not count as a finished assembly ---
            fc = write_fasta(os.path.join(wd, "c.fasta"))
            open(fc + ".cap.aln", "w").close()          # what an old crash left
            os.environ["PATH"] = stub_cap3(wd, "#!/bin/sh\necho REALOUTPUT\nexit 0\n")
            res = dt.cap3assembly(fc)
            check(res == fc + ".cap.aln" and
                  open(fc + ".cap.aln").read().strip() == "REALOUTPUT",
                  "an empty .cap.aln from a crashed run is retried, not trusted")

            # --- a non-empty .cap.aln is still reused --------------------------
            os.environ["PATH"] = stub_cap3(wd, "#!/bin/sh\necho SHOULD_NOT_RUN\nexit 0\n")
            res = dt.cap3assembly(fc)
            check(open(fc + ".cap.aln").read().strip() == "REALOUTPUT",
                  "an existing non-empty assembly is reused, not recomputed")
        finally:
            os.environ["PATH"] = old_path

        # --- the guard sits below the measured crash point --------------------
        # Measured on cap3 10.2011: 1.060 Gbp assembles, 1.073 Gbp segfaults.
        check(dt.CAP3_MAX_BASES <= 1_060_000_000,
              "CAP3_MAX_BASES (%.3f Gbp) stays at or below the largest input measured to work"
              % (dt.CAP3_MAX_BASES / 1e9))
        # 6300 bp regions fragment to ~8993 bases each (run-000129: 168,012
        # regions -> 1.51 Gbp), so the guard corresponds to ~118k regions/part.
        check(int(dt.CAP3_MAX_BASES / 8993) > 100_000,
              "the guard leaves room for ~%d regions per part"
              % int(dt.CAP3_MAX_BASES / 8993))
    finally:
        shutil.rmtree(wd, ignore_errors=True)

    # --- the memory budget that bounds CAP3 concurrency --------------------
    wd2 = tempfile.mkdtemp()
    try:
        # Predictions must track the measured points (docs/round3_scaling.md
        # 8.2.1): 8.98 Mbp -> 0.98 GB, 35.9 -> 5.30, 71.9 -> 12.45.
        for mbp, measured in ((8.98, 0.98), (17.95, 2.25), (35.95, 5.30),
                              (71.9, 12.45)):
            pred = dt.cap3_expected_rss_gb(mbp * 1e6)
            check(abs(pred - measured) / measured < 0.05,
                  "RSS prediction at %.1f Mbp is %.2f GB, measured %.2f"
                  % (mbp, pred, measured))

        # A budget that fits three of these but not four.
        files = []
        for i in range(6):
            p = os.path.join(wd2, "part%d.fasta" % i)
            with open(p, "w") as f:                     # ~1 Mbp each
                for j in range(100):
                    f.write(">r%d\n%s\n" % (j, "ACGT" * 2500))
            files.append(p)
        one = dt.cap3_expected_rss_gb(dt.fasta_total_bases(files[0]))
        n, msg = dt.cap3_pool_size(files, cpu=64, memory_budget_gb=one * 3.5)
        check(n == 3, "a budget of 3.5x one assembly permits 3 at a time (got %d)" % n)
        check("GB budget" in msg, "the decision is reported: %s" % msg.split("\n")[0])

        n, _ = dt.cap3_pool_size(files, cpu=2, memory_budget_gb=one * 100)
        check(n == 2, "a generous budget still respects --cpu (got %d)" % n)

        n, msg = dt.cap3_pool_size(files, cpu=64, memory_budget_gb=one / 10)
        check(n == 1, "an impossible budget still runs one at a time (got %d)" % n)
        check("WARNING" in msg, "and warns that one assembly exceeds the budget")

        check(dt.system_memory_gb() > 0.5,
              "system memory is detectable (%.0f GB)" % dt.system_memory_gb())
    finally:
        shutil.rmtree(wd2, ignore_errors=True)

    if FAILS:
        print("test_cap3_guard: %d FAILURES" % FAILS)
        return 1
    print("test_cap3_guard: PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
