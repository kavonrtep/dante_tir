# Plan: seed the fragmentation jitter per sequence

Implementation plan for making a region's fragments depend only on that region,
so that `--max_class_size` stops changing results when it splits nothing.

Status: **implemented in 0.3.0.** All six invariants hold and the acceptance
criterion is met — see §9 for what it measured.

---

## 1. The problem

`dict_fasta_to_dict_fragments()` (`dt_utils.py:270`) cuts each flanking region
into 100 bp fragments every 70 bp, offsetting each cut by a random jitter:

```python
for seq_id, seq in fasta_dict.items():
    for i in range(0, seq_len, step):
        if i > 0:
            ii = i + random.randint(-jitter, jitter)   # global RNG
```

The draws come from one global stream, seeded once by `random.seed(42)`
(`dante_tir.py:109`). A sequence's fragments therefore depend on **how many
draws happened before it** — that is, on how many sequences were processed
first and how long they were.

Three consequences, in increasing order of nuisance:

1. **Grouping changes fragments.** With `--max_class_size`, sequences are
   regrouped into parts and fragmented per part, so every sequence sees a
   different position in the draw sequence. Measured on `tests/data/short`:
   enabling `--max_class_size 10000` moved the result from 14 records to 9
   **with one part per class** — no split at all, only reordering.
2. **The flag can never be a no-op**, so it cannot be enabled by default, so
   the CAP3 crash at ~1.07 Gbp stays reachable (docs/round3_scaling.md §8.2).
3. **Fragment identity is a function of the whole run.** Re-running one class
   in isolation cannot reproduce the fragments it had inside a full run, which
   makes debugging a specific class awkward.

For context on magnitude: on that same dataset, `--seed` alone moves the result
between 14, 10 and 10 records, so the differences above sit inside the
pipeline's own noise band. This is a correctness-of-design problem, not
evidence that current results are wrong.

## 2. The change

Give every sequence its own generator, seeded from the pipeline seed and the
sequence's identity:

```python
def dict_fasta_to_dict_fragments(fasta_dict, step=70, length=100, jitter=10,
                                 seed=42, salt=""):
    for seq_id, seq in fasta_dict.items():
        rng = random.Random(_fragment_seed(seed, salt, seq_id))
        ...
            ii = i + rng.randint(-jitter, jitter)
```

with a **stable** hash — not Python's `hash()`, which is randomised per process
and would reintroduce exactly the hazard `tests/test_aa_fasta_order.py` guards
against:

```python
def _fragment_seed(seed, salt, seq_id):
    key = f"{seed}\t{salt}\t{seq_id}".encode("utf-8")
    return int.from_bytes(hashlib.blake2b(key, digest_size=8).digest(), "big")
```

After this, a region's fragments are a pure function of
`(seq_id, sequence, step, length, jitter, seed, salt)`. Order, grouping and
neighbouring sequences drop out.

`salt` separates the two directions so that the upstream and downstream regions
of one element do not receive an identical jitter pattern. Harmless either way
(the sequences differ), but free to do properly.

## 3. Files to change

| File | Change |
|---|---|
| `dt_utils.py:270` | `dict_fasta_to_dict_fragments()`: add `seed`/`salt`, per-sequence `random.Random`, add `_fragment_seed()` helper; `import hashlib` |
| `dt_utils.py:302` | `make_fragment_files()`: accept `seed`, pass `salt="upstream"` / `"downstream"` to the two calls (lines 324, 327) |
| `dante_tir.py:217-218` | split path: pass `seed=args.seed` and the matching salt |
| `dante_tir.py:109` | `random.seed(42)` becomes vestigial — fragmentation was its only consumer (verified by grep: the only other `random.` use in the Python stage). Keep it, with a comment, in case a future caller expects a seeded global RNG |
| `dante_tir.py` | thread `args.seed` into `make_fragment_files()` |

Note that `--seed` currently reaches only the R stage; after this it also
controls fragmentation, which is what a user would expect it to mean.

## 4. Tests

New `tests/test_fragmentation.py`, wired into `tests/unit.sh` (pure Python, no
data, instant). Each asserts an invariant the current code violates:

| # | Invariant | Why |
|---|---|---|
| I1 | **Order invariance** — fragmenting `{a,b,c}` and `{c,a,b}` gives identical fragments per sequence | the direct cause of the 14 → 9 shift |
| I2 | **Subset invariance** — `frag({a,b})` equals the `a`,`b` entries of `frag({a,b,c,d})` | this is what makes `--max_class_size` a no-op when it does not split |
| I3 | **Hash-seed invariance** — identical output across three `PYTHONHASHSEED` values, in separate interpreters | catches a future switch to `hash()` |
| I4 | **Seed sensitivity** — a different `seed` gives different fragments | the knob must still do something |
| I5 | **Salt separation** — same id, different salt, different jitter | upstream/downstream independence |
| I6 | **Shape** — every fragment is `length` long, lies inside the sequence, starts within `jitter` of a multiple of `step`, and the count matches the current algorithm's | the rewrite must not change the fragmentation *scheme*, only its seeding |

Existing suites: `smoke` and `short` will change output (expected, see §5).
`tests/test_extract_flanking_regions.py` is unaffected — it tests the stage
before fragmentation.

## 5. Acceptance criteria

1. **The point of the exercise:** `dante_tir.py --max_class_size N` with `N`
   larger than every class produces **byte-identical** `DANTE_TIR_final.gff3`
   and `.fasta` to a run without the flag. Today it does not.
2. All six invariants above hold.
3. `tests.sh short` changes exactly once, and the new output is stable across
   repeated runs at a fixed seed.
4. The magnitude of that one-time change is reported in the changelog, measured
   as record counts before/after on `tests/data/short` and on `long` if it is
   available.

## 6. Rollout

1. Land after the run-000129 exact run completes, so no comparison spans the
   change.
2. Changelog entry stating plainly that Round-1 results shift once, why, and
   that the shift is within the seed-sensitivity band already documented.
3. **Then, separately**, evaluate turning `--max_class_size` on by default
   (CARP ships its equivalent at 1000). That is a second decision with its own
   evidence — docs/round3_scaling.md §8.2.1 puts the useful range at
   8,000–16,000 copies per part — and it must not ride along silently with this
   change.

## 7. Risks

- **Round-1 output shifts once.** Unavoidable; the whole point is that the
  current fragments are an artefact of processing order. Bounded by the
  measured seed-sensitivity band.
- **Stale `working_dir` reuse.** CAP3 skips assembly when `<fasta>.cap.aln`
  exists, so re-running into an existing working directory would pair *new*
  fragment FASTAs with *old* assemblies. This hazard exists today; the change
  makes it likelier to bite. Suggested guard, cheap and independent: skip only
  when the `.cap.aln` is newer than its FASTA, otherwise reassemble.
- **Slight cost.** One `blake2b` and one `random.Random()` per sequence —
  microseconds against fragmenting 6,300 bp, and it removes a shared mutable
  RNG that would otherwise block parallelising fragmentation later.

## 8. Alternatives considered

- **Sort sequences before fragmenting.** Makes order deterministic but not
  grouping-independent: a part still starts its draws fresh, so I2 still fails
  and the flag still changes results.
- **Keep the old behaviour behind `--legacy_fragmentation`.** Rejected: the old
  behaviour is not reproducible in the first place (it depends on grouping), so
  the flag would preserve an artefact rather than a result, and it would double
  the paths every future test has to cover.
- **Drop the jitter entirely.** Would be the simplest determinism fix, but the
  jitter exists so that fragment boundaries do not align across copies, which
  is what lets CAP3 find staggered overlaps. Removing it is a scientific
  change, not a plumbing one, and would need its own evaluation.


## 9. Outcome

Implemented as designed: `_fragment_seed()` (blake2b over
`seed \t salt \t seq_id`) plus a per-sequence `random.Random` in
`dict_fasta_to_dict_fragments()`, with `salt='upstream'` / `'downstream'` at
the two call sites and `seed` threaded from `--seed`. Both the split and
unsplit paths run through the same loop in `dante_tir.py`, so both are covered.

**Acceptance criterion met.** On `tests/data/short`, `--max_class_size 10000`
(above every class, so it splits nothing) now produces `DANTE_TIR_final.gff3`,
`DANTE_TIR_final.fasta` and `TIR_classification_summary.txt` **byte-identical**
to a run without the flag. Before the change the same comparison was 14 records
against 9.

**One-time result shift, as predicted.** `tests/data/short` moves from 15
records to 11. That is the fragments changing once because they are no longer a
function of processing order, and it sits inside the pipeline's own
seed-sensitivity band on this dataset (14 / 10 / 10 records for seeds 42 / 1 / 7
before the change). Repeat runs at a fixed seed are byte-identical.

All six invariants are asserted by `tests/test_fragmentation.py`, wired into
`tests/unit.sh`; it needs no data and runs instantly.

Not done here, deliberately: turning `--max_class_size` on by default. That is
now *safe* to consider, but it is a separate decision with its own evidence
(§8.2.1 of docs/round3_scaling.md), and it would want the CAP3 concurrency
memory budget alongside it.
