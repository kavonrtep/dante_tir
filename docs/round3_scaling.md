# Round 3 scaling: streaming BLAST → coverage profile → switch point

Design document for fixing the Round-3 self-BLAST blow-up on very large,
high-copy genomes. Written against `0.2.8`; supersedes the `awk` prefilter
introduced in `0.2.7` (commit `dcddcca`).

Status: **P0 and P1 implemented** (see §10); P2 still open. Keep this document
updated as steps are completed.

---

## 1. The failure

Reference run: `/nfsroot/projects/darwin/runs/run-000129/scratch/out/`
(pipeline 1.1.5, genome `GCA_963277665.1` / `drVisAlbu1.1`, **89 Gb**,
host `storage-budejovice2a`, 96 threads).

DANTE domain counts for this genome:

| Superfamily | Domains |
|---|---|
| Class_II\|Subclass_1\|TIR\|**EnSpm/CACTA** | **168,012** |
| Class_II\|Subclass_1\|TIR\|MuDR/Mutator | 11,260 |
| Class_II\|Subclass_1\|TIR\|hAT | 1,548 |
| Class_II\|Subclass_1\|TIR\|PIF/Harbinger | 333 |
| Class_II\|Subclass_1\|TIR\|Sola2 | 1 |

`DANTE_TIR/log/stderr.txt` ends at:

```
---- Identification of elements - Round 3 ----
An error occurred during the pipeline
Error message:
long vectors not supported yet: ../../../include/Rinlinedfuns.h:537
```

Measured evidence:

| Item | Value |
|---|---|
| `working_dir/blastn/…EnSpm_CACTA_upstream3.filtered.tsv` | **41.5 GB — already post-`awk`** |
| rows in that file | **≈ 2.54 × 10⁹** (16.4 bytes/row, measured) |
| the matching `…downstream3.filtered.tsv` | never produced — the run died on the upstream file |
| `benchmarks/dante_tir.tsv` | 13 h 08 m wall, **max_rss 477 GB**, 128 GB written |
| mean load (same benchmark) | **24.7 of 96 cores** |
| CAP3 for this class (Round 1) | segfaulted — see §8 |
| elements found in Round 1 | 75 (genome-wide) |

The `0.2.7` `awk` prefilter (`blast3_reduce_awk()`, `dt_utils.R:850`) worked
exactly as designed and was still about four orders of magnitude short of what
`read.table()` can ingest.

## 2. Why prefiltering cannot fix this

The `0.2.7` filter is **per row**: it shrinks each record (12 columns → 3) and
drops records failing `filter_blast3()`'s predicates. The problem is the
**number of rows**, which is intrinsic to the analysis:

Round 3 self-BLASTs each class's flanking regions against themselves
(`round3()`, `dt_utils.R:1627` — `query_db == blast_db`). For CACTA that is
168,012 × 6,300 bp regions against a 1.06 GB database. After filtering, each
query still hits ≈ 15,100 subjects. No per-row filter can change that; only
reducing the number of query–subject comparisons can.

### 2.1 The data is ~180× oversaturated for its purpose

From a 0.45 % sample (nine 20 MB chunks spread across the 41.5 GB file), the
per-subject **mean coverage depth** over the 6,300 bp region, extrapolated to
the full file:

```
min 5   p5 93   q25 1365   median 3600   q75 5292   p95 6837   p99.9 8575
```

The switch-point criteria only ever test `mcov2 > 20` / `mcov1 < 3`
(`find_switch_point_from_blast_coverage3()`, `dt_utils.R:499`) or
`mcov1 < 3 & mcov2 > 8` (`…coverage2()`, `dt_utils.R:478`). The median subject
carries ~180× more depth than any decision needs.

Alignment lengths in the same sample: mean 1,341 bp, with 59 % of hits in the
1,250–1,500 bp bin.

**This oversaturation is the lever** — it is what makes sampling (§6) a
statistically safe way to bound the work rather than a quality compromise.

### 2.2 A second, independent bottleneck

`benchmarks/dante_tir.tsv` records a mean load of 24.7 cores against
`-num_threads 96`. The Round-3 stage was a 96-thread `blastn` formatting a
~300 GB 12-column text stream into a **single-threaded `awk`** that then wrote
41.5 GB. The serial text path — not the alignment — plausibly capped the run
at ~25 effective cores for most of its 13 hours.

## 3. What Round 3 actually needs

`run_blast_tir_analysis()` (`dt_utils.R:858`) returns
`list(blast_df, blast_cov, cp_vals)`. Tracing the callers
(`dt_utils.R:1688–1689`), `round3()` consumes **only `cp_vals`** — a named
integer vector mapping `saccver → switch point`. `blast_df` and `blast_cov`
exist purely as intermediates and are never read.

So the entire useful output of that 41.5 GB file is a **168k-line, two-column
table**. Everything between the BLAST stream and that table can be replaced.

Also relevant: `swt_pt_function()` (`dt_utils.R:528`) only ever dispatches to
`…coverage2` or `…coverage3` in Round 3 (unknown classes fall back to
`…coverage2`). **Rbeast is not involved in Round 3 at all** — the
BEAST-based `find_switch_point_from_blast_coverage()` is used elsewhere.
Nothing in this path needs to stay in R.

## 4. Proposed architecture

Three independent stages, each shippable on its own.

### P0 — streaming Python parser: BLAST → coverage profile → switch point

Replace `awk` + `read.table()` + `get_coverage_from_blast()` +
`mclapply(swt_pt_fun)` with a single tool consuming outfmt-6 on **stdin** and
emitting the small cp table:

```
blastn -query … -db … -outfmt "6 qaccver saccver pident length sstart send" … | \
  blast_cp.py --method cumsum \
              --columns qaccver,saccver,pident,length,sstart,send \
              --min-length 150 --min-identity 80 \
              --subjects <db.fasta> [--keep-hits hits.tsv.gz] \
  > cp_upstream.tsv
```

(`--subjects` is the searched database: its headers define the subject set and
its longest sequence the profile width. `pident` is requested even though
`-perc_identity` already filters, so the `> 80` boundary stays exactly where
`filter_blast3()` had it.)

Output: `saccver<TAB>cp`, one line per subject that produced a switch point
(`NA` allowed, dropped downstream by `complete.cases`).

`run_blast_tir_analysis()` becomes a thin wrapper: build the command,
`system2()`, read the small TSV, return `cp_vals` with identical names and
semantics. `blast_df` / `blast_cov` are dropped from the return value —
nothing reads them.

`--method` is chosen by R (`swt_pt_function()` stays the single source of
truth for class → method): `win200` = `…coverage2`, `cumsum` = `…coverage3`.

#### Internals

- **Parse** in large blocks: `np.array(buf.split(), dtype=np.int32)` —
  measured **4.1 M rows/s per core**. Shard across a few reader processes.
- **Accumulate a difference array** (`+1` at `sstart`, `−1` at `send+1`), not
  per-position increments. Partition each block by subject-block and use
  `np.bincount` per block: measured **82 M updates/s**.
  For comparison, `np.add.at` on the same data measured **5.6 M/s** — the
  obvious route is 15× slower and should be avoided.
- **Memory**: dense `n_subjects × max_len` int32 = **4.6 GB** for this class
  (181,155 × 6,301 — subject IDs are not contiguous, so size by max ID or map
  IDs to indices). Trivial on a node that just peaked at 477 GB.
  A `--reservoir C` mode (keep ≤ C random intervals per subject, rescale by
  nₛ/C) brings it to ~400 MB if a low-memory mode is ever wanted.
- **Switch point** vectorised in numpy: measured **61 µs/subject → ~10 s for
  all 168k on one core**. Today `…coverage2` runs a 5,600-iteration `apply()`
  per subject inside `mclapply`.

Expected result for this class-direction: from *cannot complete* to ~10–20 min
of parsing at the current BLAST volume, in ~5 GB of RAM.

#### Semantics that must be replicated exactly

1. `coverage(gr)[[1]]` (`get_coverage_from_blast()`, `dt_utils.R:416`) returns
   a vector of length **max(send) for that subject**, not the sequence length.
   `L` inside both switch-point functions is therefore data-dependent, and
   positions before the first `sstart` are 0.
2. Subjects with no surviving hit never enter the list → no cp → dropped later.
3. The self-hit rule is `gsub("_+", "", qaccver) != saccver`
   (`filter_blast3()`, `dt_utils.R:362`). Round-3 IDs are bare integers with no
   underscores, so it removes exactly the diagonal. The same rule applied to
   `id_start_end`-style names elsewhere can never match. The Python tool must
   encode a deliberate choice here, not inherit the ambiguity by accident.
4. `-evalue 1e-10` is passed to blastn while the filter tests `evalue < 1e-5`
   — the evalue predicate is already unreachable (§5).

#### Dependency

numpy must be added to `requirements.txt` and `conda/…/meta.yaml`. Python is
currently stdlib-only by policy (`CLAUDE.md`), but pure-Python accumulation is
~100× too slow at this scale; this is the one place where that policy costs
more than it buys.

### P1 — stop generating the volume in the first place

Cheap constant-factor wins, no change to results:

| Change | Effect |
|---|---|
| `-outfmt "6 qaccver saccver length sstart send"` | blastn currently formats 12 columns and `awk` discards 9. ~75 → ~35 bytes/hit; ~300 GB → ~90 GB of text |
| `-perc_identity 80` | moves the pident test into BLAST (note `>=` vs the current `>`) |
| drop the evalue column and predicate | dead code: `-evalue 1e-10` is already stricter than the `1e-5` filter |
| pipe directly, no intermediate file | removes the 41.5 GB write; `--keep-hits` (gzipped) behind `--debug` |
| test `-mt_mode 1` | query-based threading may scale past the ~25-core ceiling with 168k queries against a 1 GB db |

Measure before adopting as defaults (both bias which hits survive):

- **`-max_hsps 1`** — `filter_blast3()` keeps every HSP of a pair, so
  repeat-rich pairs contribute several intervals each.
- **`-max_target_seqs`** — currently 500,000, i.e. effectively unlimited.
  Capping biases coverage toward the best-scoring subjects.

### P2 — bound the compute, not just the parsing

P0+P1 still leave ~6–10 h of BLAST for this class. §2.1 says we can cut that
honestly.

**Query subsampling with per-subject rescaling.** The *subject* set must stay
complete (every element needs its own cp), but the *query* set is a sample and
BLAST cost is linear in it. Draw a random Q of N queries, accumulate as usual,
rescale coverage by N/Q before applying the absolute thresholds.
At Q = 20k (12 %): ~1.2 h of BLAST, ~3 × 10⁸ hits, median depth still ~430 —
a few percent error on window means.

**The caveat, stated plainly:** the thresholds are absolute copy counts, so
subjects whose *true* depth sits near the threshold get noisy under sampling
(true depth 20 → sampled ≈ 2.4, Poisson sd 1.55 → rescaled sd ≈ 13 against a
threshold of 20). Uniform subsampling alone would silently lose marginal
low-copy subfamilies.

**Adaptive second pass** fixes that: after pass 1, collect the subjects that
are unresolved or near-threshold, build a small BLAST db from **only those**,
and run the **full** query set against it. BLAST cost scales with db size, so a
10k-subject db costs ~6 % of a full run.

```
pass 1: 12 % queries × full db      ≈ 1.2 h
pass 2: 100 % queries × ~6 % db     ≈ 0.6 h
                                    ---------
                                    ≈ 1.8 h   (vs ~10 h today)
```

Exposed as a parameter (e.g. `--max_round3_queries`) that **defaults to
unlimited**, so small genomes stay bit-identical to today and only genomes
that currently fail change behaviour.

## 5. Validation plan

The 41.5 GB filtered file **still exists on disk**, so a new parser can be
validated against real data at real scale without repeating the 10 h BLAST.

1. **Unit equivalence.** Extend `tests/test_blast_reduce.R` (which already
   proves `awk` ≡ `filter_blast3()`) with a third leg: R `…coverage2` /
   `…coverage3` vs the Python implementation on the bundled `scaffold_1` data
   — assert identical cp vectors, not merely similar ones.
2. **Scale check.** Run `blast_cp.py` over the existing 41.5 GB file; confirm
   a sane cp table, bounded RSS, and wall time in the expected range.
3. **End-to-end.** `tests.sh` output must be unchanged with P0+P1 and with P2
   at its default (unlimited) setting.
4. **P2 quality check.** On a genome where the full run is feasible (e.g. the
   MuDR/Mutator class here, 11,260 copies), compare cp vectors and final GFF3
   between full and sampled runs; report how many elements are gained/lost.

## 6. Expected end state (CACTA upstream, this genome)

| Stage | BLAST text | On-disk intermediate | R/Python memory | Wall (this class-direction) |
|---|---|---|---|---|
| today (0.2.8) | ~300 GB | 41.5 GB | 477 GB peak → **crash** | ~10 h then fails |
| P0 | ~300 GB | none | ~5 GB | ~10 h + ~15 min |
| P0+P1 | ~90 GB | none | ~5 GB | ~6–8 h |
| P0+P1+P2 | ~11 GB | none | ~5 GB | ~1.5–2 h |

## 7. Staged plan

1. ~~**P0** — Python streaming parser + cp; R reduced to a thin caller; numpy
   added to the run deps. Fixes the crash.~~ **done, §10**
2. ~~**P1** — column pushdown, `-perc_identity`, drop the dead evalue filter,
   direct piping. Cheap, no semantic change.~~ **done, §10**
3. **P2** — query subsampling + adaptive second pass, default off.
4. Separately: the `…coverage3` index bug and the CAP3 splitting default (§8).

## 8. Adjacent findings (separate from the BLAST work)

### 8.1 `find_switch_point_from_blast_coverage3()` has an off-by-W index bug

`dt_utils.R:499` — and this is the function used for **CACTA**, the class that
fails here.

```r
swp <- seq(W, L-200, by = 1)
swp <- seq(1, L, by = 1)      # silently overwrites the line above
...
m1 <- sumsum_left/swp
m2 <- sumsum_right/(L - swp)
...
cp <- which.max(m12)
mcov1 <- m1[cp - W]           # m1/m2 are indexed by POSITION here
mcov2 <- m2[cp - W]
```

`m1`/`m2` are indexed by position, so the QC thresholds are evaluated 200 bp
before the detected switch point. In `…coverage2` the identical expression is
correct, because there `m1`/`m2` are indexed by window offset and
`cp = which.max(...) + W` — it looks copy-pasted. Additionally
`m2 <- sumsum_right/(L - swp)` is off by one (the mean of `cvrg[i..L]` needs
`L - i + 1`) and yields `Inf` at `i = L`, masked by the edge zeroing.

A Python reimplementation forces a decision: **replicate bug-for-bug first**
so P0 is provably equivalent, then fix in a separate change with a visible
before/after on real data.

### 8.2 CAP3 segfaults on unsplit large classes

```
working_dir/Class_II_Subclass_1_TIR_EnSpm_CACTA_upstream.part_001.fasta.cap.err:
  Segmentation fault (core dumped)
```

Both CACTA `.cap.aln` files are 0 bytes; the inputs were 1.78 GB each.
`--max_class_size` was not passed, so no splitting happened. That is why
Round 1 found only 75 elements genome-wide and why this genome depends
entirely on Round 3. A finite default for `--max_class_size`, or a size-based
auto-split, is warranted independently of everything above.

### 8.3 Other stages not yet stress-tested at this scale

Round 3's per-element TIR detection, Round 4, and `dante_tir_summary.R`
(mmseqs clustering) have never run on ~168k elements of one class in this
pipeline — no evidence of a problem, but no evidence of safety either. Re-run
run-000129 to completion before assuming Round 3 was the only wall.

## 9. Open questions

- Is per-subject **reservoir sampling** (bounded memory, §4/P0) wanted as a
  default, or only as a `--low-mem` escape hatch? The dense array is cheap on
  the target hardware.
- Should P2's near-threshold second-pass criterion be expressed in coverage
  units or as "cp not resolved"? The former is more principled, the latter is
  simpler to implement and audit.
- Does `-max_hsps 1` change any cp on the small test data? If not, it is a
  free ~2–3× reduction and belongs in P1 rather than "measure first".

## 10. Implementation notes (P0 + P1, landed)

What changed:

| File | Change |
|---|---|
| `blast_cp.py` | new: streams outfmt-6 → per-subject coverage profile → switch points → `id`/`cp` table |
| `dt_utils.R` | `run_blast_tir_analysis()` rewritten as a thin caller (`blastn \| blast_cp.py`, read the small table); `blast3_reduce_awk()` deleted; `swt_pt_method()` added; `DT_UTILS_DIR` / `blast_cp_script()` locate the helper in both a source checkout and `share/dante_tir/` |
| `dt_utils.R` (`round3`) | Round-3 outputs are now `*_upstream3.cp.tsv` / `*_downstream3.cp.tsv` |
| `tests/test_round3_cp.R` | replaces `tests/test_blast_reduce.R`; asserts R ≡ python switch points |
| `requirements.txt`, `conda/dante_tir/meta.yaml` | numpy added |

Decisions taken on the §9 open questions:

- **Reservoir sampling: not implemented.** The dense profile array is 4.2 GB
  for the worst class the pipeline has ever seen, on a machine that peaked at
  477 GB; a second sampling code path would have to be maintained and tested
  for equivalence to buy memory nobody is short of. `--reservoir` can be added
  if a genuinely memory-constrained run ever appears.
- **`-max_hsps 1`: not adopted.** It is a real ~2–3× reduction but it changes
  which intervals reach the profile, so it belongs with P2's sampling work
  where the quality effect is measured on a genome that can be run both ways —
  not smuggled into a change whose whole claim is "results are identical".
- **P2's second-pass criterion** stays open; it is a P2 decision.

### Measured on the run-000129 data

`tests/test_round3_cp.R` proves equivalence; `tests.sh short` was run at HEAD
and on the new code and produced byte-identical `DANTE_TIR_final.gff3`,
`DANTE_TIR_final.fasta`, `TIR_classification_summary.txt` and
`dante_tir_summary.R` outputs (14 TIR records both ways).

The stage that killed run-000129 — the 41.5 GB CACTA upstream table — was then
replayed through `blast_cp.py` (`--columns saccver,sstart,send`, the layout
that file already has):

| | |
|---|---|
| rows read | **2,538,889,291** |
| subjects with coverage | 167,284 of 168,012 |
| switch points found | **92,141** |
| peak RSS | **5.94 GB** |
| wall / CPU | 51 min / **23 min** — the rest is NFS read wait at ~14 MB/s; in production the rows arrive on a pipe and are never re-read |
| output | a **1.7 MB** table |

Sustained parse rate on that 3-column layout: ~1.8 M rows/s single-process.

One observation, not a regression: 13,753 of the 92,141 switch points land on
position 201, the leftmost position `m12[1:W] <- 0` admits. That is the R
algorithm's own edge behaviour (it is what `…coverage3` would have returned),
and it usually means the whole upstream window is inside a repeat with no flank
to step down from. Worth a look when the `…coverage3` bug in §8.1 is revisited.

Notes for whoever picks up P2:

- `filter_blast3()`, `get_coverage_from_blast()`, `swt_pt_function()` and
  `find_switch_point_from_blast_coverage2/3()` stay in `dt_utils.R` as the
  **reference implementations** the test asserts against. Production no longer
  calls them. Keep `swt_pt_function()` and `swt_pt_method()` in sync.
- `blast_cp.py --keep-hits` writes the surviving `saccver/sstart/send` rows
  (gzipped if the name ends in `.gz`) and `--input` reads a table back, so a
  profile can be recomputed without re-running BLAST. `run_blast_tir_analysis()`
  exposes it as `keep_hits =` but Round 3 does not pass it: on the CACTA class
  that file is the 41.5 GB the whole change exists to avoid.
- Restart behaviour: the cp table is the cache. A run that dies after Round 3's
  BLAST reuses it; a run that dies during it repeats it (the `.tmp` file is
  renamed only on success, so a partial table is never mistaken for a cached
  one).
- The parser is single-process by design: it sustains ~1.5 M rows/s on the
  6-column layout, roughly 20× faster than `blastn` produced hits in
  run-000129, so threading it would only add complexity.
