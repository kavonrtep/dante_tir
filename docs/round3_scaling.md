# Round 3 scaling: streaming BLAST → coverage profile → switch point

Design document for fixing the Round-3 self-BLAST blow-up on very large,
high-copy genomes. Written against `0.2.8`; supersedes the `awk` prefilter
introduced in `0.2.7` (commit `dcddcca`).

Status: **P0, P1 and P2 implemented** (§10, §11) and validated on run-000129
(§12). Keep this document updated as steps are completed.

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
3. ~~**P2** — query subsampling + adaptive second pass, default off.~~
   **done, §11**
4. Separately: the `…coverage3` index bug and the CAP3 splitting default (§8).

## 8. Adjacent findings (separate from the BLAST work)

### 8.1 `find_switch_point_from_blast_coverage3()` off-by-W bug — FIXED in 0.2.9

`dt_utils.R` — and this is the function used for **CACTA**, the class that
motivated all of the above. It read its quality thresholds at `cp - W`:

```r
m1 <- sumsum_left/swp          # indexed by POSITION
m2 <- sumsum_right/(L - swp)
cp <- which.max(m12)
mcov1 <- m1[cp - W]            # 200 bp before the switch point
mcov2 <- m2[cp - W]
```

`m1`/`m2` are indexed by position, so those thresholds tested a stretch 200 bp
away from the switch point they had just found. The expression was copied from
`…coverage2()`, where it is correct because there `m1`/`m2` are indexed by
window offset and `cp = which.max(...) + W`. Two further defects sat alongside
it: `swp <- seq(W, L-200, by = 1)` was overwritten on the next line, and
`m2 <- sumsum_right/(L - swp)` divided by one element too few (`Inf` at
`i = L`, masked by the edge zeroing).

**Why it mattered.** At the smallest admissible cp of 201, `m1[cp - W]` is
`m1[1]` — the coverage of a *single base* — tested against "flank mean < 3".
Almost anything passes that, which is why 13,753 of run-000129's CACTA calls
piled onto position 201.

**Measured effect of the fix**, recomputing the exact cp table from that run's
own 41.5 GB Round-3 output (all 168,012 queries, 2.54e9 hits):

| | buggy | fixed |
|---|---|---|
| switch points | 92,141 | **75,135** (−18.5 %) |
| calls at position 201 | 13,753 (14.9 % of calls) | **1,572 (2.1 %)** |
| median cp | 1,859 | 2,043 |
| lost / gained | — | 17,323 lost, 317 gained |
| positions of calls kept by both | — | **99.9 % identical, 100 % within 10 bp** |

So the fix is almost purely a *rejection* change: the boundaries it keeps sit
where they always did, and the ~19 % it drops are calls the flank test should
never have accepted — 12,181 of them the position-201 artefact. On
`tests/data/short` the element count moved 14 → 15, all 14 originals unchanged.

`blast_cp.py`'s `cp_cumsum()` carries the same fix, and
`tests/test_round3_cp.R` Part 6 pins the semantics: a profile whose flank is
dirty *at* the switch point but clean 200 bp earlier must be rejected, a clean
one must still be found, and the thresholds must hold at the returned position.

### 8.2 CAP3 segfaults on unsplit large classes — DIAGNOSED AND GUARDED

```
working_dir/Class_II_Subclass_1_TIR_EnSpm_CACTA_upstream.part_001.fasta.cap.err:
  Segmentation fault (core dumped)
```

Both CACTA `.cap.aln` files are 0 bytes; the inputs were 1.78 GB each
(15,112,884 reads, 1.51 Gbp). `--max_class_size` was not passed, so no
splitting happened, and Round 1 found only 75 elements genome-wide.

**The limit, measured** (cap3 10.2011, synthetic random reads, so it reproduces
away from the real data):

| input | result |
|---|---|
| 10.60 M × 100 bp = 1.060 Gbp | assembles (still running at 200 s) |
| 10.73 M × 100 bp = 1.073 Gbp | **SEGFAULT after 26 s** |
| 10.80 M × 100 bp = 1.080 Gbp | **SEGFAULT after 27 s** |
| 12.00 M × 100 bp = 1.200 Gbp | **SEGFAULT after 30 s** |
| 15.20 M × 100 bp = 1.520 Gbp | **SEGFAULT after 54 s** |
| 5.40 M × **200 bp** = 1.080 Gbp | **SEGFAULT after 25 s** |

The last row is the control: the same base count with *half* the reads fails
identically, so the limit is **total bases, not read count** — and content is
irrelevant, since random sequence fails exactly like real fragments. The
threshold sits between 1.060 and 1.073 Gbp, matching a signed 32-bit index over
the concatenated forward + reverse-complement sequence:
2 × (1.073e9 + 10.7e6) = 2.17e9 > 2³¹, while 2 × (1.060e9 + 10.6e6) = 2.14e9 <
2³¹. It always dies within a minute — an index overflow, not slow exhaustion.

For the pipeline: fragmentation yields ~8,993 bases per 6,300 bp region
(168,012 regions → 1.51 Gbp), so CAP3 breaks above **~118,000 regions per
part**.

**What was fixed.** `cap3assembly()` never checked CAP3's exit status: it wrote
the empty stdout to `.cap.aln`, returned that path as a success, and the
`os.path.exists` guard then treated the zero-byte file as finished work on
every re-run. It now refuses oversized inputs up front, checks the exit status,
leaves no zero-byte `.cap.aln`, and `dante_tir.py` reports which classes lost
their assembly. `tests/test_cap3_guard.py` stubs `cap3` so CI needs neither an
assembler nor a gigabase of input.

**The default was deliberately left off.** Enabling `--max_class_size` changes
results even when it splits nothing — on `tests/data/short` (31 CACTA copies)
it moved the count from 14 to 9 with a single part per class, because the
grouping path perturbs the random stream that jitters fragmentation. That is
not evidence that splitting is harmful: the same dataset gives **14, 10 and 10
records for `--seed` 42, 1 and 7 with no splitting at all**, Round 1's BEAST
being an MCMC. On a 31-copy class the pipeline's own seed sensitivity is
±30 %, so the flag stays opt-in rather than silently perturbing every existing
run.

### 8.2.1 What CAP3 actually costs, and how to split (measured)

Real CACTA fragments from run-000129, one CAP3 per part, pipeline options
(`-o 40 -p 80 -x cap -w`), single core:

| part | copies | fragments | wall | peak RSS | contigs | fragments assembled |
|---|---|---|---|---|---|---|
| clustered | 1,000 | 89,843 | 571 s | 0.98 GB | 4,235 | 63.8 % |
| clustered | 2,000 | 179,524 | 1,497 s | 2.25 GB | 8,232 | 64.7 % |
| clustered | 4,000 | 359,524 | 3,853 s | 5.30 GB | 14,862 | 66.0 % |
| clustered | 8,000 | 719,217 | 10,359 s | 12.45 GB | 24,149 | 67.4 % |
| positional | 4,000 | 352,167 | 1,136 s | 1.89 GB | 23,797 | 59.3 % |
| random | 4,000 | 359,784 | 1,146 s | 1.86 GB | 21,043 | 53.3 % |

**Cost grows as ~n^1.4, not n².** Doubling a part multiplies wall time by 2.62,
2.57 and 2.69 across the three doublings measured, and memory by ~2.3. Projected: 8,000 copies ≈ 2.8 h / 12 GB, 16,000 ≈ 7 h /
29 GB, and the 118,000-copy crash limit ≈ 4.6 days / ~300 GB — so **runtime,
not the 32-bit overflow, is what really caps part size**, but not as brutally
as a quadratic would.

**Cluster-coherent parts are worth their cost.** At the same 4,000 copies,
grouping by mmseqs cluster costs 3.4× the wall time and 2.8× the memory of a
positional or random split — and that is the assembly actually happening:
it yields *half* as many contigs (14,862 vs 23,797) holding *twice* as many
fragments each (16.0 vs 8.8), and assembles 66 % of the input against 53 % for
a random split. Round 1 reads element boundaries off contig coverage, so deeper
contigs from genuinely related copies are exactly the signal it needs. A random
split of the same size does the cheap thing and learns the least.

That is what `--max_class_size` already does — `group_sequences_by_clusters()`
keeps mmseqs clusters together and only splits a cluster that exceeds the
threshold on its own. The cluster structure of this class makes that
worthwhile: 168,012 domains in 85,606 clusters, with the largest at 14,825,
12,471 and 10,565 members and 17 clusters ≥ 1,000 covering 42 % of all copies —
but also 82,745 singleton clusters covering 49 % of copies, which will be
merged into mixed groups where little will assemble.

**Sizing guidance.** A part of 8,000–16,000 copies keeps each CAP3 to ~3–7 h
and 12–29 GB while leaving the three largest clusters nearly intact. For CACTA
that is 11–21 parts per direction.

**But watch concurrency memory.** `dante_tir.py` runs CAP3 through
`Pool(processes=args.cpu)` over every (class, part, direction) FASTA, so on a
96-core node it would start dozens of these at once: 21 parts × 2 directions at
12 GB each is ~500 GB resident. That fits 768 GB but not much else, and nothing
currently bounds it. Splitting a class therefore needs a matching cap on how
many CAP3 jobs run concurrently — a memory budget rather than a core count.

**For run-000129** the practical setting is `--max_class_size` in the 8k–16k
range with bounded CAP3 concurrency, not the ~118,000 the crash limit allows.

**For run-000129**: `--max_class_size 100000` splits CACTA into two parts of
~0.76 Gbp, under the limit. Runtime is the open question — MuDR's 1.01 M
fragments took ~1.9 h and CAP3 is superlinear, so a 7.6 M-fragment part may be
impractical. The trade-off between assembly context and runtime is unmeasured.

### 8.3 Other stages not yet stress-tested at this scale — ANSWERED

Round 3's per-element TIR detection, Round 4, and the mmseqs clustering had
never run on ~168k elements of one class. The re-run in §12 took them all the
way through in three minutes, so Round 3 was indeed the only wall in this
chain — but only after one more fix: both places that extract TIR sequences
did `readDNAStringSet(genome)`, which on this 94.3 Gbp assembly asks for
~94 GB of RAM. They now use indexed access (`genome_fa_handle`,
`tests/test_genome_access.R`). Neither call had ever been reached before,
because Round 3 failed first.

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

## 11. Implementation notes (P2, landed)

`--max_round3_queries N` (`dante_tir.py` → `detect_tirs.R` → `round3()` →
`run_blast_tir_analysis()`) caps how many sequences are used as queries in the
Round-3 self-BLAST. **Default 0 = no cap**, i.e. the exact path; a class with
fewer copies than the cap is unaffected, so enabling it changes nothing except
for the classes that are actually too big.

Above the cap:

1. **Pass 1** BLASTs a random Q-of-N sample of the queries (seeded from
   `--seed`, so it is reproducible) against the full database, and scales the
   coverage profiles by N/Q — an unbiased estimate of the full coverage, so
   the absolute thresholds keep their meaning.
2. `blast_cp.py` reports, next to each cp, the **support**: the element-side
   coverage behind the call in *observed* (unscaled) units, i.e. how many
   sampled relatives actually back it.
3. **Pass 2** re-resolves every subject with `support < min_support` (default
   30), plus every subject that produced no profile at all, by BLASTing the
   **full** query set against a database holding only those subjects. For them
   the answer is not an estimate: it is what the unsampled path would produce.
   It is cheap because BLAST cost scales with database size.

### Measured cost

BLAST cost for the CACTA class is linear in the query count (q500/q2000
samples against the full 168k database, 14 threads):

| queries | hits | CPU-min | CPU-min per query |
|---|---|---|---|
| 500 | 7,960,733 | 40.1 | 0.0801 |
| 2,000 | 32,100,791 | ~153 | 0.0765 |

Extrapolating: 168,012 queries ≈ 216 CPU-h, which matches the 13 h × ~25
effective cores the production run actually spent. Hit volume scales exactly
(4.03× hits for 4× queries).

Combined with the measured support distribution of the CACTA class (median
element-side support **5,881**, p5 = 30 — the oversaturation from §2.1 seen
directly), the predicted pass-2 fraction and total cost are:

| Q | scale | pass 2 (at min_support 30) | total cost vs exact |
|---|---|---|---|
| 5,000 | 33.6 | 19.6 % | ~23 % |
| 10,000 | 16.8 | 14.1 % | ~20 % |
| **20,000** | **8.4** | **11.2 %** | **~23 %** |
| 40,000 | 4.2 | 8.6 % | ~33 % |

The predicted cost column turned out to be **pessimistic** — see §12: pass 2
is far cheaper than this model assumes, because the subjects it re-resolves
are precisely the sparse ones.

### Measured quality

Measured on MuDR/Mutator upstream from the same run (11,260 copies — small
enough to run exactly, 19.3 min), sampling 2,000 queries (17.8 %):

| | exact | sampled |
|---|---|---|
| switch points | 8,400 | 8,323 |
| called in both | — | 8,175 |
| lost / gained | — | 225 / 148 |
| \|Δcp\| = 0 | — | 61.5 % |
| \|Δcp\| ≤ 10 bp | — | 96.4 % |
| \|Δcp\| ≤ 50 bp | — | 99.3 % |

Positions barely move (p95 = 9 bp, and Round 3 searches ±200 bp around cp), but
**~4 % of subjects change between called and not-called**.

That was the only quality figure available before the class that motivates
sampling could be run. It is *not* representative — see §12, where CACTA turns
out to be about three times worse, not better.

**The gate has a ceiling.** Calibrating `min_support` on the MuDR data
(pass-1 estimate vs the exact answer, by observed support):

| support | n | called↔NA flip | \|Δcp\| ≤ 10 bp |
|---|---|---|---|
| 1–2 | 595 | 37 % | 52 % |
| 2–5 | 798 | 29 % | 69 % |
| 5–10 | 1,092 | 15 % | 84 % |
| 10–20 | 1,464 | 10 % | 89 % |
| 20–30 | 1,002 | 8 % | 92 % |
| 30–60 | 1,589 | 8 % | 91 % |
| 60–150 | 3,255 | 4 % | 94 % |
| ≥150 | 1,088 | 9 % | 95 % |

The flip rate falls steeply up to support ≈ 10 and then **plateaus at 4–9 %
instead of going to zero** — raising `min_support` does not buy the rest.
That residue is not element-side noise, which is what the gate measures. It is
almost certainly the *flank*-side test: `mcov1 < 3` in scaled units means the
sampled flank coverage must be under 3·Q/N ≈ 0.5, a quantity estimated from a
handful of alignments, so subjects whose true flank coverage sits near the
threshold are close to a coin flip however deep the element side is.

**The obvious next refinement** is therefore to have `blast_cp.py` report
`mcov1` as well and send a subject to pass 2 when *either* threshold is within
a noise band of its boundary, rather than gating on element-side support
alone. That was not done here: it adds a second calibration nobody has data
for on the target class, and the current gate is already conservative.

### Reproducibility

Sampling introduces the only randomness in Round 3, so it is pinned end to end.
A sampled run is a deterministic function of (input, `--seed`,
`--max_round3_queries`, `min_support`):

- **The draw depends on `--seed` and the input alone.** `with_seed()` fixes the
  generator explicitly (`Mersenne-Twister` / `Inversion` / `Rejection`) rather
  than inheriting it. Both leaks this closes are real: R changed the default
  sampler in 3.6, and `parallel` code commonly switches the process to
  `L'Ecuyer-CMRG` — under which `set.seed(42)` demonstrably draws a *different*
  sample (verified, `tests/test_round3_cp.R` Part 5). The caller's RNG stream
  and generator are restored afterwards, so nothing downstream shifts either.
- **The sample is recorded in the log.** Each sampled class logs its seed and a
  fingerprint of the drawn indices (`n`, sum, first, last) to `log/stderr.txt`,
  which survives without `--debug` — the sampled query FASTA itself lives in
  `working_dir` and does not. Two runs can be confirmed to have used the same
  queries from the logs alone.
- **Thread count does not change results.** The profile is an accumulation of
  integer counts, which is order-independent, so however `blastn` interleaves
  its output the profile and the switch points are the same. Verified on real
  data: 200 MuDR queries against the full 11,260-sequence database at
  `-num_threads 1` and `-num_threads 8` produce byte-identical cp tables
  (172 s vs 34 s).
- **Pass 2 membership is a deterministic function of pass 1**, and each pass
  caches its own cp table, so a restart resumes rather than redraws.

Verified end to end in Part 5: re-running with the same seed reproduces the cp
table byte for byte; a different seed draws different queries but still
describes the same subjects.

### What P2 does not do

On a mid-size class the two-pass scheme is close to a wash: the same
experiment on MuDR took 17.4 min sampled vs 19.3 min exact (10 %), because at
a 71 MB database BLAST's fixed costs dominate and 47 % of subjects fell below
`min_support` anyway. The win is specific to genuinely oversaturated
mega-families, which is exactly the population the cap selects.

## 12. What the production re-run showed (run-000129, 2026-07-26)

Rounds 1–4 re-run on run-000129's own working_dir (`dev_scripts/rerun_from_working_dir.sh`,
96 threads, `--max_round3_queries 20000`). **It completed**, producing 7,831
elements where the original run produced none: 7,442 EnSpm/CACTA, 293 hAT,
89 MuDR/Mutator, 7 PIF/Harbinger. Round 3 contributed 5,781 of them, Round 4
another 1,970, Rounds 1–2 just 80.

### The model was right about volume, wrong about cost

| | predicted | actual |
|---|---|---|
| pass-1 hits (CACTA upstream) | 302.2 M | **302.5 M** |
| pass-2 fraction | 11.2 % | **11.6 %** (19,445 subjects) |

But the cost model assumed pass 2 costs about as much as pass 1 (the full query
set against an 11.6 % database). Measured: **pass 1 63 min, pass 2 5 min**.
Pass 2 is ~12× cheaper than modelled because its subjects are the *sparse*
ones — 8 M hits against pass 1's 302 M — and blastn time follows hit count,
not database size. Sampling is therefore **~7.8× faster than exact, not the
~4.3× predicted above**.

Wall clock, Round-3 start to finish (4 h 11 m):

| phase | wall |
|---|---|
| CACTA upstream (pass 1 / pass 2) | 68 min (63 + 5) |
| CACTA downstream (pass 1 / pass 2) | 141 min (133 + 8) |
| hAT + MuDR + PIF, exact (below the cap) | 39 min |
| Round 4 + mmseqs clustering + final extraction | 3 min |

Downstream cost twice upstream because it genuinely carried 39 % more hits
(419 M kept vs 302 M).

### §8.3 answered: no further scaling walls

Round 4, the mmseqs clustering of TIR sequences and the final sequence
extraction all completed at 168k-element scale, in three minutes. That required
the FaFile change (`genome_fa_handle`) — on this 94.3 Gbp assembly the previous
`readDNAStringSet(genome)` would have asked for ~94 GB twice.

### Quality on CACTA — worse than the MuDR stand-in, not better

Sampled vs the exact upstream switch points (recomputed from the original
run's own 41.5 GB Round-3 output, all 168,012 queries):

| | |
|---|---|
| called in both | 86,492 |
| lost / gained | 5,649 / 7,078 |
| **churn** | **6.1 % of exact calls lost, 7.7 % gained** |
| \|Δcp\| = 0 / ≤10 bp / ≤50 bp / ≤200 bp | 60.1 % / 90.1 % / 94.8 % / 98.3 % |
| p95 / p99 / max \|Δcp\| | 53 bp / 307 bp / 4,954 bp |

§11 predicted the opposite — that CACTA's deeper support would leave it
*better* determined than MuDR. It is roughly three times worse. The likely
reason is the method, not the depth: CACTA uses `cp_cumsum`, whose argmax over
global cumulative means is far more sensitive to profile shape than MuDR's
local 200 bp windows, and it already piles 13,753 exact calls onto position
201 (§10).

**Element level is much better than cp level.** Of the 7,442 CACTA elements,
6,109 join to an exact upstream switch point:

| | elements | |
|---|---|---|
| cp identical to exact | 2,753 | 37.0 % |
| within the ±200 bp search window | 3,337 | 44.8 % |
| moved beyond that window | 19 | **0.3 %** |
| rests on a cp the exact run never produced | 301 | 4.0 % |

Median \|Δcp\| = 1 bp, p90 = 3 bp, p99 = 78 bp. Elements are much less affected
than switch points because an element only survives if TIR *and* TSD detection
succeed, which needs a well-defined boundary; the marginal boundary-sitters
that flip in the cp table rarely become elements. What this cannot measure is
the other direction — elements the exact run would have found and sampling
missed — which needs the exact run.

### Sampling has a churn floor

Flip rate against how much data a subject actually had (CACTA, pass-1-resolved):

| observed support | n | called↔NA flip | moves >200 bp |
|---|---|---|---|
| 30–60 | 4,798 | 10.9 % | 10.7 % |
| 60–120 | 9,218 | 13.6 % | 4.6 % |
| 120–250 | 10,704 | 14.2 % | 2.3 % |
| 250–500 | 17,397 | 10.8 % | 1.2 % |
| 500–1,000 | 53,099 | 9.5 % | 0.1 % |
| 1,000–2,000 | 49,425 | 4.6 % | 0.0 % |
| ≥2,000 | 3,926 | 5.3 % | 0.0 % |

Two different behaviours, and the distinction drives the recommendation:

- **Position errors are cured by more data** — "moves >200 bp" falls from
  10.7 % to ~0 % as support grows. A larger sample fixes these.
- **Called/not-called flips are not** — they fall from ~14 % to ~5 % and then
  flatten, even for subjects with thousands of supporting relatives. That is
  the signature of subjects sitting *on* the decision boundary, where any
  perturbation tips the result. A bigger sample shrinks the perturbation but
  never removes it.

(An attempt to extrapolate churn to other Q from this curve was discarded: it
bottoms out near 6 % even at Q = N, which is definitionally 0, because observed
support alone does not capture the large-Q regime. Only the two qualitative
behaviours above are supported by the data.)

### Pass 2 is exact to 99.93 %, not 100 %

Of the 18,661 pass-2 subjects that appear in both tables, **18,648 (99.93 %)
match the exact answer**. The 13 that differ are the small-database e-value
effect flagged when the two-pass scheme was designed: pass 2 searches a
database ~11 % the size of the full one, so identical alignments get slightly
different e-values and a marginal hit can cross `-evalue 1e-10`. The synthetic
test in `tests/test_round3_cp.R` Part 4 asserts exact equality and passes; at
production scale the claim needs this 0.07 % qualifier.

### Recommendation, revised

Sampling is a screening tool, not the default:

- **For a production library on an extreme genome, run exact**
  (`--max_round3_queries 0`). That is now possible — it is what P0 bought — and
  on this genome costs ~28 h for Round 3 (extrapolating the measured pass-1
  times: ~8.8 h upstream, ~18.6 h downstream), against 4 h sampled.
- **Use sampling to iterate**: parameter sweeps, pipeline debugging, or a first
  look at a new genome, where ~7 % switch-point churn and ~4 % of elements
  resting on boundaries the exact path would not produce are acceptable.
- **A larger sample does not buy away the churn** — it buys accurate positions.
  If churn is what matters, exact is the answer, not `--max_round3_queries
  40000`.


## 13. The exact run (run-000129, 2026-07-27)

Rounds 1–4 re-run with `--max_round3_queries 0` — no sampling — on the code
including the §8.1 fix, with direction-level and chunk-level concurrency.

**It completed in 10 h 02 m** (05:36 → 15:38) against ~28 h for the unchunked
path, producing 7,051 elements: 6,653 EnSpm/CACTA, 293 hAT, 98 MuDR/Mutator,
7 PIF/Harbinger. Rounds contributed 81 / 8 / 5,140 / 1,822.

Round 3 ran 05:49 → 15:34. CACTA upstream took 4 h 16 m and downstream 9 h 29 m
— a 2.2× ratio, matching the 2.1× seen in the sampled run, and confirming that
downstream carries disproportionately more work on this genome. The log shows
the intended layout: `2 directions x 4 query chunks, 12 threads each
(8 searches)`, followed by `merged 4 profiles`.

### The switch points reproduce exactly

The CACTA upstream cp table from this run is **identical, line for line, to the
table computed independently** from the *original* run's 41.5 GB Round-3 output
(`blast_cp.py` over that file, same method): all 167,284 subjects, all 75,135
switch points, same `cp` on every row.

Those two numbers come from different BLAST runs — different day, different
invocation (six columns plus `-perc_identity` versus the old twelve-column awk
path), chunked-and-merged versus a single stream, different machine. 433 of the
167,284 `support` values differ in the third decimal (e.g. 7272.184 vs
7272.237), which is the two BLAST runs reporting marginally different hit sets;
not one of them moved a switch point.

That is the end-to-end check for everything in §10–§11 at once: the new BLAST
invocation, the chunk merge, the streaming parser and the fixed `cumsum`
method.

### Exact versus sampled

| class | exact | sampled (§12) |
|---|---|---|
| EnSpm/CACTA | 6,653 | 7,442 |
| hAT | 293 | 293 |
| MuDR/Mutator | 98 | 89 |
| PIF/Harbinger | 7 | 7 |
| **total** | **7,051** | 7,831 |

The comparison is **confounded** and should not be read as "sampling found 789
more CACTA elements": the sampled run predates the §8.1 fix, which by itself
removed 18.5 % of upstream switch points. hAT and PIF are identical across the
two runs — both were below the sampling cap and therefore exact in both, so
they also demonstrate reproducibility across two independent runs.

Round 1 contributed 81 of 7,051 elements, unchanged, because this rerun reuses
the original CAP3 output — which is empty for CACTA (§8.2). Closing that gap
needs a full `dante_tir.py` run, which is what the `--max_class_size` default
introduced in 0.3.0 makes safe.
