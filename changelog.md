## 0.3.0 — 2026-07-27

The release that makes DANTE_TIR finish on very large, high-copy genomes. On an
89 Gb assembly with 168,012 EnSpm/CACTA copies, Round 3 previously aborted after
13 h with `long vectors not supported yet`, having produced nothing; it now
completes in about 10 h and the pipeline yields a full element set. Most of the
work below came out of that one run, and every number quoted was measured on it.

Two deliberate changes to results, both explained in their entries: the
EnSpm/CACTA switch-point method had its quality thresholds read at the wrong
position, and fragmentation is now seeded per sequence rather than from one
global stream. `--max_class_size` is also on by default now. Small genomes are
largely unaffected — on the 30 Mb `long` test set, every element found by 0.2.8
is found again at identical coordinates, plus one more.

### Round 3 no longer materialises the self-BLAST table

- The `awk` reduction added in 0.2.7 kept the Round-3 table off the heap but
  could not shrink it enough: on an 89 Gb genome with 168k EnSpm/CACTA domains
  the *filtered* table was still 41.5 GB / 2.5e9 rows and `read.table` aborted
  with `long vectors not supported yet`. The row count is intrinsic to an
  all-vs-all self-BLAST, so no per-row filter can fix it.

- `run_blast_tir_analysis` now streams `blastn` straight into `blast_cp.py`,
  which applies `filter_blast3`'s predicates, folds each surviving hit into a
  per-subject coverage profile and computes the switch points, returning only
  the `id`/`cp` table Round 3 actually consumes. Nothing proportional to the
  hit count is ever stored: memory is a dense int32 profile array
  (`n_subjects x max_length`, 4.2 GB for the class above) and the BLAST stream
  itself is discarded as it is read. Round-3 outputs are now
  `*_upstream3.cp.tsv` / `*_downstream3.cp.tsv`.

- Ask `blastn` for only the six columns the predicates need (it was formatting
  twelve and `awk` discarded nine) and push the identity cut-off down into
  `-perc_identity`. `filter_blast3`'s evalue predicate was already unreachable
  behind `-evalue 1e-10` and is no longer evaluated.

- Results are unchanged: `tests/test_round3_cp.R` (replacing
  `tests/test_blast_reduce.R`) asserts the new path reproduces the R
  implementation's switch points exactly, across predicate edge cases, random
  coverage profiles for both switch-point methods, and a real self-BLAST.

### Round 3 uses the machine properly

- Round 3 runs the two directions concurrently. A single `blastn` saturates at
  roughly a quarter of the threads it is given (24.7 of 96 on run-000129, 6.1 of
  14 locally) and `-mt_mode 1` does not change that -- measured identical wall
  time for 7x the memory -- so the way to use a large machine is more searches,
  not more threads per search. Measured locally: two 7-thread searches sustain
  1.81 queries/s against 1.28 for one 14-thread search, 1.42x on the same
  cores. The upstream and downstream searches now run side by side with the
  threads split between them. Output is unaffected (`tests.sh short` is
  byte-identical) because the searches write separate files and the coverage
  profiles are order-independent integer accumulations.

- Each Round-3 search is split into concurrent query chunks whose coverage
  profiles are summed. `blast_cp.py` gained `--dump-profile` and `--merge`; the
  profile is a difference array, so summing chunk profiles is exactly the
  arithmetic of one search over all the queries, and profile files carry a
  header pinning the subject set so profiles from different databases cannot be
  added together. Measured on the CACTA database, four 3-thread searches
  processed the same queries **3.0x faster than one 14-thread search, on fewer
  threads**. `blast_chunk_count()` sizes the split so no search is asked to
  scale past the point where it stops paying (12 threads), capped at 8 chunks
  because each costs one coverage profile in RAM and on disk.

### Optional query sampling, off by default

- New `--max_round3_queries N` (default 0 = off, i.e. the exact path). Above
  that many copies, a class runs its Round-3 self-BLAST on a random N-query
  sample with the coverage scaled back up, and every subject whose sampled
  support falls below `min_support` (default 30) is then re-resolved exactly
  against the full query set in a cheap second pass over just those subjects.
  Validated on run-000129: at `--max_round3_queries 20000` that genome's
  Round 3 finished in 4 h and the pipeline produced 7,831 elements where it had
  previously produced none, ~7.8x faster than the exact path. Sampling is not
  free, and on the CACTA class it costs more than the earlier MuDR measurement
  suggested: 6.1% of switch points lost and 7.7% gained, though at the element
  level 99.7% of calls land within the +/-200 bp window Round 3 searches and
  4% of elements rest on a boundary the exact run would not have produced.
  More sampling buys accurate positions but does not remove the churn, so the
  exact path (the default) remains the recommendation for a production library.
  Numbers, calibration and the limits of the support gate are in
  `docs/round3_scaling.md` 11-12.

- Sampled runs are deterministic: the draw is a function of `--seed` and the
  input alone. The generator is pinned explicitly rather than inherited, which
  closes a real hole -- under `RNGkind("L'Ecuyer-CMRG")`, which `parallel` code
  commonly sets, `set.seed()` draws a different sample -- and the caller's RNG
  stream is left untouched. Each sampled class logs its seed and a fingerprint
  of the draw to `log/stderr.txt` (which survives without `--debug`, unlike the
  sampled query FASTA). Thread count does not affect results either: the
  coverage profile is an order-independent integer accumulation, verified
  byte-identical at `-num_threads 1` vs `8` on real data.

### Round 1: CAP3 is bounded, and its failures are visible

- `--max_class_size` is on by default (10,000 sequences), now that seeding
  fragmentation per sequence makes splitting a genuine no-op below the
  threshold -- verified byte-identical on `tests/data/short`, whose classes sit
  far below it. Splitting a large class is not only what keeps CAP3 inside the
  ~1.07 Gbp it can index; it is also *faster*, because CAP3 costs ~n^1.4, so
  two halves are cheaper than the whole (measured: 1,000 copies 571 s, 8,000
  copies 10,359 s). Grouping follows mmseqs2 clusters, so related copies stay
  together -- at 4,000 copies that yields half as many contigs holding twice as
  many fragments each, and 66% of the input assembled against 53% for a random
  split. `--max_class_size 0` restores the old unsplit behaviour.

- New `--cap3_max_memory` bounds concurrent CAP3 assemblies by memory rather
  than core count. CAP3 needs ~`0.067 * Mbp^1.22` GB (a fit reproducing four
  measured points within 1%), so a class split into 21 parts of ~12.5 GB would
  have asked for ~1.2 TB when every core started one. The pool is now sized to
  fit a budget defaulting to 60% of detected memory (cgroup-aware), and the
  decision is logged.

- CAP3 failures are no longer silent. `cap3assembly` ignored CAP3's exit
  status, so a crash produced a zero-byte `.cap.aln` that the pipeline accepted
  as a finished assembly -- and that the `os.path.exists` guard then reused on
  every re-run. On run-000129 this cost the 168,012-copy EnSpm/CACTA class its
  entire Round 1 without a word in the log. CAP3's limit is now measured: it
  segfaults once one input exceeds ~1.07 Gbp (between 1.060 and 1.073 Gbp,
  driven by total bases rather than read count or content), which for 6300 bp
  flanking regions is ~118,000 sequences per part. `cap3assembly` refuses such
  inputs up front, checks the exit status, leaves no bogus `.cap.aln`, and
  `dante_tir.py` reports which classes lost their assembly.
  `--max_class_size` remains off by default: enabling it perturbs the random
  fragmentation and so changes results even when it splits nothing, and on
  small inputs that shift is within the pipeline's own seed sensitivity
  (`tests/data/short` gives 14/10/10 records for seeds 42/1/7).

- Fragmentation is seeded per sequence. `dict_fasta_to_dict_fragments` drew its
  jitter from one global RNG, so a region's fragments depended on how many draws
  had happened before it -- that is, on how many other sequences were processed
  first. Regrouping the input therefore changed the fragments even when it split
  nothing: `--max_class_size 10000` on `tests/data/short` moved the result from
  14 records to 9 with one part per class, which is why the flag could not be
  enabled by default. Each sequence now gets its own generator seeded from
  (`--seed`, direction, sequence id) via blake2b -- not `hash()`, which Python
  randomises per process. `--max_class_size` above every class size is now
  byte-identical to omitting it. **Round-1 results shift once** as a
  consequence (`tests/data/short`: 15 records to 11), within the band the seed
  alone already spans on that dataset (14/10/10 for seeds 42/1/7).
  `--seed` now also controls the Python stage, which it did not before.
  `tests/test_fragmentation.py` asserts order-, subset- and hash-seed
  invariance plus seed sensitivity, salt separation and the unchanged
  fragmentation scheme.

### Correctness fixes

- Fix `find_switch_point_from_blast_coverage3` (the EnSpm/CACTA method), which
  read its quality thresholds at `cp - W` -- 200 bp away from the switch point
  it had just found. `m1`/`m2` there are indexed by position, so the test
  applied to the wrong stretch of sequence; the expression had been copied from
  `...coverage2`, where it is correct because those vectors are indexed by
  window offset. At the smallest admissible switch point, 201, it compared the
  coverage of a *single base* against "flank mean < 3", so almost anything
  passed. Recomputed on run-000129's own 41.5 GB Round-3 output, the fix drops
  switch points from 92,141 to 75,135 (-18.5%) and the pile-up on position 201
  from 13,753 calls to 1,572, while 99.9% of the calls kept by both land on
  exactly the same position -- it rejects, it does not move boundaries. Two
  neighbouring defects went with it: a dead `swp` assignment and an
  off-by-one in `m2`'s divisor. `blast_cp.py` carries the same fix and
  `tests/test_round3_cp.R` Part 6 pins the semantics. **This changes results
  for EnSpm/CACTA**; other superfamilies use `...coverage2` and are unaffected.

- Fix a hang on large inputs: `cluster_aa_sequences_mmseqs2` and
  `make_blast_db` ran their child processes with `subprocess.check_call(...,
  stdout=PIPE, stderr=PIPE)`. `check_call` waits without ever reading those
  pipes, so the child blocks forever once it writes past the 64 KB pipe buffer
  -- `mmseqs easy-cluster` reached that on real data and deadlocked the whole
  run in `pipe_write`. mmseqs now goes through `subprocess.run` (which drains
  while waiting) and reports its stderr when clustering fails; `makeblastdb`
  discards its output via `DEVNULL`. Added `tests/test_subprocess_pipes.py`
  (in `tests/unit.sh`) to keep the pattern from coming back.

### Memory

- Stop loading the whole genome into R. `detect_tirs.R` and
  `cluster_tir_sequences` both did `readDNAStringSet(genome)` to pull out a few
  thousand short TIR ranges; a `DNAStringSet` costs ~1 byte per base, so for
  run-000129's 94.3 Gbp assembly that is ~94 GB of RAM, twice. Neither call had
  ever been reached on that genome because Round 3 failed first. Both now use
  `Rsamtools::FaFile` random access through the `.fai`, which reads only the
  requested ranges (`genome_fa_handle`, which indexes the genome if needed and
  fails loudly if it cannot). Output is unchanged: `tests/test_genome_access.R`
  asserts the two paths return identical sequences including
  reverse-complementing of minus-strand ranges, and `tests.sh short` produces a
  byte-identical `DANTE_TIR_final.fasta`. `bioconductor-rsamtools` is now a
  declared dependency (it arrived via BSgenome before). **The genome must be
  indexable**: the `.fai` is created on first use, so this matters only when
  the directory holding the genome is not writable, in which case the run stops
  with an error naming the file and the `samtools faidx` command that fixes it.
  0.2.8 would have loaded such a genome whole. See the README.

### Reproducibility

- Pin the record order of the amino-acid FASTA that mmseqs2 clusters, since
  `--max_class_size` splits classes along those clusters. Measured on
  run-000129's MuDR domains with the pipeline's own parameters, mmseqs2 is
  stable for a fixed input -- rerunning is identical and 1 thread matches 4 --
  but *order-sensitive*: shuffling the input alone moved 258 of 5,983 clusters
  and changed the cluster count to 6,016. The order the pipeline produces is
  the GFF3's, carried through dicts and lists, and no `set` sits in that chain;
  `tests/test_aa_fasta_order.py` builds the FASTA in separate interpreters
  under different `PYTHONHASHSEED` values and asserts the bytes match, so a
  future `set` cannot silently make split runs irreproducible.

### Dependencies

- `bioconductor-genomeinfodbdata` is declared explicitly. It is a transitive
  requirement of `GenomeInfoDb`, which every Bioconductor package here loads,
  and a fresh solve pulls it in — but the release workflow installs into an env
  that already holds `conda-build`, and that constrained solve left it out. The
  R stage then died at startup with *there is no package called
  'GenomeInfoDbData'*. Naming it removes the dependence on solver order.

- New runtime dependency: `numpy`.

## 0.2.8 — 2026-07-21

Round-3 memory fix plus a release guard.

- Round 3 no longer reads the full self-BLAST table into R. On large,
  high-copy genomes the all-vs-all self-BLAST can exceed hundreds of GB, and
  `read.table` aborted with `long vectors not supported yet` (R's 2^31-element
  limit) before any element was called. `run_blast_tir_analysis` now streams
  `blastn` through an `awk` filter (`blast3_reduce_awk`) that reproduces
  `filter_blast3`'s predicates and keeps only the `saccver/sstart/send` columns
  coverage needs, so the full table is never written to disk nor loaded into R.
  Results are unchanged — verified identical at the filter, coverage, switch-
  point, and final-GFF3 level. Round-3 BLAST outputs are now named
  `*_upstream3.filtered.tsv` / `*_downstream3.filtered.tsv`.
- Add `tests/test_blast_reduce.R` (run from `tests/smoke.sh`) proving the
  streamed reduction is identical to the previous read.table/filter path across
  predicate edge cases and a real self-BLAST.
- Release plumbing: add `dev_scripts/check_release_version.sh`, a guard that
  refuses to (re)release a version already published on the conda channel (and,
  with `--require-untagged`, one that already has a local git tag). Wired into
  `conda-release.yml` as a fail-fast preflight so a duplicate version bump can
  no longer waste a full build and then fail at upload.

## 0.2.7 — 2026-07-15

Memory-usage fix plus CI plumbing.

- `extract_flanking_regions` no longer loads the whole genome into RAM
  (previously ~genome_bp, causing OOM on large assemblies such as a 90 Gbp
  genome). It now takes the FASTA path and streams one sequence at a time
  via a new `fasta_record_generator`, so peak memory is the largest single
  sequence rather than the whole assembly. Output is byte-identical to the
  previous dict-based logic.
- Add `tests/test_extract_flanking_regions.py` (plus `tests/unit.sh` and a
  `unit` level in `tests.sh`) verifying the streaming output matches the
  original logic across both strands, window clamping, and line-wrapped
  multi-sequence FASTA.
- CI: switch the conda-release workflow from the removed `conda mambabuild`
  to the standalone `conda-build` executable.

## 0.2.6 — 2026-06-30

- `dante_tir_summary.R`: reorder the TIR consensus logo figure in the HTML
  report so the 5' TIR is shown first (top panel) and the 3' TIR second
  (bottom panel).

## 0.2.5 — 2026-04-28

CI / release-plumbing migration with a small batch of robustness fixes.

CI / release plumbing:

- Add in-repo conda recipe under `conda/dante_tir/`.
- Add tag-driven GitHub Actions release pipeline under
  `.github/workflows/conda-release.yml`; tests on every push/PR via
  `.github/workflows/tests.yml`.
- Refactor `tests.sh` into a `{smoke|short|long|all}` dispatcher; add
  tiered scripts under `tests/`.
- Commit tiered datasets under `tests/data/{smoke,short,long}/` derived
  from the first 2 / 20 / 30 Mb of `tiny_pea` Chr1.
- Move old developer harnesses (`tests2.sh`, `tests3.sh`, `tests2.py`,
  `run_parameter_tests.sh`) under `dev_scripts/` — they reference HPC
  paths and are not run in CI.
- Add `requirements.txt` for use by both the recipe and the CI workflows.

Robustness:

- When no TIRs are detected, write valid empty `DANTE_TIR_final.gff3`
  (with `##gff-version 3` + `##DANTE_TIR version` banner) and
  `DANTE_TIR_final.fasta` instead of leaving the output dir bare. Print
  a final summary line so users can tell apart "ran cleanly, found 0"
  from "wrote no output / crashed".
- `dante_tir_summary.R`: capture per-class processing errors and emit
  a "Per-class summary unavailable: <reason>" placeholder in the HTML
  report instead of blowing up the whole report.
- `dt_utils.R`: keep the consensus matrix in matrix shape via
  `drop = FALSE`, return an empty `DNAStringSet` when every column has
  zero counts.

## 0.2.4

(Pre-migration release. See git log for details.)
