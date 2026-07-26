#!/usr/bin/env Rscript
# tests/test_round3_cp.R
#
# Proves that the Round-3 streaming path (blastn | blast_cp.py) produces
# switch points IDENTICAL to the R implementation it replaces:
#
#   read.table -> filter_blast3 -> get_coverage_from_blast
#              -> find_switch_point_from_blast_coverage{2,3}
#
# The R functions above stay in dt_utils.R purely as this reference; production
# never calls them (see docs/round3_scaling.md).
#
# Three layers:
#   Part 1 — filter equivalence on a hand-crafted outfmt-6 table that hits every
#            predicate boundary (evalue, length, pident, sstart, self-hit,
#            underscore names, sstart >= send). No BLAST needed.
#   Part 2 — switch-point equivalence on randomly generated hit sets whose
#            coverage profiles have real element/flank steps, for both methods.
#            This is the part that exercises the numpy reimplementations.
#   Part 3 — a real self-BLAST on a small FASTA: OLD (full table -> filter ->
#            coverage -> cp in R) vs the actual production
#            run_blast_tir_analysis() stream.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(Biostrings)
  library(parallel)
})

ROOT <- normalizePath(file.path(dirname(sub("--file=", "",
        grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
source(file.path(ROOT, "dt_utils.R"))

BLAST_CP <- file.path(ROOT, "blast_cp.py")

fail <- function(msg) { message("FAIL: ", msg); quit(status = 1) }
ok   <- function(msg) message("  ok: ", msg)

# outfmt-6 column order used throughout.
COLS <- c("qaccver", "saccver", "pident", "length", "mismatch", "gapopen",
          "qstart", "qend", "sstart", "send", "evalue", "bitscore")

# --- helpers ---------------------------------------------------------------

# Reference: the pre-0.2.9 R path, from a full outfmt-6 table to named cps.
r_cp <- function(full_file, fn, min_length = 150, min_identity = 80) {
  df <- read.table(full_file, header = FALSE, sep = "\t", col.names = COLS)
  kept <- filter_blast3(df, min_length = min_length, min_identity = min_identity)
  cov <- get_coverage_from_blast(kept)
  cp <- unlist(lapply(cov, fn))
  if (is.null(cp)) return(setNames(numeric(0), character(0)))
  sort_cp(cp)
}

# Production: blast_cp.py over the same table.
py_cp <- function(full_file, subjects_fa, method, columns = "std",
                  min_length = 150, min_identity = 80, max_evalue = 1e-5) {
  out <- tempfile(fileext = ".tsv")
  args <- c("--subjects", shQuote(subjects_fa), "--method", method,
            "--columns", columns, "--input", shQuote(full_file),
            "--out", shQuote(out),
            "--min-length", min_length, "--min-identity", min_identity)
  if (!is.null(max_evalue)) args <- c(args, "--max-evalue", format(max_evalue))
  st <- system2(BLAST_CP, args, stdout = FALSE, stderr = FALSE)
  if (!identical(as.integer(st), 0L)) fail(paste("blast_cp.py exited", st))
  df <- read_cp_table(out)
  if (nrow(df) == 0) return(setNames(numeric(0), character(0)))
  sort_cp(setNames(df$cp, df$id))
}

sort_cp <- function(x) {
  if (length(x) == 0) return(setNames(numeric(0), character(0)))
  nm <- names(x)
  if (is.null(nm)) fail("sort_cp: cp vector has no names")
  out <- as.numeric(x)[order(nm)]      # as.numeric() drops names, so re-attach
  names(out) <- nm[order(nm)]          # them from the same ordering
  out
}

# Named numeric vectors are compared by name set and value, NA == NA.
cp_equal <- function(a, b) {
  identical(names(a), names(b)) && identical(as.numeric(a), as.numeric(b))
}

cp_report <- function(a, b, label) {
  if (!identical(names(a), names(b))) {
    fail(sprintf("%s: subject sets differ (R %d, python %d; symmetric diff: %s)",
                 label, length(a), length(b),
                 paste(setdiff(union(names(a), names(b)),
                               intersect(names(a), names(b))), collapse = ",")))
  }
  if (!cp_equal(a, b)) {
    d <- which(!(is.na(a) & is.na(b)) & (is.na(a) | is.na(b) | a != b))
    fail(sprintf("%s: %d/%d switch points differ, e.g. %s: R=%s python=%s",
                 label, length(d), length(a), names(a)[d[1]],
                 a[d[1]], b[d[1]]))
  }
  ok(sprintf("%s: %d subjects identical (%d with a switch point)",
             label, length(a), sum(!is.na(a))))
}

# Minimal FASTA covering every subject accession, long enough for the profile.
write_subjects <- function(names_vec, len, path) {
  set.seed(7)
  seqs <- DNAStringSet(setNames(
    vapply(seq_along(names_vec),
           function(i) paste(sample(c("A", "C", "G", "T"), len, replace = TRUE),
                             collapse = ""),
           character(1)),
    names_vec))
  writeXStringSet(seqs, path)
  path
}

# Write hits as a full outfmt-6 table (unused columns filled plausibly).
write_full_table <- function(hits, path) {
  df <- data.frame(
    qaccver = hits$qaccver, saccver = hits$saccver,
    pident = hits$pident, length = hits$length,
    mismatch = 0, gapopen = 0, qstart = 1, qend = hits$length,
    sstart = hits$sstart, send = hits$send,
    evalue = hits$evalue, bitscore = hits$length,
    stringsAsFactors = FALSE)
  write.table(df, path, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  path
}

if (!file.exists(BLAST_CP)) fail(paste("blast_cp.py not found at", BLAST_CP))
if (!identical(as.integer(system2(BLAST_CP, "--help",
                                  stdout = FALSE, stderr = FALSE)), 0L)) {
  fail("blast_cp.py --help failed (is numpy installed?)")
}

# ---------------------------------------------------------------------------
message("=== Part 1: filter equivalence on crafted edge cases ===")

rows <- list(
  # keep: clean hit
  c(1, 2, 90,   200, 0,0, 1,200,  20, 400, "1e-20", 300),
  # drop: evalue == threshold (not strictly <)
  c(1, 3, 90,   200, 0,0, 1,200,  20, 400, "1e-5",  300),
  # drop: evalue above threshold
  c(1, 4, 90,   200, 0,0, 1,200,  20, 400, "1e-4",  300),
  # drop: length == min_length (not strictly >)
  c(1, 6, 90,   150, 0,0, 1,150,  20, 400, "1e-20", 300),
  # keep: length just over
  c(1, 7, 90,   151, 0,0, 1,151,  20, 400, "1e-20", 300),
  # drop: pident == min_identity (not strictly >)
  c(1, 8, 80,   200, 0,0, 1,200,  20, 400, "1e-20", 300),
  # keep: pident just over
  c(1, 9, 80.5, 200, 0,0, 1,200,  20, 400, "1e-20", 300),
  # drop: sstart == 12 (not strictly > 12)
  c(1, 10, 90,  200, 0,0, 1,200,  12, 400, "1e-20", 300),
  # keep: sstart == 13
  c(1, 11, 90,  200, 0,0, 1,200,  13, 400, "1e-20", 300),
  # drop: sstart > send
  c(1, 12, 90,  200, 0,0, 1,200, 400,  50, "1e-20", 300),
  # drop: sstart == send
  c(1, 13, 90,  200, 0,0, 1,200, 100, 100, "1e-20", 300),
  # drop: self-hit (q == s)
  c(5, 5, 99,   300, 0,0, 1,300,  20, 500, "1e-30", 500),
  # drop: self after underscore strip ("a_b" -> "ab" == "ab")
  c("a_b", "ab", 95, 200, 0,0, 1,200, 20, 400, "1e-20", 300),
  # keep: underscore strip, not self ("a_b" -> "ab" != "ac")
  c("a_b", "ac", 95, 200, 0,0, 1,200, 20, 400, "1e-20", 300),
  # keep: second overlapping hit on subject 2 (exercises coverage union)
  c(7, 2, 88,   250, 0,0, 1,250, 100, 600, "1e-15", 280)
)
full <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
names(full) <- COLS
for (nm in setdiff(COLS, c("qaccver", "saccver")))
  full[[nm]] <- as.numeric(full[[nm]])

wd1 <- tempfile(); dir.create(wd1)
full_file <- file.path(wd1, "full.tsv")
write.table(full, full_file, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
subj1 <- write_subjects(unique(as.character(full$saccver)), 800,
                        file.path(wd1, "subjects.fasta"))

# The surviving subject set is hand-computed: 7 (length 151), 9 (pident),
# 11 (sstart 13), ac (underscore strip), 2 (two overlapping hits).
kept1 <- filter_blast3(read.table(full_file, header = FALSE, sep = "\t",
                                  col.names = COLS),
                       min_length = 150, min_identity = 80)
if (!identical(sort(as.character(kept1$saccver)), sort(c("2", "7", "9", "11", "ac", "2"))))
  fail(paste("Part 1: unexpected surviving subjects:",
             paste(sort(kept1$saccver), collapse = ",")))
ok("filter_blast3 surviving subjects match hand-computed expectation")

for (m in c("win200", "cumsum")) {
  fn <- if (m == "win200") find_switch_point_from_blast_coverage2
        else find_switch_point_from_blast_coverage3
  cp_report(r_cp(full_file, fn), py_cp(full_file, subj1, m),
            paste("Part 1", m))
}

# ---------------------------------------------------------------------------
message("=== Part 2: switch points on random coverage profiles ===")

# Build hit sets whose per-subject coverage looks like the real thing: a
# low-coverage flank, a step at a random boundary, then deep element coverage
# running to the end. Some subjects get shallow or absent steps so that the
# NA / threshold branches are exercised too.
set.seed(11)
n_subj <- 40
seq_len_bp <- 3000
hits <- list(qaccver = character(), saccver = character(), pident = numeric(),
             length = numeric(), sstart = numeric(), send = numeric(),
             evalue = character())
add_hit <- function(s, a, b) {
  if (b <= a) return(invisible(NULL))
  hits$qaccver[[length(hits$qaccver) + 1]] <<- paste0("q", length(hits$qaccver))
  hits$saccver[[length(hits$saccver) + 1]] <<- s
  hits$pident[[length(hits$pident) + 1]] <<- 85 + runif(1) * 14
  hits$length[[length(hits$length) + 1]] <<- b - a + 1
  hits$sstart[[length(hits$sstart) + 1]] <<- a
  hits$send[[length(hits$send) + 1]] <<- b
  hits$evalue[[length(hits$evalue) + 1]] <<- "1e-20"
  invisible(NULL)
}
for (i in seq_len(n_subj)) {
  s <- paste0("s", i)
  boundary <- sample(600:2200, 1)
  depth <- sample(c(0, 2, 5, 15, 40, 90), 1)         # element-side depth
  flank <- sample(c(0, 1, 2, 4), 1)                  # flank-side depth
  for (k in seq_len(depth)) {
    start <- max(13, boundary + sample(-40:40, 1))
    add_hit(s, start, min(seq_len_bp, start + sample(400:(seq_len_bp), 1)))
  }
  for (k in seq_len(flank)) {
    start <- max(13, sample(13:max(14, boundary - 300), 1))
    add_hit(s, start, min(boundary, start + sample(151:600, 1)))
  }
  # a couple of hits that the filter must drop, mixed in
  add_hit(s, 5, 300)                                  # sstart <= 12
  add_hit(s, 400, 380)                                # sstart > send
}
hits <- as.data.frame(hits, stringsAsFactors = FALSE)
# short hits must be dropped by min_length; make a few explicitly short
hits$length[seq(1, nrow(hits), by = 17)] <- 100

wd2 <- tempfile(); dir.create(wd2)
full2 <- write_full_table(hits, file.path(wd2, "full.tsv"))
subj2 <- write_subjects(paste0("s", seq_len(n_subj)), seq_len_bp,
                        file.path(wd2, "subjects.fasta"))

for (m in c("win200", "cumsum")) {
  fn <- if (m == "win200") find_switch_point_from_blast_coverage2
        else find_switch_point_from_blast_coverage3
  a <- r_cp(full2, fn)
  b <- py_cp(full2, subj2, m)
  if (sum(!is.na(a)) < 3)
    fail(paste("Part 2", m, "- fixture produced too few switch points to be a real test"))
  cp_report(a, b, paste("Part 2", m))
}

# ---------------------------------------------------------------------------
message("=== Part 3: real self-BLAST end-to-end ===")

blastn_bin <- Sys.which("blastn")
mkdb_bin   <- Sys.which("makeblastdb")
if (blastn_bin == "" || mkdb_bin == "") {
  message("  skip: blastn/makeblastdb not on PATH")
} else {
  set.seed(1)
  rnd <- function(n) paste(sample(c("A","C","G","T"), n, replace = TRUE), collapse = "")
  # A shared 500bp core makes several sequences BLAST-similar; unique flanks vary.
  core <- rnd(500)
  seqs <- character()
  for (i in 1:30) {
    if (i %% 2 == 0) {
      s <- paste0(rnd(sample(80:200, 1)), core, rnd(sample(80:200, 1)))
    } else {
      s <- rnd(sample(700:1100, 1))
    }
    seqs[as.character(i)] <- s
  }
  wd  <- tempfile(); dir.create(wd)
  fa  <- file.path(wd, "regions.fasta")
  writeXStringSet(DNAStringSet(seqs), fa)

  system2(mkdb_bin, c("-in", fa, "-dbtype", "nucl"),
          stdout = FALSE, stderr = FALSE)

  # OLD path derived from a plain full table.
  full3 <- file.path(wd, "full.tsv")
  st <- system2(blastn_bin, c("-query", fa, "-db", fa, "-outfmt", "6",
                              "-evalue", "1e-10", "-strand", "plus",
                              "-num_threads", "1", "-out", full3))
  if (!identical(as.integer(st), 0L)) fail("Part 3: blastn failed")
  cp_old <- r_cp(full3, find_switch_point_from_blast_coverage2)

  # NEW path: the actual production function (blastn | blast_cp.py).
  res <- run_blast_tir_analysis(
    query_db = fa, out_cp_file = file.path(wd, "cp.tsv"), blast_db = fa,
    method = "win200", evalue = "1e-10", strand = "plus",
    min_length = 150, min_identity = 80, mc.cores = 1)
  cp_report(cp_old, sort_cp(res$cp_vals), "Part 3 run_blast_tir_analysis")

  # Cached result must be reused rather than re-running BLAST.
  res2 <- run_blast_tir_analysis(
    query_db = "/nonexistent", out_cp_file = file.path(wd, "cp.tsv"),
    blast_db = "/nonexistent", method = "win200", mc.cores = 1)
  if (!cp_equal(sort_cp(res$cp_vals), sort_cp(res2$cp_vals)))
    fail("Part 3: cached cp table was not reused")
  ok("existing cp table is reused without re-running BLAST")
}

# ---------------------------------------------------------------------------
message("=== Part 4: query sampling + exact second pass (P2) ===")

if (blastn_bin == "" || mkdb_bin == "") {
  message("  skip: blastn/makeblastdb not on PATH")
} else {
  # Two families plus singletons, so subjects differ in how many relatives
  # support them: family A is well covered even in a small query sample, while
  # family B and the singletons are exactly the under-supported subjects that
  # sampling would get wrong and the second pass has to rescue.
  set.seed(3)
  rnd <- function(n) paste(sample(c("A","C","G","T"), n, replace = TRUE), collapse = "")
  mutate <- function(s, rate = 0.03) {
    ch <- strsplit(s, "")[[1]]
    hit <- which(runif(length(ch)) < rate)
    ch[hit] <- sample(c("A","C","G","T"), length(hit), replace = TRUE)
    paste(ch, collapse = "")
  }
  element_a <- rnd(1200)
  element_b <- rnd(1200)
  n_a <- 50; n_b <- 6; n_single <- 4
  n_copies <- n_a + n_b + n_single
  seqs <- character()
  for (i in seq_len(n_a))
    seqs[paste0("a", i)] <- paste0(rnd(sample(900:1400, 1)), mutate(element_a))
  for (i in seq_len(n_b))
    seqs[paste0("b", i)] <- paste0(rnd(sample(900:1400, 1)), mutate(element_b))
  for (i in seq_len(n_single))
    seqs[paste0("s", i)] <- rnd(sample(2100:2600, 1))
  wd4 <- tempfile(); dir.create(wd4)
  fa4 <- file.path(wd4, "regions.fasta")
  writeXStringSet(DNAStringSet(seqs), fa4)
  system2(mkdb_bin, c("-in", fa4, "-dbtype", "nucl"), stdout = FALSE, stderr = FALSE)

  run_p2 <- function(tag, max_queries, min_support) {
    run_blast_tir_analysis(
      query_db = fa4, out_cp_file = file.path(wd4, paste0(tag, ".tsv")),
      blast_db = fa4, method = "win200", evalue = "1e-10", strand = "plus",
      min_length = 150, min_identity = 80, mc.cores = 1,
      max_queries = max_queries, min_support = min_support, seed = 42)
  }

  # Reference: every sequence used as a query (max_queries = 0).
  exact <- run_p2("exact", 0, 30)
  cp_exact <- sort_cp(exact$cp_vals)
  if (sum(!is.na(cp_exact)) < 10)
    fail("Part 4: fixture produced too few switch points to be a real test")
  if (length(exact$pass2_ids) != 0)
    fail("Part 4: max_queries = 0 must not trigger sampling")
  ok(sprintf("max_queries = 0 is the exact path (%d subjects, %d switch points)",
             length(cp_exact), sum(!is.na(cp_exact))))

  # min_support = Inf sends every subject to pass 2, which BLASTs the FULL
  # query set against them -- so the merged result must equal the exact run.
  # This is the load-bearing invariant: pass 2 is not an approximation.
  all_p2 <- run_p2("allpass2", 15, Inf)
  if (length(all_p2$pass2_ids) != length(cp_exact) &&
      length(all_p2$pass2_ids) < n_copies - 5)
    fail(sprintf("Part 4: expected ~every subject in pass 2, got %d",
                 length(all_p2$pass2_ids)))
  cp_report(cp_exact, sort_cp(all_p2$cp_vals),
            "Part 4 pass-2-only (min_support = Inf)")

  # Realistic setting: sample 15 of 60 queries, resolve the under-supported
  # ones exactly. min_support is 4 here rather than the production default of
  # 30 because observed support cannot exceed the 15 sampled queries; the
  # point is the mix of pass-1 and pass-2 subjects, not the exact cut-off.
  # Subject sets must match the exact run, and every subject that went through
  # pass 2 must carry the exact run's answer.
  sampled <- run_p2("sampled", 15, 4)
  cp_s <- sort_cp(sampled$cp_vals)
  if (!identical(names(cp_exact), names(cp_s)))
    fail("Part 4: sampled run covers a different subject set than the exact run")
  ok(sprintf("sampled run covers the same %d subjects", length(cp_s)))

  if (length(sampled$pass2_ids) == 0)
    fail("Part 4: fixture never exercised pass 2 - lower min_support or widen the copy-number spread")
  p2 <- intersect(sampled$pass2_ids, names(cp_exact))
  if (!identical(as.numeric(cp_exact[p2]), as.numeric(cp_s[p2])))
    fail("Part 4: pass-2 subjects disagree with the exact run")
  ok(sprintf("%d pass-2 subjects match the exact run exactly", length(p2)))

  # Sampling does move some switch points (the argmax is estimated from fewer
  # relatives). Quantify it here rather than assert a bound: on a 60-sequence
  # fixture at a 1-in-4 sample the noise is far worse than in production. The
  # measurement that matters is on real data -- docs/round3_scaling.md P2.
  both <- !is.na(cp_exact) & !is.na(cp_s)
  d <- abs(cp_exact[both] - cp_s[both])
  message(sprintf(
    "  note: %d subjects called in both; |delta cp| median %.0f max %.0f; %d flipped to/from NA",
    sum(both), if (any(both)) median(d) else 0, if (any(both)) max(d) else 0,
    sum(xor(is.na(cp_exact), is.na(cp_s)))))
}

# ---------------------------------------------------------------------------
message("=== Part 5: the query sample is deterministic ===")

# The sample must depend on the seed and the input alone. These are the two
# ways it could silently stop doing so.
draw <- function(seed) with_seed(seed, sort(sample.int(10000, 50)))

if (!identical(draw(42), draw(42)))
  fail("Part 5: same seed drew a different sample")
ok("same seed draws the same queries")

if (identical(draw(42), draw(43)))
  fail("Part 5: different seeds drew the same sample (seed is being ignored)")
ok("a different seed draws a different sample")

# 1. The caller's position in the RNG stream must not leak in.
set.seed(1); invisible(runif(1)); a <- draw(42)
set.seed(9); invisible(runif(1000)); b <- draw(42)
if (!identical(a, b))
  fail("Part 5: the sample depends on preceding RNG use")
ok("preceding RNG use does not change the sample")

# 2. Neither must the generator in force -- parallel code often switches the
#    process to L'Ecuyer-CMRG, and R changed the default sampler in 3.6.
old_kind <- RNGkind()
RNGkind("L'Ecuyer-CMRG")
c_ecuyer <- draw(42)
suppressWarnings(RNGkind(kind = old_kind[1], normal.kind = old_kind[2],
                         sample.kind = old_kind[3]))
if (!identical(a, c_ecuyer))
  fail("Part 5: the sample depends on the RNGkind in force")
ok("the sample is independent of RNGkind()")

# 3. with_seed must leave the caller's RNG exactly as it found it.
set.seed(7)
before <- get(".Random.seed", envir = .GlobalEnv)
invisible(draw(42))
if (!identical(before, get(".Random.seed", envir = .GlobalEnv)))
  fail("Part 5: with_seed disturbed the caller's RNG stream")
if (!identical(old_kind, RNGkind()))
  fail("Part 5: with_seed left RNGkind() changed")
ok("with_seed restores the caller's RNG stream and generator")

# 4. End to end: the same seed must reproduce the cp table byte for byte, and
#    a different seed must still describe the same subjects.
if (blastn_bin == "" || mkdb_bin == "") {
  message("  skip: blastn/makeblastdb not on PATH")
} else {
  again <- run_p2("sampled_again", 15, 4)
  first <- readLines(file.path(wd4, "sampled.tsv"))
  if (!identical(first, readLines(file.path(wd4, "sampled_again.tsv"))))
    fail("Part 5: re-running with the same seed produced a different cp table")
  ok("a re-run with the same seed reproduces the cp table byte for byte")

  other <- run_blast_tir_analysis(
    query_db = fa4, out_cp_file = file.path(wd4, "seed99.tsv"), blast_db = fa4,
    method = "win200", evalue = "1e-10", strand = "plus", min_length = 150,
    min_identity = 80, mc.cores = 1, max_queries = 15, min_support = 4,
    seed = 99)
  if (!identical(sort(names(other$cp_vals)), sort(names(exact$cp_vals))))
    fail("Part 5: a different seed changed which subjects are described")
  ok("a different seed still covers the same subjects")
}

# ---------------------------------------------------------------------------
message("=== Part 6: coverage3 reads its thresholds at the switch point ===")

# Until 0.2.9, find_switch_point_from_blast_coverage3() read its quality
# thresholds at cp - W: 200 bp before the switch point it had just found. m1/m2
# there are indexed by position, so that tested the wrong stretch of sequence.
# At the smallest admissible cp of 201 it compared m1[1] -- the coverage of a
# single base -- against "flank mean < 3", which let almost anything through
# and piled 13,753 of run-000129's CACTA calls onto position 201.

# A profile whose flank is genuinely dirty at the switch point (mean 46) but
# looks clean 200 bp earlier (mean 2). The old code accepted cp = 201 here.
cvrg_dirty <- c(rep(2, 121), rep(113, 192), rep(246, 298))
cp_dirty <- find_switch_point_from_blast_coverage3(cvrg_dirty)
if (!is.na(cp_dirty))
  fail(sprintf("Part 6: accepted cp=%s although the flank mean at it is %.1f (needs < 3)",
               cp_dirty, mean(cvrg_dirty[1:cp_dirty])))
ok("a switch point whose flank is dirty AT it is rejected")

# The same shape with a genuinely clean flank must still be found.
cvrg_clean <- c(rep(0, 700), rep(150, 1300))
cp_clean <- find_switch_point_from_blast_coverage3(cvrg_clean)
if (is.na(cp_clean) || abs(cp_clean - 700) > 5)
  fail(sprintf("Part 6: expected a switch point near 700, got %s", cp_clean))
ok(sprintf("a clean switch point is still found (cp = %d)", cp_clean))

# And the thresholds must hold at cp itself, not 200 bp before it.
m1_at_cp <- mean(cvrg_clean[1:cp_clean])
m2_at_cp <- mean(cvrg_clean[cp_clean:length(cvrg_clean)])
if (!(m1_at_cp < 3 && m2_at_cp > 20))
  fail("Part 6: the returned cp does not satisfy the thresholds at its own position")
ok(sprintf("thresholds hold at cp: flank mean %.2f < 3, element mean %.1f > 20",
           m1_at_cp, m2_at_cp))

message("ALL ROUND-3 SWITCH-POINT IDENTITY TESTS PASSED")
