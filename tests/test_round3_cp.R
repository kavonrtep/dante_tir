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
  df <- read.table(out, header = TRUE, sep = "\t",
                   colClasses = c("character", "numeric"), na.strings = "NA")
  if (nrow(df) == 0) return(setNames(numeric(0), character(0)))
  sort_cp(setNames(df$cp, df$id))
}

sort_cp <- function(x) {
  x <- as.numeric(x)[order(names(x))]
  names(x) <- sort(names(x))
  x
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

message("ALL ROUND-3 SWITCH-POINT IDENTITY TESTS PASSED")
