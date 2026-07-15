#!/usr/bin/env Rscript
# tests/test_blast_reduce.R
#
# Proves that the Round-3 streaming awk reduction (blast3_reduce_awk +
# run_blast_tir_analysis) produces results IDENTICAL to the previous
# read.table() + filter_blast3() path, so the OOM/long-vector fix does not
# change any output.
#
# Two layers:
#   Part 1 — filter equivalence on a hand-crafted outfmt-6 table that hits
#            every predicate boundary (evalue, length, pident, sstart, self-hit,
#            underscore names, sstart>=send). Also checks coverage + switch
#            points end-to-end. No BLAST needed.
#   Part 2 — a real self-BLAST on a small FASTA: derive the OLD (full table ->
#            filter_blast3) and NEW (awk reduction) results from the SAME blastn
#            output, and additionally run the real run_blast_tir_analysis()
#            stream, asserting identical switch points throughout.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(Biostrings)
  library(parallel)
})

ROOT <- normalizePath(file.path(dirname(sub("--file=", "",
        grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
source(file.path(ROOT, "dt_utils.R"))

fail <- function(msg) { message("FAIL: ", msg); quit(status = 1) }
ok   <- function(msg) message("  ok: ", msg)

# outfmt-6 column order used throughout.
COLS <- c("qaccver", "saccver", "pident", "length", "mismatch", "gapopen",
          "qstart", "qend", "sstart", "send", "evalue", "bitscore")

# Reduce a data frame to the canonical (saccver, sstart, send) row multiset,
# sorted, so two paths can be compared regardless of row order.
canon <- function(df) {
  m <- data.frame(saccver = as.character(df$saccver),
                  sstart  = as.integer(df$sstart),
                  send    = as.integer(df$send),
                  stringsAsFactors = FALSE)
  m <- m[order(m$saccver, m$sstart, m$send), , drop = FALSE]
  rownames(m) <- NULL      # compare content, not the source row indices
  m
}

# Run the production awk program over a full outfmt-6 file, returning the
# reduced 3-column data frame (empty-safe).
awk_reduce <- function(full_file, min_length, min_identity, max_evalue = 1e-5) {
  reduced <- tempfile(fileext = ".tsv")
  prog <- blast3_reduce_awk(min_length, min_identity, max_evalue)
  cmd  <- paste0("awk ", shQuote(prog), " ", shQuote(full_file),
                 " > ", shQuote(reduced))
  st <- system2("bash", c("-c", shQuote(cmd)))
  if (!identical(as.integer(st), 0L)) fail(paste("awk exited", st))
  if (file.info(reduced)$size == 0)
    return(data.frame(saccver = character(), sstart = integer(),
                      send = integer(), stringsAsFactors = FALSE))
  read.table(reduced, header = FALSE, col.names = c("saccver", "sstart", "send"))
}

cov_equal <- function(a, b) {
  a <- a[order(names(a))]; b <- b[order(names(b))]
  identical(names(a), names(b)) && identical(lapply(a, as.numeric),
                                             lapply(b, as.numeric))
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
# type the numeric columns as read.table would
for (nm in setdiff(COLS, c("qaccver", "saccver")))
  full[[nm]] <- as.numeric(full[[nm]])

full_file <- tempfile(fileext = ".tsv")
write.table(full, full_file, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

# OLD path
old_df  <- read.table(full_file, header = FALSE, sep = "\t", col.names = COLS)
old_kept <- filter_blast3(old_df, min_length = 150, min_identity = 80)
# NEW path
new_kept <- awk_reduce(full_file, min_length = 150, min_identity = 80)

if (!identical(canon(old_kept), canon(new_kept)))
  fail("Part 1: kept row set differs between filter_blast3 and awk")
ok(sprintf("kept-row set identical (%d rows)", nrow(canon(new_kept))))

# expected: subjects 7(len151),9(pident),11(sstart13),ac,2(two hits) survive
expected_subjects <- sort(c("7", "9", "11", "ac", "2", "2"))
if (!identical(sort(as.character(canon(new_kept)$saccver)), expected_subjects))
  fail(paste("Part 1: unexpected surviving subjects:",
             paste(canon(new_kept)$saccver, collapse = ",")))
ok("surviving subjects match hand-computed expectation")

# coverage + switch points identical
cov_old <- get_coverage_from_blast(old_kept)
cov_new <- get_coverage_from_blast(new_kept)
if (!cov_equal(cov_old, cov_new)) fail("Part 1: coverage vectors differ")
ok("coverage vectors identical")

for (fn_name in c("find_switch_point_from_blast_coverage2",
                  "find_switch_point_from_blast_coverage3")) {
  fn <- get(fn_name)
  cp_old <- unlist(lapply(cov_old, fn))
  cp_new <- unlist(lapply(cov_new, fn))
  cp_old <- cp_old[order(names(cp_old))]
  cp_new <- cp_new[order(names(cp_new))]
  if (!identical(cp_old, cp_new)) fail(paste("Part 1: cp differ for", fn_name))
  ok(paste("switch points identical for", fn_name))
}

# ---------------------------------------------------------------------------
message("=== Part 2: real self-BLAST end-to-end ===")

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
  dss <- DNAStringSet(seqs)
  wd  <- tempfile(); dir.create(wd)
  fa  <- file.path(wd, "regions.fasta")
  writeXStringSet(dss, fa)

  system2(mkdb_bin, c("-in", fa, "-dbtype", "nucl"),
          stdout = FALSE, stderr = FALSE)

  # Single deterministic full table shared by OLD and NEW derivations.
  full2 <- file.path(wd, "full.tsv")
  st <- system2(blastn_bin, c("-query", fa, "-db", fa, "-outfmt", "6",
                              "-evalue", "1e-10", "-strand", "plus",
                              "-num_threads", "1", "-out", full2))
  if (!identical(as.integer(st), 0L)) fail("Part 2: blastn failed")

  old2 <- read.table(full2, header = FALSE, sep = "\t", col.names = COLS)
  old_kept2 <- filter_blast3(old2, min_length = 150, min_identity = 80)
  new_kept2 <- awk_reduce(full2, min_length = 150, min_identity = 80)
  if (!identical(canon(old_kept2), canon(new_kept2)))
    fail("Part 2: kept row set differs on real BLAST output")
  ok(sprintf("real-BLAST kept-row set identical (%d rows)",
             nrow(canon(new_kept2))))

  cov_old2 <- get_coverage_from_blast(old_kept2)
  cov_new2 <- get_coverage_from_blast(new_kept2)
  if (!cov_equal(cov_old2, cov_new2)) fail("Part 2: coverage differs")
  ok("real-BLAST coverage identical")

  fn <- find_switch_point_from_blast_coverage2
  cp_old2 <- unlist(lapply(cov_old2, fn)); cp_old2 <- cp_old2[order(names(cp_old2))]
  cp_new2 <- unlist(lapply(cov_new2, fn)); cp_new2 <- cp_new2[order(names(cp_new2))]
  if (!identical(cp_old2, cp_new2)) fail("Part 2: switch points differ")
  ok("real-BLAST switch points identical")

  # Exercise the actual production function (blastn | awk stream) and compare
  # its switch points to the OLD full-table computation.
  reduced_out <- file.path(wd, "prod_reduced.tsv")
  res <- run_blast_tir_analysis(
    query_db = fa, out_blast_file = reduced_out, blast_db = fa,
    swt_pt_fun = find_switch_point_from_blast_coverage2,
    filter_fun = filter_blast3, coverage_fun = get_coverage_from_blast,
    evalue = "1e-10", strand = "plus",
    min_length = 150, min_identity = 80, mc.cores = 1)
  cp_prod <- res$cp_vals; cp_prod <- cp_prod[order(names(cp_prod))]
  if (!identical(cp_prod, cp_old2))
    fail("Part 2: run_blast_tir_analysis stream differs from OLD path")
  ok("run_blast_tir_analysis() stream matches OLD path end-to-end")
}

message("ALL BLAST-REDUCE IDENTITY TESTS PASSED")
