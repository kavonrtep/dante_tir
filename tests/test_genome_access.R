#!/usr/bin/env Rscript
# tests/test_genome_access.R
#
# genome_fa_handle() replaced readDNAStringSet(genome) at the two places that
# extract TIR sequences (detect_tirs.R and cluster_tir_sequences), because
# holding a ~94 Gbp assembly in RAM at 1 byte/base is not viable. This asserts
# the substitution is exact: the same GRanges must yield the same sequences
# through an FaFile as through a fully loaded DNAStringSet.
#
# Strand is the trap. Both paths must reverse-complement minus-strand ranges;
# if only one did, TIR sequences would come out silently wrong.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(Biostrings)
  library(BSgenome)     # provides getSeq() for DNAStringSet, as detect_tirs.R does
  library(Rsamtools)
})

ROOT <- normalizePath(file.path(dirname(sub("--file=", "",
        grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
source(file.path(ROOT, "dt_utils.R"))

fail <- function(msg) { message("FAIL: ", msg); quit(status = 1) }
ok   <- function(msg) message("  ok: ", msg)

set.seed(5)
rnd <- function(n) paste(sample(c("A", "C", "G", "T"), n, replace = TRUE),
                         collapse = "")

wd <- tempfile(); dir.create(wd)
fa <- file.path(wd, "genome.fasta")
# Headers carry descriptions, as real assemblies do: the DNAStringSet path had
# to strip them by hand, the .fai keeps only the first token.
seqs <- DNAStringSet(c("chr1 length=4000 assembled" = rnd(4000),
                       "chr2 some description here" = rnd(2500),
                       "scaffold_3" = rnd(3200)))
writeXStringSet(seqs, fa)

gr <- GRanges(
  seqnames = c("chr1", "chr1", "chr2", "scaffold_3", "chr2", "chr1", "scaffold_3"),
  ranges = IRanges(start = c(1, 1500, 10, 900, 2400, 3990, 1),
                   end   = c(120, 1800, 610, 1500, 2500, 4000, 3200)),
  strand = c("+", "-", "-", "+", "-", "+", "-"))

# OLD: whole genome in RAM, descriptions stripped by hand.
g_old <- readDNAStringSet(fa)
names(g_old) <- sub(" .*", "", names(g_old))
old <- getSeq(g_old, gr)

# NEW: indexed random access.
handle <- genome_fa_handle(fa)
if (!file.exists(paste0(fa, ".fai"))) fail("genome_fa_handle did not create the .fai")
ok("genome_fa_handle indexes a genome that has no .fai yet")

new <- getSeq(handle, gr)

if (length(old) != length(new)) fail("different number of sequences returned")
# Compare the sequences themselves: the two paths label the result differently
# (DNAStringSet returns them unnamed, FaFile names them after the seqnames),
# which is inert because both call sites assign names on the next line. That
# assumption is asserted below.
if (!identical(unname(as.character(old)), unname(as.character(new)))) {
  d <- which(unname(as.character(old)) != unname(as.character(new)))
  fail(sprintf("%d/%d ranges differ, first at %s:%d-%d strand %s",
               length(d), length(old), as.character(seqnames(gr))[d[1]],
               start(gr)[d[1]], end(gr)[d[1]], as.character(strand(gr))[d[1]]))
}
ok(sprintf("%d ranges identical, both strands (%d minus)",
           length(gr), sum(as.character(strand(gr)) == "-")))

# Guard the trap explicitly rather than trusting the mixed set above: a
# minus-strand range must be the reverse complement of the plus-strand one.
gr_p <- GRanges("chr1", IRanges(100, 200), strand = "+")
gr_m <- GRanges("chr1", IRanges(100, 200), strand = "-")
if (!identical(as.character(reverseComplement(getSeq(handle, gr_p))),
               as.character(getSeq(handle, gr_m))))
  fail("FaFile getSeq does not reverse-complement minus-strand ranges")
ok("minus strand is reverse-complemented")

# A pre-existing index must be used as is, not rebuilt or ignored.
h2 <- genome_fa_handle(fa)
if (!identical(as.character(getSeq(h2, gr)), as.character(new)))
  fail("second handle over the same genome returned different sequences")
ok("an existing .fai is reused")

# An unindexable genome must fail loudly rather than silently returning nothing.
ro <- file.path(wd, "readonly"); dir.create(ro)
fa_ro <- file.path(ro, "genome.fasta")
writeXStringSet(seqs, fa_ro)
Sys.chmod(ro, "555")
res <- tryCatch({ genome_fa_handle(fa_ro); "no error" },
                error = function(e) "error")
Sys.chmod(ro, "755")
if (!identical(res, "error")) {
  message("  note: .fai creation succeeded despite a read-only directory ",
          "(running as root?) - skipping the failure-mode check")
} else {
  ok("an unindexable genome raises a clear error")
}

# Both call sites must overwrite the names right after getSeq(), which is what
# makes the labelling difference above harmless. Guard that with a source
# check, so moving the assignment away from the call cannot pass unnoticed.
for (f in c("detect_tirs.R", "dt_utils.R")) {
  src <- suppressWarnings(readLines(file.path(ROOT, f)))
  hits <- grep("tir_seqs <- getSeq\\(genome", src)
  for (h in hits) {
    if (!any(grepl("^\\s*names\\(tir_seqs\\) <-", src[(h + 1):(h + 2)])))
      fail(sprintf("%s:%d - getSeq(genome, ...) is no longer followed by names(tir_seqs) <-, so the FaFile labelling difference stops being inert",
                   f, h))
  }
  if (length(hits)) ok(sprintf("%s: getSeq result is renamed immediately (%d site%s)",
                               f, length(hits), if (length(hits) > 1) "s" else ""))
}

message("ALL GENOME-ACCESS TESTS PASSED")
