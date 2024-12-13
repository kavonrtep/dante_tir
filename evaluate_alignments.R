#!/usr/bin/env Rscript
library(optparse)

options_list <- list(
    make_option(c("-c", "--contig_dir"), action = "store", type = "character",
                help = "directory with contigs"),
    make_option(c("-o", "--output"), action = "store", type = "character",
                help = "output prefix", default = NA)
)

parser <- OptionParser(option_list = options_list)
opt <- parse_args(parser, args = commandArgs(TRUE))

suppressPackageStartupMessages({
  library(Biostrings)
  library(BSgenome)
  library(Rbeast)
})

get_consensus <- function (CM, perc=70){
  totals <- colSums(CM)
  cm_perc <- t(t(CM)/totals) * 100
  cm_pass <- cm_perc >= perc
  # add N column
  N_row <- colSums(cm_pass) == 0

  cm_pass <- rbind(cm_pass, N_row)
  nucleotides <- c("A", "C", "G", "T", "N")

  idx <- apply(cm_pass,2, which, arr.ind = TRUE)
  cons <- nucleotides[idx] |> paste(collapse = "")
}

extract_TSD <- function(s1, s2){
  s1 <- strsplit(s1, "")[[1]]
  s2 <- strsplit(s2, "")[[1]]
  L <- length(s1) # length of the sequences is the same
  is_tsd <- logical(L)
  for (i in 1:L){
    tds1 <- s1[(L-i+1):L]
    tsd2 <- s2[1:i]
    if (all(tds1 == tsd2)){
      is_tsd[i] <- TRUE
    }
  }
  if (all(!is_tsd)){
    return(NA)
  }
  TSD_length <- max(which(is_tsd))
  TSD <- paste(s1[(L-TSD_length+1):L], collapse = "")
}



analyze_consensus_matrix <- function(CM) {
    # set order of rows
  nucleotides <- c("A", "C", "G", "T")
  CM <- CM[nucleotides,]
  # consensusMatrix is a numeric matrix with rows = A,C,G,T and columns = positions
  totals <- colSums(CM)
  # Compute frequencies per position, add pseudocounts
  # freqs <- t(t(CM + 0.1)/(totals+0.4))
  # allow zero frequencies and define  0 * log2(0) = 0
  freqs <- t(t(CM)/(totals))
  p_log2p <- freqs * log2(freqs)
  p_log2p[is.na(p_log2p)] <- 0 # set NA to 0
  # Information content per position:
  IC <- 2 + colSums(p_log2p)
  pssm <- log2(freqs / 0.25)
  cons <- get_consensus(CM, perc = 60)
  # Return results as a list
  aln_coverage <- colSums(CM)
  list(
    aln_coverage = aln_coverage,
    cons = cons,
    PSSM = pssm,
    InformationContent = IC
  )
}

analyze_consensus_matrix <- function(CM, min_coverage = 10, base_pseudocount = 1, scaling = TRUE) {
  # Set order of rows
  nucleotides <- c("A", "C", "G", "T")
  CM <- CM[nucleotides, , drop=FALSE]

  totals <- colSums(CM)

  # Define pseudocount approach
  # If scaling is TRUE, reduce pseudocount where coverage is high.
  # For example, pseudocount = base_pseudocount * (min_coverage / (totals + min_coverage))
  if (scaling) {
    pseudocounts <- sapply(totals, function(tot) {
      pc <- base_pseudocount * (min_coverage / (tot + min_coverage))
      pc
    })
  } else {
    pseudocounts <- rep(base_pseudocount, length(totals))
  }

  # Create a matrix of pseudocounts for each column
  pseudo_matrix <- matrix(rep(pseudocounts, each=nrow(CM)), nrow=nrow(CM), byrow=FALSE)

  # Compute frequencies with pseudocount
  freqs <- (CM + pseudo_matrix) / (matrix(totals + 4 * pseudocounts, nrow=nrow(CM), ncol=ncol(CM), byrow=TRUE))

  # Compute Information content
  p_log2p <- freqs * log2(freqs)
  IC <- 2 + colSums(p_log2p, na.rm = TRUE)

  # Compute PSSM = log2(freq / 0.25)
  pssm <- log2(freqs / 0.25)

  # Compute consensus (assuming get_consensus is defined similarly as before)
  cons <- get_consensus(CM, perc = 60)

  list(
    aln_coverage = totals,
    cons = cons,
    PSSM = pssm,
    InformationContent = IC
  )
}





find_switch_point <- function (aln_info, plotit = FALSE) {
  switch_info <- beast(aln_info$InformationContent, season = 'none', tcp.minmax = c(0,5), torder.minmax = c(0,0), quiet = TRUE)
  cp <- switch_info$trend$cp
  cpPr <- switch_info$trend$cpPr
  abrupt_change <- switch_info$trend$cpAbruptChange
  # condition to pass switch point detection
  # the strongest change point is the first one
  # position of switch point must be between 11-90 nucleotide of the alignment
  # abrubt change must be greater than 0.5
  # posterior probability of change point must be greater than 0.5
  # proportion of unambiguous nucleotides must be greater than 0.7
  c1 <- which.min(cp) == 1
  c2 <- cp[1] > 10 && cp[1] < 90
  c3 <- abrupt_change[1] > 0.7
  c4 <- cpPr[1] > 0.5
  if (c1 & c2 & c3 & c4) {
    Ncount <- sum(unlist(strsplit(substr(aln_info$cons,cp[1] - 10, cp[1]-1),""))=="N")
    if (Ncount < 6){
      return(NA)
    }
    if (plotit){
      plot(switch_info)
    }
    return(cp[1])
  } else {
    return(NA)
  }
}

extract_info_from_switch_points <- function(aln_info,ctg, cp){
  # extract ID of sequences which does not have gaps in the switch point
  # extract potential consensus of TIR sequence
  # For each positive ID, extract the sequence from the alignment - two sequences
  # TIR, TSD part
  cp_occupied <- subseq(ctg,cp, cp) != "-"
  # min length befor co must be 8 nucleotides
  ctg_upstream_cp <- subseq(ctg, cp - 8, cp - 1)
  ctg_upstream_cp_n_bases <- nchar(gsub("-", "", ctg_upstream_cp))
  ctg_donwstream_cp <- subseq(ctg, cp , cp + 50)
  ctg_donwstream_cp_n_bases <- nchar(gsub("-", "", ctg_donwstream_cp))
  pass <- ctg_upstream_cp_n_bases >= 8 & cp_occupied & ctg_donwstream_cp_n_bases >= 45
  ctg_pass <- ctg[pass]
  N <- names(ctg)[pass]
  # N hase syntax ID_start_end
  IDinfo <- sapply(strsplit(N, "_"),
                   function(x)c(ID = as.integer(x[1]),
                                start = as.integer(x[2]),
                                end = as.integer(x[3]))) |> t() |> as.data.frame()
  # count number of gaps before cp
  gap_count <- alphabetFrequency(subseq(ctg_pass, 1, cp))[, "-"]
  # recalculate cp position for each sequence
  CP <- IDinfo$start - gap_count + cp
  IDinfo$CP <- CP
  IDinfo$coverage <- length(ctg_pass)
  TIR <- subseq(ctg_pass, cp, cp + 50)
  TSD <- subseq(ctg_pass, cp - 12, cp-1)
  IDinfo$TIR <- as.character(TIR)
  IDinfo$TSD <- as.character(TSD)
  IDinfo
}


upstream_contigs <- dir(opt$contig_dir, pattern = "upstream_Contig", full.names = TRUE)
downstream_contigs <- dir(opt$contig_dir, pattern = "downstream_Contig", full.names = TRUE)

upstream_info <- list()
for (f in upstream_contigs) {
  i <- i + 1
  ctg <- readDNAStringSet(f)
  CM <- consensusMatrix(ctg)[1:4,]
  aln_info <- analyze_consensus_matrix(CM)
  cp <- find_switch_point(aln_info)
  if (!is.na(cp)) {
    # extract classification and Contig ID
    contig_id <- gsub("[.]fasta","", gsub(".+_Contig", "", f))
    classification <- gsub("_upstream.+","",basename(f))
    info <- extract_info_from_switch_points(aln_info, ctg, cp)
    info$ContigID <- contig_id
    info$Classification <- classification
    info$Filename <- f  # remove later - for debugging
    upstream_info[[contig_id]] <- info

  }
}

downstream_info <- list()
for (f in downstream_contigs) {
  i <- i + 1
  ctg <- readDNAStringSet(f)
  CM <- consensusMatrix(ctg)[1:4,]
  aln_info <- analyze_consensus_matrix(CM)
  cp <- find_switch_point(aln_info)
  if (!is.na(cp)) {
    # extract classification and Contig ID
    contig_id <- gsub("[.]fasta","", gsub(".+_Contig", "", f))
    classification <- gsub("_downstream.+","",basename(f))
    info <- extract_info_from_switch_points(aln_info, ctg, cp)
    info$ContigID <- contig_id
    info$Classification <- classification
    info$Filename <- f  # remove later - for debugging
    downstream_info[[contig_id]] <- info
  }
}

downstream_info <- do.call(rbind, downstream_info)
upstream_info <- do.call(rbind, upstream_info)


both_side_id <-  intersect(downstream_info$ID, upstream_info$ID)

for (i in seq_along(both_side_id)) {
  id <- both_side_id[i]
  downstream <- downstream_info[downstream_info$ID == id,]
  upstream <- upstream_info[upstream_info$ID == id,]
  message("ID: ", id)
  print(downstream)
  print(upstream)
  tir_aln <- pairwiseAlignment(substr(upstream$TIR,1,9), substr(downstream$TIR,1,9), type = "global")
  message("TIR Alignment partial:")
  print(tir_aln)
  message("TIR Alignment full:")
  print(pairwiseAlignment(upstream$TIR, downstream$TIR, type = "global-local"))
  TSD <- extract_TSD(upstream$TSD, as.character(reverseComplement(DNAString(downstream$TSD))))
  message("TSD: ", TSD)
  message("----------------------")

}
