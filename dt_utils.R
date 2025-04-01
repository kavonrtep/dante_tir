get_consensus_from_aln <- function (aln, perc=70){
    nucleotides <- c("A", "C", "G", "T", "N")
    CM <- consensusMatrix(aln)[nucleotides,]
    out <- DNAStringSet(get_consensus(CM, perc))
    out
}

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
  TSD
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
  switch_info <- beast(aln_info$InformationContent, season = 'none',
                       tcp.minmax = c(1,2),
                       torder.minmax = c(0,0),
                       tseg.leftmargin = 10,
                       tseg.rightmargin = 10,
                       quiet         = TRUE)
  cp <- switch_info$trend$pos_cp
  cpPr <- switch_info$trend$pos_cpPr
  abrupt_change <- switch_info$trend$pos_cpAbruptChange
  # condition to pass switch point detection
  # the strongest change point is the first one
  # position of switch point must be between 11-90 nucleotide of the alignment
  # abrubt change must be greater than 0.5
  # posterior probability of change point must be greater than 0.5
  # proportion of unambiguous nucleotides must be greater than 0.7
  c1 <- which.min(cp) == 1
  if (length(c1) == 0){
    return(NA)
  }
  c2 <- cp[1] > 12 && cp[1] < 120
  c3 <- abrupt_change[1] > 0.2
  c4 <- cpPr[1] > 0.3
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

eval_aln_length <- function(s1,s2,L0=4){
  L <- nchar(s1)
  sc<-numeric()
  for (i in L0:L){
    mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -1, baseOnly = TRUE)
    sc[i] <- pairwiseAlignment(substr(s1,1,i), substr(s2,1,i),
                               type = "global", scoreOnly = TRUE,
                               substitutionMatrix = mat)



  }
  L <- which.max(sc)
  aln <- pairwiseAlignment(substr(s1,1,L), substr(s2,1,L), type = "global")
  aln
}

eval_aln_length_alt <- function(s1,s2, threshold = 0.74, plot_it = FALSE){
  L <- nchar(s1)
  sc<-numeric()
  s1L <- utf8ToInt(s1)
  s2L <- utf8ToInt(s2)
    # use lapply insted of for loop
  sc <- sapply(1:L, function(i) sum(s1L[1:i] == s2L[1:i])/i)
  p0 <- 0.5  # was 0.25
  pvals <- sapply(1:L, function(i) {
    x <- sum(s1L[1:i] == s2L[1:i])
    binom.test(x, i, p=p0, alternative="greater")$p.value
  })
  L1 <- which.min(pvals * sc)
  mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -1, baseOnly = TRUE)
  s1g <- gsub("-", "", substr(s1,1,L1))
  s2g <- gsub("-", "", substr(s2,1,L1))
  aln <- pairwiseAlignment(s1g, s2g, type = "global", substitutionMatrix = mat)
  if (plot_it){
    plot(sc, type = "l", ylim=c(0,1))
    abline(v = L1, col = "red", )
    plot(-log10(pvals), type = "l")
    plot(-log10(pvals) * sc, type = "l")
  }
  aln
}

get_positions <- function(ctg_names){
  s <- strsplit(ctg_names, "_")
  s1 <- sapply(s, function(x) as.integer(x[2]))
  s2 <- sapply(s, function(x) as.integer(x[3]))
  cbind(s1, s2)
}

filter_blast <- function(blast_df, min_length = 30, max_evalue = 1e-5, min_identity = 80){
  f1 <- blast_df$qstart == 1  # alignment starts at the beginning of the sequence
  f2 <- blast_df$evalue < max_evalue
  f3 <- blast_df$length > min_length
  f4 <- blast_df$pident > min_identity
  f5 <- blast_df$sstart > 12  # to extract TIR
  blast_df <- blast_df[f1 & f2 & f3 & f4 & f5,]
  # get best sstart for each saccver
  # sort by bitscore
  blast_df <- blast_df[order(blast_df$bitscore, decreasing = TRUE),]
  # get the first row for each saccver
  blast_df <- blast_df[!duplicated(blast_df$saccver),]
  blast_df
}


## Helper to count leading gaps:
count_leading_gaps <- function(seq_char) {
  # Locate the first non-'-' character
  # regexpr("[^-]", seq_char) gives the position of the first character that is not '-'
  # If everything is '-', regexpr returns -1
  pos <- regexpr("[^-]", seq_char)
  if (pos < 0) {
    # If there is no non-gap, the sequence is all '-'
    nchar(seq_char)
  } else {
    pos - 1  # Everything before this position is '-'
  }
}

## Helper to count trailing gaps:
count_trailing_gaps <- function(seq_char) {
  # A simple trick is to reverse the string and do the same counting as leading:
  rev_seq <- paste0(rev(strsplit(seq_char, "")[[1]]), collapse = "")
  pos <- regexpr("[^-]", rev_seq)
  if (pos < 0) {
    nchar(seq_char)
  } else {
    pos - 1
  }
}

## Main function:
count_gaps_at_ends <- function(dna_set) {
  seq_char <- as.character(dna_set)  # Convert DNAStringSet to character vectors
  leading_counts  <- sapply(seq_char, count_leading_gaps)
  trailing_counts <- sapply(seq_char, count_trailing_gaps)

  data.frame(
    sequence = names(dna_set),
    leading  = leading_counts,
    trailing = trailing_counts
  )
}

filter_blast2 <- function (blast_df, min_length = 30, max_evalue = 1e-5, min_identity = 80){
  f2 <- blast_df$evalue < max_evalue
  f3 <- blast_df$length > min_length
  f4 <- blast_df$pident > min_identity
  f5 <- blast_df$sstart > 12  # to extract TIR
  blast_df <- blast_df[f2 & f3 & f4 & f5,]
  # sstart must be less than send
  blast_df <- blast_df[blast_df$sstart < blast_df$send,]
  # remove self hits - id of query is in format "id_start_end", keep only id
  id1 <- gsub("_+","", blast_df$qaccver)
  blast_df <- blast_df[id1 != blast_df$saccver,]
  # keep each query-subject pair only once
  # sort by bitscore
  blast_df <- blast_df[order(blast_df$bitscore, decreasing = TRUE),]
  blast_df <- blast_df[!duplicated(paste(blast_df$qaccver, blast_df$saccver)),]
  blast_df
}

filter_blast3 <- function (blast_df, min_length = 30, max_evalue = 1e-5, min_identity = 80){
  f2 <- blast_df$evalue < max_evalue
  f3 <- blast_df$length > min_length
  f4 <- blast_df$pident > min_identity
  f5 <- blast_df$sstart > 12  # to extract TIR
  blast_df <- blast_df[f2 & f3 & f4 & f5,]
  # sstart must be less than send
  blast_df <- blast_df[blast_df$sstart < blast_df$send,]
  # remove self hits - id of query is in format "id_start_end", keep only id
  id1 <- gsub("_+","", blast_df$qaccver)
  blast_df <- blast_df[id1 != blast_df$saccver,]
  # keep each query-subject pair only once
  # sort by bitscore
  blast_df <- blast_df[order(blast_df$bitscore, decreasing = TRUE),]
  blast_df
}

filter_blast4 <- function(blast_df, upstream=TRUE, min_length = 100, max_offset=10){
  if (upstream){
    c1 <- blast_df$qstart < max_offset
  }else{
    c1 <- blast_df$qend > (blast_df$qlen - max_offset)
  }
  c2 <- blast_df$length > min_length
  blast_df <- blast_df[c1 & c2,]
  blast_df
}

get_cp_from_blast4 <- function(blast_df, upstream=TRUE){
  # split by saccver
  blast_parts <- split(blast_df, blast_df$saccver)
  # calculate coverage for each saccver
  if (upstream){
    cp_list <- lapply(blast_parts, function(x) {
      cl <- table(x$sstart - x$qstart + 1) |> sort(decreasing = TRUE)
      cl[1]
    })
    }else{
    cp_list <- lapply(blast_parts, function(x) {
      cl <- table(x$send -(x$qlen -x$qend)) |> sort(decreasing = TRUE)
      cl[1]
    })
  }
  id <- names(cp_list)
  cp_multiplicity <- as.numeric(unlist(cp_list))
  cp_pos <- as.numeric(sapply(cp_list, names))
  if (upstream){
    cp_df <- data.frame(ID=id, up=cp_pos, mult_up=cp_multiplicity, stringsAsFactors = FALSE)
  }else{
    cp_df <- data.frame(ID=id, down=cp_pos, mult_down=cp_multiplicity, stringsAsFactors = FALSE)
  }
  return(cp_df)
}

get_coverage_from_blast <- function(bl) {
  # spit by saccver
  bl_parts <- split(bl, bl$saccver)
  # calculate coverage for each saccver
  cov_list <- lapply(bl_parts, function(x) {
    gr <- GRanges(seqnames = x$saccver,
                  ranges = IRanges(start = x$sstart, end = x$send),
                  strand = "*")
    gr_cov <- as.vector(coverage(gr)[[1]])
  })
  cov_list
}

find_switch_point_from_blast_coverage <- function(cvrg) {
  switch_info <- beast(cvrg, season = 'none',
                       tcp.minmax = c(1, 2),
                       torder.minmax = c(0, 0),
                       tseg.leftmargin = 500,
                       tseg.rightmargin = 300,
                       quiet = TRUE)
  cp <- switch_info$trend$pos_cp[1]
  cp_neg <- switch_info$trend$pos_cpNeg
  if (!is.null(cp_neg)) {
    if (any(cp_neg > cp)) {
      return(NA)
    }
  }
  if (is.null(cp) | is.na(cp)) {
    return(NA)
  }
  cpPr <- switch_info$trend$pos_cpPr[1]
  abrupt_change <- switch_info$trend$pos_cpAbruptChange[1]
  # coverage before change point
  cvrg0 <- mean(cvrg[1:cp[1]])
  cvrg1 <- mean(cvrg[(cp[1] + 1):length(cvrg)])

  c1 <- cpPr > 0.5
  c4 <- cvrg1 > 9
  c3 <- (cvrg1 + 1) / (cvrg0 + 1) > 5
  if (c1 & c3 & c4) {
    return(cp)
  } else {
    return(NA)
  }
}


find_switch_point_from_blast_coverage2 <- function (cvrg){
  W <- 200
  L <- length(cvrg)
  if (L < 500){
    return(NA)
  }
  swp <- seq(W, L-200, by = 1)
  w1 <- cbind(swp - W + 1, swp)
  w2 <- cbind(swp + 1 , swp + W)
  m1 <- apply(w1, 1, function(x) mean(cvrg[x[1]:x[2]]))
  m2 <- apply(w2, 1, function(x) mean(cvrg[x[1]:x[2]]))
  cp <- which.max((m2 + 1)/(m1 + 2)) + W
  mcov1 <- m1[cp - W]
  mcov2 <- m2[cp - W]
  if (mcov1 < 3 & mcov2 > 8 | mcov1 < 2 & mcov2 > 6){
    return(cp)
  }else{
    return(NA)
  }
}

find_switch_point_from_blast_coverage3 <- function (cvrg){
  W <- 200
  L <- length(cvrg)
  if (L < 500){
    return(NA)
  }
  swp <- seq(W, L-200, by = 1)
  swp <- seq(1, L, by = 1)
  sumsum_left <- cumsum(cvrg)
  sumsum_right <- rev(cumsum(rev(cvrg)))
  # means of coverage 1 to swp and swp to end
  m1 <- sumsum_left/swp
  m2 <- sumsum_right/(L - swp)
  #m12 <- (m2 + 1) / (m1 + 1)
  m12 <- (m2 + mean(cvrg)*0.04) / (m1 + mean(cvrg)*0.04)
  m12[1:W] <- 0
  m12[(L-W):L] <- 0
  cp <- which.max(m12)
  mcov1 <- m1[cp - W]
  mcov2 <- m2[cp - W]
  if (mcov1 < 3 & mcov2 > 20 | mcov1 < 2 & mcov2 > 10){
    return(cp)
  }else{
    return(NA)
  }

  return(cp)
}

swt_pt_function <- function(tir_class){
  defined_functions <- list(
  "Class_II_Subclass_1_TIR_Tc1_Mariner" = find_switch_point_from_blast_coverage2,
  "Class_II_Subclass_1_TIR_PIF_Harbinger" = find_switch_point_from_blast_coverage2,
  "Class_II_Subclass_1_TIR_MuDR_Mutator" = find_switch_point_from_blast_coverage2,
  "Class_II_Subclass_1_TIR_hAT" = find_switch_point_from_blast_coverage2,
  "Class_II_Subclass_1_TIR_EnSpm_CACTA" = find_switch_point_from_blast_coverage3
  )
  out_fn <- defined_functions[[tir_class]]
  if (is.null(out_fn)){
    out_fn <- find_switch_point_from_blast_coverage2
  }
  return(out_fn)
}


detect_tir_pif_harbinger <- function(seq_up, seq_down, cp_up, cp_down,
                                    L = 100, min_aln_score = 8) {
  W <- 10
  tsd_length <- 3
  seq_up_w <- subseq(seq_up, cp_up - L / 2, cp_up + L / 2)
  seq_down_w <- subseq(seq_down, cp_down - L / 2, cp_down + L / 2)
  ID <- names(seq_up)
  tir_up <- getSeq(seq_up_w, GRanges(seqnames = ID,
                                     ranges = IRanges((tsd_length + 1):(L - W),
                                                      width = W),
                                     strand = "*"))
  tir_down <- getSeq(seq_down_w, GRanges(seqnames = ID,
                                         ranges = IRanges((tsd_length + 1):(L - W),
                                                          width = W),
                                         strand = "*"))

  tsd_up3 <- getSeq(seq_up_w, GRanges(seqnames = ID,
                                      ranges = IRanges(1:(L - W - tsd_length),
                                                       width = tsd_length),
                                      strand = "*"))
  tsd_down3 <- getSeq(seq_down_w, GRanges(seqnames = ID,
                                          ranges = IRanges(1:(L - W - tsd_length),
                                                           width = tsd_length),
                                          strand = "*")) |> reverseComplement()
  tir_match <- outer(as.character(tir_up), as.character(tir_down), aln_score)
  tsd_match <- outer(as.character(tsd_up3), as.character(tsd_down3), "==")
  tir_tsd_match <- tir_match * tsd_match
  max_score <- max(tir_tsd_match)
  if (max_score < min_aln_score) {
      return(NA)
  }
  ind <- which(tir_tsd_match == max_score, arr.ind = TRUE)
  N_best <- nrow(ind)
  if (N_best > 1){
    message("more than one max score")
    return(NA)
  }
  cp_left <- ind[1, 1] + cp_up - L / 2 - 1 + tsd_length
  cp_right <- ind[1, 2] + cp_down - L / 2 - 1 + tsd_length
  TIR_full <- eval_aln_length_alt(substring(seq_up, cp_left, cp_left + 100),
                                  substring(seq_down, cp_right, cp_right + 100))
  # test TSD again in longer region
  TSD_test1 <- subseq(seq_up, cp_left - 12, cp_left - 1)
  TSD_test2 <- subseq(seq_down, cp_right - 12, cp_right - 1) |>
    reverseComplement() |>
    as.character()
  TSD_full <- extract_TSD(as.character(TSD_test1), as.character(TSD_test2))
  if (is.na(TSD_full)) {
    return(NA)
  }
  return(list(
    aln_upstream_start = cp_left,
    aln_downstream_start = cp_right,
    TIR_up = as.character(pattern(TIR_full)),
    TIR_down = as.character(subject(TIR_full)),
    TIR_score = score(TIR_full),
    TSD = TSD_full
  ))
}


detect_tir_CACTA <- function(seq_up, seq_down, cp_up, cp_down,
                             L = 400, min_aln_score = 15
) {
  W <- 20
  tsd_length <- 3
  min_aln_score <- 8
  seq_up_w <- subseq(seq_up, cp_up - L / 2, cp_up + L / 2)
  seq_down_w <- subseq(seq_down, cp_down - L / 2, cp_down + L / 2)
  ID <- names(seq_up)
  tir_up <- getSeq(seq_up_w, GRanges(seqnames = ID,
                                     ranges = IRanges((tsd_length + 1):(L - W),
                                                      width = W),
                                     strand = "*"))
  tir_down <- getSeq(seq_down_w, GRanges(seqnames = ID,
                                         ranges = IRanges((tsd_length + 1):(L - W),
                                                          width = W),
                                         strand = "*"))

  tsd_up3 <- getSeq(seq_up_w, GRanges(seqnames = ID,
                                      ranges = IRanges(1:(L - W - tsd_length),
                                                       width = tsd_length),
                                      strand = "*"))
  tsd_down3 <- getSeq(seq_down_w, GRanges(seqnames = ID,
                                          ranges = IRanges(1:(L - W - tsd_length),
                                                           width = tsd_length),
                                          strand = "*")) |> reverseComplement()
  tsd_up2 <- subseq(tsd_up3, 2)
  tsd_down2 <- subseq(tsd_down3, 2)

  cacta_up <- subseq(tir_up, 1, 5) == "CACTA"
  cacta_down <- subseq(tir_down, 1, 5) == "CACTA"
  cacta_match <- outer(cacta_up, cacta_down, "&")
  if (!any(cacta_match)) {
    return(NA)
  }
  tir_match <- outer(as.character(tir_up), as.character(tir_down), aln_score)
  tsd_match_3 <- outer(as.character(tsd_up3), as.character(tsd_down3), "==")
  tsd_match_2 <- outer(as.character(tsd_up2), as.character(tsd_down2), "==")
  tir_tsd_match <- tir_match * tsd_match_3 * cacta_match
  tir_tsd_match <- tir_match * (tsd_match_3 | tsd_match_2) * cacta_match

  max_score <- max(tir_tsd_match)
  if (max_score < min_aln_score) {
    return(NA)
  }

  ind <- which(tir_tsd_match == max_score, arr.ind = TRUE)
  if (length(ind) == 0) {
    message("more than one max score")
    return(NA)
  }
  cp_left <- ind[1, 1] + cp_up - L / 2 - 1 + tsd_length
  cp_right <- ind[1, 2] + cp_down - L / 2 - 1 + tsd_length

  TIR_full <- eval_aln_length_alt(substring(seq_up, cp_left, cp_left + 100),
                                  substring(seq_down, cp_right, cp_right + 100))
  # test TSD again in longer region
  TSD_test1 <- subseq(seq_up, cp_left - 12, cp_left - 1)
  TSD_test2 <- subseq(seq_down, cp_right - 12, cp_right - 1) |>
    reverseComplement() |>
    as.character()
  TSD_full <- extract_TSD(as.character(TSD_test1), as.character(TSD_test2))
  if (is.na(TSD_full)) {
    return(NA)
  }


  return(list(
    aln_upstream_start = cp_left,
    aln_downstream_start = cp_right,
    TIR_up = as.character(pattern(TIR_full)),
    TIR_down = as.character(subject(TIR_full)),
    TIR_score = score(TIR_full),
    TSD = TSD_full
  ))

}

detect_tir_tc1 <- function(seq_up, seq_down, cp_up, cp_down,
                           L = 100, min_aln_score = 30){
  seq_up_w <- subseq(seq_up, cp_up - L / 2, cp_up + L / 2)
  seq_down_w <- subseq(seq_down, cp_down - L / 2, cp_down + L / 2)
  ID <- names(seq_up)
  aln <- pairwiseAlignment(seq_up_w, seq_down_w, type = "local")
  if (score(aln) < min_aln_score) {
    return(NA)
  }
  aln_up_start <- start(pattern(aln))
  aln_down_start <- start(subject(aln))
  TIR_up <- as.character(pattern(aln))
  TIR_down <- as.character(subject(aln))

  # TA TSD is part of the TIR, find it
  # it could be withing the first 3 nt of the TIR
  TSD_window <- 10
  putative_TSD_up <- substring(TIR_up, 1:(TSD_window), 2:(TSD_window + 1))
  putative_TSD_down <- substring(TIR_down, 1:(TSD_window), 2:(TSD_window + 1))
  TIRu_2nt <- substring(TIR_up, 3:(TSD_window + 2), 4:(TSD_window + 3))
  TIRd_2nt <- substring(TIR_down, 3:(TSD_window + 2), 4:(TSD_window + 3))

  TA <- putative_TSD_up == "TA" & putative_TSD_down == "TA" & TIRu_2nt == TIRd_2nt
  if (!any(TA)){
    return(NA)
  }
  first_TA <- which(TA)[1]
  aln_up_start <- aln_up_start + first_TA + tsd_length
  aln_down_start <- aln_down_start + first_TA + tsd_length
  TIR_up <- substring(TIR_up, first_TA + tsd_length)
  TIR_down <- substring(TIR_down, first_TA + tsd_length)
  cp_left <- aln_up_start + cp_up - L / 2 - 2
  cp_right <- aln_down_start + cp_down - L / 2 - 2
  TIR_full <- eval_aln_length_alt(substring(seq_up, cp_left, cp_left + 100),
                                  substring(seq_down, cp_right, cp_right + 100))
  if (score(TIR_full) <= (min_aln_score * 0.5)) {
    return(NA)
  }
  return(list(
    aln_upstream_start = cp_left,
    aln_downstream_start = cp_right,
    TIR_up = as.character(pattern(TIR_full)),
    TIR_down = as.character(subject(TIR_full)),
    TIR_score = score(TIR_full),
    TSD = "TA"
  ))
}

detect_tir <- function(seq_up, seq_down, cp_up, cp_down,
                       max_window = 500, min_aln_score = 25,
                       min_tsd_length = 7, min_window = 50) {
  for (w in seq(min_window, max_window, by = 50)) {
    x_part <- substr(seq_up, cp_up - w, cp_up + w)
    y_part <- substr(seq_down, cp_down - w, cp_down + w)
    # check length of sequences is as expected
    if (nchar(x_part) != 2 * w + 1 | nchar(y_part) != 2 * w + 1) {
      return(NA)
    }

    aln <- pairwiseAlignment(x_part, y_part, type = "local")
    if (score(aln) < min_aln_score) {
      next
    }
    cp_left <- start(pattern(aln)) + cp_up - w - 1
    cp_right <- start(subject(aln)) + cp_down - w - 1
    TSD_test1 <- substr(seq_up, cp_left - 12, cp_left - 1)
    TSD_test2 <- substr(seq_down, cp_right - 12, cp_right - 1) |>
      DNAString() |>
      reverseComplement() |>
      as.character()
    if (nchar(TSD_test1) != 12 | nchar(TSD_test2) != 12) {
      return(NA)
    }
    TSD_full <- extract_TSD(TSD_test1, TSD_test2)
    if (is.na(TSD_full) | nchar(TSD_full) < min_tsd_length) {
      next
    }else {
      TIR_full <- eval_aln_length_alt(substring(seq_up, cp_left, cp_left + 200),
                                      substring(seq_down, cp_right, cp_right + 200))

      if (nchar(pattern(TIR_full)) > nchar(pattern(aln)) & score(TIR_full) > score(aln)) {
        TIR_up <- pattern(TIR_full)
        TIR_down <- subject(TIR_full)
        sc <- score(TIR_full)
      }else {
        TIR_up <- subject(aln)
        TIR_down <- pattern(aln)
        sc <- score(aln)
      }

      return(list(
        aln_upstream_start = cp_left,
        aln_downstream_start = cp_right,
        TIR_up = as.character(TIR_up),
        TIR_down = as.character(TIR_down),
        TIR_score = score(aln),
        TSD = TSD_full,
        w = w
      ))
    }
  }
  return(NA)
}



aln_score <- function(x,y){
  if (length(x) == 0 | length(y) == 0){
    return(0)
  }
  mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -1, baseOnly = FALSE)
  aln <- pairwiseAlignment(x,y, type = "global", substitutionMatrix = mat, scoreOnly = TRUE)
  return(aln)
}

run_blast_tir_analysis <- function(
    query_db,                   # e.g., upstream_db or downstream_db
    out_blast_file,             # e.g., blast_upstream or blast_downstream
    blast_db,                   # same as query_db in your example
    swt_pt_fun,                 # swt_pt_fun[[cls]]
    filter_fun,                 # e.g., filter_blast3
    coverage_fun,               # e.g., get_coverage_from_blast
    evalue          = "1e-10",
    max_target_seqs = 500000L,
    strand          = "plus",
    min_length      = 150,
    min_identity    = 80,
    mc.cores        = 1
) {
  # Run BLAST if output file doesn't exist
  if (!file.exists(out_blast_file)) {
    system(
      paste(
        "blastn -query", query_db,
        "-db", blast_db,
        "-out", out_blast_file,
        "-outfmt 6",
        paste0("-evalue ", evalue),
        paste0("-max_target_seqs ", max_target_seqs),
        paste0("-strand ", strand),
        paste0("-num_threads ", mc.cores)
#        ," -word_size 7 -gapextend 1 -gapopen 2 -reward 1 -penalty -1"
      )
    )
  }
  # Read in the BLAST results
  blast_df <- read.table(out_blast_file, header = FALSE, sep = "\t",
                         col.names = c("qaccver", "saccver", "pident", "length",
                                       "mismatch", "gapopen",
                                       "qstart", "qend", "sstart", "send",
                                       "evalue", "bitscore"))

  # Filter the BLAST results
  blast_df <- filter_fun(blast_df, min_length = min_length, min_identity = min_identity)
  # Compute coverage
  blast_cov <- coverage_fun(blast_df)
  # Compute the switch points in parallel
  cp_vals <- unlist(mclapply(blast_cov, FUN = swt_pt_fun, mc.cores = mc.cores))
  return(list (blast_df = blast_df, blast_cov = blast_cov, cp_vals = cp_vals))
}

process_region_files <- function(file_list, side) {
  info_list <- list()
  CM_list <- list()
  ctg_list <- list()

  if (side == "upstream") {
    pattern <- "_upstream.+"
  } else if (side == "downstream") {
    pattern <- "_downstream.+"
  } else {
    stop("side must be either 'upstream' or 'downstream'")
  }

  for (f in file_list) {
    ctg <- readDNAStringSet(f)
    CM <- consensusMatrix(ctg)[1:4, ]
    if (any(colSums(CM) == 0)) {
      next
    }
    aln_info <- analyze_consensus_matrix(CM)
    cp <- find_switch_point(aln_info)
    classification <- gsub(pattern, "", basename(f))
    mean_coverage <- mean(colSums(CM))
    # (dpos is calculated in your original code, but it is not used later)

    if (!is.na(cp)) {
      contig_id <- gsub("[.]fasta", "", gsub(".+_Contig", "", f))
      info <- extract_info_from_switch_points(aln_info, ctg, cp)
      info$ContigID <- contig_id
      info$Classification <- classification
      info$Filename <- f  # for debugging; remove if not needed
      info_list[[contig_id]] <- info

      # Store consensus matrix starting at the switch point
      CM_list[[contig_id]] <- CM[, cp:ncol(CM)]

      # Extract and filter the contig part
      ctg_part <- subseq(ctg, cp)
      ctg_part <- ctg_part[nchar(gsub("-", "", ctg_part)) > 10]
      ctg_list[[contig_id]] <- ctg_part
    }
  }

  return(list(info = info_list, CM_list = CM_list, ctg_list = ctg_list))
}

# (1) Compute element boundaries from a merged data frame
compute_element_boundaries <- function(df) {
  df$element_start <- ifelse(df$Strand == "+",
                             df$Upstream_start + df$CP_left,
                             df$Upstream_end - df$CP_left)
  df$element_end <- ifelse(df$Strand == "+",
                           df$Downstream_end - df$CP_right,
                           df$Downstream_start + df$CP_right)
  return(df)
}

# (2) Create a GFF-like data frame and return GRanges (assuming boundaries are computed)
create_granges_from_df <- function(df, iter) {
  df_gff <- data.frame(
    seqid           = df$SeqID,
    source          = "DANTE_TIR",
    start           = ifelse(df$element_start < df$element_end, df$element_start, df$element_end),
    end             = ifelse(df$element_start < df$element_end, df$element_end, df$element_start) + 1,
    strand          = df$Strand,
    tir_seq5        = df$TIR_left_aln,
    tir_seq3        = df$TIR_right_aln,
    tsd             = df$TSD,
    ID              = paste(df$Classification, df$ID, sep = "_"),
    Classification  = df$Classification,
    Iter            = iter,
    stringsAsFactors = FALSE
  )
  gr <- makeGRangesFromDataFrame(df_gff, keep.extra.columns = TRUE)
  return(gr)
}

# (3) Merge with flank coordinates, compute boundaries, and generate GRanges
prepare_granges <- function(res_df, tir_flank_coordinates, iter) {
  merged_df <- merge(res_df, tir_flank_coordinates, by = "ID")
  merged_df <- compute_element_boundaries(merged_df)
  gr <- create_granges_from_df(merged_df, iter)
  return(gr)
}


cluster_tir_sequences <- function(genome_file, gr_fin, output, threads) {
  # Read the genome sequences.

  genome <- readDNAStringSet(genome_file)
  names(genome) <- gsub(" .*", "", names(genome))
  # Extract TIR sequences from gr_fin.
  tir_seqs <- getSeq(genome, gr_fin)
  names(tir_seqs) <- paste0(
    gsub("Class_II_Subclass_1_TIR_", "", gr_fin$ID),
    "#",
    gsub("Class_II_Subclass_1_", "Class_II/Subclass_1/", gr_fin$Classification)
  )
  # Split TIR sequences by classification.
  tir_seqs_part <- split(tir_seqs, gsub("Class_II_Subclass_1_TIR_", "", gr_fin$Classification))
  # Create output FASTA filenames for each classification.
  tir_seqs_parts_files <- paste0(output, "/DANTE_TIR_", names(tir_seqs_part), ".fasta")
  names(tir_seqs_parts_files) <- names(tir_seqs_part)
  # Write each classification-specific FASTA file.
  for (cls in names(tir_seqs_part)) {
    writeXStringSet(tir_seqs_part[[cls]], tir_seqs_parts_files[cls])
  }
  # Ensure the mmseqs2 output directory exists.
  mmseqs2_dir <- file.path(output, "mmseqs2")
  dir.create(mmseqs2_dir, showWarnings = FALSE, recursive = TRUE)
  tmp_dir <- tempdir()
  tir_cls <- list()
  cls_size_list <- list()
  # Loop over each classification file and run mmseqs2 clustering.
  for (i in seq_along(tir_seqs_parts_files)) {
    class_name <- names(tir_seqs_parts_files)[i]  # classification name (already stripped)
    out_prefix <- file.path(mmseqs2_dir, class_name)
    cluster_file <- paste0(out_prefix, "_cluster.tsv")
    rep_seq_file <- paste0(out_prefix, ".rep_seq.fasta")

    # Build the mmseqs2 command.
    cmd <- paste(
      "mmseqs easy-cluster",
      tir_seqs_parts_files[i],
      out_prefix,
      tmp_dir,
      "--cluster-mode 0 -v 1 --min-seq-id 0.8",
      "--cov-mode 0 --mask 0 -s 6 --threads", threads
    )
    print(cmd)
    system(cmd, intern = TRUE, ignore.stderr=FALSE)

    # Read clustering output.
    df <- read.table(cluster_file, header = FALSE, sep = "\t",
                     col.names = c("Representative_ID", "ID"),
                     comment.char = "")
    df$Cluster_ID <- as.numeric(factor(df$Representative_ID))
    cls_size <- table(df$Cluster_ID)
    cls_size_list[[class_name]] <- sort(cls_size, decreasing = TRUE)
    df$Multiplicity <- cls_size[df$Cluster_ID]
    df$ID2 <- paste("Class_II_Subclass_1_TIR", gsub("#.+", "", df$ID), sep = "_")
    tir_cls[[i]] <- df
  }

  tir_cls_df <- do.call(rbind, tir_cls)
  rownames(tir_cls_df) <- tir_cls_df$ID2

  # Update the GRanges object with clustering information.
  gr_fin$Multiplicity <- tir_cls_df[gr_fin$ID, "Multiplicity"]
  gr_fin$Cluster_ID   <- tir_cls_df[gr_fin$ID, "Cluster_ID"]

  # Return the updated GRanges, clustering results, and tir_seqs.
  return(list(gr_fin = gr_fin, tir_cls_df = tir_cls_df, tir_seqs = tir_seqs))
}

blastn_db_exists <- function(fasta_file) {
  db_file <- paste0(fasta_file, ".nsq")
  return(file.exists(paste0(db_file, ".nsq")))
}

# Helper function that creates BLAST databases, runs BLAST, reads, filters, and synchronizes the outputs.
blast_helper <- function(query_up, query_down, db_up, db_down,
                         blast_up_out, blast_down_out,
                         blast_args_up, blast_args_down,
                         col_names, filter_fun_up, filter_fun_down) {
  # check if the BLAST databases already exist
  if (!blastn_db_exists(db_up)) {
    system(paste("makeblastdb -in", db_up, "-dbtype nucl -out", db_up), intern = TRUE)
  }
  # Create BLAST databases for upstream and downstream regions.
  if (!blastn_db_exists(db_down)){
    system(paste("makeblastdb -in", db_down, "-dbtype nucl -out", db_down), intern = TRUE)
  }
  # Run BLAST searches for upstream and downstream queries.
  system(paste("blastn -query", query_up, "-db", db_up, "-out", blast_up_out, blast_args_up), intern = TRUE)
  system(paste("blastn -query", query_down, "-db", db_down, "-out", blast_down_out, blast_args_down), intern = TRUE)

  # Read BLAST outputs.
  blast_up_df <- read.table(blast_up_out, header = FALSE, sep = "\t",
                            col.names = col_names, comment.char = "")
  blast_down_df <- read.table(blast_down_out, header = FALSE, sep = "\t",
                              col.names = col_names, comment.char = "")

  # Apply filtering functions.
  blast_up_df <- filter_fun_up(blast_up_df)
  blast_down_df <- filter_fun_down(blast_down_df)

  # Keep only common subject IDs.
  common_ids <- intersect(blast_up_df$saccver, blast_down_df$saccver)
  blast_up_df <- blast_up_df[blast_up_df$saccver %in% common_ids, ]
  blast_down_df <- blast_down_df[blast_down_df$saccver %in% common_ids, ]

  # Order data frames by subject ID to synchronize rows.
  blast_up_df <- blast_up_df[order(blast_up_df$saccver), ]
  blast_down_df <- blast_down_df[order(blast_down_df$saccver), ]

  return(list(blast_up_df = blast_up_df, blast_down_df = blast_down_df))
}

round1 <- function(contig_dir, tir_flank_file) {
  # Read TIR flank coordinates
  tir_flank_coordinates <- read.table(tir_flank_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  # Identify upstream and downstream contig files
  upstream_contigs <- dir(contig_dir, pattern = "upstream_Contig", full.names = TRUE)
  downstream_contigs <- dir(contig_dir, pattern = "downstream_Contig", full.names = TRUE)

  message("Processing upstream regions")
  upstream_results <- process_region_files(upstream_contigs, "upstream")
  # Merge upstream info (assuming process_region_files() returns a list with element 'info')
  upstream_info <- do.call(rbind, upstream_results$info)
  ctg_list_upstream <- upstream_results$ctg_list

  message("Processing downstream regions")
  downstream_results <- process_region_files(downstream_contigs, "downstream")
  downstream_info <- do.call(rbind, downstream_results$info)
  ctg_list_downstream <- downstream_results$ctg_list

  # Find common IDs between upstream and downstream results
  both_side_id <- intersect(downstream_info$ID, upstream_info$ID)

  res <- list()
  for (id in both_side_id) {
    # Subset rows for the current id
    downstream <- downstream_info[downstream_info$ID == id, ]
    upstream   <- upstream_info[upstream_info$ID == id, ]
    # Ensure exactly one match on each side
    if (nrow(downstream) != 1 || nrow(upstream) != 1) next

    # Evaluate TIR alignment and extract TSD
    tir_aln <- eval_aln_length_alt(upstream$TIR, downstream$TIR)
    tsd <- extract_TSD(upstream$TSD, as.character(reverseComplement(DNAString(downstream$TSD))))
    if (!is.na(tsd) && score(tir_aln) > 0) {
      res[[length(res) + 1]] <- list(
        ID = id,
        Classification = downstream$Classification,
        start_left = upstream$start,
        end_left = upstream$end,
        CP_left = upstream$CP,
        start_right = downstream$start,
        end_right = downstream$end,
        CP_right = downstream$CP,
        TIR_left_aln = as.character(pattern(tir_aln)),
        TIR_right_aln = as.character(subject(tir_aln)),
        TIR_score = score(tir_aln),
        TSD = tsd,
        upstream_contig = upstream$ContigID,
        downstream_contig = downstream$ContigID
      )
    }
  }

  # Combine the list of results into a data frame
  res_df <- do.call(rbind, lapply(res, as.data.frame))

  # Create GRanges for round 1 using the helper function prepare_granges()
  gr1 <- prepare_granges(res_df, tir_flank_coordinates, iter = 1)

  message("------------------------------------------------------------------")
  message("Number of elements found in the first round: ", length(gr1))

  # Return outputs needed by subsequent rounds
  return(list(
    gr1 = gr1,
    res_df = res_df,
    ctg_list_upstream = ctg_list_upstream,
    ctg_list_downstream = ctg_list_downstream,
    tir_flank_coordinates = tir_flank_coordinates
  ))
}

### Updated round2() using the new blast_helper() ###
round2 <- function(res_df, ctg_list_upstream, ctg_list_downstream, gr1,
                   tir_flank_coordinates, output, threads) {
  message("Identification of elements - Round 2")
  ## Retrieve passed contigs from Round 1
  ctg_upstream_pass <- ctg_list_upstream[unique(res_df$upstream_contig)]
  ctg_upstream_pass_class <- res_df$Classification[match(names(ctg_upstream_pass), res_df$upstream_contig)]
  ctg_downstream_pass <- ctg_list_downstream[unique(res_df$downstream_contig)]
  ctg_downstream_pass_class <- res_df$Classification[match(names(ctg_downstream_pass), res_df$downstream_contig)]

  ## Build consensus sequences for upstream and downstream contigs.
  cons_seqs <- list()
  for (i in seq_along(ctg_upstream_pass)) {
    label <- paste0(ctg_upstream_pass_class[i], "_", names(ctg_upstream_pass)[i], "_upstream")
    cons_seqs[[label]] <- get_consensus_from_aln(ctg_upstream_pass[[i]], 60)
  }
  for (i in seq_along(ctg_downstream_pass)) {
    label <- paste0(ctg_downstream_pass_class[i], "_", names(ctg_downstream_pass)[i], "_downstream")
    cons_seqs[[label]] <- get_consensus_from_aln(ctg_downstream_pass[[i]], 60)
  }
  cons_seq_out <- unlist(DNAStringSetList(cons_seqs))

  ## Extract metadata from consensus sequence names.
  side <- sapply(strsplit(names(cons_seqs), "_"), function(x) tail(x, 1))
  cons_class <- sapply(strsplit(names(cons_seqs), "_"), function(x) paste(head(x, -2), collapse = "_"))

  ## Create a data frame of output file names (one row per classification).
  unique_classes <- unique(cons_class)
  cons_file_df <- data.frame(
    classification = unique_classes,
    upstream_cons = NA,
    downstream_cons = NA,
    upstream_db = NA,
    downstream_db = NA,
    blast_downstream = NA,
    blast_upstream = NA,
    stringsAsFactors = FALSE
  )

  dir.create(file.path(output, "TIR_plus_consensus"), showWarnings = FALSE)
  for (i in seq_along(unique_classes)) {
    cls <- unique_classes[i]
    cons_file_df$upstream_cons[i] <- file.path(output, "TIR_plus_consensus", paste0(cls, "_upstream_consensus.fasta"))
    cons_file_df$downstream_cons[i] <- file.path(output, "TIR_plus_consensus", paste0(cls, "_downstream_consensus.fasta"))
    cons_file_df$upstream_db[i] <- file.path(output, paste0(cls, "_upstream_regions.fasta"))
    cons_file_df$downstream_db[i] <- file.path(output, paste0(cls, "_downstream_regions.fasta"))

    upstream_cons <- cons_seq_out[cons_class == cls & side == "upstream"]
    downstream_cons <- cons_seq_out[cons_class == cls & side == "downstream"]
    writeXStringSet(upstream_cons, cons_file_df$upstream_cons[i])
    writeXStringSet(downstream_cons, cons_file_df$downstream_cons[i])
  }

  dir.create(file.path(output, "blastn"), showWarnings = FALSE)
  # Define BLAST output column names.
  col_names <- c("qaccver", "saccver", "pident", "length", "mismatch",
                 "gapopen", "qstart", "qend", "sstart", "send", "evalue",
                 "bitscore")
  blast_args <- "-outfmt 6 -evalue 1e-5"
  res2_list <- list()
  gff2_list <- list()
  for (i in seq_along(unique_classes)) {
    cls <- unique_classes[i]
    upstream_cons <- cons_file_df$upstream_cons[i]
    downstream_cons <- cons_file_df$downstream_cons[i]
    upstream_db <- cons_file_df$upstream_db[i]
    downstream_db <- cons_file_df$downstream_db[i]
    blast_upstream <- file.path(output, "blastn", paste0(cls, "_upstream.blastn"))
    blast_downstream <- file.path(output, "blastn", paste0(cls, "_downstream.blastn"))

    # Call blast_helper() to create databases, run BLAST, and process outputs.
    blast_res <- blast_helper(query_up = upstream_cons,
                              query_down = downstream_cons,
                              db_up = upstream_db,
                              db_down = downstream_db,
                              blast_up_out = blast_upstream,
                              blast_down_out = blast_downstream,
                              blast_args_up = blast_args,
                              blast_args_down = blast_args,
                              col_names = col_names,
                              filter_fun_up = filter_blast,
                              filter_fun_down = filter_blast)

    blast_up_df <- blast_res$blast_up_df
    blast_down_df <- blast_res$blast_down_df

    upstream_regions <- readDNAStringSet(upstream_db)
    downstream_regions <- readDNAStringSet(downstream_db)

    gr_up_TIR <- GRanges(seqnames = blast_up_df$saccver,
                         ranges = IRanges(start = blast_up_df$sstart, end = blast_up_df$sstart + 100),
                         strand = "*")
    gr_down_TIR <- GRanges(seqnames = blast_down_df$saccver,
                           ranges = IRanges(start = blast_down_df$sstart, end = blast_down_df$sstart + 100),
                           strand = "*")
    gr_up_TSD <- GRanges(seqnames = blast_up_df$saccver,
                         ranges = IRanges(start = blast_up_df$sstart - 12, end = blast_up_df$sstart - 1),
                         strand = "*")
    gr_down_TSD <- GRanges(seqnames = blast_down_df$saccver,
                           ranges = IRanges(start = blast_down_df$sstart - 12, end = blast_down_df$sstart - 1),
                           strand = "*")

    blast_up_df$TIR <- as.character(getSeq(upstream_regions, gr_up_TIR))
    blast_down_df$TIR <- as.character(getSeq(downstream_regions, gr_down_TIR))
    blast_up_df$TSD <- as.character(getSeq(upstream_regions, gr_up_TSD))
    blast_down_df$TSD <- as.character(getSeq(downstream_regions, gr_down_TSD))

    res2 <- list()
    for (j in seq_along(blast_up_df$TIR)) {
      tir_aln <- eval_aln_length_alt(blast_up_df$TIR[j], blast_down_df$TIR[j])
      tsd <- extract_TSD(blast_up_df$TSD[j], as.character(reverseComplement(DNAString(blast_down_df$TSD[j]))))
      if (!is.na(tsd) && score(tir_aln) > 0) {
        res2[[j]] <- list(
          ID = blast_up_df$saccver[j],
          Classification = cons_file_df$classification[i],
          CP_left = blast_up_df$sstart[j],
          CP_right = blast_down_df$sstart[j],
          TIR_left_aln = as.character(pattern(tir_aln)),
          TIR_right_aln = as.character(subject(tir_aln)),
          TIR_score = score(tir_aln),
          TSD = tsd
        )
      }
    }
    res2_df <- do.call(rbind, lapply(res2, as.data.frame))
    res2_df <- merge(res2_df, tir_flank_coordinates, by = "ID")
    res2_df <- compute_element_boundaries(res2_df)

    df_gff <- data.frame(
      seqid = res2_df$SeqID,
      source = "DANTE_TIR",
      start = ifelse(res2_df$element_start < res2_df$element_end, res2_df$element_start, res2_df$element_end),
      end = ifelse(res2_df$element_start < res2_df$element_end, res2_df$element_end, res2_df$element_start) + 1,
      strand = res2_df$Strand,
      tir_seq5 = res2_df$TIR_left_aln,
      tir_seq3 = res2_df$TIR_right_aln,
      tsd = res2_df$TSD,
      ID = paste(res2_df$Classification, res2_df$ID, sep = "_"),
      Classification = res2_df$Classification,
      Iter = 2,
      stringsAsFactors = FALSE
    )
    gff2_list[[cons_file_df$classification[i]]] <- makeGRangesFromDataFrame(df_gff, keep.extra.columns = TRUE)
    res2_list[[cons_file_df$classification[i]]] <- res2_df
  } # end for each classification

  ## Combine the GRanges from all classifications into gr2,
  ## remove any entries already in gr1, and then combine with gr1 to obtain gr_fin.
  gr2 <- unlist(GRangesList(gff2_list))
  gr2_unique <- gr2[!gr2$ID %in% gr1$ID, ]
  gr_fin <- sort(c(gr1, gr2_unique), by = ~ seqnames * start)
  message("-----------------------------------------------------------------")
  message("Number of elements found in the second round: ", length(gr2_unique))
  return(list(gr2 = gr2, gr_fin = gr_fin))
}

make_detection_worker <- function(detection_fun, cls, seq_upstream, seq_downstream, cp_boundaries,
                                    boundary_mode = c("byName", "byIndex"), ...) {
  boundary_mode <- match.arg(boundary_mode)
  extra_params <- list(...)

  if (boundary_mode == "byName") {
    # Assumes cp_boundaries is indexed by ID matching names(seq_upstream)
    return(function(i) {
      ID <- names(seq_upstream)[i]
      args_list <- c(list(seq_upstream[i],
                          seq_downstream[i],
                          cp_boundaries[ID, 1],
                          cp_boundaries[ID, 2]),
                     extra_params)
      tir_info <- do.call(detection_fun, args_list)
      if (is.na(tir_info[[1]])) return(NULL)
      list(ID = ID,
           Classification = cls,
           CP_left = tir_info$aln_upstream_start,
           CP_right = tir_info$aln_downstream_start,
           TIR_left_aln = tir_info$TIR_up,
           TIR_right_aln = tir_info$TIR_down,
           TIR_score = tir_info$TIR_score,
           TSD = tir_info$TSD)
    })
  } else {
    # Assumes cp_boundaries is a data frame with columns 'up', 'down', and 'ID'
    return(function(i) {
      args_list <- c(list(seq_upstream[i],
                          seq_downstream[i],
                          cp_boundaries$up[i],
                          cp_boundaries$down[i]),
                     extra_params)
      tir_info <- do.call(detection_fun, args_list)
      if (is.na(tir_info[[1]])) return(NULL)
      list(ID = cp_boundaries$ID[i],
           Classification = cls,
           CP_left = tir_info$aln_upstream_start,
           CP_right = tir_info$aln_downstream_start,
           TIR_left_aln = tir_info$TIR_up,
           TIR_right_aln = tir_info$TIR_down,
           TIR_score = tir_info$TIR_score,
           TSD = tir_info$TSD)
    })
  }
}

round3 <- function(contig_dir, output, tir_flank_coordinates, gr_fin, threads) {
  message("Identification of elements - Round 3")

  upstream_regions_file <- dir(contig_dir, pattern = "upstream_regions.fasta$", full.names = TRUE)
  downstream_regions_file <- dir(contig_dir, pattern = "downstream_regions.fasta$", full.names = TRUE)
  names(upstream_regions_file) <- gsub("_upstream_regions.fasta", "", basename(upstream_regions_file))
  names(downstream_regions_file) <- gsub("_downstream_regions.fasta", "", basename(downstream_regions_file))

  tir_class_defined <- c("Class_II_Subclass_1_TIR_Tc1_Mariner",
                         "Class_II_Subclass_1_TIR_PIF_Harbinger",
                         "Class_II_Subclass_1_TIR_MuDR_Mutator",
                         "Class_II_Subclass_1_TIR_hAT",
                         "Class_II_Subclass_1_TIR_EnSpm_CACTA")

  res3_list <- list()

  for (cls in names(upstream_regions_file)) {
    upstream_db   <- upstream_regions_file[[cls]]
    downstream_db <- downstream_regions_file[[cls]]
    blast_upstream   <- file.path(output, "blastn", paste0(cls, "_upstream3.blastn"))
    blast_downstream <- file.path(output, "blastn", paste0(cls, "_downstream3.blastn"))

    cp_detect_info_up <- run_blast_tir_analysis(
      query_db       = upstream_db,
      out_blast_file = blast_upstream,
      blast_db       = upstream_db,
      swt_pt_fun     = swt_pt_function(cls),
      filter_fun     = filter_blast3,
      coverage_fun   = get_coverage_from_blast,
      mc.cores       = threads
    )
    cp_detect_info_down <- run_blast_tir_analysis(
      query_db       = downstream_db,
      out_blast_file = blast_downstream,
      blast_db       = downstream_db,
      swt_pt_fun     = swt_pt_function(cls),
      filter_fun     = filter_blast3,
      coverage_fun   = get_coverage_from_blast,
      mc.cores       = threads
    )
    cp_up <- cp_detect_info_up$cp_vals
    cp_down <- cp_detect_info_down$cp_vals
    valid_names <- intersect(names(cp_up), names(cp_down))

    cp_up_down <- data.frame(
      up = cp_up[valid_names],
      down = cp_down[valid_names],
      row.names = valid_names,
      stringsAsFactors = FALSE
    )
    cp_up_down <- cp_up_down[complete.cases(cp_up_down), , drop = FALSE]
    if (nrow(cp_up_down) == 0) {
      message("No valid boundaries found for ", cls)
      next
    }

    seq_upstream <- readDNAStringSet(upstream_db)
    seq_downstream <- readDNAStringSet(downstream_db)
    # Reorder sequences based on the rownames of cp_up_down:
    seq_upstream <- seq_upstream[match(rownames(cp_up_down), names(seq_upstream))]
    seq_downstream <- seq_downstream[match(rownames(cp_up_down), names(seq_downstream))]

    # Use the helper to generate a worker function. Here we use boundary_mode = "byName"
    if (cls == "Class_II_Subclass_1_TIR_PIF_Harbinger") {
      worker_fun <- make_detection_worker(detection_fun = detect_tir_pif_harbinger,
                                            cls = cls,
                                            seq_upstream = seq_upstream,
                                            seq_downstream = seq_downstream,
                                            cp_boundaries = cp_up_down,
                                            boundary_mode = "byName",
                                            L = 100, min_aln_score = 8)
    } else if (cls == "Class_II_Subclass_1_TIR_Tc1_Mariner") {
      worker_fun <- make_detection_worker(detection_fun = detect_tir_tc1,
                                            cls = cls,
                                            seq_upstream = seq_upstream,
                                            seq_downstream = seq_downstream,
                                            cp_boundaries = cp_up_down,
                                            boundary_mode = "byName",
                                            L = 100, min_aln_score = 20)
    } else if (cls == "Class_II_Subclass_1_TIR_EnSpm_CACTA") {
      worker_fun <- make_detection_worker(detection_fun = detect_tir_CACTA,
                                            cls = cls,
                                            seq_upstream = seq_upstream,
                                            seq_downstream = seq_downstream,
                                            cp_boundaries = cp_up_down,
                                            boundary_mode = "byName",
                                            L = 400, min_aln_score = 15)
    } else if (cls == "Class_II_Subclass_1_TIR_MuDR_Mutator") {
      worker_fun <- make_detection_worker(detection_fun = detect_tir,
                                            cls = cls,
                                            seq_upstream = seq_upstream,
                                            seq_downstream = seq_downstream,
                                            cp_boundaries = cp_up_down,
                                            boundary_mode = "byName",
                                            min_tsd_length = 7, min_aln_score = 25, max_window = 500)
    } else if (cls == "Class_II_Subclass_1_TIR_hAT") {
      worker_fun <- make_detection_worker(detection_fun = detect_tir,
                                            cls = cls,
                                            seq_upstream = seq_upstream,
                                            seq_downstream = seq_downstream,
                                            cp_boundaries = cp_up_down,
                                            boundary_mode = "byName",
                                            min_tsd_length = 6, min_aln_score = 20, max_window = 500)
    } else if (!cls %in% tir_class_defined) {
      message("Processing ", cls)
      worker_fun <- make_detection_worker(detection_fun = detect_tir,
                                            cls = cls,
                                            seq_upstream = seq_upstream,
                                            seq_downstream = seq_downstream,
                                            cp_boundaries = cp_up_down,
                                            boundary_mode = "byName",
                                            min_tsd_length = 8, min_aln_score = 15, max_window = 500)
    }

    parallel_res <- mclapply(seq_along(seq_upstream), worker_fun, mc.cores = threads)
    res3_list <- c(res3_list, Filter(Negate(is.null), parallel_res))
  }

  res3_df <- do.call(rbind, lapply(res3_list, as.data.frame))
  gr3 <- prepare_granges(res3_df, tir_flank_coordinates, iter = 3)
  gr3_unique <- gr3[!gr3$ID %in% gr_fin$ID, ]
  gr_fin <- c(gr_fin, gr3_unique)
  message("-------------------------------------------------------------------")
  message("Number of elements found in the third round: ", length(gr3_unique))
  return(list(gr3 = gr3, gr3_unique = gr3_unique, gr_fin = gr_fin))
}

round4 <- function(gr_fin, tir_cls_df, tir_seqs, tir_flank_coordinates, output, threads, multiplicity_threshold = 5) {
  message("Identification of elements - Round 4")

  repre_ID <- unique(tir_cls_df$Representative_ID[tir_cls_df$Multiplicity >= multiplicity_threshold])
  repre_seqs <- tir_seqs[names(tir_seqs) %in% repre_ID]
  repre_cls <- paste0("Class_II_Subclass_1_", gsub(".+/", "", names(repre_seqs)))
  repre_seqs_groups <- split(repre_seqs, repre_cls)

  repre_seqs_files <- paste0(output, "/mmseqs2/", names(repre_seqs_groups), "_repre.fasta")
  names(repre_seqs_files) <- names(repre_seqs_groups)

  for (cls in names(repre_seqs_groups)) {
    writeXStringSet(repre_seqs_groups[[cls]], repre_seqs_files[cls])
    system(paste("makeblastdb -in", repre_seqs_files[cls], "-dbtype nucl"), intern = TRUE)
  }

  blast_cols <- c("qaccver", "saccver", "pident", "length", "mismatch",
                  "gapopen", "qstart", "qend", "sstart", "send", "evalue",
                  "bitscore", "slen", "qlen")
  outfmt_str <- paste0(' -outfmt "6 ', paste(blast_cols, collapse = " "), '"')

  res4_list <- list()

  for (cls in names(repre_seqs_groups)) {
    upstream_db <- paste0(output, "/", cls, "_upstream_regions.fasta")
    downstream_db <- paste0(output, "/", cls, "_downstream_regions.fasta")
    blast_upstream <- paste0(output, "/blastn/", cls, "_upstream_round4.blastn")
    blast_downstream <- paste0(output, "/blastn/", cls, "_downstream_round4.blastn")

    blast_args_up <- paste0(outfmt_str, " -evalue 1e-5 -strand plus")
    blast_args_down <- paste0(outfmt_str, " -evalue 1e-5 -strand minus")

    blast_res <- blast_helper(query_up = repre_seqs_files[cls],
                              query_down = repre_seqs_files[cls],
                              db_up = upstream_db,
                              db_down = downstream_db,
                              blast_up_out = blast_upstream,
                              blast_down_out = blast_downstream,
                              blast_args_up = blast_args_up,
                              blast_args_down = blast_args_down,
                              col_names = blast_cols,
                              filter_fun_up = function(x) filter_blast4(x, upstream = TRUE, max_offset = 120),
                              filter_fun_down = function(x) filter_blast4(x, upstream = FALSE, max_offset = 120)
    )
    blast_up_df <- blast_res$blast_up_df
    blast_down_df <- blast_res$blast_down_df

    cp_up <- get_cp_from_blast4(blast_up_df, upstream = TRUE)
    cp_down <- get_cp_from_blast4(blast_down_df, upstream = FALSE)
    cp_up_down <- merge(cp_up, cp_down, by = "ID")
    already_annotated <- cp_up_down$ID %in% gsub(".+_", "", gr_fin$ID)
    cp_up_down <- cp_up_down[!already_annotated, ]
    rownames(cp_up_down) <- cp_up_down$ID
    if (nrow(cp_up_down) == 0) {
      message("No additional elements found for ", cls, " in round 4")
      next
    }

    seq_upstream <- readDNAStringSet(upstream_db)
    seq_downstream <- readDNAStringSet(downstream_db)
    seq_upstream <- seq_upstream[match(rownames(cp_up_down), names(seq_upstream))]
    seq_downstream <- seq_downstream[match(rownames(cp_up_down), names(seq_downstream))]

    # For round4 we use boundary_mode = "byIndex" since cp_up_down has columns "up", "down", and "ID"
    if (cls == "Class_II_Subclass_1_TIR_hAT") {
      worker_fun <- make_detection_worker(detection_fun = detect_tir,
                                            cls = cls,
                                            seq_upstream = seq_upstream,
                                            seq_downstream = seq_downstream,
                                            cp_boundaries = cp_up_down,
                                            boundary_mode = "byIndex",
                                            min_tsd_length = 6, min_aln_score = 13,
                                            max_window = 50, min_window = 20)
    } else if (cls == "Class_II_Subclass_1_TIR_EnSpm_CACTA") {
      worker_fun <- make_detection_worker(detection_fun = detect_tir_CACTA,
                                            cls = cls,
                                            seq_upstream = seq_upstream,
                                            seq_downstream = seq_downstream,
                                            cp_boundaries = cp_up_down,
                                            boundary_mode = "byIndex",
                                            L = 50, min_aln_score = 15)
    } else if (cls == "Class_II_Subclass_1_TIR_MuDR_Mutator") {
      worker_fun <- make_detection_worker(detection_fun = detect_tir,
                                            cls = cls,
                                            seq_upstream = seq_upstream,
                                            seq_downstream = seq_downstream,
                                            cp_boundaries = cp_up_down,
                                            boundary_mode = "byIndex",
                                            min_tsd_length = 7, min_aln_score = 25,
                                            max_window = 50)
    } else if (cls == "Class_II_Subclass_1_TIR_Tc1_Mariner") {
      worker_fun <- make_detection_worker(detection_fun = detect_tir_tc1,
                                            cls = cls,
                                            seq_upstream = seq_upstream,
                                            seq_downstream = seq_downstream,
                                            cp_boundaries = cp_up_down,
                                            boundary_mode = "byIndex",
                                            L = 50, min_aln_score = 20)
    } else if (cls == "Class_II_Subclass_1_TIR_PIF_Harbinger") {
      worker_fun <- make_detection_worker(detection_fun = detect_tir_pif_harbinger,
                                            cls = cls,
                                            seq_upstream = seq_upstream,
                                            seq_downstream = seq_downstream,
                                            cp_boundaries = cp_up_down,
                                            boundary_mode = "byIndex",
                                            L = 50, min_aln_score = 8)
    }

    parallel_res <- mclapply(seq_along(seq_upstream), worker_fun, mc.cores = threads)
    res4_list <- c(res4_list, Filter(Negate(is.null), parallel_res))
  }

  res4_df <- do.call(rbind, lapply(res4_list, as.data.frame))
  gr4 <- prepare_granges(res4_df, tir_flank_coordinates, iter = 4)
  gr_fin <- c(gr_fin, gr4)
  message("-------------------------------------------------------------------")
  message("Number of elements found in the fourth round: ", length(gr4))
  return(list(gr4 = gr4, gr_fin = gr_fin, res4_df = res4_df))
}
