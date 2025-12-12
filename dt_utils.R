
concatGRanges <- function(gr1, gr2) {
  # 1. Build superset of seqlevels
  sup <- union(seqlevels(gr1), seqlevels(gr2))
  sup_si <- Seqinfo(seqnames = sup)

  # 2. Map superset back to each object:
  #    positions in 'sup' matching existing levels, NA for new ones
  map1 <- match(sup, seqlevels(gr1))
  map2 <- match(sup, seqlevels(gr2))

  # 3. Assign Seqinfo with mapping (adds missing levels)
  seqinfo(gr1, new2old = map1) <- sup_si
  seqinfo(gr2, new2old = map2) <- sup_si
  c(gr1, gr2)
}

get_consensus_from_aln <- function (aln, perc=70, clean = FALSE){
  nucleotides <- c("A", "C", "G", "T", "N")
  CM <- consensusMatrix(aln)[nucleotides,]
  if (clean){
     CM <- CM[, colSums(CM) > 0]
  }
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



find_switch_point <- function (aln_info, plotit = FALSE, mcmc_seed = 42, n_beast_iter = 1) {
  valid_cps <- c()

  for (iter in 1:n_beast_iter) {
    # Generate a unique seed for each iteration
    current_seed <- mcmc_seed + iter - 1

    invisible(capture.output(
    switch_info <- beast(aln_info$InformationContent, season = 'none',
                         tcp.minmax = c(1,2),
                         torder.minmax = c(0,0),
                         tseg.leftmargin = 10,
                         tseg.rightmargin = 10,
                         mcmc.seed = current_seed,
                         quiet         = TRUE)
    ))
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
      next
    }
    c2 <- cp[1] > 12 && cp[1] < 120
    c3 <- abrupt_change[1] > 0.2
    c4 <- cpPr[1] > 0.3
    if (c1 & c2 & c3 & c4) {
      Ncount <- sum(unlist(strsplit(substr(aln_info$cons,cp[1] - 10, cp[1]-1),""))=="N")
      if (Ncount < 6){
        next
      }
      if (plotit){
        plot(switch_info)
      }
      valid_cps <- c(valid_cps, cp[1])
    }
  }

  if (length(valid_cps) == 0) {
    return(NA)
  } else {
    # Return the most frequent switch point, or the first one if all are unique
    return(as.numeric(names(sort(table(valid_cps), decreasing = TRUE)[1])))
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
  # check if cp + 50 is not out of bounds
  if ((cp + 50) > width(ctg)[1]){
    ctg_donwstream_cp <- subseq(ctg, cp , width(ctg)[1])
  }else{
    ctg_donwstream_cp <- subseq(ctg, cp , cp + 50)
  }
  ctg_donwstream_cp_n_bases <- nchar(gsub("-", "", ctg_donwstream_cp))
  pass <- ctg_upstream_cp_n_bases >= 8 & cp_occupied & ctg_donwstream_cp_n_bases >= 45
  if (sum(pass) == 0){
    return(NULL)
  }

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
  # check if cp + 50 is not out of bounds
  if (cp + 50 > width(ctg_pass)[1]){
    TIR <- subseq(ctg_pass, cp, width(ctg_pass)[1])

    # add "-" to extend it ot 51
    TIR2 <- DNAStringSet(paste0(TIR, paste(rep("-", 51 - width(TIR)[1]), collapse = "")))
    names(TIR2) <- names(TIR)
    TIR <- TIR2
  }else{
    TIR <- subseq(ctg_pass, cp, cp + 50)
  }
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
  # mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -1, baseOnly = TRUE)
  mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -1, baseOnly = FALSE)
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

find_switch_point_from_blast_coverage <- function(cvrg, mcmc_seed = 42, n_beast_iter = 1) {
  valid_cps <- c()

  for (iter in 1:n_beast_iter) {
    # Generate a unique seed for each iteration
    current_seed <- mcmc_seed + iter - 1

    invisible(capture.output(
    switch_info <- beast(cvrg, season = 'none',
                         tcp.minmax = c(1, 2),
                         torder.minmax = c(0, 0),
                         tseg.leftmargin = 500,
                         tseg.rightmargin = 300,
                         mcmc.seed = current_seed,
                         quiet = TRUE)
    ))
    cp <- switch_info$trend$pos_cp[1]
    cp_neg <- switch_info$trend$pos_cpNeg
    if (!is.null(cp_neg)) {
      if (any(cp_neg > cp)) {
        next
      }
    }
    if (is.null(cp) | is.na(cp)) {
      next
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
      valid_cps <- c(valid_cps, cp)
    }
  }

  if (length(valid_cps) == 0) {
    return(NA)
  } else {
    # Return the most frequent switch point
    return(as.numeric(names(sort(table(valid_cps), decreasing = TRUE)[1])))
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
  # verify then the extracted sequence will not be out of bounds
  if (cp_up - L / 2 < 1 | cp_down - L / 2 < 1) {
    return(NA)
  }
  if (cp_up + L / 2 > nchar(seq_up) | cp_down + L / 2 > nchar(seq_down)) {
    return(NA)
  }

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

  if (cp_left - 12 < 1 | cp_right - 12 < 1) {
    return(NA)
  }

  TSD_test1 <- subseq(seq_up, cp_left - 12, cp_left - 1)
  TSD_test2 <- subseq(seq_down, cp_right - 12, cp_right - 1) |>
    reverseComplement() |>
    as.character()
  if (nchar(TSD_test1) != 12 | nchar(TSD_test2) != 12) {
    return(NA)
  }

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
  # verify then the extracted sequence will not be out of bounds
  if (cp_up - L / 2 < 1 | cp_down - L / 2 < 1) {
    return(NA)
  }
  if (cp_up + L / 2 > nchar(seq_up) | cp_down + L / 2 > nchar(seq_down)) {
    return(NA)
  }
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
  if (cp_left -12 < 1 | cp_right - 12 < 1) {
    return(NA)
  }
  TSD_test1 <- subseq(seq_up, cp_left - 12, cp_left - 1)
  TSD_test2 <- subseq(seq_down, cp_right - 12, cp_right - 1) |>
    reverseComplement() |>
    as.character()
  if (nchar(TSD_test1) != 12 | nchar(TSD_test2) != 12) {
    return(NA)
  }
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

  # verify then the extracted sequence will not be out of bounds
  if (cp_up - L / 2 < 1 | cp_down - L / 2 < 1) {
    return(NA)
  }
  if (cp_up + L / 2 > nchar(seq_up) | cp_down + L / 2 > nchar(seq_down)) {
    return(NA)
  }

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
  tsd_length <- 2
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
    if (cp_left -12 < 1 | cp_right - 12 < 1) {
      return(NA)
    }
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

process_region_files <- function(file_list, side, mcmc_seed = 42, n_beast_iter = 1) {
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
    cp <- find_switch_point(aln_info, mcmc_seed = mcmc_seed, n_beast_iter = n_beast_iter)
    classification <- gsub(pattern, "", basename(f))
    mean_coverage <- mean(colSums(CM))
    # (dpos is calculated in your original code, but it is not used later)

    if (!is.na(cp)) {
      contig_id <- gsub("[.]fasta", "", gsub(".+_Contig", "", f))
      # Remove _part00X suffix if present (from split contig assemblies)
      contig_id <- gsub("_part\\d+$", "", contig_id)
      info <- extract_info_from_switch_points(aln_info, ctg, cp)
      if (is.null(info)) {
        next
      }
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
    cmd <- "mmseqs"


    args_list <- c(
      "easy-cluster",
      tir_seqs_parts_files[i],
      out_prefix,
      tmp_dir,
      "--cluster-mode",  "0",  "-v", "1", "--min-seq-id", "0.8",
      "--cov-mode", "0", "--mask", "0",  "-s", "6", "--threads", threads
    )
    res <- system2(cmd, args = args_list, stdout = TRUE, stderr = TRUE)
    # export mmeqs2 cluster output for debugging
    log_file <- paste0(out_prefix, "_log.txt")
    # write command output to log file
    cat(paste("Command:", cmd, paste(args_list, collapse = " "),"\n"), file = log_file)
    cat(res, file = log_file, append = TRUE)


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
  system(paste("blastn -task blastn -query", query_up, "-db", db_up,
               "-out", blast_up_out, blast_args_up),
         intern = TRUE)
  system(paste("blastn -task blastn -query", query_down, "-db", db_down,
               "-out", blast_down_out, blast_args_down),
         intern = TRUE)

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


round1 <- function(contig_dir, tir_flank_file, mcmc_seed = 42, n_beast_iter = 1) {
  message("\n ---- Identification of elements - Round 1 ----")
  # Read TIR flank coordinates
  tir_flank_coordinates <- read.table(tir_flank_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  if (nrow(tir_flank_coordinates) == 0) {
    message("TIR flank coordinates file is empty.")
  }

  # Identify upstream and downstream contig files
  upstream_contigs <- dir(contig_dir, pattern = "upstream_Contig", full.names = TRUE)
  if (length(upstream_contigs) == 0) {
    message("No upstream contig files found in the specified directory.")
  }
  downstream_contigs <- dir(contig_dir, pattern = "downstream_Contig", full.names = TRUE)
  if (length(downstream_contigs) == 0) {
    message("No downstream contig files found in the specified directory.")
  }

  message("Processing upstream regions")
  upstream_results <- process_region_files(upstream_contigs, "upstream", mcmc_seed = mcmc_seed,
                                           n_beast_iter = n_beast_iter)
  # Ensure process_region_files returned non-empty 'info'
  if (length(upstream_results$info) == 0) {
    message("No upstream information returned from process_region_files.")
    upstream_info <- data.frame()
  } else {
    upstream_info <- do.call(rbind, upstream_results$info)
  }
  ctg_list_upstream <- upstream_results$ctg_list

  message("Processing downstream regions")
  downstream_results <- process_region_files(downstream_contigs, "downstream", mcmc_seed = mcmc_seed,
                                             n_beast_iter = n_beast_iter)
  if (length(downstream_results$info) == 0) {
    message("No downstream information returned from process_region_files.")
    downstream_info <- data.frame()
  } else {
    downstream_info <- do.call(rbind, downstream_results$info)
  }
  ctg_list_downstream <- downstream_results$ctg_list

  # If either upstream_info or downstream_info is empty, return early with empty results.
  if (nrow(upstream_info) == 0 || nrow(downstream_info) == 0) {
    message("Warning: Upstream or downstream information is empty. Returning empty results.")
    return(list(
      gr1 = GRanges(),
      res_df = data.frame(),
      ctg_list_upstream = ctg_list_upstream,
      ctg_list_downstream = ctg_list_downstream,
      tir_flank_coordinates = tir_flank_coordinates
    ))
  }

  # Find common IDs between upstream and downstream results
  both_side_id <- intersect(downstream_info$ID, upstream_info$ID)
  if (length(both_side_id) == 0) {
    message("Warning: No common IDs found between upstream and downstream information. Returning empty results.")
    return(list(
      gr1 = GRanges(),
      res_df = data.frame(),
      ctg_list_upstream = ctg_list_upstream,
      ctg_list_downstream = ctg_list_downstream,
      tir_flank_coordinates = tir_flank_coordinates
    ))
  }

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
  if (length(res) == 0) {
    message("No valid elements found in the loop. Creating empty result data frame.")
    res_df <- data.frame()
  } else {
    res_df <- do.call(rbind, lapply(res, as.data.frame))
  }

  # Create GRanges for round 1 using the helper function prepare_granges()
  if (nrow(res_df) == 0) {
    message("Result data frame is empty, generating empty GRanges object.")
    gr1 <- GRanges()
  } else {
    gr1 <- prepare_granges(res_df, tir_flank_coordinates, iter = 1)
  }

  message("Number of elements found in the first round: ", length(gr1))
  message("------------------------------------------------------------------")

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
  message("\n---- Identification of elements - Round 2 ----")

  # 1. Early check: if the input result data frame is empty, return empty GRanges.
  if (nrow(res_df) == 0) {
    message("Input result data frame is empty. Returning empty GRanges for Round 2.")
    return(list(gr2 = GRanges(), gr_fin = GRanges()))
  }

  ## Retrieve passed contigs from Round 1
  ctg_upstream_pass <- ctg_list_upstream[unique(res_df$upstream_contig)]
  ctg_upstream_pass_class <- res_df$Classification[match(names(ctg_upstream_pass), res_df$upstream_contig)]
  ctg_downstream_pass <- ctg_list_downstream[unique(res_df$downstream_contig)]
  ctg_downstream_pass_class <- res_df$Classification[match(names(ctg_downstream_pass), res_df$downstream_contig)]

  # 2. Check if any contig sequences were retrieved.
  if (length(ctg_upstream_pass) == 0 || length(ctg_downstream_pass) == 0) {
    message("No upstream or downstream contig sequences passed from Round 1. Returning empty results.")
    return(list(gr2 = GRanges(), gr_fin = GRanges()))
  }

  ## Build consensus sequences for upstream and downstream contigs.
  cons_seqs <- list()
  for (i in seq_along(ctg_upstream_pass)) {
    label <- paste0(ctg_upstream_pass_class[i], "_", names(ctg_upstream_pass)[i], "_upstream")
    cons_seq <- get_consensus_from_aln(ctg_upstream_pass[[i]], 60)
    # 3. Check if the consensus sequence is empty and skip it if so.
    if (length(cons_seq) == 0) {
      message("Consensus sequence for upstream contig ", names(ctg_upstream_pass)[i], " is empty. Skipping.")
      next
    }
    cons_seqs[[label]] <- cons_seq
  }
  for (i in seq_along(ctg_downstream_pass)) {
    label <- paste0(ctg_downstream_pass_class[i], "_", names(ctg_downstream_pass)[i], "_downstream")
    cons_seq <- get_consensus_from_aln(ctg_downstream_pass[[i]], 60)
    if (length(cons_seq) == 0) {
      message("Consensus sequence for downstream contig ", names(ctg_downstream_pass)[i], " is empty. Skipping.")
      next
    }
    cons_seqs[[label]] <- cons_seq
  }
  # 4. If no consensus sequences were generated, exit early.
  if (length(cons_seqs) == 0) {
    message("No consensus sequences were generated. Returning empty results.")
    return(list(gr2 = GRanges(), gr_fin = GRanges()))
  }
  cons_seq_out <- unlist(DNAStringSetList(cons_seqs))
  # 5. Verify that the unlisted consensus sequences are non-empty.
  if (length(cons_seq_out) == 0) {
    message("No consensus sequences  generated, Returning empty results.")
    return(list(gr2 = GRanges(), gr_fin = GRanges()))
  }

  ## Extract metadata from consensus sequence names.
  side <- sapply(strsplit(names(cons_seqs), "_"), function(x) tail(x, 1))
  cons_class <- sapply(strsplit(names(cons_seqs), "_"), function(x) paste(head(x, -2), collapse = "_"))

  ## Create a data frame of output file names (one row per classification).
  unique_classes <- unique(cons_class)
  # 6. Check if any unique classifications were identified.
  if (length(unique_classes) == 0) {
    message("No unique classifications found. Returning empty results.")
    return(list(gr2 = GRanges(), gr_fin = GRanges()))
  }
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
    # 7. Ensure consensus sequences exist for the current classification.
    if (length(upstream_cons) == 0 || length(downstream_cons) == 0) {
      message("Consensus sequences for classification ", cls, " are empty. Skipping this classification.")
      next
    }
    writeXStringSet(upstream_cons, cons_file_df$upstream_cons[i])
    writeXStringSet(downstream_cons, cons_file_df$downstream_cons[i])
  }

  dir.create(file.path(output, "blastn"), showWarnings = FALSE)
  # Define BLAST output column names.
  col_names <- c("qaccver", "saccver", "pident", "length", "mismatch",
                 "gapopen", "qstart", "qend", "sstart", "send", "evalue",
                 "bitscore")
  blast_args <- paste0("-outfmt 6 -evalue 1e-5 -num_threads ", threads)
  res2_list <- list()
  gff2_list <- list()

  for (i in seq_along(unique_classes)) {
    cls <- unique_classes[i]
    message("DEBUG: Processing classification: ", cls)

    # 8. Check that the consensus file names for the current classification are defined.
    if (is.na(cons_file_df$upstream_cons[i]) || is.na(cons_file_df$downstream_cons[i])) {
      message("Consensus files for ", cls, " are not defined. Skipping BLAST for this class.")
      next
    }
    upstream_cons <- cons_file_df$upstream_cons[i]
    downstream_cons <- cons_file_df$downstream_cons[i]
    upstream_db <- cons_file_df$upstream_db[i]
    downstream_db <- cons_file_df$downstream_db[i]
    blast_upstream <- file.path(output, "blastn", paste0(cls, "_upstream.blastn"))
    blast_downstream <- file.path(output, "blastn", paste0(cls, "_downstream.blastn"))

    message("DEBUG: Running BLAST for ", cls)
    # 9. Call blast_helper() to create databases, run BLAST, and process outputs.
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

    message("DEBUG: BLAST completed. Upstream hits: ", nrow(blast_up_df), ", Downstream hits: ", nrow(blast_down_df))

    # 10. Validate blast results; if empty, skip this classification.
    if (nrow(blast_up_df) == 0 || nrow(blast_down_df) == 0) {
      message("BLAST returned empty results for ", cls, ". Skipping this class.")
      next
    }

    message("DEBUG: Reading sequence databases for ", cls)
    # 11. Read sequence databases and verify that they are non-empty.
    upstream_regions <- readDNAStringSet(upstream_db)
    downstream_regions <- readDNAStringSet(downstream_db)
    message("DEBUG: Upstream regions: ", length(upstream_regions), ", Downstream regions: ", length(downstream_regions))

    if (length(upstream_regions) == 0 || length(downstream_regions) == 0) {
      message("No sequences in upstream and donwstream db for ", cls, ". Skipping..")
      next
    }

    message("DEBUG: Getting sequence lengths for boundary checking")
    # Get sequence lengths to filter out-of-range BLAST hits
    upstream_lengths <- width(upstream_regions)
    names(upstream_lengths) <- names(upstream_regions)
    downstream_lengths <- width(downstream_regions)
    names(downstream_lengths) <- names(downstream_regions)

    message("DEBUG: Checking for valid sequence names in BLAST results")
    # Filter BLAST results to remove hits too close to sequence boundaries
    # Also filter out hits where sequence name doesn't exist in the database
    valid_up_seq <- blast_up_df$saccver %in% names(upstream_regions)
    valid_down_seq <- blast_down_df$saccver %in% names(downstream_regions)

    message("DEBUG: Valid upstream sequences: ", sum(valid_up_seq), "/", length(valid_up_seq))
    message("DEBUG: Valid downstream sequences: ", sum(valid_down_seq), "/", length(valid_down_seq))

    message("DEBUG: Checking coordinate boundaries")
    # For valid sequences, check if coordinates are within bounds
    # Need: sstart - 12 >= 1 (for TSD) and sstart + 100 <= seq_length (for TIR)
    # Only check coordinates for sequences that exist in the database
    valid_up_coords <- rep(FALSE, nrow(blast_up_df))
    valid_down_coords <- rep(FALSE, nrow(blast_down_df))

    # Check coordinates only for valid sequence names
    if (sum(valid_up_seq) > 0) {
      message("DEBUG: Checking upstream coordinate boundaries for valid sequences")
      valid_idx <- which(valid_up_seq)
      for (idx in valid_idx) {
        seq_name <- blast_up_df$saccver[idx]
        seq_len <- upstream_lengths[seq_name]
        sstart <- blast_up_df$sstart[idx]
        if (!is.na(seq_len) && !is.na(sstart)) {
          valid_up_coords[idx] <- (sstart >= 13) && ((sstart + 100) <= seq_len)
        }
      }
    }

    if (sum(valid_down_seq) > 0) {
      message("DEBUG: Checking downstream coordinate boundaries for valid sequences")
      valid_idx <- which(valid_down_seq)
      for (idx in valid_idx) {
        seq_name <- blast_down_df$saccver[idx]
        seq_len <- downstream_lengths[seq_name]
        sstart <- blast_down_df$sstart[idx]
        if (!is.na(seq_len) && !is.na(sstart)) {
          valid_down_coords[idx] <- (sstart >= 13) && ((sstart + 100) <= seq_len)
        }
      }
    }

    message("DEBUG: Valid upstream coordinates: ", sum(valid_up_coords), "/", nrow(blast_up_df))
    message("DEBUG: Valid downstream coordinates: ", sum(valid_down_coords), "/", nrow(blast_down_df))

    # Report filtering statistics
    if (sum(!valid_up_seq) > 0) {
      message("Filtered out ", sum(!valid_up_seq), " upstream BLAST hits with invalid sequence names")
    }
    if (sum(!valid_down_seq) > 0) {
      message("Filtered out ", sum(!valid_down_seq), " downstream BLAST hits with invalid sequence names")
    }
    if (sum(valid_up_seq & !valid_up_coords) > 0) {
      message("Filtered out ", sum(valid_up_seq & !valid_up_coords), " upstream BLAST hits too close to sequence boundaries")
    }
    if (sum(valid_down_seq & !valid_down_coords) > 0) {
      message("Filtered out ", sum(valid_down_seq & !valid_down_coords), " downstream BLAST hits too close to sequence boundaries")
    }

    # Apply filters
    blast_up_df <- blast_up_df[valid_up_coords, ]
    blast_down_df <- blast_down_df[valid_down_coords, ]

    # Check if any valid hits remain
    if (nrow(blast_up_df) == 0 || nrow(blast_down_df) == 0) {
      message("No valid BLAST hits remain after boundary filtering for ", cls, ". Skipping.")
      next
    }

    # Now create ranges with the filtered data (no boundary issues possible)
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
      tsd <- extract_TSD(blast_up_df$TSD[j],
                         as.character(reverseComplement(DNAString(blast_down_df$TSD[j]))))
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
    # 12. If no valid BLAST alignment was found, skip this classification.
    if (length(res2) == 0) {
      message("No valid BLAST alignment found for ", cls, ". Skipping this class.")
      next
    }
    res2_df <- do.call(rbind, lapply(res2, as.data.frame))
    res2_df <- merge(res2_df, tir_flank_coordinates, by = "ID")
    res2_df <- compute_element_boundaries(res2_df)

    df_gff <- data.frame(
      seqid = res2_df$SeqID,
      source = "DANTE_TIR",
      start = ifelse(res2_df$element_start < res2_df$element_end,
                     res2_df$element_start, res2_df$element_end),
      end = ifelse(res2_df$element_start < res2_df$element_end,
                   res2_df$element_end, res2_df$element_start) + 1,
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
  if (length(gff2_list) == 0) {
    gr2 <- GRanges()
  } else {
    gr2 <- unlist(GRangesList(gff2_list))
  }
  gr2_unique <- gr2[!gr2$ID %in% gr1$ID, ]
  gr_fin <- sort(concatGRanges(gr1, gr2_unique), by = ~ seqnames * start)
  message("Number of elements found in the second round: ", length(gr2_unique))
  message("-----------------------------------------------------------------")
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
  message("\n---- Identification of elements - Round 3 ----")

  # Retrieve fasta files for upstream and downstream regions.
  upstream_regions_file <- dir(contig_dir, pattern = "upstream_regions.fasta$", full.names = TRUE)
  downstream_regions_file <- dir(contig_dir, pattern = "downstream_regions.fasta$", full.names = TRUE)

  # Check for file existence before proceeding.
  if (length(upstream_regions_file) == 0 || length(downstream_regions_file) == 0) {
    message("Upstream or Downstream regions files not found in the directory. Returning empty GRanges for Round 3.")
    return(list(gr3 = GRanges(), gr3_unique = GRanges(), gr_fin = gr_fin))
  }

  names(upstream_regions_file) <- gsub("_upstream_regions.fasta", "", basename(upstream_regions_file))
  names(downstream_regions_file) <- gsub("_downstream_regions.fasta", "", basename(downstream_regions_file))

  tir_class_defined <- c("Class_II_Subclass_1_TIR_Tc1_Mariner",
                         "Class_II_Subclass_1_TIR_PIF_Harbinger",
                         "Class_II_Subclass_1_TIR_MuDR_Mutator",
                         "Class_II_Subclass_1_TIR_hAT",
                         "Class_II_Subclass_1_TIR_EnSpm_CACTA")

  res3_list <- list()

  # Iterate over classification names present in the upstream files.
  for (cls in names(upstream_regions_file)) {

    # Check that a corresponding downstream file is available.
    if (!(cls %in% names(downstream_regions_file))) {
      message("No downstream file found for ", cls, ". Skipping.")
      next
    }

    upstream_db   <- upstream_regions_file[[cls]]
    downstream_db <- downstream_regions_file[[cls]]
    blast_upstream   <- file.path(output, "blastn", paste0(cls, "_upstream3.blastn"))
    blast_downstream <- file.path(output, "blastn", paste0(cls, "_downstream3.blastn"))
    dir.create(file.path(output, "blastn"), showWarnings = FALSE)

    # Run BLAST analysis for upstream and downstream.
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

    # Check that cp_vals are not NULL.
    if (is.null(cp_up) || is.null(cp_down)) {
      message("CP detection returned NULL for ", cls, ". Skipping.")
      next
    }

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

    # Read sequences from the upstream and downstream databases.
    seq_upstream <- readDNAStringSet(upstream_db)
    seq_downstream <- readDNAStringSet(downstream_db)

    # Confirm that sequences were retrieved.
    if (length(seq_upstream) == 0 || length(seq_downstream) == 0) {
      message("No sequences found in either upstream or downstream for ", cls, ". Skipping.")
      next
    }

    # Reorder sequences based on the rownames of cp_up_down.
    seq_upstream <- seq_upstream[match(rownames(cp_up_down), names(seq_upstream))]
    seq_downstream <- seq_downstream[match(rownames(cp_up_down), names(seq_downstream))]

    # Verify that reordering did not produce any missing values.
    if (any(is.na(names(seq_upstream))) || any(is.na(names(seq_downstream)))) {
      message("Sequence reordering resulted in NA names for ", cls, ". Skipping.")
      next
    }

    # Generate the appropriate worker function based on the classification.
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
    } else {
      message("No detection worker defined for class ", cls, ". Skipping.")
      next
    }

    # Apply the worker function in parallel.
    parallel_res <- mclapply(seq_along(seq_upstream), worker_fun, mc.cores = threads)
    filtered_res <- Filter(Negate(is.null), parallel_res)
    if (length(filtered_res) > 0) {
      res3_list <- c(res3_list, filtered_res)
    } else {
      message("No valid results returned for ", cls)
    }
  } # end for each classification

  # Check that some valid results were obtained.
  if (length(res3_list) == 0) {
    message("No valid elements detected in Round 3. Returning empty GRanges.")
    return(list(gr3 = GRanges(), gr3_unique = GRanges(), gr_fin = gr_fin))
  }

  res3_df <- do.call(rbind, lapply(res3_list, as.data.frame))
  if (nrow(res3_df) == 0) {
    message("Resulting data frame from Round 3 is empty. Returning empty GRanges.")
    return(list(gr3 = GRanges(), gr3_unique = GRanges(), gr_fin = gr_fin))
  }

  gr3 <- prepare_granges(res3_df, tir_flank_coordinates, iter = 3)
  gr3_unique <- gr3[!gr3$ID %in% gr_fin$ID, ]
  # gr_fin <- c(gr_fin, gr3_unique)
  gr_fin <- sort(concatGRanges(gr_fin, gr3_unique), by = ~ seqnames * start)

  message("Number of elements found in the third round: ", length(gr3_unique))
  message("-------------------------------------------------------------------")
  return(list(gr3 = gr3, gr3_unique = gr3_unique, gr_fin = gr_fin))
}

round4 <- function(gr_fin, tir_cls_df, tir_seqs, tir_flank_coordinates, output, threads, multiplicity_threshold = 5) {
  message("\n---- Identification of elements - Round 4 ----")

  # 1. Check for valid input in tir_cls_df.
  if (nrow(tir_cls_df) == 0) {
    message("tir_cls_df is empty. No representative IDs available. Returning original gr_fin.")
    return(list(gr4 = GRanges(), gr_fin = gr_fin, res4_df = data.frame()))
  }

  repre_ID <- unique(tir_cls_df$Representative_ID[tir_cls_df$Multiplicity >= multiplicity_threshold])
  if (length(repre_ID) == 0) {
    message("No representative IDs pass the multiplicity threshold. Skipping round 4 and returning original gr_fin.")
    return(list(gr4 = GRanges(), gr_fin = gr_fin, res4_df = data.frame()))
  }

  # 2. Subset representative sequences.
  repre_seqs <- tir_seqs[names(tir_seqs) %in% repre_ID]
  if (length(repre_seqs) == 0) {
    message("No representative sequences found for the selected representative IDs. Skipping round 4.")
    return(list(gr4 = GRanges(), gr_fin = gr_fin, res4_df = data.frame()))
  }

  repre_cls <- paste0("Class_II_Subclass_1_", gsub(".+/", "", names(repre_seqs)))
  repre_seqs_groups <- split(repre_seqs, repre_cls)
  if (length(repre_seqs_groups) == 0) {
    message("No representative sequence groups formed. Skipping round 4.")
    return(list(gr4 = GRanges(), gr_fin = gr_fin, res4_df = data.frame()))
  }

  repre_seqs_files <- paste0(output, "/mmseqs2/", names(repre_seqs_groups), "_repre.fasta")
  names(repre_seqs_files) <- names(repre_seqs_groups)

  # 3. Ensure output subdirectory exists.
  mmseqs2_dir <- file.path(output, "mmseqs2")
  if (!dir.exists(mmseqs2_dir)) {
    dir.create(mmseqs2_dir, recursive = TRUE)
  }

  # Write FASTA files and create BLAST databases.
  for (cls in names(repre_seqs_groups)) {
    writeXStringSet(repre_seqs_groups[[cls]], repre_seqs_files[cls])
    system(paste("makeblastdb -in", repre_seqs_files[cls], "-dbtype nucl"), intern = TRUE)
  }

  blast_cols <- c("qaccver", "saccver", "pident", "length", "mismatch",
                  "gapopen", "qstart", "qend", "sstart", "send", "evalue",
                  "bitscore", "slen", "qlen")
  outfmt_str <- paste0(' -outfmt "6 ', paste(blast_cols, collapse = " "), '"')

  res4_list <- list()

  # 4. Loop through each representative sequence group (classification).
  for (cls in names(repre_seqs_groups)) {
    upstream_db <- paste0(output, "/", cls, "_upstream_regions.fasta")
    downstream_db <- paste0(output, "/", cls, "_downstream_regions.fasta")

    # 5. Verify that the upstream and downstream database files exist.
    if (!file.exists(upstream_db) || !file.exists(downstream_db)) {
      message("Database files for ", cls, " do not exist. Skipping.")
      next
    }

    blast_upstream <- paste0(output, "/blastn/", cls, "_upstream_round4.blastn")
    blast_downstream <- paste0(output, "/blastn/", cls, "_downstream_round4.blastn")

    blast_args_up <- paste0(outfmt_str, " -evalue 1e-5 -strand plus -num_threads ", threads)
    blast_args_down <- paste0(outfmt_str, " -evalue 1e-5 -strand minus -num_threads ", threads)

    # 6. Execute BLAST helper and ensure valid output.
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
                              filter_fun_down = function(x) filter_blast4(x, upstream = FALSE, max_offset = 120))

    if (is.null(blast_res) || is.null(blast_res$blast_up_df) || is.null(blast_res$blast_down_df)) {
      message("BLAST did not return valid results for ", cls, ". Skipping.")
      next
    }

    blast_up_df <- blast_res$blast_up_df
    blast_down_df <- blast_res$blast_down_df

    cp_up <- get_cp_from_blast4(blast_up_df, upstream = TRUE)
    cp_down <- get_cp_from_blast4(blast_down_df, upstream = FALSE)
    if (is.null(cp_up) || is.null(cp_down) || nrow(cp_up) == 0 || nrow(cp_down) == 0) {
      message("No boundaties found for ", cls, ". Skipping.")
      next
    }

    cp_up_down <- merge(cp_up, cp_down, by = "ID")
    already_annotated <- cp_up_down$ID %in% gsub(".+_", "", gr_fin$ID)
    cp_up_down <- cp_up_down[!already_annotated, ]
    rownames(cp_up_down) <- cp_up_down$ID
    if (nrow(cp_up_down) == 0) {
      message("No additional elements found for ", cls, " in round 4")
      next
    }

    # 7. Read upstream and downstream sequences, ensuring they are available.
    seq_upstream <- readDNAStringSet(upstream_db)
    seq_downstream <- readDNAStringSet(downstream_db)
    if (length(seq_upstream) == 0 || length(seq_downstream) == 0) {
      message("No sequences read from databases for ", cls, ". Skipping.")
      next
    }
    seq_upstream <- seq_upstream[match(rownames(cp_up_down), names(seq_upstream))]
    seq_downstream <- seq_downstream[match(rownames(cp_up_down), names(seq_downstream))]

    # 8. Define a detection worker based on classification.
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
    } else {
      Message("No valid detection worker defined for ", cls, ". Skipping.")
      next
    }

    # 9. Run detection in parallel and add non-null results.
    parallel_res <- mclapply(seq_along(seq_upstream), worker_fun, mc.cores = threads)
    filtered_res <- Filter(Negate(is.null), parallel_res)
    if (length(filtered_res) == 0) {
      message("No detection results for ", cls, " in round 4")
      next
    }
    res4_list <- c(res4_list, filtered_res)
  }

  # 10. Check that valid detection results were obtained.
  if (length(res4_list) == 0) {
    message("No valid elements found in round 4. Returning original gr_fin.")
    return(list(gr4 = GRanges(), gr_fin = gr_fin, res4_df = data.frame()))
  }

  res4_df <- do.call(rbind, lapply(res4_list, as.data.frame))
  if (nrow(res4_df) == 0) {
    message("Resulting data frame from round 4 is empty. Returning original gr_fin.")
    return(list(gr4 = GRanges(), gr_fin = gr_fin, res4_df = data.frame()))
  }

  # 11. Prepare GRanges for round 4 and update overall results.
  gr4 <- prepare_granges(res4_df, tir_flank_coordinates, iter = 4)
  # gr_fin <- c(gr_fin, gr4)
  gr_fin <- sort(concatGRanges(gr_fin, gr4), by = ~ seqnames * start)

  message("Number of elements found in the fourth round: ", length(gr4))
  message("-------------------------------------------------------------------\n")
  return(list(gr4 = gr4, gr_fin = gr_fin, res4_df = res4_df))
}




##### Functions fo TIR summary:



## FUNCTIONS:
# Helper: write the HTML header (with basic CSS)
create_html_header <- function(htmlfile, title = "TIR Summary Report") {
  cat(
    "<!DOCTYPE html>\n<html>\n<head>\n",
    "<meta charset='utf-8'>\n",
    "<title>", title, "</title>\n",
    "<style>\n",
    "  body { font-family: sans-serif; margin:0; padding:0; display: flex; }\n",
    "  #sidebar { width: 300px; background: #eee; padding: 1em; overflow-y: auto; height:100vh; }\n",
    "  #main    { flex: 1; padding: 1em; overflow-y: auto; height:100vh; }\n",
    "  #sidebar a { display:block; margin:0.5em 0; text-decoration:none; color:#0077cc; font-size:0.85em; }\n",
    "  img { max-width: 100%; margin:0.5em 0; }\n",
    "  table { border-collapse: collapse; width: auto; font-size:0.85em; }\n",
    "  th, td { border: 1px solid #ddd; padding: 4px; }\n",
    "  th { text-align: left; background: #f9f9f9; }
",
    "  td { text-align: left; }
",
    "</style>\n",
    "</head>\n<body>\n",
    file = htmlfile, sep = ""
  )
}

# Helper: close the HTML document
close_html <- function(htmlfile) {
  cat("</body>\n</html>\n", file = htmlfile, append = TRUE)
}

# Main: generate the report
# - file_info: named list of lists, each with png paths and CSV path
# - outdir: directory containing those files
# - htmlfile: output HTML file path
# - max_fig_height: CSS max-height for each figure, e.g. "400px"
generate_html_report <- function(
  file_info,
  outdir,
  htmlfile = file.path(outdir, "report.html"),
  max_fig_height = "400px",
  element_counts
) {
  # 1) header
  create_html_header(htmlfile, title = "DANTE_TIR Summary")

  # 2) sidebar
  cat(
    "<div id='sidebar'>\n",
    "<h2>Contents</h2>\n",
    file = htmlfile,
    append = TRUE
  )
  for (name in names(file_info)) {
    safe_id <- gsub("\\W+", "_", name)
    cat(
      sprintf("<a href='#%s'>%s</a>\n", safe_id, name),
      file = htmlfile,
      append = TRUE
    )
  }


  cat("</div>\n", file = htmlfile, append = TRUE)
  # 3) main content
  cat("<div id='main'>\n", file = htmlfile, append = TRUE)
  # table with element counts

  # add title
  cat("<h1>DANTE_TIR summary</h1>\n", file = htmlfile, append = TRUE)

  cat("<table>\n", file = htmlfile, append = TRUE)
  # 2) Write header row
  header_cells <- paste0("<th>", colnames(element_counts), "</th>", collapse = "")
  cat("  <tr>", header_cells, "</tr>\n", file = htmlfile, append = TRUE)
  # 3) Write each data row
  for(i in seq_len(nrow(element_counts))) {
    # Coerce each cell to character (to avoid factors, etc.)
    row_vals  <- as.character(element_counts[i, , drop = TRUE])
    row_cells <- paste0("<td>", row_vals, "</td>", collapse = "")
    cat("  <tr>", row_cells, "</tr>\n", file = htmlfile, append = TRUE)
  }

  # 4) Close table
  cat("</table>\n", file = htmlfile, append = TRUE)
  cat("<hr/>\n", file = htmlfile, append = TRUE)
  for (name in names(file_info)) {
    info <- file_info[[name]]
    safe_id <- gsub("\\W+", "_", name)

    # Section header
    cat(
      sprintf("<h2 id='%s'>%s</h2>\n", safe_id, name),
      file = htmlfile,
      append = TRUE
    )

    # Display each image with max height
    imgs <- c(
      "Detected TIR Lengths" = info$hist_length_name,
      "TIR Consensus" = info$logo_name,
      "Complete Element Lengths" = info$element_length_hist_name
    )
    for (i in seq_along(imgs)) {
      img <- imgs[i]
      title <- names(imgs)[i]
      # get relative path from outdir!
      rel <- rel_path(img, outdir)
      # make header
      cat(
        sprintf("<h3>%s</h3>\n", title),
        file = htmlfile,
        append = TRUE
      )
      cat(
        sprintf("<img src='%s' alt='%s' style='max-height:%s;'>\n", rel, rel, max_fig_height),
        file = htmlfile,
        append = TRUE
      )
      # add break
      cat("<br>\n", file = htmlfile, append = TRUE)
    }

    # Inline table for clusters of size >= 3
    df <- read.table(info$df_name, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    df_sub <- subset(df, Cluster_Size >= 3)
    if (nrow(df_sub) > 0) {
      cat(
        "<table>\n<tr>",
        paste0(sprintf("<th>%s</th>", names(df_sub)), collapse = ""),
        "</tr>\n",
        file = htmlfile,
        append = TRUE
      )
      for (i in seq_len(nrow(df_sub))) {
        cat(
          "<tr>",
          paste0(sprintf("<td>%s</td>", df_sub[i, ]), collapse = ""),
          "</tr>\n",
          file = htmlfile,
          append = TRUE
        )
      }
      cat("</table>\n<hr/>\n", file = htmlfile, append = TRUE)
    } else {
      cat(
        "<p><em>No clusters of size &ge; 3.</em></p>\n<hr/>\n",
        file = htmlfile,
        append = TRUE
      )
    }
  }
  cat("</div>\n", file = htmlfile, append = TRUE)

  # 4) footer
  close_html(htmlfile)
}


letterA <- function(x.pos, y.pos, ht, wt, id = NULL) {

  x <- c(0, 4, 6, 10, 8, 6.8, 3.2, 2, 0, 3.6, 5, 6.4, 3.6)
  y <- c(0, 10, 10, 0, 0, 3, 3, 0, 0, 4, 7.5, 4, 4)
  x <- 0.1 * x
  y <- 0.1 * y

  x <- x.pos + wt * x
  y <- y.pos + ht * y

  if (is.null(id)) {
    id <- c(rep(1, 9), rep(2, 4))
  }else {
    id <- c(rep(id, 9), rep(id + 1, 4))
  }

  fill <- c("green", "white")

  list(x = x, y = y, id = id, fill = fill)
}

## T
letterT <- function(x.pos, y.pos, ht, wt, id = NULL) {

  x <- c(0, 10, 10, 6, 6, 4, 4, 0)
  y <- c(10, 10, 9, 9, 0, 0, 9, 9)
  x <- 0.1 * x
  y <- 0.1 * y

  x <- x.pos + wt * x
  y <- y.pos + ht * y

  if (is.null(id)) {
    id <- rep(1, 8)
  }else {
    id <- rep(id, 8)
  }

  fill <- "red"

  list(x = x, y = y, id = id, fill = fill)
}

## C
letterC <- function(x.pos, y.pos, ht, wt, id = NULL) {
  angle1 <- seq(0.3 + pi / 2, pi, length = 100)
  angle2 <- seq(pi, 1.5 * pi, length = 100)
  x.l1 <- 0.5 + 0.5 * sin(angle1)
  y.l1 <- 0.5 + 0.5 * cos(angle1)
  x.l2 <- 0.5 + 0.5 * sin(angle2)
  y.l2 <- 0.5 + 0.5 * cos(angle2)

  x.l <- c(x.l1, x.l2)
  y.l <- c(y.l1, y.l2)

  x <- c(x.l, rev(x.l))
  y <- c(y.l, 1 - rev(y.l))

  x.i1 <- 0.5 + 0.35 * sin(angle1)
  y.i1 <- 0.5 + 0.35 * cos(angle1)
  x.i1 <- x.i1[y.i1 <= max(y.l1)]
  y.i1 <- y.i1[y.i1 <= max(y.l1)]
  y.i1[1] <- max(y.l1)

  x.i2 <- 0.5 + 0.35 * sin(angle2)
  y.i2 <- 0.5 + 0.35 * cos(angle2)

  x.i <- c(x.i1, x.i2)
  y.i <- c(y.i1, y.i2)

  x1 <- c(x.i, rev(x.i))
  y1 <- c(y.i, 1 - rev(y.i))

  x <- c(x, rev(x1))
  y <- c(y, rev(y1))

  x <- x.pos + wt * x
  y <- y.pos + ht * y

  if (is.null(id)) {
    id <- rep(1, length(x))
  }else {
    id <- rep(id, length(x))
  }

  fill <- "blue"

  list(x = x, y = y, id = id, fill = fill)
}

## G
letterG <- function(x.pos, y.pos, ht, wt, id = NULL) {
  angle1 <- seq(0.3 + pi / 2, pi, length = 100)
  angle2 <- seq(pi, 1.5 * pi, length = 100)
  x.l1 <- 0.5 + 0.5 * sin(angle1)
  y.l1 <- 0.5 + 0.5 * cos(angle1)
  x.l2 <- 0.5 + 0.5 * sin(angle2)
  y.l2 <- 0.5 + 0.5 * cos(angle2)

  x.l <- c(x.l1, x.l2)
  y.l <- c(y.l1, y.l2)

  x <- c(x.l, rev(x.l))
  y <- c(y.l, 1 - rev(y.l))

  x.i1 <- 0.5 + 0.35 * sin(angle1)
  y.i1 <- 0.5 + 0.35 * cos(angle1)
  x.i1 <- x.i1[y.i1 <= max(y.l1)]
  y.i1 <- y.i1[y.i1 <= max(y.l1)]
  y.i1[1] <- max(y.l1)

  x.i2 <- 0.5 + 0.35 * sin(angle2)
  y.i2 <- 0.5 + 0.35 * cos(angle2)

  x.i <- c(x.i1, x.i2)
  y.i <- c(y.i1, y.i2)

  x1 <- c(x.i, rev(x.i))
  y1 <- c(y.i, 1 - rev(y.i))

  x <- c(x, rev(x1))
  y <- c(y, rev(y1))

  r1 <- max(x.l1)

  h1 <- 0.4
  x.add <- c(r1, 0.5, 0.5, r1 - 0.2, r1 - 0.2, r1, r1)
  y.add <- c(h1, h1, h1 - 0.1, h1 - 0.1, 0, 0, h1)


  if (is.null(id)) {
    id <- c(rep(1, length(x)), rep(2, length(x.add)))
  }else {
    id <- c(rep(id, length(x)), rep(id + 1, length(x.add)))
  }

  x <- c(rev(x), x.add)
  y <- c(rev(y), y.add)

  x <- x.pos + wt * x
  y <- y.pos + ht * y


  fill <- c("orange", "orange")

  list(x = x, y = y, id = id, fill = fill)

}

Letter <- function(which, x.pos, y.pos, ht, wt) {

  if (which == "A") {
    letter <- letterA(x.pos, y.pos, ht, wt)
  }else if (which == "C") {
    letter <- letterC(x.pos, y.pos, ht, wt)
  }else if (which == "G") {
    letter <- letterG(x.pos, y.pos, ht, wt)
  }else if (which == "T") {
    letter <- letterT(x.pos, y.pos, ht, wt)
  }else {
    stop("which must be one of A,C,G,T")
  }

  letter
}


plot_multiline_logo <- function(cons.logo, read = NULL, W = 50, setpar = TRUE, gaps =
  NULL) {
  ## logo - base order  - A C G T
  if (ncol(cons.logo) == 5) {
    gaps_prob <- cons.logo[, 5]
  }else {
    gaps_prob <- NULL
  }
  tm <- 4
  pwm <- as.matrix(cons.logo[, 1:4])
  N <- nrow(pwm)
  Nori <- N
  if (N < W) {
    W <- N
  }
  s1 <- seq(1, N, by = W)
  s2 <- seq(W, N, by = W)
  if (length(s2) < length(s1)) {
    pwm <- rbind(pwm, matrix(0, nrow = W * length(s1) - N, ncol = 4, dimnames = list
     (NULL, c('A', 'C', 'G', 'T'))))
    if (!is.null(read)) {
      pwm_read <- rbind(read, matrix(0, nrow = W * length(s1) - N, ncol = 4, dimnames =
        list(NULL, c('A', 'C', 'G', 'T'))))
    }
    N <- nrow(pwm)
    s2 <- seq(W, N, by = W)
  }
  if (setpar) {
    par(mfrow = c(ceiling(N / W), 1), mar = c(1, 4, 1, 0))
  }
  for (i in seq_along(s1)) {
    if (!is.null(read)) {
      plot.logo(pwm_read[s1[i]:s2[i],], maxh = 2)
    }
    plot.logo(pwm[s1[i]:s2[i],], maxh = max(rowSums(cons.logo)))
    if (!is.null(gaps)) {
      ## transparent rectangles
      rect((gaps[, 'start'] - s1[i] + 1), 0, (gaps[, 'end'] - s1[i] + 2), max(pwm), col
        = "#00000005")

    }
    if (!is.null(gaps_prob)) {
      rect(seq_along(s1[i]:s2[i]),
           max(rowSums(cons.logo)),
           seq_along(s1[i]:s2[i]) + 1,
           max(rowSums(cons.logo)) - gaps_prob[s1[i]:s2[i]],
           col = "#00000030")


    }
    ticks <- intersect(intersect(pretty(pretty(s1[i]:s2[i]) + 1), s1[i]:s2[i]), 1:Nori)
    axis(1, at = ticks + 1.5 - s1[i], label = ticks, tick = FALSE)
    y <- pretty(c(0, max(pwm)), n = tm)
    axis(2, at = y, label = y, las = 2, cex.axis = .7)
  }
}

plot.logo <- function(pwm, maxh = NULL) {
  acgt <- c("A", "C", "G", "T")
  pwm <- pwm[, acgt]
  nbp <- dim(pwm)[1]
  if (is.null(maxh)) { maxh <- max(rowSums(pwm)) }

  plot(0, 0, xlim = c(0, nbp + 1), ylim = c(0, maxh), type = "n", axes = F, xlab = "",
       ylab = "")
  for (i in 1:nbp) {
    S <- order(pwm[i,])
    hgts <- pwm[i, S]
    nts <- acgt[S]
    ypos <- c(0, cumsum(hgts)[1:3])
    for (j in 1:4) {
      if (hgts[j] == 0) next
      L <- Letter(which = nts[j], x.pos = i, y.pos = ypos[j], ht = hgts[j], wt = 1)
      Id <- L$id == 1
      polygon(L$x[Id], L$y[Id], lty = 0, col = L$fill[1])
      if (sum(L$id == 2) > 0) {
        polygon(L$x[!Id], L$y[!Id], lty = 0, col = L$fill[2])
      }
    }
  }
}



make_seq_equal_width <- function(x) {
    nmax <- max(width(x))
    gap_size <- nmax - width(x)
    gaps <- sapply(gap_size, function(x) {
        if (x > 0) {
            return(paste(rep("-", x), collapse = ""))
        } else {
            return("")
        }
    })
    snames <- names(x)
    new_seq <- DNAStringSet(paste0(as.character(x), gaps))
    names(new_seq) <- snames
    return(new_seq)
}
cluster_mmseqs2 <- function(
  seqSet,
  output_dir   = tempfile("mmseqs2_out_"),
  threads      = 1,
  min_seq_id   = 0.8,
  cov_mode     = 0,
  cluster_mode = 0,
  sensitivity  = 6,
  mmseqs_bin   = "mmseqs",
  tmp_dir      = tempfile("mmseqs2_tmp_")
) {
  #— sanity checks ----------------------------
  if (!inherits(seqSet, "XStringSet"))
    stop("`seqSet` must be a DNAStringSet/RNAStringSet/AAStringSet (i.e. XStringSet).")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tmp_dir,    recursive = TRUE, showWarnings = FALSE)

  #— write input FASTA ------------------------
  in_fa <- file.path(output_dir, "input_seqs.fasta")
  Biostrings::writeXStringSet(seqSet, filepath = in_fa)

  #— run mmseqs2 easy‐cluster -----------------
  out_pref     <- file.path(output_dir, "cluster")
  cluster_tsv  <- paste0(out_pref, "_cluster.tsv")
  args <- c(
    "easy-cluster",
    in_fa,
    out_pref,
    tmp_dir,
    "--threads",      threads,
    "--min-seq-id",   min_seq_id,
    "--cov-mode",     cov_mode,
    "--cluster-mode", cluster_mode,
    "-s",             sensitivity
  )
  # capture stdout/stderr in case you want to debug
  res <- system2(mmseqs_bin, args = args, stdout = TRUE, stderr = TRUE)

  #— parse clustering result ------------------
  if (!file.exists(cluster_tsv))
    stop("MMseqs2 clustering output not found. Check `mmseqs_bin`, parameters or logs:\n", paste(res, collapse = "\n"))

  df <- read.table(
    cluster_tsv,
    header       = FALSE,
    sep          = "\t",
    col.names    = c("Representative_ID", "Member_ID"),
    stringsAsFactors = FALSE
  )
  return(df)
}
calc_info_content <- function(cnt,
                               occ_thresh = 0.5,     # hard‐filter threshold (0–1)
                               weight_by_occ = TRUE, # if FALSE, hard‐filter only
                               bias_correct = TRUE) { # apply Miller‐Madow
  nuc <- intersect(c("A","C","G","T"), rownames(cnt))
  m   <- cnt[nuc, , drop=FALSE]
  Nmax <- max(colSums(m))    # full alignment depth (assuming gaps excluded)
  S    <- nrow(m)            # 4

  ic <- apply(m, 2, function(x) {
    N <- sum(x)
    f <- N / Nmax

    # if below occupancy threshold
    if (N == 0 || (!weight_by_occ && f < occ_thresh))
      return(NA_real_)

    p <- x / N
    p <- p[p > 0]
    H <- -sum(p * log2(p))    # raw entropy in bits

    # Miller–Madow bias correction?
    if (bias_correct) {
      K <- length(p)
      H <- H + (K - 1) / (2 * N * log(2))
    }

    raw_ic <- log2(S) - H

    # if soft weighting, scale by occupancy fraction
    if (weight_by_occ) raw_ic <- raw_ic * f

    # if hard‐filter but below threshold, mask
    if (!weight_by_occ && f < occ_thresh) raw_ic <- NA_real_

    raw_ic
  })

  return(ic)
}
rel_path <- function(file, outdir) {
  # normalize without requiring existence
  norm <- function(x) normalizePath(x, winslash = "/", mustWork = FALSE)
  file_abs   <- norm(file)
  outdir_abs <- norm(outdir)
  # ensure outdir ends with a slash
  prefix <- paste0(outdir_abs, "/")
  if (!startsWith(file_abs, prefix)) {
    stop(sprintf("File '%s' is not inside directory '%s'", file, outdir))
  }
  # strip off the common prefix
  substring(file_abs, nchar(prefix) + 1)
}




### END OF FUNCTIONS ###