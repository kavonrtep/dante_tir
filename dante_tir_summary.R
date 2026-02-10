#!/usr/bin/env Rscript

library(optparse)
# parse command line arguments - input - gff3 , genome file, output prefix
option_list <- list(
  make_option(c("-g", "--gff"), type = "character", default = NULL, help = "GFF3 file"),
  make_option(c("-o", "--output"), type = "character", default = NULL, help =
    "Output dir"),
  make_option(c("-f", "--fasta"), type = "character", default = NULL, help = "Fasta file with sequences used in GFF3"),
  make_option(c("-t", "--threads"), type = "integer", default = 1, help = "Number of threads to use"),
  make_option(c("-d", "--debug"), action = "store_true", default = FALSE, help =
    "save intermediate data for debugging"),
  make_option(c("-m", "--min_cluster_size"), type = "integer", default = 3, help =
    "Minimum cluster size to report representative elements")
)

description <-
  "This script summarize DANTE_TIR output.
    - Generate graphical summary
    - Report representative elements for each classification in FASTA format
    - Report cluster sizes
  "

initial_options <- commandArgs(trailingOnly = FALSE)
opt_parser <- OptionParser(option_list = option_list,
                           description = description)
opt <- parse_args(opt_parser )
if (is.null(opt$gff) ||
  is.null(opt$output) ||
  is.null(opt$fasta)) {
  stop("Please provide gff3 file, output prefix and fasta file")
}

## — check that input files are present
if (!file.exists(opt$gff))  stop(sprintf("GFF3 file not found: '%s'", opt$gff))
if (!file.exists(opt$fasta)) stop(sprintf("Fasta file not found: '%s'", opt$fasta))

# loading take which so it is after options are parsed wso help works without loading these packages
suppressPackageStartupMessages({
  library(rtracklayer)
  library(Biostrings)
  library(BSgenome)
})


file_arg_name <- "--file="
script_name <- normalizePath(sub(file_arg_name, "",
                                  initial_options[grep(file_arg_name, initial_options)]))
script_dir <- dirname(script_name)
utils_file <- file.path(script_dir, "dt_utils.R")
if (!file.exists(utils_file)) {
  stop(sprintf("Required helper script not found: '%s'", utils_file))
}
source(utils_file)


mmseqs_bin <- 'mmseqs' # it should be in PATH!
## — check that mmseqs binary is available
if (Sys.which(mmseqs_bin) == "") {
  stop(sprintf("mmseqs binary not found: please install it or add to PATH ('%s')", mmseqs_bin))
}

# for testing
if (FALSE) {
  opt <- list()
  opt$gff <- "/mnt/raid/454_data/workshop/2025/example_analysis/DANTE_TIR
  /DANTE_TIR_final.gff"
  opt$output <- "./tmp/test"
  opt$fasta <- "/mnt/raid/454_data/workshop/2025/example_analysis/input_data
  /Pisum_assembly_cameor_ver_2.fasta"
  mmseqs_bin <- "/home/petr/data/miniforge3/envs/dante_tir/bin/mmseqs"
  source("dt_utils.R")
}

dir.create(opt$output, showWarnings = FALSE)
dir.create(paste0(opt$output, "/mmseqs2"), showWarnings = FALSE)
dir.create(paste0(opt$output, "/img"), showWarnings = FALSE)


g <- import(opt$gff, format = "gff3")
# Verify that the GFF3 file has the required columns
required_cols <- c("ID", "Classification", "tir_seq3", "tir_seq5")
missing_cols  <- setdiff(required_cols, colnames(mcols(g)))
if (length(missing_cols) > 0) {
  stop(sprintf(
    "GFF3 is missing required attribute columns: %s",
    paste(missing_cols, collapse = ", ")
  ))
}



s <- readDNAStringSet(opt$fasta)
# remove strings after the first space
names(s) <- gsub(" .*", "", names(s))
#
## — validate that FASTA contains all seqnames used in the GFF
seqlevels_gff <- unique(as.character(seqnames(g)))
missing_seqs  <- setdiff(seqlevels_gff, names(s))
if (length(missing_seqs) > 0) {
  stop(sprintf(
    "Sequences referenced in GFF but not found in FASTA: %s",
    paste(missing_seqs, collapse = ", ")
  ))
}

# split by classification
g_split <- split(g, g$Classification)

element_counts <- sapply(g_split, length)
element_counts_df <- data.frame("Classification" = names(element_counts),
                             "Number of elements"= element_counts,
                             check.names = FALSE)



if (opt$debug) {
  save.image(file = paste0(opt$output, "/debug.RData"))
}

file_info <- list()
repr_seq_min <- DNAStringSetList()
for (i in names(g_split)) {
  message("Processing ", i)
  file_info[[i]] <- list()
  tryCatch({
  # get the name of the TIR
    tir3 <- DNAStringSet(g_split[[i]]$tir_seq3)
    names(tir3) <- paste0(g_split[[i]]$ID)

    tir5 <- DNAStringSet(g_split[[i]]$tir_seq5)
    names(tir5) <- paste0(g_split[[i]]$ID)
    # convert to DNAStringSet of the same length

    # make hist/table of sequence lengths
    hist_length_name <- paste0(opt$output, "/img/", i, "_tir_lengths.png")
    file_info[[i]]$hist_length_name <- hist_length_name

    png(hist_length_name, width = 1500, height = 800, pointsize = 30)
    par(mfrow = c(1, 1))
    l3 <- tabulate(width(tir3))
    l5 <- tabulate(width(tir5))
    lims <- c(0, max(width(tir3), width(tir5)))
    plot(l3, main = "", xlab = "Length [bp]", ylab = "Count", xlim = lims, type = 'h',
         lwd = 3)
    points(seq_along(l5) + 0.15, l5, main = paste0("TIR 5'"), xlab = "Length", ylab =
      "Count", xlim = lims, type = 'h', lwd = 3, col = "red")
    legend("topright", legend = c("3' TIR", "5' TIR"), col = c("black", "red"),
           lty = 1, lwd = 3)
    box()
    dev.off()

    logo_name <- paste0(opt$output, "/img/", i, "_logo.png")
    file_info[[i]]$logo_name <- logo_name

    t3 <- consensusMatrix(tir3, as.prob = FALSE, baseOnly = TRUE)
    t5 <- consensusMatrix(tir5, as.prob = FALSE, baseOnly = TRUE)
    # Limit logo to first 200 positions to avoid cairo size errors
    max_logo_pos <- 100
    t3_plot <- t3[, seq_len(min(ncol(t3), max_logo_pos)), drop = FALSE]
    t5_plot <- t5[, seq_len(min(ncol(t5), max_logo_pos)), drop = FALSE]
    n_pos <- ncol(t5_plot)
    logo_truncated <- ncol(t5) > max_logo_pos || ncol(t3) > max_logo_pos
    png(logo_name, height = 800, width = 100 + n_pos * 50)
    par(mfrow = c(2, 1))
    plot.logo(t(t3_plot))
    axis(side = 1, at = 1:n_pos + 0.5, labels = 1:n_pos, tick = FALSE)
    label3 <- if (logo_truncated) paste0("3' TIR (first ", ncol(t3_plot), " of ", ncol(t3), " bp)") else "3' TIR"
    mtext(label3, cex = 2, adj = 1)
    plot.logo(t(t5_plot))
    axis(side = 1, at = 1:n_pos + 0.5, labels = 1:n_pos, tick = FALSE)
    label5 <- if (logo_truncated) paste0("5' TIR (first ", ncol(t5_plot), " of ", ncol(t5), " bp)") else "5' TIR"
    mtext(label5, cex = 2, adj = 1)
    dev.off()

    t3aln <- DNAMultipleAlignment(make_seq_equal_width(tir3))
    t5aln <- DNAMultipleAlignment(make_seq_equal_width(tir5))

    max_elements_len <- max(width(g_split[[i]]))
    element_length_hist_name <- paste0(opt$output, "/img/", i, "_element_length_hist.png")
    file_info[[i]]$element_length_hist_name <- element_length_hist_name

    png(element_length_hist_name, width = 1500, height = 800, pointsize = 30)
    hist(width(g_split[[i]]), breaks = seq(0, max_elements_len + 50, by = 50),
         main = "", xlab = "Length [bp]", ylab = "Count", col = "#365070",
         border = "#365070", xaxs = "i",yaxs = "i")
    box()
    dev.off()
    # clustering with mmseqs2

    tir_complete <- getSeq(s, g_split[[i]])
    names(tir_complete) <- paste0(g_split[[i]]$ID)

    cls <- cluster_mmseqs2(
      tir_complete,
      output_dir = paste(opt$output,"/mmseqs2/mmseqs2_", i, sep = ""),
      threads = opt$threads,
      min_seq_id = 0.8,
      cov_mode = 0,
      cluster_mode = 0,
      sensitivity = 6,
      mmseqs_bin = mmseqs_bin,
      tmp_dir = tempfile("mmseqs2_tmp_")
    )

    cls$tir5 <- as.character(getSeq(tir5, cls$Member_ID))
    cls$tir3 <- as.character(getSeq(tir3, cls$Member_ID))
    cls$element_length <- width(g_split[[i]])[match(cls$Member_ID, g_split[[i]]$ID)]

    cluster_sizes <- sort(table(cls$Representative_ID), decreasing = TRUE)
      # make it data.frame
    cluster_sizes_df <- data.frame(Index = seq_along(cluster_sizes),
                                   Representative_Element = names(cluster_sizes),
                                   Cluster_Size = as.numeric(cluster_sizes))

    mean_length <- tapply(cls$element_length, cls$Representative_ID, mean)
    upper_quartile <- tapply(cls$element_length, cls$Representative_ID, function(x) {
      quantile(x, 0.75)
    })
    lower_quartile <- tapply(cls$element_length, cls$Representative_ID, function(x) {
      quantile(x, 0.25)
    })

    upper_10 <- tapply(cls$element_length, cls$Representative_ID, function(x) {
      quantile(x, 0.90)
    })
    lower_10 <- tapply(cls$element_length, cls$Representative_ID, function(x) {
      quantile(x, 0.10)
    })

    tir3_consensus <- tapply(cls$tir3, cls$Representative_ID, function(x) {
      x <- DNAStringSet(x)
      if (length(x) ==1){
        return(as.character(x))
      }
      x <- make_seq_equal_width(x)
      cons <- as.character(get_consensus_from_aln(x, clean = TRUE))
      cons
    })

    tir5_consensus <- tapply(cls$tir5, cls$Representative_ID, function(x) {
      x <- DNAStringSet(x)
      if (length(x) ==1){
        return(as.character(x))
      }
      x <- make_seq_equal_width(x)
      cons <- as.character(get_consensus_from_aln(x, clean = TRUE))
      cons
    })

    cluster_sizes_df$Cluster_Mean_Length <- round(mean_length[match(cluster_sizes_df$Representative_Element,
                                                      names(mean_length))], 0)

    q10 <- lower_10[match(cluster_sizes_df$Representative_Element,
                                                            names(lower_10))]
    q90 <- upper_10[match(cluster_sizes_df$Representative_Element,
                                                              names(upper_10))]
    cluster_sizes_df$Cluster_Length_IDR <- paste0(round(q10), " - ", round(q90))


    cluster_sizes_df$Representative_Element_Length <- width(tir_complete)[match(cluster_sizes_df$Representative_Element,
                                                                        names(tir_complete))]

    cluster_sizes_df$TIR3 <- cls$tir3[match(cluster_sizes_df$Representative_Element,
                                                            cls$Member_ID)]
    cluster_sizes_df$TIR5 <- cls$tir5[match(cluster_sizes_df$Representative_Element,
                                                            cls$Member_ID)]

    # export table with cluster sizes
    df_name <- paste0(opt$output, "/", i, "_representative_elements.csv")
    file_info[[i]]$df_name <- df_name

    write.table(cluster_sizes_df, file = df_name, sep = "\t", quote = FALSE,
                row.names = FALSE)

    # export representative sequences with cluster sizes at least 3
    repr_seq <- getSeq(tir_complete, cluster_sizes_df$Representative_Element)
    # add cluster size to the names

    multiplicity <- cluster_sizes_df$Cluster_Size[match(names(repr_seq),
                                                          cluster_sizes_df$Representative_Element)]


    # make classification string for repeatmasker:
    cls_string <- paste0("Class_II/Subclass_1/TIR/", gsub("^.+_TIR_","", names(repr_seq)))
    # remove numerical suffixes _1, __2, etc.
    cls_string <- gsub("_\\d+$", "", cls_string)
    names(repr_seq) <- paste0(names(repr_seq), "#", cls_string)


    # export only those with cluster size at least 3
    repr_seq_min[[i]] <- DNAStringSet(repr_seq[multiplicity >= opt$min_cluster_size])



    writeXStringSet(repr_seq[multiplicity >= opt$min_cluster_size],
                  file = paste0(opt$output, "/", i, "_representative_elements_min",
                                opt$min_cluster_size, ".fasta"))
    # export all representative sequences
    writeXStringSet(repr_seq,
                  file = paste0(opt$output, "/", i, "_representative_elements_all.fasta"))


  }, error = function(e) {
    message("Error processing ", i, ": ", e$message)
  })
}

# combine repr_seq_min to one FASTA file
writeXStringSet(unlist(repr_seq_min),
                file = paste0(opt$output, "/all_representative_elements_min",
                              opt$min_cluster_size, ".fasta"))

generate_html_report(file_info, opt$output, element_counts = element_counts_df)

