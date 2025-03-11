#!/usr/bin/env Rscript
library(optparse)
initial_options <- commandArgs(trailingOnly = FALSE)
options_list <- list(
  make_option(c("-c", "--contig_dir"), action = "store", type = "character",
              help = "directory with contigs"),
  make_option(c("-o", "--output"), action = "store", type = "character",
              help = "output prefix", default = NA),
  make_option(c("-t", "--threads"), action = "store", type = "integer",
              help = "number of threads to use", default = 1),
  make_option(c("-g", "--genome"), action = "store", type = "character",
              help = "genome fasta file", default = NA)
)

parser <- OptionParser(option_list = options_list)
opt <- parse_args(parser, args = commandArgs(TRUE))

suppressPackageStartupMessages({
  library(Biostrings)
  library(BSgenome)
  library(Rbeast)
  library(rtracklayer)
  library(parallel)
})
file_arg_name <- "--file="
script_name <- normalizePath(sub(file_arg_name, "",
                                 initial_options[grep(file_arg_name, initial_options)]))
script_dir <- dirname(script_name)
# get functions
source(file.path(script_dir, "dt_utils.R"))


# Define the path to the TIR flank coordinates file
tir_flank_file <- file.path(opt$contig_dir, "tir_flank_coords.txt")

tryCatch({

  ########################################################################################
  ####                           ROUND 1: TIR EXTRACTION                              ####
  ########################################################################################

  # Call the round1 function
  round1_results <- round1(opt$contig_dir, tir_flank_file)

  # Extract outputs for use in later rounds
  gr1 <- round1_results$gr1
  res_df <- round1_results$res_df
  ctg_list_upstream <- round1_results$ctg_list_upstream
  ctg_list_downstream <- round1_results$ctg_list_downstream
  tir_flank_coordinates <- round1_results$tir_flank_coordinates


  ########################################################################################
  ####                           ROUND 2: CONSENSUS & BLAST                           ####
  ########################################################################################
  round2_results <- round2(res_df, ctg_list_upstream, ctg_list_downstream, gr1,
                           tir_flank_coordinates, opt$output, opt$threads)
  gr2 <- round2_results$gr2
  gr_fin <- round2_results$gr_fin

  ########################################################################################
  ####                            THIRD ROUND OF DETECTION                            ####
  ########################################################################################
  round3_results <- round3(opt$contig_dir, opt$output, tir_flank_coordinates, gr_fin,
                           opt$threads)
  gr3 <- round3_results$gr3
  gr3_unique <- round3_results$gr3_unique
  gr_fin <- round3_results$gr_fin

  ########################################################################################
  ####                       CLUSTERING STEP (mmseqs2)                                ####
  ########################################################################################
  # This block clusters the TIR sequences (from gr_fin) and creates tir_cls_df.
  # It is required by round 4.

  message("Clustering TIR sequences with mmseqs2")
  save.image(paste0(opt$output, "/DANTE_TIR.RData"))
  clustering_results <- cluster_tir_sequences(opt$genome, gr_fin, opt$output, opt$threads)
  gr_fin <- clustering_results$gr_fin
  tir_cls_df <- clustering_results$tir_cls_df
  tir_seqs <- clustering_results$tir_seqs  # This is now available for use in Round 4.


  ########################################################################################
  # FOURTH ROUND OF DETECTION OF ADDITIONAL TIR ELEMENTS
  ########################################################################################

  round4_results <- round4(gr_fin, tir_cls_df, tir_seqs, tir_flank_coordinates,
                           opt$output, opt$threads)
  gr4 <- round4_results$gr4
  gr_fin <- round4_results$gr_fin
  res4_df <- round4_results$res4_df


  class_table <- as.data.frame.array(table(gr_fin$Classification))
  colnames(class_table) <- "Number of Elements"
  message("\nTIR classification summary:")
  print(class_table)
  class_table$Classification <- rownames(class_table)
  write.table(class_table[, 2:1], file = paste0(opt$output,
                                                "/TIR_classification_summary.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  export(gr_fin, paste0(opt$output, "/DANTE_TIR_final.gff"), format = "gff3")
  save.image(paste0(opt$output, "/DANTE_TIR.RData"))


  genome <- readDNAStringSet(opt$genome)
  names(genome) <- sub(" .*", "", names(genome))
  tir_seqs <- getSeq(genome, gr_fin)
  names(tir_seqs) <- paste0(gsub("Class_II_Subclass_1_TIR_", "", gr_fin$ID),
                            "#",
                            gsub("Class_II_Subclass_1_", "Class_II/Subclass_1/",
                                 gr_fin$Classification))
  writeXStringSet(tir_seqs, paste0(opt$output, "/DANTE_TIR_final.fasta"))
  tir_seqs_part <- split(tir_seqs, gsub("Class_II_Subclass_1_TIR_", "",
                                        gr_fin$Classification))
  tir_seqs_parts_files <- paste0(opt$output, "/DANTE_TIR_", names(tir_seqs_part), "
  .fasta")
  names(tir_seqs_parts_files) <- names(tir_seqs_part)
  for (i in names(tir_seqs_part)) {
    writeXStringSet(tir_seqs_part[[i]], tir_seqs_parts_files[i])
  }

}, error = function(e) {
  message("An error occurred during the pipeline")
  message("Error message:")
  message(e)
  message("Error traceback:")
  message(traceback())
  save.image(paste0(opt$output, "/DANTE_TIR.RData"))
  message("Intermediate results saved to DANTE_TIR.RData")
  q(status = 1)
})
