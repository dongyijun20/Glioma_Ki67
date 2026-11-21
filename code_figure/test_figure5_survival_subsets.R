#!/usr/bin/env Rscript
# Figure 5C must subset IDH-mutant TCGA cases. filter(c) passes the base
# concatenate function into dplyr::filter and aborts the script before the
# IDH-stratified Kaplan-Meier curves and the Cox model.
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
fig5 <- file.path(dirname(dirname(normalizePath(file_arg))), "code_figure", "Figure5.R")
if (!file.exists(fig5)) {
  fig5 <- file.path(dirname(normalizePath(file_arg)), "Figure5.R")
}
src <- paste(readLines(fig5, warn = FALSE), collapse = "\n")

if (grepl("filter\\(\\s*c\\s*\\)", src)) {
  stop("Figure5.R still calls filter(c), which aborts at the IDH-mutant panel")
}
if (!grepl('filter\\(IDH\\.status == "Mutant"\\)', src)) {
  stop('Figure5.R is missing filter(IDH.status == "Mutant")')
}
if (!grepl('by\\.x = "patient", by\\.y = "sample_id"', src)) {
  stop("Figure5 lineage merge still keys on sample_id, which meta_combined does not have")
}

cat("OK figure5 survival subsets\n")
