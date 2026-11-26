#!/usr/bin/env Rscript
# Locks the Figure 4 crashes: duplicate donor rows, and the ecotype palette
# being used before it exists.
args_cmd <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", args_cmd[grep("^--file=", args_cmd)])
root <- if (length(this_file)) dirname(dirname(normalizePath(this_file))) else getwd()
setwd(root)

source(file.path(root, "code_figure", "collapse_patient_meta.R"))

meta <- data.frame(
  Patient = c("LGG-04", "LGG-04", "GBM057"),
  Age = c(45, 45, 60),
  Sex = c("F", "F", "M"),
  Grade = c("LGG (II-III)", "LGG (II-III)", "GBM (IV)"),
  Stage = c("Primary", "Primary", "Primary"),
  Location = c("Frontal", "Temporal", "Frontal"),
  IDH = c("Mutant", "Mutant", "Wildtype"),
  MGMT.status = c("methylated", "methylated", "unmethylated"),
  Diagnosis = c("Astrocytoma", "Astrocytoma", "Glioblastoma"),
  batch = c("b1", "b2", "b1"),
  stringsAsFactors = FALSE
)

# distinct() keeps both LGG-04 rows because location and batch differ.
# factor() refuses duplicated levels and aborts Figure4 before any panel.
distinct_rows <- meta[!duplicated(meta), , drop = FALSE]
stopifnot(sum(distinct_rows$Patient == "LGG-04") == 2L)
dup_err <- tryCatch(
  factor(distinct_rows$Patient, levels = distinct_rows$Patient),
  error = function(e) e
)
stopifnot(inherits(dup_err, "error"))
stopifnot(grepl("duplicat", dup_err$message, ignore.case = TRUE))

one <- collapse_patient_meta(meta)
stopifnot(nrow(one) == 2L)
stopifnot(!anyDuplicated(one$Patient))
stopifnot(one$batch[one$Patient == "LGG-04"] == "b1")
stopifnot(one$Location[one$Patient == "LGG-04"] == "Frontal")
ok <- factor(one$Patient, levels = one$Patient)
stopifnot(length(ok) == 2L)
stopifnot(!anyNA(ok))

# A blank clinical field must not block a later real value on the same donor.
blank <- data.frame(
  Patient = c("P1", "P1"),
  Location = c("", "Parietal"),
  stringsAsFactors = FALSE
)
picked <- collapse_patient_meta(blank)
stopifnot(picked$Location == "Parietal")

fig4 <- readLines(file.path(root, "code_figure", "Figure4.R"))
code_lines <- fig4
code_lines[grepl("^\\s*#", code_lines)] <- ""
line_of <- function(pattern) {
  hits <- grep(pattern, code_lines)
  if (!length(hits)) stop("missing pattern: ", pattern)
  hits[[1]]
}
stopifnot(line_of("collapse_patient_meta\\(meta\\)") < line_of("factor\\(Patient, levels = Patient\\)"))
stopifnot(line_of("^cluster_colors <-") < line_of("scale_color_manual\\(values = cluster_colors\\)"))
stopifnot(line_of("library\\(circlize\\)") < line_of("colorRamp2\\("))
stopifnot(line_of("library\\(ComplexHeatmap\\)") < line_of("ht <- Heatmap\\("))
stopifnot(line_of("library\\(grid\\)") < line_of('unit\\(0\\.5, "cm"\\)'))
stopifnot(!any(grepl("distinct\\(\\)", fig4)))

cat("FIGURE4_ECOTYPE_SETUP_OK\n")
