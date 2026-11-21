#!/usr/bin/env Rscript
# Locks donor-level metadata to one row per Patient after sample consolidation.
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
src <- file.path(dirname(normalizePath(file_arg)), "collapse_patient_meta.R")
source(src)

meta <- data.frame(
  Patient = c("LGG-04", "LGG-04", "LGG-04", "SF11488"),
  Age = c(NA, "40", "40", "55"),
  Sex = c("F", "F", "F", "M"),
  Grade = c("LGG (II-III)", "LGG (II-III)", "LGG (II-III)", "GBM (IV)"),
  Stage = c("Primary", "Primary", "Primary", "Recurrent"),
  Location = c("Frontal", "Frontal", "Temporal", "Temporal"),
  IDH = c("Mutant", "Mutant", "Mutant", "Wildtype"),
  MGMT.status = c("methylated", "methylated", "methylated", "unmethylated"),
  Diagnosis = c("Astrocytoma", "Astrocytoma", "Astrocytoma", "Glioblastoma"),
  batch = c("b1", "b2", "b3", "b1"),
  stringsAsFactors = FALSE
)

# The pre-fix distinct() key keeps one row per batch, so LGG-04 is repeated
# and the Figure 3A factor() call aborts.
distinct_like <- meta[!duplicated(meta), , drop = FALSE]
stopifnot(sum(distinct_like$Patient == "LGG-04") == 3L)
factor_err <- tryCatch(
  factor(distinct_like$Patient, levels = distinct_like$Patient),
  error = function(e) e
)
stopifnot(inherits(factor_err, "error"))
stopifnot(grepl("duplicated", conditionMessage(factor_err)))

out <- collapse_patient_meta(meta)
stopifnot(nrow(out) == 2L)
stopifnot(sum(out$Patient == "LGG-04") == 1L)
stopifnot(out$Age[out$Patient == "LGG-04"] == "40")
ids <- factor(out$Patient, levels = out$Patient)
stopifnot(nlevels(ids) == 2L)

cat("OK collapse_patient_meta\n")
