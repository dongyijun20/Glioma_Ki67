# Collapse cell-level clinical metadata to one row per Patient.
#
# Donor consolidation (LGG-04-1/2/3 -> LGG-04) leaves sample-level fields such
# as batch on multiple rows that now share a Patient id. distinct() keeps
# those rows, and factor(Patient, levels = Patient) then fails because levels
# are duplicated. Later Wilcoxon / Pearson calls would also count that donor
# once per leftover row.
collapse_patient_meta <- function(meta) {
  meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  if (!"Patient" %in% names(meta)) {
    stop("collapse_patient_meta requires a Patient column")
  }
  meta$Patient <- as.character(meta$Patient)
  patients <- unique(meta$Patient[!is.na(meta$Patient)])
  picked <- lapply(patients, function(p) {
    sub <- meta[meta$Patient == p, , drop = FALSE]
    keep <- sub[1, , drop = FALSE]
    for (col in names(sub)) {
      vals <- sub[[col]]
      usable <- !is.na(vals) &
        trimws(as.character(vals)) != "" &
        as.character(vals) != "NA"
      if (any(usable)) {
        keep[[col]] <- vals[which(usable)[1]]
      }
    }
    keep
  })
  out <- do.call(rbind, picked)
  rownames(out) <- NULL
  out
}
