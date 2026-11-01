## Standalone regeneration of revision-locked Fig3 clinical panels/stats.
## Uses script_composite_union donor metrics (NOT manuscript_text_union).
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})
root <- "/publicfs10/fs10-m9/home/m9s004832/TTLL7"
setwd(root)
metrics_csv <- "revision/20261023_reproduction/exact_cellset/definition_sensitivity_20261029/donor_lineage_metrics.csv"
donor_meta_csv <- "revision/20261021/results/donor_metadata.csv"
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

metrics <- read.csv(metrics_csv, check.names = FALSE, stringsAsFactors = FALSE)
meta_d <- read.csv(donor_meta_csv, check.names = FALSE, stringsAsFactors = FALSE)
meta_d$donor_key <- as.character(meta_d$donor_id)
d <- subset(metrics, definition == "script_composite_union")
agg <- aggregate(cbind(N, P) ~ donor + IDH, d, sum, na.rm = TRUE)
agg$global_fraction <- agg$P / agg$N
agg$donor_key <- as.character(agg$donor)
z <- merge(agg, meta_d, by = "donor_key", all = FALSE)
z$Age_num <- suppressWarnings(as.numeric(as.character(z$Age)))
sex <- z[!is.na(z$Sex) & is.finite(z$global_fraction), ]
stage <- z[!is.na(z$Stage) & is.finite(z$global_fraction), ]
age <- z[is.finite(z$Age_num) & is.finite(z$global_fraction), ]
psex <- wilcox.test(global_fraction ~ Sex, data = sex, exact = FALSE)$p.value
pstage <- wilcox.test(global_fraction ~ Stage, data = stage, exact = FALSE)$p.value
ct <- cor.test(age$Age_num, age$global_fraction, method = "pearson")
stats_lines <- c(
  sprintf("Definition=script_composite_union; donors=%d; age complete=%d; sex complete=%d; recurrence complete=%d",
          nrow(z), nrow(age), nrow(sex), nrow(stage)),
  sprintf("Sex Wilcoxon P=%.8g; group counts=%s", psex, paste(names(table(sex$Sex)), as.integer(table(sex$Sex)), collapse = "/")),
  sprintf("Stage Wilcoxon P=%.8g; group counts=%s", pstage, paste(names(table(stage$Stage)), as.integer(table(stage$Stage)), collapse = "/")),
  sprintf("Age Pearson r=%.8g; P=%.8g; n=%d", unname(ct$estimate), ct$p.value, nrow(age)),
  sprintf("Union positives (sum P)=%.0f / cells (sum N)=%.0f (%.4f%%)", sum(agg$P), sum(agg$N), 100 * sum(agg$P) / sum(agg$N))
)
writeLines(stats_lines, "figures/Fig3_clinical_locked_stats.txt")
cat(paste(stats_lines, collapse = "\n"), "\n")
stopifnot(abs(psex - 0.15360542) < 1e-5)
stopifnot(abs(pstage - 0.83825649) < 1e-5)
stopifnot(abs(unname(ct$estimate) - 0.36543776) < 1e-5)
stopifnot(nrow(age) == 34L)

sex$group_label <- ifelse(sex$Sex == "F", sprintf("Female\n(n=%d)", sum(sex$Sex == "F")), sprintf("Male\n(n=%d)", sum(sex$Sex == "M")))
stage$group_label <- ifelse(stage$Stage == "Primary", sprintf("Primary\n(n=%d)", sum(stage$Stage == "Primary")), sprintf("Recurrent\n(n=%d)", sum(stage$Stage == "Recurrent")))
long <- bind_rows(
  data.frame(comparison = "Sex", group = sex$group_label, value = sex$global_fraction * 100),
  data.frame(comparison = "Recurrence status", group = stage$group_label, value = stage$global_fraction * 100)
)
long$comparison <- factor(long$comparison, levels = c("Sex", "Recurrence status"))
pvals <- data.frame(comparison = factor(c("Sex", "Recurrence status"), levels = levels(long$comparison)),
                    label = sprintf("Wilcoxon P = %.3f", c(psex, pstage)))
pB <- ggplot(long, aes(x = group, y = value, color = group, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = .18, width = .55) + geom_jitter(width = .13, size = 1.5, alpha = .75) +
  facet_wrap(~comparison, scales = "free_x") +
  geom_text(data = pvals, aes(x = 1.5, y = Inf, label = label), inherit.aes = FALSE, vjust = 1.5, size = 3.2, color = "black") +
  labs(title = "Ki-67+ cells across sex and recurrence status", x = NULL, y = "Global Ki-67+ cells / all cells (%)") +
  theme_classic(base_size = 12) + theme(legend.position = "none", strip.background = element_blank())
pearson_text <- sprintf("Pearson r = %.2f, P = %.3g, n = %d", unname(ct$estimate), ct$p.value, nrow(age))
idh_col <- if ("IDH.x" %in% names(age)) age$IDH.x else age$IDH
pC <- ggplot(age, aes(x = Age_num, y = global_fraction * 100, color = idh_col)) +
  geom_point(size = 2, alpha = .8) +
  geom_smooth(aes(x = Age_num, y = global_fraction * 100), inherit.aes = FALSE, method = "loess", formula = y ~ x, se = TRUE, color = "#555555") +
  annotate("text", x = Inf, y = Inf, label = pearson_text, hjust = 1.05, vjust = 1.3, size = 3.2) +
  labs(title = "Trend of Ki-67 along increasing age", x = "Age (years)", y = "Global Ki-67+ cells / all cells (%)", color = "IDH") +
  theme_classic(base_size = 12)
pdf("figures/Fig3_clinical_locked.pdf", height = 4.5, width = 10)
print(pB + pC)
dev.off()
message("OK locked Fig3 clinical panels verified")
