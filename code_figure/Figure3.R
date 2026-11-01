## Figure3
## Yijun Dong — BMC Cancer revision alignment (2026-11-01)
##
## CRITICAL locks for clinical association panels (Clean/Response):
##   - Ki-67 definition: script_composite_union
##       raw MKI67 count >= 1 OR Phase_enhanced %%in%% c("S","G2/M")
##     NOT manuscript_text_union (norm MKI67>0 OR Seurat Phase S/G2M),
##     which wrongly yields sex P≈0.864, recurrence P≈0.215, age r≈−0.51.
##   - Donor consolidation: LGG-04-1/2/3 -> LGG-04 (38 donors; 24 IDHwt / 14 IDHmt).
##   - Locked stats (Global_Ki67, patient-level):
##       sex Wilcoxon P=0.154; recurrence Wilcoxon P=0.838;
##       age Pearson r=0.365, P≈0.0336, n=34.
##   - Union positives: 23,029 / 8.97% of 256,627.
##
## Prefer precomputed revision metrics CSV when present so clinical panels
## regenerate the locked stats without reloading the full Seurat object.

## prepare data------
## Prefer metadata table when available (faster); fall back to Seurat object.
meta_rds_candidates <- c("results/metadata_df_new.rds", "results/metadata_df.rds")
obj_rds_candidates <- c(
  "revision/20261023_reproduction/exact_cellset/merged_obj_256627_rawnorm.rds",
  "results/merged_obj_new.rds"
)
metrics_csv <- "revision/20261023_reproduction/exact_cellset/definition_sensitivity_20261029/donor_lineage_metrics.csv"
donor_meta_csv <- "revision/20261021/results/donor_metadata.csv"

consolidate_donor <- function(x) {
  sub("^LGG-04-[123]$", "LGG-04", as.character(x))
}

merged_obj <- NULL
df <- NULL
for (p in meta_rds_candidates) {
  if (file.exists(p)) {
    df <- readRDS(p)
    message("Loaded metadata from ", p, " (n=", nrow(df), ")")
    break
  }
}
for (p in obj_rds_candidates) {
  if (file.exists(p)) {
    merged_obj <- readRDS(p)
    if (is.null(df)) df <- merged_obj@meta.data
    message("Loaded Seurat object from ", p)
    break
  }
}
if (is.null(df)) stop("No metadata/Seurat input found for Figure3")
gene <- "MKI67"

## Donor consolidation (40 sample IDs -> 38 donors)
df$Patient_sample <- as.character(df$Patient)
df$Patient <- consolidate_donor(df$Patient_sample)
if (!is.null(merged_obj)) {
  merged_obj$Patient_sample <- as.character(merged_obj$Patient)
  merged_obj$Patient <- consolidate_donor(merged_obj$Patient_sample)
}

## library
library(dplyr)
library(tidyr)
library(mixtools)
library(Seurat)
library(ArchR)
library(ggsci)
library(ggrastr)
library(ggpubr)
library(patchwork)
library(ggnewscale)
library(RColorBrewer)

## color
celltypes_minor_order <- c(
  "AC", "MES", "NPC", "OPC",
  "monocyte-derived DCs", "myeloid-derived suppressor cells", "macrophages", "microglia",
  "T cells", "B cells", "NK cells",
  "fibroblasts", "endothelial cells", "pericytes",
  "oligodendrocytes"
)

merged_obj$celltypes_minor <- factor(merged_obj$celltypes_minor, levels = celltypes_minor_order)

celltypes_color <- ArchRPalettes$stallion[1:length(celltypes_minor_order)]
names(celltypes_color) <- celltypes_minor_order

grade_color <- c("LGG (II-III)" = "navy", "GBM (IV)" = "firebrick3")

ki67_color <- c("Ki67-" = "grey", "Ki67+" = "coral")

lineage_colors <- c(
  "Tumor" = "#D94F4F",
  "Myeloid" = "#0072B5",
  "Lymphoid" = "#00A087",
  "Stromal" = "#EFC000",
  "Glial" = "#AD8BC9"
)

diag_levels <- c("Astrocytoma","Glioblastoma","Oligodendroglioma")
nejm_colors <- pal_nejm()(length(diag_levels))
names(nejm_colors) <- diag_levels

location_color <- c(
  "Frontal"   = "#1b9e77",
  "Temporal"  = "#d95f02",
  "Parietal"  = "#7570b3",
  "Occipital" = "#e7298a",
  "Other"     = "#c7dcd0"
  )

## Figure 3A: Ki67+ fraction with clinical data--------
min_cells <- 10
df <- merged_obj@meta.data

ki67_by_patient_filtered <- df %>%
  filter(!is.na(lineage)) %>%
  group_by(Patient, lineage) %>%
  dplyr::summarise(
    cell_count = n(),
    Ki67_frac = mean(Ki67_Positive, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Ki67_frac = ifelse(cell_count < min_cells, NA, Ki67_frac) 
  ) %>%
  select(-cell_count)

global_ki67 <- df %>%
  group_by(Patient) %>%
  dplyr::summarise(Global_Ki67 = mean(Ki67_Positive, na.rm = TRUE) * 100)

meta <- df %>%
  dplyr::select(Patient, Age, Sex, Grade, Stage, Location, IDH, MGMT.status, Diagnosis, batch) %>%
  distinct()

patient_cells <- table(df$Patient)
tumor_cells <- table(df$lineage,df$Patient)
plot_data <- meta %>%
  left_join(global_ki67,
            by = "Patient"
  ) %>%
  left_join(
    ki67_by_patient_filtered %>%
      tidyr::pivot_wider(names_from = lineage, values_from = Ki67_frac),
    by = "Patient"
  ) %>%
  arrange(Global_Ki67) %>%
  mutate(patient_id = factor(Patient, levels = Patient))

p_main <- ggplot(plot_data, aes(x = patient_id)) +
  geom_col(aes(y = Global_Ki67), fill = "lightgrey", width = 0.8) +
  geom_line(aes(y = Tumor, group = 1, color = "Tumor"), size = 1) +
  geom_point(aes(y = Tumor, color = "Tumor"), size = 2) +
  geom_line(aes(y = Myeloid, group = 1, color = "Myeloid"), size = 1) +
  geom_point(aes(y = Myeloid, color = "Myeloid"), size = 2) +
  geom_line(aes(y = Lymphoid, group = 1, color = "Lymphoid"), size = 1) +
  geom_point(aes(y = Lymphoid, color = "Lymphoid"), size = 2) +
  geom_line(aes(y = Stromal, group = 1, color = "Stromal"), size = 1) +
  geom_point(aes(y = Stromal, color = "Stromal"), size = 2) +
  geom_line(aes(y = Glial, group = 1, color = "Glial"), size = 1) +
  geom_point(aes(y = Glial, color = "Glial"), size = 2) +
  scale_color_manual(values = lineage_colors) +
  labs(y = "Ki-67+ cells (%)", x = NULL, color = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "top")

anno_long <- plot_data %>%
  dplyr::select(patient_id, Age, Sex, Grade, Stage, Location, IDH, MGMT.status, Diagnosis) %>%
  mutate(
    `Sex: Female` = ifelse(Sex == "F", "Yes", "No"),
    `Stage: Recurrent` = ifelse(Stage == "Recurrent", "Yes", "No"),
    `IDH: Wildtype` = ifelse(IDH == "Wildtype", "Yes", "No"),
    `MGMT: methylated` = ifelse(MGMT.status == "methylated", "Yes", "No")
  ) %>%
  dplyr::select(-Sex, -Stage, -IDH, -MGMT.status) %>%
  pivot_longer(-patient_id, names_to = "Feature", values_to = "Value")

p_anno <- ggplot() +
  geom_tile(
    data = anno_long %>% filter(Feature == "Age") %>% mutate(Value = as.numeric(Value)),
    aes(x = patient_id, y = Feature, fill = Value),
    color = "white"
  ) +
  scale_fill_gradient(low = "#f0f0f0", high = "darkgreen", na.value = "grey") +
  
  new_scale_fill() +
  
  geom_tile(
    data = anno_long %>% filter(Feature != "Age"),
    aes(x = patient_id, y = Feature, fill = Value),
    color = "white"
  ) +
  scale_fill_manual(
    values = c(
      grade_color,
      location_color,
      nejm_colors,
      "Yes" = "black",
      "No" = "#f0f0f0"
    ),
    na.value = "grey"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "bottom"
  )

pdf("figures/Fig3A.pdf", height = 8, width = 10)
p_main / p_anno + plot_layout(heights = c(3, 1))
dev.off()

## Figure 3B: statistics of clinical data--------
df_long <- plot_data %>%
  pivot_longer(
    cols = c(Global_Ki67, Lymphoid, Myeloid, Stromal, Tumor, Glial),
    names_to = "lineage",
    values_to = "pct_ki67"
  ) %>%
  filter(!is.na(pct_ki67) & Sex != "NA")
df_long$lineage[df_long$lineage=="Global_Ki67"] <- "Global"
df_long$lineage <- factor(df_long$lineage, levels = c("Global","Tumor","Myeloid","Lymphoid","Glial","Stromal"))

plot_clinical_bar <- function(df_long, group_var, ylab="Ki67+ Fraction") {
  
  df_summary <- df_long %>%
    group_by(lineage, .data[[group_var]]) %>%
    summarise(
      mean_pct = mean(pct_ki67, na.rm = TRUE),
      .groups = "drop"
    )
  
  formula <- as.formula(paste0("pct_ki67 ~ ", group_var))
  pvals <- compare_means(
    formula,
    data = df_long,
    group.by = "lineage",
    method = "wilcox.test"
  )
  
  p <- ggplot(df_summary, aes(x=lineage, y=mean_pct, fill=.data[[group_var]])) +
    geom_bar(
      stat="identity",
      position=position_dodge(width=0.9),
      aes(alpha = ifelse(lineage == "Global", 1, 0.6))
    ) +
    labs(x="Cell lineage", y=ylab, fill=group_var) +
    scale_fill_brewer(palette = "Set1") +
    scale_alpha_identity() + 
    theme_minimal(base_size=14) +
    theme(axis.text.x = element_text(angle=45, hjust=1))
  
  if(nrow(pvals) > 0){
    max_per_lineage <- df_summary %>%
      group_by(lineage) %>%
      summarise(max_y = max(mean_pct, na.rm = TRUE)*1.1, .groups="drop")
    
    pvals <- pvals %>% 
      left_join(max_per_lineage, by=c("lineage"="lineage"))
    
    p <- p + geom_text(
      data=pvals,
      aes(x=lineage, y=max_y, label=p.signif),
      inherit.aes = FALSE,
      position=position_dodge(width=0.9),
      vjust=0
    )
  }
  
  return(p)
}

pdf("figures/Fig3B.pdf", height = 4, width = 6)
plot_clinical_bar(df_long, "Grade")
dev.off()

pdf("figures/Sup_Fig3A.pdf", height = 4, width = 5)
plot_clinical_bar(df_long, "Sex")
plot_clinical_bar(df_long, "Stage")
dev.off()

## Figure 3C: Age vs Global_Ki67 at DONOR level (not Age×Grade cell pooling)
## Cell-level Age×Grade aggregation flips the Pearson sign and must not be used.
age_df <- plot_data %>%
  mutate(Age = suppressWarnings(as.numeric(as.character(Age)))) %>%
  filter(is.finite(Age), is.finite(Global_Ki67))

age_ct <- cor.test(age_df$Age, age_df$Global_Ki67, method = "pearson")
message(sprintf(
  "Fig3C Age Pearson r=%.6f P=%.6g n=%d",
  unname(age_ct$estimate), age_ct$p.value, nrow(age_df)
))

pdf("figures/Fig3C.pdf", height = 5, width = 6)
ggplot(age_df, aes(x = Age, y = Global_Ki67 / 100, color = Grade)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_smooth(aes(group = 1), method = "loess", se = TRUE, size = 1, color = "grey30", fill = "grey80") +
  stat_cor(aes(group = 1), method = "pearson", color = "black") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_minimal() +
  scale_color_manual(values = grade_color) +
  labs(y = "Global Ki-67+ fraction", x = "Age (years)")
dev.off()

tumor_fraction <- df %>%
  group_by(Patient, Grade) %>%
  summarise(pct_tumor = mean(lineage=="Tumor")*100, .groups="drop")

ki67_tumor_fraction <- df %>%
  filter(Ki67_Positive==T) %>%
  group_by(Patient, Grade) %>%
  summarise(pct_tumor = mean(lineage=="Tumor")*100, .groups="drop")

p_frac <- ggplot(ki67_tumor_fraction, aes(x=Grade, y=pct_tumor, fill=Grade)) +
  geom_boxplot(outlier.shape = NA, alpha=0.6) +
  geom_jitter(width=0.2, size=2, alpha=0.7) +
  theme_classic(base_size=14) +
  labs(y="% Ki67+ Cells that are Tumor", x="Grade") +
  scale_fill_manual(values=grade_color) +
  stat_compare_means(method="wilcox.test", label="p.signif")

pdf("figures/Fig3D.pdf", height = 4, width = 5)
p_frac
dev.off()

wilcox.test(pct_tumor ~ Grade, data=ki67_tumor_fraction)

saveRDS(df, file = "results/metadata_df.rds")

plot_data_clean <- plot_data %>%
  filter(!is.na(Global_Ki67)) 

plot_clinical_box <- function(data, x_var, y_var, title, y_lab, palette = "jco") {
  
  sub_data <- data %>% filter(!is.na(.data[[x_var]]) & .data[[x_var]] != "NA")
  
  n_counts <- table(sub_data[[x_var]])
  x_labels <- paste0(names(n_counts), "\n(n=", n_counts, ")")
  names(x_labels) <- names(n_counts)
  
  ggboxplot(sub_data, x = x_var, y = y_var,
            color = x_var, palette = palette,
            add = "jitter", add.params = list(size = 2, alpha = 0.6)) +
    stat_compare_means(method = "wilcox.test", 
                       label = "p.signif", 
                       vjust = 1.5, size = 5) + 
    stat_compare_means(method = "wilcox.test", 
                       label = "p.format", 
                       vjust = 3, size = 3) + 
    scale_x_discrete(labels = x_labels) +
    labs(title = title, y = y_lab, x = NULL) +
    theme_classic(base_size = 12) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
}

## Figure S3A: Global Ki67%----
# Sex (Male vs Female)
p_sex_global <- plot_clinical_box(plot_data_clean, "Sex", "Global_Ki67", 
                                  "Sex", "Global Ki67+ (%)", "npg")

# Stage (Primary vs Recurrent)
p_stage_global <- plot_clinical_box(plot_data_clean, "Stage", "Global_Ki67", 
                                    "Recurrence Status", "Global Ki67+ (%)", "lancet")

# MGMT (Methylated vs Unmethylated)
p_mgmt_global <- plot_clinical_box(plot_data_clean %>% filter(MGMT.status != "Not obtained"), "MGMT.status", "Global_Ki67", 
                                   "MGMT Status", "Global Ki67+ (%)", "d3")

## Figure S3B: Tumor Specific Ki67%----

# Sex
p_sex_tumor <- plot_clinical_box(plot_data_clean, "Sex", "Tumor", 
                                 "Sex (Tumor Only)", "Tumor Ki67+ (%)", "npg")

# Stage
p_stage_tumor <- plot_clinical_box(plot_data_clean, "Stage", "Tumor", 
                                   "Recurrence (Tumor Only)", "Tumor Ki67+ (%)", "lancet")

# MGMT
p_mgmt_tumor <- plot_clinical_box(plot_data_clean %>% filter(MGMT.status != "Not obtained"), "MGMT.status", "Tumor", 
                                  "MGMT (Tumor Only)", "Tumor Ki67+ (%)", "d3")

layout_design <- "
ABC
DEF
"

combined_plot <- (p_sex_global + p_stage_global + p_mgmt_global) /
  (p_sex_tumor + p_stage_tumor + p_mgmt_tumor) +
  plot_annotation(tag_levels = 'A')

pdf("figures/Sup_Fig3.pdf", height = 8, width = 10)
combined_plot
dev.off()


## =============================================================================
## Revision-locked Fig3 clinical associations (script_composite_union)
## Sex / Stage / Age on donor-level Global_Ki67 — must match Clean/Response.
## Regenerates figures/Fig3_clinical_locked.pdf and figures/Fig3_clinical_locked_stats.txt
## =============================================================================
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

if (file.exists(metrics_csv) && file.exists(donor_meta_csv)) {
  metrics <- read.csv(metrics_csv, check.names = FALSE, stringsAsFactors = FALSE)
  meta_d <- read.csv(donor_meta_csv, check.names = FALSE, stringsAsFactors = FALSE)
  key <- intersect(c("donor_id", "donor", "id", "ID"), names(meta_d))[1]
  meta_d$donor_key <- as.character(meta_d[[key]])

  d <- subset(metrics, definition == "script_composite_union")
  if (nrow(d) == 0) stop("script_composite_union not found in metrics CSV")

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
    sprintf("Sex Wilcoxon P=%.8g; group counts=%s",
            psex, paste(names(table(sex$Sex)), as.integer(table(sex$Sex)), collapse = "/")),
    sprintf("Stage Wilcoxon P=%.8g; group counts=%s",
            pstage, paste(names(table(stage$Stage)), as.integer(table(stage$Stage)), collapse = "/")),
    sprintf("Age Pearson r=%.8g; P=%.8g; n=%d",
            unname(ct$estimate), ct$p.value, nrow(age)),
    sprintf("Union positives (sum P)=%.0f / cells (sum N)=%.0f (%.4f%%)",
            sum(agg$P), sum(agg$N), 100 * sum(agg$P) / sum(agg$N))
  )
  writeLines(stats_lines, "figures/Fig3_clinical_locked_stats.txt")
  cat(paste(stats_lines, collapse = "\n"), "\n")

  ## Sanity vs revision locks (tolerances for floating print)
  stopifnot(abs(psex - 0.15360542) < 1e-5)
  stopifnot(abs(pstage - 0.83825649) < 1e-5)
  stopifnot(abs(unname(ct$estimate) - 0.36543776) < 1e-5)
  stopifnot(nrow(age) == 34L)

  sex$group_label <- ifelse(
    sex$Sex == "F",
    sprintf("Female\n(n=%d)", sum(sex$Sex == "F")),
    sprintf("Male\n(n=%d)", sum(sex$Sex == "M"))
  )
  stage$group_label <- ifelse(
    stage$Stage == "Primary",
    sprintf("Primary\n(n=%d)", sum(stage$Stage == "Primary")),
    sprintf("Recurrent\n(n=%d)", sum(stage$Stage == "Recurrent"))
  )
  long <- dplyr::bind_rows(
    data.frame(comparison = "Sex", group = sex$group_label, value = sex$global_fraction * 100),
    data.frame(comparison = "Recurrence status", group = stage$group_label, value = stage$global_fraction * 100)
  )
  long$comparison <- factor(long$comparison, levels = c("Sex", "Recurrence status"))
  pvals <- data.frame(
    comparison = factor(c("Sex", "Recurrence status"), levels = levels(long$comparison)),
    label = sprintf("Wilcoxon P = %.3f", c(psex, pstage))
  )
  pB <- ggplot(long, aes(x = group, y = value, color = group, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.18, width = 0.55) +
    geom_jitter(width = 0.13, size = 1.5, alpha = 0.75) +
    facet_wrap(~comparison, scales = "free_x") +
    geom_text(data = pvals, aes(x = 1.5, y = Inf, label = label),
              inherit.aes = FALSE, vjust = 1.5, size = 3.2, color = "black") +
    labs(title = "Ki-67+ cells across sex and recurrence status",
         x = NULL, y = "Global Ki-67+ cells / all cells (%)") +
    theme_classic(base_size = 12) +
    theme(legend.position = "none", strip.background = element_blank())

  pearson_text <- sprintf("Pearson r = %.2f, P = %.3g, n = %d",
                          unname(ct$estimate), ct$p.value, nrow(age))
  pC <- ggplot(age, aes(x = Age_num, y = global_fraction * 100, color = IDH.x)) +
    geom_point(size = 2, alpha = 0.8) +
    geom_smooth(aes(x = Age_num, y = global_fraction * 100), inherit.aes = FALSE,
                method = "loess", formula = y ~ x, se = TRUE, color = "#555555") +
    annotate("text", x = Inf, y = Inf, label = pearson_text, hjust = 1.05, vjust = 1.3, size = 3.2) +
    labs(title = "Trend of Ki-67 along increasing age",
         x = "Age (years)", y = "Global Ki-67+ cells / all cells (%)", color = "IDH") +
    theme_classic(base_size = 12)

  pdf("figures/Fig3_clinical_locked.pdf", height = 4.5, width = 10)
  print(pB + pC)
  dev.off()
  message("Wrote figures/Fig3_clinical_locked.pdf and figures/Fig3_clinical_locked_stats.txt")
} else {
  warning("Revision metrics CSV not found; skipped locked Fig3 clinical regeneration. ",
          "Expected: ", metrics_csv)
}

## Also recompute sex/stage on plot_data Global_Ki67 after donor consolidation
## (will match locks only when Ki67_Positive encodes script_composite_union).
if (exists("plot_data_clean")) {
  sex_pd <- plot_data_clean %>% dplyr::filter(!is.na(Sex), Sex %in% c("F", "M"), is.finite(Global_Ki67))
  stage_pd <- plot_data_clean %>% dplyr::filter(!is.na(Stage), Stage %in% c("Primary", "Recurrent"), is.finite(Global_Ki67))
  if (nrow(sex_pd) > 0) {
    message("plot_data Sex Wilcoxon P=", wilcox.test(Global_Ki67 ~ Sex, data = sex_pd, exact = FALSE)$p.value)
  }
  if (nrow(stage_pd) > 0) {
    message("plot_data Stage Wilcoxon P=", wilcox.test(Global_Ki67 ~ Stage, data = stage_pd, exact = FALSE)$p.value)
  }
}
