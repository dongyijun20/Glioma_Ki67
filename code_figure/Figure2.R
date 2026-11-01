## Figure2
## Yijun Dong, 15/10/2025
## prepare data------
merged_obj <- readRDS("~/neurosurgery/TTLL7/results/idhmut_strict21_plus_snu4_core_verified_merged.rds")
gene <- "MKI67"

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
library(scales)
library(ggrepel)

## color
celltypes_minor_order <- c(
  "AC", "MES", "NPC", "OPC",
  "macrophages", "microglia",
  "T cells", "B cells",
  "oligodendrocytes"
)

merged_obj$celltypes_minor <- factor(merged_obj$celltypes_minor, levels = celltypes_minor_order)

celltypes_color <- ArchRPalettes$stallion[1:length(celltypes_minor_order)]
names(celltypes_color) <- celltypes_minor_order

grade_color <- c("G2" = "lightblue", "G3" = "orange", "G4" = "firebrick3")

ki67_color <- c("Ki67-" = "grey", "Ki67+" = "coral")

lineage_colors <- c(
  "Tumor" = "#D94F4F",
  "Myeloid" = "#0072B5",
  "Lymphoid" = "#00A087",
  "Stromal" = "#EFC000",
  "Glial" = "#AD8BC9"
)

## Figure 2A-B: definition of Ki67+, phase is S or G2/M, or MKI67 expression > 0---------
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
merged_obj <- CellCycleScoring(merged_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

mki67_expr <- FetchData(merged_obj, vars = "MKI67")$MKI67

# sample_list <- SplitObject(merged_obj, split.by = "Patient")
# 
# run_sample_specific_cc <- function(obj) {
#   DefaultAssay(obj) <- "RNA"
#   obj <- CellCycleScoring(
#     object = obj,
#     s.features = s.genes,
#     g2m.features = g2m.genes,
#     set.ident = T
#   )
#   obj$CC_Difference <- obj$S.Score - obj$G2M.Score
#   return(obj)
# }
# 
# sample_list <- lapply(sample_list, run_sample_specific_cc)
# 
# merged_obj$Phase[unlist(sapply(sample_list, function(x) colnames(x)))] <- unlist(sapply(sample_list, function(x) x$Phase))
# merged_obj$G2M.Score[unlist(sapply(sample_list, function(x) colnames(x)))] <- unlist(sapply(sample_list, function(x) x$G2M.Score))
# merged_obj$S.Score[unlist(sapply(sample_list, function(x) colnames(x)))] <- unlist(sapply(sample_list, function(x) x$S.Score))
# merged_obj$CC_Difference <- NA
# merged_obj$CC_Difference[unlist(sapply(sample_list, function(x) colnames(x)))] <- unlist(sapply(sample_list, function(x) x$CC_Difference))
# 
# merged_obj$Phase_enhanced <- case_when(
#   merged_obj$G2M.Score > 0.1 & merged_obj$G2M.Score > merged_obj$S.Score ~ "G2/M",
#   merged_obj$S.Score > 0.1 & merged_obj$S.Score > merged_obj$G2M.Score ~ "S",
#   abs(merged_obj$CC_Difference) < 0.1 ~ "G1/S",
#   TRUE ~ "G1"
# )

pdf("figure_v2/panels/Fig2A.pdf", height = 3, width = 5)
# RidgePlot(merged_obj, features = gene,
#           assay = "RNA", group.by = "Phase_enhanced") +
#   scale_x_continuous(limits = c(0.01, NA))
RidgePlot(merged_obj, features = gene,
          assay = "RNA", group.by = "Phase") +
  scale_x_continuous(limits = c(0.01, NA))
DimPlot(merged_obj, group.by = "Phase")
dev.off()

table(merged_obj$Phase_enhanced)
table(merged_obj$Phase)

merged_obj$Ki67_status <- ifelse(
  merged_obj$Phase == c("G2M") | mki67_expr > 0, 
  "Ki67+",
  "Ki67-"
)
table(merged_obj$Ki67_status)

df <- merged_obj@meta.data
df$MKI67 <- mki67_expr
df <- df %>%
  mutate(
    Ki67_Positive = Ki67_status == "Ki67+"
  )

saveRDS(df, file = "results/metadata_df_IDHmt_only.rds")

df <- readRDS("results/metadata_df_new.rds")

ki67_high_pct <- round(mean(df$Ki67_status == "Ki67+") * 100, 2)
ki67_high_pct

df_density <- df[df$Ki67_status == "Ki67+", ]

p_density <- ggplot(df_density, aes(x = MKI67)) +
  geom_density(fill = "coral", alpha = 0.5) +
  labs(y = "Density") +
  theme_minimal(base_size = 12) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

p_points <- ggplot(df, aes(x = MKI67, y = 1, color = Ki67_status)) +
  geom_point_rast(position = position_jitter(height = 0.1), size = 0.4, alpha = 0.4) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = ki67_color) +
  labs(x = "MKI67 Expression", y = "") +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank())

combined_plot <- p_density / p_points + plot_layout(heights = c(2, 1))

pdf("figure_v2/panels/Fig2B.pdf", height = 2.5, width = 5)
combined_plot
dev.off()

merged_obj_subset <- merged_obj[,merged_obj$Ki67_status == "Ki67+"]

## Figure 2C-F: MKI67 expression in Ki67+ cells--------
expr_patient <- df[,c(gene, "lineage", "Grade", "Diagnosis", "Patient", "Tumor.Grade")] %>%
  dplyr::group_by(Patient, lineage, Grade, Diagnosis, Tumor.Grade) %>%
  dplyr::summarise(
    mean_expr = mean(get(gene), na.rm = TRUE), 
    n_cell = n(), 
    .groups = 'drop'
  ) %>%
  dplyr::mutate(mean_expr_plot = mean_expr + 0.01)

plot_compact <- function(df, x_var, y_var, fill_var, fill_colors = NULL, is_nejm = FALSE) {
  lvls <- levels(factor(df[[x_var]]))
  comp_list <- combn(lvls, 2, simplify = FALSE)
  
  p <- ggboxplot(df, x = x_var, y = y_var, fill = fill_var,
                 add = "jitter", add.params = list(size = 1.2, alpha = 0.5),
                 alpha = 0.6, size = 0.3) +
    stat_compare_means(comparisons = comp_list,
                       method = "wilcox.test",
                       label = "p.signif", 
                       hide.ns = TRUE,
                       tip.length = 0.01,
                       step.increase = 0.05) +
    theme_classic(base_size = 14) +
    scale_y_log10() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    ylab(paste("Mean", gene))
  
  if(is_nejm) {
    p <- p + scale_fill_nejm()
  } else if (!is.null(fill_colors)) {
    p <- p + scale_fill_manual(values = fill_colors)
  }
  return(p)
}

p1 <- plot_compact(expr_patient, "lineage", "mean_expr", "lineage", lineage_colors)
p2 <- plot_compact(expr_patient[expr_patient$Grade == "LGG (II-III)",], "Tumor.Grade", "mean_expr", "Tumor.Grade")

expr_patient$Diagnosis <- factor(expr_patient$Diagnosis, 
                                 levels = c("Astrocytoma","Glioblastoma","Oligodendroglioma"))
p3 <- plot_compact(expr_patient, "Diagnosis", "mean_expr", "Diagnosis", is_nejm = TRUE)

pdf("figures/Fig2C.pdf", height = 5, width = 6)
p1
p2
p3
dev.off()

pdf("figure_v2/panels/Fig2D.pdf", height = 5, width = 9)
DimPlot(merged_obj, group.by = "Ki67_status", cols = ki67_color, split.by = "Grade")
dev.off()

## Figure 2G-H: compare Ki67+ fraction-------
ki67_by_grade_celltype <- df %>%
  dplyr::count(Grade, lineage, Ki67_status) %>%
  group_by(Grade, lineage) %>%
  mutate(Ki67_frac = n / sum(n)) %>%
  ungroup() %>%
  na.omit()

ki67_by_grade_global <- df %>%
  dplyr::count(Grade, Ki67_status) %>%
  group_by(Grade) %>%
  mutate(Ki67_frac = n / sum(n)) %>%
  ungroup() %>%
  na.omit()

library(ggbreak)
p_global <- ggplot(ki67_by_grade_global, aes(x = Grade, y = Ki67_frac*100, fill = Ki67_status)) +
  geom_bar(stat="identity", position="stack") +
  scale_fill_manual(values = ki67_color) +
  scale_y_break(c(20,100), space=0.05) +
  labs(y="Ki-67 cell %", x=NULL, fill="Ki-67 group") +
  theme_classic(base_size=14) +
  theme(
    axis.text.x = element_text(angle=45, hjust=1),
    legend.position="right"
  )

p_by_celltype <- ggplot(ki67_by_grade_celltype, aes(x=lineage, y=Ki67_frac*100, fill=Ki67_status)) +
  geom_bar(stat="identity", position="stack") +
  scale_fill_manual(values = ki67_color) +
  scale_y_break(c(20,100), scales=0.05, ticklabels=c(seq(0,20,10),100)) +
  labs(y=NULL, x="Cell type", fill="Ki-67 group") +
  theme_classic(base_size=14) +
  theme(
    axis.text.x = element_text(angle=45, hjust=1),
    axis.ticks.y = element_blank()
  ) +
  facet_wrap(~Grade, nrow=1)

pdf("figure_v2/panels/Fig2G.pdf", height = 3, width = 8)
(p_global | p_by_celltype) + plot_layout(widths = c(3, 4.5))
dev.off()

tumor_rate <- df %>%
  filter(lineage == "Tumor") %>%
  group_by(Patient, Grade) %>%
  dplyr::summarise(pct_ki67_tumor = mean(Ki67_Positive)*100, .groups="drop")

levels_grp <- unique(tumor_rate$Grade)
comparisons_list <- combn(levels_grp, 2, simplify = FALSE)

p_tumor <- ggboxplot(tumor_rate, x = "Grade", y = "pct_ki67_tumor", fill = "Grade",
                     alpha = 0.6, size = 0.3, add = "jitter", add.params = list(width = 0.2, size = 2, alpha=0.7)) +
  stat_compare_means(comparisons = comparisons_list,
                     method = "wilcox.test",
                     label = "p.signif", hide.ns = TRUE) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = grade_color) +
  labs(y = "% Ki67+ Tumor Cells", x = "Grade")

nontumor_rate <- df %>%
  filter(lineage != "Tumor") %>%
  group_by(Patient, Grade) %>%
  dplyr::summarise(pct_ki67_nontumor = mean(Ki67_Positive)*100, .groups="drop")

p_nontumor <- ggboxplot(nontumor_rate, x = "Grade", y = "pct_ki67_nontumor", fill = "Grade",
                        alpha = 0.6, size = 0.3, add = "jitter", add.params = list(width = 0.2, size = 2, alpha=0.7)) +
  stat_compare_means(comparisons = comparisons_list,
                     method = "wilcox.test",
                     label = "p.signif", hide.ns = TRUE) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = grade_color) +
  labs(y = "% Ki67+ Non-tumor Cells", x = "Grade")

pdf("figure_v2/panels/Fig2H.pdf", height = 5, width = 8)
p_tumor + p_nontumor
dev.off()

## Figure 2I: proportion of Ki67+ cells-------
plot_double_donut <- function(meta_df, title = "", all_colors) {
  
  outer_df <- meta_df %>%
    dplyr::count(lineage, celltypes_minor) %>%
    dplyr::arrange(lineage, dplyr::desc(n)) %>%
    dplyr::mutate(
      prop = n / sum(n),
      ymax = cumsum(prop),
      ymin = c(0, head(ymax, n = -1)),
      ymid = (ymin + ymax) / 2
    )
  
  inner_df <- outer_df %>%
    dplyr::group_by(lineage) %>%
    dplyr::summarise(
      ymin = min(ymin),
      ymax = max(ymax),
      n = sum(n),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      ymid = (ymin + ymax) / 2
    )
  

  ggplot() +
    geom_rect(data = inner_df, 
              aes(ymax = ymax, ymin = ymin, xmax = 2, xmin = 1, fill = lineage),
              color = "white", linewidth = 0.5) +
    
    geom_rect(data = outer_df, 
              aes(ymax = ymax, ymin = ymin, xmax = 3, xmin = 2, fill = celltypes_minor),
              color = "white", linewidth = 0.5) +
    
    geom_text(data = inner_df, 
              aes(x = 1.5, y = ymid, label = lineage),
              color = "black", size = 3.5, fontface = "bold") +
    
    geom_text_repel(data = outer_df, 
                    aes(x = 3, y = ymid, 
                        label = paste0(celltypes_minor, "\n(", percent(prop, 0.1), ")")),
                    size = 3.5, 
                    nudge_x = 1,   
                    direction = "y", 
                    segment.size = 0.4, 
                    segment.color = "grey50",
                    show.legend = FALSE) +
    
    coord_polar(theta = "y") +
    xlim(c(0, 4.5)) +
    
    theme_void() +
    scale_fill_manual(values = all_colors) +
    ggtitle(title) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      legend.position = "none"
    )
}


all_colors <- c(lineage_colors, celltypes_color)

meta_df <- df %>% filter(Ki67_Positive == TRUE)
grades <- unique(meta_df$Grade)
plots <- list()

for (g in grades) {
  message("Plotting Grade ", g)
  meta_g <- meta_df %>% filter(Grade == g)
  
  p <- plot_double_donut(meta_g, title = paste("Grade", g), all_colors = all_colors)
  plots[[g]] <- p
  
  safe_filename <- gsub("[^A-Za-z0-9_.-]", "_", g)
  ggsave(filename = paste0("figure_v2/panels/Fig2_", safe_filename, "_DoubleDonut.pdf"), 
         plot = p, width = 7, height = 7)
}

## Supp Fig 1C: Pie chart of all cells----
grades <- unique(df$Grade)
plots <- list()

for (g in grades) {
  message("Plotting Grade ", g)
  meta_g <- df %>% filter(Grade == g)
  
  p <- plot_double_donut(meta_g, title = paste("Grade", g), all_colors = all_colors)
  plots[[g]] <- p
  
  safe_filename <- gsub("[^A-Za-z0-9_.-]", "_", g)
  ggsave(filename = paste0("figure_v2/supplementary/panels/Sup_Fig1C_", safe_filename, "_DoubleDonut.pdf"), 
         plot = p, width = 7, height = 7)
}

## Supp Fig 2A: Technical Validation of Cell Cycle Scoring----
# Scatter plot to show the relationship between S and G2M scores
p_cc_scatter <- FeatureScatter(merged_obj, feature1 = "S.Score", feature2 = "G2M.Score", group.by = "Phase_enhanced") +
  theme_classic() +
  labs(title = "S Score vs G2M Score Correlation", subtitle = "Validation of Phase Assignment")

# UMAP showing cell cycle scores
p_s_score <- FeaturePlot(merged_obj, features = "S.Score", order = TRUE) + 
  scale_color_viridis_c() + theme_void() + ggtitle("S Phase Score")
p_g2m_score <- FeaturePlot(merged_obj, features = "G2M.Score", order = TRUE) + 
  scale_color_viridis_c() + theme_void() + ggtitle("G2M Phase Score")

pdf("figure_v2/supplementary/panels/Sup_Fig2A.pdf", height = 8, width = 10)
p_cc_scatter | (p_s_score / p_g2m_score)
dev.off()

## Supp Fig 2B: Patient Heterogeneity----
# Calculate Ki67+ percentage per patient to rule out outlier effects
patient_variability <- df %>%
  group_by(Patient, Grade) %>%
  dplyr::summarise(Ki67_Pct = mean(Ki67_Positive) * 100, .groups = "drop") %>%
  arrange(Grade, Ki67_Pct)

patient_variability$Patient <- factor(patient_variability$Patient, levels = patient_variability$Patient)

p_patient <- ggplot(patient_variability, aes(x = Patient, y = Ki67_Pct, fill = Grade)) +
  geom_bar(stat = "identity", width = 0.8, color="black", size=0.2) +
  scale_fill_manual(values = grade_color) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10)) +
  labs(y = "% Ki67+ Cells (per Patient)", x = "Patient ID", title = "Inter-patient Heterogeneity of Proliferation")

pdf("figure_v2/supplementary/panels/Sup_Fig2B.pdf", height = 5, width = 8)
p_patient
dev.off()

## Supp Fig 2C: Detailed Lineage Analysis----
# DotPlot for MKI67 in minor cell types
merged_obj$celltypes_minor <- factor(merged_obj$celltypes_minor, levels = celltypes_minor_order)

p_dot <- DotPlot(merged_obj, features = "MKI67", group.by = "celltypes_minor", cols = c("lightgrey", "firebrick")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = NULL, x = NULL, title = "MKI67 Expression by Cell Type Subsets")

# Barplot for Ki67 fraction in minor cell types
detailed_rate <- df %>%
  group_by(Patient, celltypes_minor, Grade) %>%
  dplyr::summarise(pct = mean(Ki67_Positive) * 100, .groups = "drop") %>%
  filter(!is.na(celltypes_minor))

p_detailed_bar <- ggplot(detailed_rate, aes(x = celltypes_minor, y = pct, fill = Grade)) +
  stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 0.8), width = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(width = 0.8), width = 0.2) +
  stat_compare_means(aes(group = Grade), label = "p.signif", method = "wilcox.test", label.y = 15) + 
  scale_fill_manual(values = grade_color) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "% Ki67+ Cells", x = NULL, title = "Proliferation Index by Cell Type")

pdf("figure_v2/supplementary/panels/Sup_Fig2C.pdf", height = 6, width = 8)
p_detailed_bar
dev.off()

