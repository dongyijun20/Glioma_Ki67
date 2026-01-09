## Figure2
## Yijun Dong, 15/10/2025
## prepare data------
merged_obj <- readRDS("results/merged_obj_new.rds")
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
library(plotly)
library(scales)

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

## Figure 2A-B: definition of Ki67+, phase is S or G2/M, or MKI67 expression > 0---------
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

mki67_expr <- FetchData(merged_obj, vars = "MKI67")$MKI67

sample_list <- SplitObject(merged_obj, split.by = "Patient")

run_sample_specific_cc <- function(obj) {
  DefaultAssay(obj) <- "RNA"
  obj <- CellCycleScoring(
    object = obj,
    s.features = s.genes,
    g2m.features = g2m.genes,
    set.ident = T
  )
  obj$CC_Difference <- obj$S.Score - obj$G2M.Score
  return(obj)
}

sample_list <- lapply(sample_list, run_sample_specific_cc)

merged_obj$Phase[unlist(sapply(sample_list, function(x) colnames(x)))] <- unlist(sapply(sample_list, function(x) x$Phase))
merged_obj$G2M.Score[unlist(sapply(sample_list, function(x) colnames(x)))] <- unlist(sapply(sample_list, function(x) x$G2M.Score))
merged_obj$S.Score[unlist(sapply(sample_list, function(x) colnames(x)))] <- unlist(sapply(sample_list, function(x) x$S.Score))
merged_obj$CC_Difference <- NA
merged_obj$CC_Difference[unlist(sapply(sample_list, function(x) colnames(x)))] <- unlist(sapply(sample_list, function(x) x$CC_Difference))

merged_obj$Phase_enhanced <- case_when(
  merged_obj$G2M.Score > 0.1 & merged_obj$G2M.Score > merged_obj$S.Score ~ "G2/M",
  merged_obj$S.Score > 0.1 & merged_obj$S.Score > merged_obj$G2M.Score ~ "S",
  abs(merged_obj$CC_Difference) < 0.1 ~ "G1/S",
  TRUE ~ "G1"
)

pdf("figures/Fig2A.pdf", height = 3, width = 5)
RidgePlot(merged_obj, features = gene,
          assay = "RNA", group.by = "Phase_enhanced") +
  scale_x_continuous(limits = c(0.01, NA))
dev.off()

table(merged_obj$Phase_enhanced)

merged_obj$Ki67_status <- ifelse(
  merged_obj$Phase_enhanced %in% c("G2/M","S") | mki67_expr > 0, 
  "Ki67+",
  "Ki67-"
)

df <- merged_obj@meta.data
df$MKI67 <- mki67_expr
df <- df %>%
  mutate(
    Ki67_Positive = Ki67_status == "Ki67+"
  ) %>% 
  group_by(Patient) %>% 
  filter(n() > 500) %>% 
  ungroup()
saveRDS(df, file = "results/metadata_df_new.rds")

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

pdf("figures/Fig2B.pdf", height = 2.5, width = 5)
combined_plot
dev.off()

merged_obj_subset <- merged_obj[,merged_obj$Ki67_status == "Ki67+"]

## Figure 2C-F: MKI67 expression in Ki67+ cells--------
expr_data <- FetchData(merged_obj_subset, vars = c(gene, "lineage", "Grade", "Diagnosis"))

levels_grp <- unique(expr_data[["lineage"]])
comparisons_list <- combn(levels_grp, 2, simplify = FALSE)

p1 <- ggboxplot(expr_data, x = "lineage", y = gene, fill = "lineage",
                alpha = 0.6, size = 0.3) +
  stat_compare_means(comparisons = comparisons_list,
                     method = "wilcox.test",
                     label = "p.signif", hide.ns = TRUE) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_fill_manual(values = lineage_colors)

levels_grp <- unique(expr_data[["Grade"]])
comparisons_list <- combn(levels_grp, 2, simplify = FALSE)

p2 <- ggboxplot(expr_data, x = "Grade", y = gene, fill = "Grade",
                alpha = 0.6, size = 0.3) +
  stat_compare_means(comparisons = comparisons_list,
                     method = "wilcox.test",
                     label = "p.signif", hide.ns = TRUE) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = grade_color)

levels_grp <- unique(expr_data[["Diagnosis"]])
comparisons_list <- combn(levels_grp, 2, simplify = FALSE)
expr_data[["Diagnosis"]] <- factor(expr_data[["Diagnosis"]], levels = c("Astrocytoma","Glioblastoma","Oligodendroglioma"))

p3 <- ggboxplot(expr_data, x = "Diagnosis", y = gene, fill = "Diagnosis",
                alpha = 0.6, size = 0.3) +
  stat_compare_means(comparisons = comparisons_list,
                     method = "wilcox.test",
                     label = "p.signif", hide.ns = TRUE) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_nejm()

pdf("figures/Fig2C.pdf", height = 5, width = 5)
p1
p2
p3
dev.off()

pdf("figures/Fig2D.pdf", height = 5, width = 9)
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

pdf("figures/Fig2G.pdf", height = 3, width = 8)
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

pdf("figures/Fig2H.pdf", height = 5, width = 8)
p_tumor + p_nontumor
dev.off()

## Figure 2I: proportion of Ki67+ cells-------
prepare_sunburst_data <- function(meta_df) {
  outer <- meta_df %>%
    dplyr::count(lineage, celltypes_minor, name = "value") %>%
    mutate(
      ids = paste(lineage, celltypes_minor, sep = "-"),
      labels = celltypes_minor,
      parents = lineage,
      color = celltypes_color[celltypes_minor]
    )
  
  inner <- meta_df %>%
    dplyr::count(lineage, name = "value") %>%
    mutate(
      ids = lineage,
      labels = lineage,
      parents = "",
      color = lineage_colors[lineage]
    )
  
  bind_rows(inner, outer)
}

plot_sunburst_chart <- function(sunburst_df, title = "") {
  plot_ly(
    type = "sunburst",
    ids = sunburst_df$ids,
    labels = sunburst_df$labels,
    parents = sunburst_df$parents,
    values = sunburst_df$value,
    marker = list(colors = sunburst_df$color),
    textinfo = "label+percent entry",
    branchvalues = "total"
  ) %>%
    layout(title = list(text = title, x = 0.5))
}

meta_df <- df %>% 
  filter(Ki67_Positive == T)
table(meta_df$Grade)
sunburst_df <- prepare_sunburst_data(meta_df)
grades <- unique(meta_df$Grade)

plots <- list()

for (g in grades) {
  message("Plotting Grade ", g)
  meta_g <- meta_df %>% filter(Grade == g)
  sunburst_df <- prepare_sunburst_data(meta_g)
  plots[[g]] <- plot_sunburst_chart(sunburst_df, title = paste("Grade", g))
}

htmlwidgets::saveWidget(plots[["LGG (II-III)"]], "figures/Fig2I.html")
htmlwidgets::saveWidget(plots[["GBM (IV)"]], "figures/Fig2J.html")



