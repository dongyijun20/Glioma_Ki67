## Figure4
## Yijun Dong, 17/10/2025
## prepare data------
merged_obj <- readRDS("results/merged_obj_new.rds")
gene <- "MKI67"

## library
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(scales)
library(tidyr)
library(pheatmap)
library(textshape)
library(RColorBrewer)
library(ggtern)
library(patchwork)
library(GeneNMF)
library(Matrix)
library(RcppML)
library(viridis)
library(clusterProfiler)
library(org.Hs.eg.db)
library(UCell)

## color
grade_color <- c("LGG (II-III)" = "navy", "GBM (IV)" = "firebrick3")

ki67_color <- c("Ki67-" = "grey", "Ki67+" = "coral")

lineage_colors <- c(
  "Tumor" = "#D94F4F",
  "Myeloid" = "#0072B5",
  "Lymphoid" = "#00A087",
  "Stromal" = "#EFC000",
  "Glial" = "#AD8BC9"
)


## Figure 4B: ternary plot of three ecotypes-----
df <- merged_obj@meta.data

ki67_by_patient_filtered <- df %>%
  dplyr::filter(!is.na(lineage)) %>%
  group_by(Patient, lineage, Ki67_status) %>%
  dplyr::summarise(
    cell_count = n(),
    .groups = "drop"
  ) %>%
  dplyr::filter(Ki67_status == "Ki67+") %>%
  dplyr::select(-Ki67_status)

global_ki67 <- ki67_by_patient_filtered %>%
  group_by(Patient) %>%
  dplyr::summarise(Global_Ki67 = sum(cell_count), .groups = "drop")

meta <- df %>%
  dplyr::select(Patient, Age, Sex, Grade, Stage, Location, IDH, MGMT.status, Diagnosis, batch) %>%
  distinct()

plot_data <- meta %>%
  left_join(global_ki67, by = "Patient") %>%
  left_join(
    ki67_by_patient_filtered %>%
      pivot_wider(names_from = lineage, values_from = cell_count, values_fill = 0),
    by = "Patient"
  ) %>%
  dplyr::arrange(Global_Ki67) %>%
  dplyr::mutate(patient_id = factor(Patient, levels = Patient))

ternary_data <- plot_data %>%
  mutate(
    Immune = Myeloid + Lymphoid,
    Immune_frac = Immune / (Immune + Tumor + Glial),
    Tumor_frac  = Tumor  / (Immune + Tumor + Glial),
    Glial_frac  = Glial  / (Immune + Tumor + Glial)
  )

ternary_data <- ternary_data %>%
  rowwise() %>%
  mutate(
    dist_Immune = sqrt((Immune_frac - 1)^2 + (Tumor_frac)^2 + (Glial_frac)^2),
    dist_Tumor  = sqrt((Immune_frac)^2 + (Tumor_frac - 1)^2 + (Glial_frac)^2),
    dist_Glial  = sqrt((Immune_frac)^2 + (Tumor_frac)^2 + (Glial_frac - 1)^2)
  )

ternary_data <- ternary_data %>%
  mutate(
    Cluster = case_when(
      dist_Immune == pmin(dist_Immune, dist_Tumor, dist_Glial) ~ "Immune-dominant",
      dist_Tumor  == pmin(dist_Immune, dist_Tumor, dist_Glial) ~ "Tumor-dominant",
      dist_Glial  == pmin(dist_Immune, dist_Tumor, dist_Glial) ~ "Glial-dominant"
    )
  )
ternary_data$Cluster[ternary_data$Patient=="SF11488"] <- "Stromal-dominant"

pdf("figures/Fig4A.pdf",height = 4, width = 5)
ggtern(data = ternary_data, aes(x = Immune_frac, y = Tumor_frac, z = Glial_frac, color = Cluster)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "Patient-level Ki67⁺ composition ternary plot",
    x = "Immune",
    y = "Tumor",
    z = "Glial"
  ) +
  scale_color_manual(values = cluster_colors)
dev.off()

## Figure 4B: four ecotype clustering--------
plot_data$Cluster <- ternary_data$Cluster
plot_data$Cluster <- factor(plot_data$Cluster, 
                            levels = c("Tumor-dominant", "Immune-dominant", "Glial-dominant","Stromal-dominant"))

plot_data$Grade <- factor(plot_data$Grade,
                          levels = c("LGG (II-III)", "GBM (IV)"))

col_order <- order(plot_data$Cluster, plot_data$Grade)

heatmap_data <- plot_data %>%
  mutate(across(c(Lymphoid, Myeloid, Tumor, Glial, Stromal), ~ . / Global_Ki67 * 100),
         Immune = Lymphoid + Myeloid) %>%
  dplyr::select(Patient, Immune, Tumor, Glial, Stromal) %>%
  column_to_rownames("Patient") %>%
  as.matrix()
heatmap_data <- heatmap_data[col_order, ]

annotation_col <- plot_data %>%
  dplyr::select(Patient, Grade, Cluster) %>%
  column_to_rownames("Patient")
annotation_col <- annotation_col[col_order, ]

gaps <- cumsum(table(plot_data$Cluster[order(plot_data$Cluster)]))

cluster_colors <- brewer.pal(length(levels(plot_data$Cluster)), "Set2")
names(cluster_colors) <- levels(plot_data$Cluster)
annotation_colors <- list(
  Grade = grade_color,
  Cluster = cluster_colors
)

pdf("figures/Fig4B.pdf", height = 2.5, width = 10)
pheatmap(
  t(heatmap_data),
  cluster_cols = FALSE, 
  cluster_rows = FALSE,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  gaps_col = gaps,
  color = colorRampPalette(c("#fff5eb", "#fdae6b", "#e6550d"))(100),
  legend = TRUE,
  show_colnames = FALSE,
  main = "Ki67⁺ lineage composition across glioma patients"
)
dev.off()

## Figure4C: UMAP of representative samples-------
ki67_cluster <- plot_data$Cluster[match(merged_obj$Patient, plot_data$Patient)]
names(ki67_cluster) <- colnames(merged_obj)
merged_obj$Ki67_cluster <- ki67_cluster
merged_obj$Ki67_status <- df$Ki67_status
saveRDS(merged_obj, file = "results/merged_obj_new.rds")

representative_samples <- c("GBM057","GBML1","GBM068")

subset_obj <- subset(merged_obj, subset = Patient %in% representative_samples)

p_celltype <- DimPlot(subset_obj, reduction = "umap", group.by = "lineage", label = TRUE, cols = lineage_colors) +
  theme(aspect.ratio = 1)

Idents(subset_obj) <- subset_obj$Ki67_cluster
p_ki67 <- DimPlot(subset_obj, reduction = "umap", group.by = "Ki67_cluster", label = TRUE, cols = cluster_colors) +
  scale_color_manual(values = c("Ki67+" = "red", "Ki67-" = "grey70")) +
  theme(aspect.ratio = 1)

plot_list <- list()
for(e in unique(subset_obj$Ki67_cluster)){
  obj_e <- subset(subset_obj, subset = Ki67_cluster == e)
  p1 <- DimPlot(obj_e, reduction = "umap", group.by = "lineage", label = TRUE) + ggtitle(paste0(e, " cell types"))
  p2 <- DimPlot(obj_e, reduction = "umap", group.by = "Ki67_status", order = "Ki67+") +
    scale_color_manual(values = c("Ki67+" = "red", "Ki67-" = "grey")) +
    ggtitle(paste0(e, " Ki67+ cells"))
  plot_list[[e]] <- p1 + p2
}

pdf("figures/Fig4C.pdf", width = 15, height = 6)
wrap_plots(plot_list, ncol = 2)
dev.off()

# Figure 4D: NMF-------
seu.list <- SplitObject(merged_obj, split.by = "Patient")

geneNMF.programs <- multiNMF(
  seu.list,
  assay = "RNA",
  slot = "data",
  k = 3:10,    
  L1 = c(0,0),    
  center = TRUE,
  scale = TRUE,
  nfeatures = 2000
)

geneNMF.metaprograms <- getMetaPrograms(
  geneNMF.programs,
  nMP = 10,      
  max.genes = 50,    
  hclust.method = "ward.D2",
  min.confidence = 0.3
)

# plotMetaPrograms(geneNMF.metaprograms, jaccard.cutoff = c(0,0.8))
mp.genes <- geneNMF.metaprograms$metaprograms.genes
saveRDS(mp.genes, file = "results/mp.genes.rds")
saveRDS(geneNMF.metaprograms, file = "results/geneNMF.metaprograms.rds")

sim_mat <- geneNMF.metaprograms$programs.similarity
module_names <- colnames(sim_mat)
mp_cluster <- geneNMF.metaprograms$programs.clusters[module_names]

patient <- sapply(strsplit(module_names, "\\."), `[`, 1)
patient <- sapply(strsplit(patient, "\\_"), `[`, 1)
Ki67_group <- plot_data$Cluster[match(patient, plot_data$Patient)]
Grade <- plot_data$Grade[match(patient, plot_data$Patient)]

annotation_col <- data.frame(
  Ki67_group = Ki67_group,
  Grade = Grade
)
rownames(annotation_col) <- module_names

valid_idx <- complete.cases(annotation_col)
annotation_col <- annotation_col[valid_idx, , drop = FALSE]
sim_mat <- sim_mat[valid_idx, valid_idx]
mp_cluster <- mp_cluster[valid_idx]
module_names <- module_names[valid_idx]

order_df <- data.frame(
  module = module_names,
  mp_cluster = mp_cluster,
  Ki67_group = annotation_col$Ki67_group
)

order_col <- order(order_df$mp_cluster, order_df$Ki67_group)
sim_mat_ordered <- sim_mat[order_col, order_col]
annotation_col_ordered <- annotation_col[order_col, , drop = FALSE]
mp_cluster_ordered <- mp_cluster[order_col]

annotation_colors <- list(
  Ki67_group = cluster_colors,
  Grade = grade_color 
)

n_colors <- 100
pal <- viridis(n_colors, option = "plasma")
breaks <- seq(0, max(sim_mat_ordered), length.out = n_colors)
colors <- c(rep("lightgrey", 8), pal[9:100])
col_fun <- colorRamp2(breaks = breaks, colors = colors)

mp_levels <- sort(unique(mp_cluster_ordered))
mp_colors <- setNames(brewer.pal(n = length(mp_levels), name = "Set3"), mp_levels)

col_ha <- HeatmapAnnotation(
  Ki67_group = annotation_col_ordered$Ki67_group,
  Grade = annotation_col_ordered$Grade,
  col = annotation_colors
)

row_ha <- rowAnnotation(
  MP_cluster = mp_cluster_ordered,
  col = list(MP_cluster = mp_colors),
  show_legend = TRUE,
  width = unit(0.5, "cm")
)

ht <- Heatmap(
  sim_mat_ordered,
  name = "Similarity",
  col = col_fun,
  top_annotation = col_ha,
  left_annotation = row_ha,
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_split = mp_cluster_ordered,
  use_raster = TRUE
)

pdf("figures/Fig4D.pdf", height = 5.5, width = 7)
draw(ht)
dev.off()

## Figure 4E-F: enrichment plot----------
mp_genes_entrez <- lapply(mp.genes, function(x) {
  ids <- mapIds(org.Hs.eg.db, keys = x, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  ids <- ids[!is.na(ids)]
  ids
})

cc <- compareCluster(geneCluster = mp_genes_entrez, fun = "enrichGO",
                     OrgDb = org.Hs.eg.db, ont = "BP",
                     pvalueCutoff = 0.05, qvalueCutoff = 0.1, readable = TRUE)

pdf("figures/Fig4E.pdf", height = 8, width = 8)
dotplot(cc, showCategory = 3) + 
  ggtitle("MP-wise GOBP enrichment")
dev.off()

merged_obj <- AddModuleScore_UCell(merged_obj, features = mp.genes, ncores=4, name = "")
VlnPlot(merged_obj, features=names(mp.genes), group.by = "Ki67_cluster",
  pt.size = 0, ncol=5)

cell_scores <- FetchData(merged_obj, vars = c("Patient", paste0("MP",1:10), "Ki67_cluster"))

sample_scores <- cell_scores %>%
  group_by(Patient, Ki67_cluster) %>%
  dplyr::summarise(across(starts_with("MP"), ~ mean(.x, na.rm=TRUE))) %>%
  ungroup()

long_scores <- sample_scores %>%
  pivot_longer(cols = starts_with("MP"),
               names_to = "MP",
               values_to = "score")
long_scores$MP <- factor(long_scores$MP, levels = paste0("MP",1:10))

subset_plot <- long_scores %>%
  filter(Ki67_cluster != "Stromal-dominant" & MP %in% paste0("MP", c(1, 2, 3, 5, 6, 7, 10))) %>% 
  mutate(score = as.numeric(as.character(score)))

groups <- as.vector(unique(subset_plot$Ki67_cluster))
comparisons_list <- combn(groups, 2, simplify = FALSE)

label_positions <- subset_plot %>%
  group_by(MP) %>%
  summarise(label_y = max(score, na.rm = TRUE) * 1.05)

p <- ggplot(subset_plot, aes(x = Ki67_cluster, y = score, fill = Ki67_cluster)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +
  facet_wrap(~ MP, scales = "free_y", ncol = 2) +
  labs(
    title = "MP activity by Ki67 cluster",
    x = "Ki67 grouping",
    y = "Mean MP score per patient"
  ) +
  scale_fill_manual(values = cluster_colors) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

for (i in seq_len(nrow(label_positions))) {
  p <- p + stat_compare_means(
    data = subset_plot %>% filter(MP == label_positions$MP[i]),
    aes(x = Ki67_cluster, y = score),
    comparisons = comparisons_list,
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = TRUE,
    label.y = label_positions$label_y[i] * c(1.05, 1.15, 1.25) # 三个比较错开
  )
}

pdf("figures/Fig4F.pdf", height = 8, width = 4.5)
print(p)
dev.off()
