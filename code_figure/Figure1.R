## Figure1
## Yijun Dong, 15/10/2025
## load library------
library(dplyr)
library(mixtools)
library(Seurat)
library(ArchR)
library(ggsci)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(readxl)
library(ggpubr)

## prepare data------
merged_obj <- readRDS("results/merged_obj_new.rds")
infercnv <- readRDS("results/infercnv_out/run.final.infercnv.obj")
gene <- "MKI67"

## prepare color------
grade_color <- c("LGG (II-III)" = "navy", "GBM (IV)" = "firebrick3")

lineage_colors <- c(
  "Tumor" = "#D94F4F",
  "Myeloid" = "#0072B5",
  "Lymphoid" = "#00A087",
  "Stromal" = "#EFC000",
  "Glial" = "#AD8BC9"
)

celltypes_to_label <- c(
  "AC", "MES", "NPC", "OPC",
  "monocyte-derived DCs", "myeloid-derived suppressor cells", "macrophages",
  "microglia", "T cells", "B cells", "NK cells", "fibroblasts",
  "endothelial cells", "pericytes", "oligodendrocytes"
)

merged_obj$celltypes_minor <- factor(merged_obj$celltypes_minor, levels = celltypes_to_label)
celltypes_color <- ArchRPalettes$stallion[1:length(celltypes_to_label)]
names(celltypes_color) <- celltypes_to_label

## Figure 1B-C, Supplementary Figure 1A: UMAP, lineage, grade, cell cycle phase, and cell types---------
gene <- "MKI67"
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

pdf("figures/Fig1B.pdf", height = 5, width = 6)
DimPlot(merged_obj, group.by = "lineage", cols = lineage_colors)
dev.off()

lineage_counts <- merged_obj@meta.data %>%
  group_by(lineage) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(desc(count))

p_lineage <- ggplot(lineage_counts, aes(x = reorder(lineage, count), y = count, fill = lineage)) +
  geom_col(width = 0.7) +
  theme_classic(base_size = 14) +
  labs(x = "Cell count", y = "Lineage", fill = "Lineage") +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))  +
  coord_flip() +
  scale_fill_manual(values = lineage_colors)

pdf("figures/Fig1B-2.pdf", height = 4, width = 4)
p_lineage
dev.off()

pdf("figures/Fig1C.pdf", height = 5, width = 6.5)
DimPlot(merged_obj, group.by = "Diagnosis") + scale_color_nejm()
dev.off()

pdf("figures/Sup_Fig1A.pdf", height = 5, width = 6)
DimPlot(merged_obj, group.by = "Phase_enhanced")
dev.off()

## Figure 1D: composition of lineage per sample-------
df <- merged_obj@meta.data

composition_df <- df %>%
  group_by(Patient, lineage, Diagnosis) %>%
  summarise(cell_count = n(), .groups = "drop")

sample_stats <- composition_df %>%
  group_by(Patient, Diagnosis) %>%
  summarise(
    total_cells = sum(cell_count),
    tumor_cells = sum(cell_count[lineage == "Tumor"]),
    tumor_ratio = tumor_cells / total_cells,
    .groups = "drop"
  )

sample_order <- sample_stats %>% arrange(desc(tumor_ratio)) %>% pull(Patient)
composition_df$Patient <- factor(composition_df$Patient, levels = sample_order)
sample_stats$Patient <- factor(sample_stats$Patient, levels = sample_order)

p_main <- ggplot(composition_df, aes(x = Patient, y = cell_count, fill = lineage)) +
  geom_bar(stat = "identity", position = "fill", width = 0.9) +
  scale_fill_manual(values = lineage_colors) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_classic() +
  labs(x = NULL, y = "Fraction of cells (%)", fill = "Lineage") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "right",
    plot.margin = margin(0, 5.5, 0, 5.5)
  )

p_diag <- ggplot(sample_stats, aes(x = Patient, y = 1, fill = Diagnosis)) +
  geom_tile() +
  scale_fill_nejm() +
  theme_void() +
  theme(legend.position = "bottom")

p_combined <- p_main / p_diag + plot_layout(heights = c(10, 1))

pdf("figures/Fig1D.pdf", width = 10, height = 5)
p_combined
dev.off()

## Supplementary Figure 1A-B: composition of overall lineage and cell types------
meta_df <- df
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

htmlwidgets::saveWidget(plots[["LGG (II-III)"]], "figures/Sup_Fig1A.html")
htmlwidgets::saveWidget(plots[["GBM (IV)"]], "figures/Sup_Fig1B.html")

## Supplementary Figure 1C: infercnv identify tumor cells----------
counts <- GetAssayData(merged_obj_subset, assay = "RNA", slot = "counts")  
if (!inherits(counts, "dgCMatrix")) {
  counts <- as(counts, "dgCMatrix")
}

annotations <- merged_obj_subset@meta.data %>%
  dplyr::select(lineage) %>%
  rownames_to_column("cell")

write.table(annotations, file = "results/cnv_annotations.txt",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

gene_order_file <- "~/neurosurgery/NES/reference/geneLocate.txt"

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = counts,
  annotations_file = "results/cnv_annotations.txt",
  delim = "\t",
  gene_order_file = gene_order_file,
  ref_group_names = c("Myeloid")
)

infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1,
  out_dir = "results/infercnv_out",
  cluster_by_groups = TRUE,
  denoise = TRUE,
  HMM = TRUE,
  resume_mode = T,
  output_format = 'pdf',
  num_threads = 10
)

## Figure 1E-F: umap and heatmap of important markers of cell types-----------
pdf("figures/Fig1F.pdf", height = 5, width = 8)
DimPlot(merged_obj, group.by = "celltypes_minor", cols = celltypes_color) +
  theme(legend.position = "right")
dev.off()

markers <- read_excel("~/neurosurgery/NES/reference/mmc2.xlsx", sheet = 1, skip = 3)
markers<-lapply(markers,function(x) na.omit(x))
markers<-markers[!names(markers)%in%c("G1.S","G2.M")]
gene_sets<-list(MES = c(markers$MES1,markers$MES2),
              NPC = c(markers$NPC1,markers$NPC2),
              OPC = markers$OPC,
              AC = markers$AC)

markers_celltype <- FindAllMarkers(
  merged_obj,
  group.by = "celltypes_minor",
  only.pos = TRUE, 
  min.pct = 0.25,  
  logfc.threshold = 0.25
)

markers_celltype_filter <- markers_celltype %>%
  filter(avg_log2FC > 1 & p_val_adj < 0.01) %>% 
  group_by(cluster) %>% 
  top_n(n = 100, wt = avg_log2FC)
dim(markers_celltype_filter)

canonical_markers <- list(
  `monocyte-derived DCs` = c("CD1C", "FCER1A", "CLEC10A", "LAMP3", "CCR7"),
  `myeloid-derived suppressor cells` = c("S100A8", "S100A9", "FCGR3A", "CD33", "IL1B"),
  macrophages = c("CD68", "CD163", "MRC1", "CSF1R", "LYZ"),
  microglia = c("CX3CR1", "P2RY12", "TMEM119", "GPR34", "C1QA"),
  `T cells` = c("CD3D", "CD3E", "CD2", "TRAC", "IL7R"),
  `B cells` = c("MS4A1", "CD79A", "CD19", "CD79B", "BANK1"),
  `NK cells` = c("NKG7", "GNLY", "KLRD1", "KLRB1", "PRF1"),
  fibroblasts = c("COL1A1", "COL1A2", "FAP", "PDGFRA", "DCN"),
  endothelial = c("CLDN5", "PECAM1", "VWF", "CDH5", "ESAM"),
  pericytes = c("RGS5", "PDGFRB", "ACTA2", "MKI67", "MCAM"),
  oligodendrocytes = c("MBP", "MOG", "OLIG1", "OLIG2", "PLP1"),
  AC = gene_sets$AC,    # assume you have gene_sets list for tumor subtypes
  MES = gene_sets$MES,
  NPC = gene_sets$NPC,
  OPC = gene_sets$OPC
)

top_markers_for_label <- list()
for(ct in celltypes_to_label){
  ct_markers <- markers_celltype_filter %>%
    dplyr::filter(cluster == ct) %>%
    dplyr::filter(gene %in% canonical_markers[[ct]]) %>%
    dplyr::arrange(desc(avg_log2FC)) %>%
    pull(gene)
  top_markers_for_label[[ct]] <- head(ct_markers, 2)
}

top_markers_for_label

select_cells <- sample(colnames(merged_obj), 10000)
merged_obj_subset <- merged_obj[,select_cells]

plot.data <- GetAssayData(merged_obj_subset, slot = "data")[markers_celltype_filter$gene, ]
plot.data <- plot.data - rowMeans(plot.data)
plot.data=na.omit(plot.data)
CellInfoS=merged_obj_subset@meta.data

column_annot <-CellInfoS[,c("celltypes_minor","Diagnosis"),drop=F]
column_annot$Diagnosis = as.factor(as.character(column_annot$Diagnosis))
column_annot=with(column_annot, column_annot[order(Diagnosis), , drop=F])
column_annot=with(column_annot, column_annot[order(celltypes_minor), , drop=F])
plot.data<-plot.data[,row.names(column_annot)]

column.colors=list()
column.colors[["Diagnosis"]]<-c("Glioblastoma" = "#0273B5","Oligodendroglioma" = "#E18728","Astrocytoma" = "#BD3D29")
column.colors[["celltypes_minor"]]<-celltypes_color
Diagnosis=as.matrix(column_annot[,c("Diagnosis"),drop=F])
celltypes_minor=as.matrix(column_annot[,c("celltypes_minor"),drop=F])

colanno <- columnAnnotation(df=column_annot,
                            show_annotation_name =T,show_legend = F,col=column.colors)
genes= unlist(top_markers_for_label)
rows=rowAnnotation(sel = anno_mark(at = match(genes,row.names(plot.data)), labels = genes,labels_gp =gpar(col = "black", fontsize = 9)))

col<- circlize::colorRamp2(breaks = c(-2, 0, 2), colors = c("#007dd1", "white", "#ab3000"))
HM=Heatmap(name="logcounts",as.matrix(plot.data),cluster_rows =F,cluster_columns = F,top_annotation = colanno, right_annotation = rows,row_names_gp = gpar(fontsize=5),
           col = col,show_column_names= F,show_row_names=F,border = F,show_heatmap_legend = F,use_raster = T)        
lgd=Legend(title = "logcounts", at=  c(-2,0, 2),col_fun = col)
lgd1=Legend(labels = levels(as.factor(column_annot$celltypes_minor)),title="celltypes_minor",legend_gp = gpar(fill=column.colors[["celltypes_minor"]][levels(as.factor(column_annot$celltypes_minor))],fontsize=5))
lgd2=Legend(labels = levels(as.factor(column_annot$Diagnosis)),title="Diagnosis",legend_gp = gpar(fill=column.colors[["Diagnosis"]][levels(as.factor(column_annot$Diagnosis))],fontsize=5))

pdf("figures/Fig1F.pdf",width=7,height=6)
draw(HM,heatmap_legend_list = list(lgd,lgd1, lgd2), heatmap_legend_side = "right")
dev.off()
