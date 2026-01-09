## Figure5
## Yijun Dong, 22/10/2025
##library
library(survival)
library(survminer)
library(tidyverse)
library(GSVA)
library(Matrix)
library(Seurat)
library(readr)
library(ggpubr)
library(reshape2)
library(stats)
library(GSVA)

## prepare data------
merged_obj <- readRDS("results/merged_obj_new.rds")
load("data/GBM-LGG-TCGA-exprSet_names_by_symbol.Rdata")
load("data/meta-TCGA-GBM-LGG.rdata")

fpkm_data1 <- read_table("~/neurosurgery/NES/reference/CGGA/CGGA.mRNAseq_693.RSEM-genes.20200506.txt")
fpkm_data2 <- read_table("~/neurosurgery/NES/reference/CGGA/CGGA.mRNAseq_325.RSEM-genes.20200506.txt")

clinical_data1 <- read.table("~/neurosurgery/NES/reference/CGGA/CGGA.mRNAseq_693_clinical.20200506.txt", sep = "\t", header = TRUE)
clinical_data2 <- read.table("~/neurosurgery/NES/reference/CGGA/CGGA.mRNAseq_325_clinical.20200506.txt", sep = "\t", header = TRUE)

gene <- "MKI67"

## Figure 5A: Survival Analysis with TGGA-------------
barcode <- colnames(exprSet)
first_three <- unlist(lapply(strsplit(barcode, "-", fixed = TRUE),
                             function(x) paste(x[1:3], collapse = "-")))
colnames(exprSet) <- first_three

meta_combined <- bind_rows(meta[[1]], meta[[2]])

meta_combined <- meta_combined %>%
  mutate(gene_expr = exprSet[gene, rownames(meta_combined)])

meta_combined <- meta_combined %>%
  mutate(group = case_when(
    gene_expr <= quantile(gene_expr, 0.5, na.rm = TRUE) ~ "Low",
    gene_expr >= quantile(gene_expr, 0.5, na.rm = TRUE) ~ "High",
    TRUE ~ NA_character_
  ))

meta_combined <- meta_combined %>%
  mutate(OS_time = ifelse(!is.na(days_to_death),
                          days_to_death,
                          days_to_last_follow_up),
         OS_status = ifelse(vital_status == "Dead", 1, 0))

surv_obj <- Surv(time = meta_combined$OS_time, event = meta_combined$OS_status)
fit <- survfit(surv_obj ~ group, data = meta_combined)

pdf("figures/Fig5A.pdf", height = 5, width = 4.5)
ggsurvplot(
  fit, data = meta_combined,
  pval = TRUE,
  risk.table = TRUE,
  legend.labs = c("High", "Low"),
  xlab = "Days",
  ylab = "Overall survival probability",
  title = "Survival Analysis with MKI67 in TCGA Glioma (GBM + LGG)",
  palette = c("#8B235A","#005787")
)
dev.off()

## Figure 5B: Survival Analysis with CGGA-------------
identical_genes <- intersect(fpkm_data1$Gene_Name, fpkm_data2$Gene_Name)
fpkm_data1 <- fpkm_data1 %>% column_to_rownames("Gene_Name")
fpkm_data2 <- fpkm_data2 %>% column_to_rownames("Gene_Name")
fpkm_data <- cbind(fpkm_data1[identical_genes, ], fpkm_data2[identical_genes, ])

clinical_data <- rbind(clinical_data1, clinical_data2)

clinical_data <- clinical_data %>%
  mutate(gene_expr = unlist(fpkm_data[gene, clinical_data$CGGA_ID]))

clinical_data <- clinical_data %>%
  mutate(group = case_when(
    gene_expr <= quantile(gene_expr, 0.5, na.rm = TRUE) ~ "Low",
    gene_expr >= quantile(gene_expr, 0.5, na.rm = TRUE) ~ "High",
    TRUE ~ NA_character_
  ))

surv_obj <- Surv(time = clinical_data$OS, event = clinical_data$Censor..alive.0..dead.1.)
fit <- survfit(surv_obj ~ group, data = clinical_data)

pdf("figures/Fig5B.pdf", height = 5, width = 4.5)
ggsurvplot(
  fit, data = clinical_data,
  pval = TRUE,
  risk.table = TRUE,
  legend.labs = c("High", "Low"),
  xlab = "Days",
  ylab = "Overall survival probability",
  title = "Survival Analysis with MKI67 in CGGA Glioma (GBM + LGG)",
  palette = c("#FDDE00", "#1068B1"))
dev.off()

## Figure 5C: cibersortx-----
fpkm_data_TCGA <- cbind(Gene = rownames(exprSet), exprSet)
fpkm_data_TCGA <- fpkm_data_TCGA[!duplicated(rownames(fpkm_data_TCGA)),] %>% 
  na.omit()
dim(fpkm_data_TCGA)

write.table(fpkm_data_TCGA, "results/cybersort/bulk_expression.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

merged_obj$lineage_correct <- merged_obj$lineage
merged_obj$lineage_correct[merged_obj$lineage %in% c("Lymphoid","Myeloid")] <- "Immune"
lineages <- unique(merged_obj@meta.data$lineage_correct)

reference_matrix <- GetAssayData(merged_obj, assay = "RNA", slot = "data")
colnames(reference_matrix) <- merged_obj$lineage_correct

selected_cells <- sapply(unique(colnames(reference_matrix)), function(x) {
  cell_indices <- which(colnames(reference_matrix) == x)
  sample_size <- 1000
  as.vector(sample(cell_indices, sample_size))
})

gene_used <- intersect(rownames(fpkm_data_TCGA), rownames(reference_matrix))
reference_matrix <- reference_matrix[gene_used, as.vector(selected_cells)]
reference_matrix <- as.matrix(reference_matrix)
reference_matrix <- cbind(GeneSymbol = rownames(reference_matrix), reference_matrix)
dim(reference_matrix)

write.table(reference_matrix, "results/cybersort/single_cell_reference_matrix.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

set.seed(123)
selected_cells <- lapply(lineages, function(lin) {
  cells_in_lineage <- colnames(merged_obj)[merged_obj@meta.data$lineage_correct == lin]
  if (length(cells_in_lineage) > 2000) {
    sample(cells_in_lineage, 2000) 
  } else {
    cells_in_lineage 
  }
}) %>% unlist()

merged_obj_subset <- merged_obj[, selected_cells]

Idents(merged_obj_subset) <- "lineage_correct"
markers <- FindAllMarkers(
  merged_obj_subset,
  only.pos = TRUE,     
  min.pct = 0.25,    
  logfc.threshold = 0.25,
  return.thresh = 0.05  
)

top_markers <- markers %>%
  filter(gene %in% gene_used & p_val_adj < 0.01) %>% 
  group_by(cluster) %>%            
  arrange(desc(avg_log2FC)) %>%   
  slice_head(n = 100) %>% 
  ungroup()

genes_interest <- c("MKI67")
gene_subset <- unique(c(top_markers$gene, genes_interest))

write.table(
  gene_subset,                      
  file = "results/cybersort/cybersort_gene_subset.txt", 
  quote = FALSE,                       
  row.names = FALSE,                 
  col.names = FALSE                  
)

## after using cibersortx website to generate the results
lineage_files <- list.files("results/cybersort", pattern = "CIBERSORTxHiRes_Job26_.*_Window16.txt", full.names = TRUE)

analyze_lineage <- function(file_path) {
  lineage_name <- gsub("CIBERSORTxHiRes_Job26_|_Window16.txt", "", basename(file_path))
  exp_data <- read.delim(file_path, row.names = 1) %>% 
    t() %>% 
    as.data.frame() %>%
    mutate(sample_id = gsub("\\.", "-", rownames(.))) %>%
    select(sample_id, MKI67)
  merged_data <- merge(meta_combined, exp_data, by = "sample_id")
  mki67_values <- merged_data$MKI67[!is.na(merged_data$MKI67) & merged_data$MKI67 != 1]
  cutoff <- median(mki67_values)
  merged_data <- merged_data %>%
    mutate(group = case_when(
      MKI67 <= quantile(MKI67, 0.5, na.rm = TRUE) ~ "Low",
      MKI67 >= quantile(MKI67, 0.5, na.rm = TRUE) ~ "High",
      TRUE ~ NA_character_ 
    ))
  
  mki67_values <- merged_data$MKI67[!is.na(merged_data$MKI67) & merged_data$MKI67 != 1]
  if (length(mki67_values) == 0) {
    message("Skipping lineage ", lineage_name, ": No valid MKI67 expression")
    return(NULL)
  }
  
  fit <- survfit(Surv(time, event) ~ group, data = merged_data)
  km_plot <- ggsurvplot(
    fit,
    data = merged_data,
    pval = TRUE,
    risk.table = TRUE,
    title = paste("Survival by MKI67 in", lineage_name, "Cells"),
    legend.title = "MKI67",
    palette = c("#D55E00","#0072B2")
  )
  
  cox_fit <- coxph(Surv(time, event) ~ group, data = merged_data)
  cox_summary <- summary(cox_fit)
  
  return(list(
    lineage = lineage_name,
    km_plot = km_plot,
    cox_hr = exp(cox_summary$coefficients[1, 1]),  # HR
    cox_lcl = exp(cox_summary$conf.int[1, 3]),     # 95% CI lower
    cox_ucl = exp(cox_summary$conf.int[1, 4]),     # 95% CI upper
    cox_p = cox_summary$coefficients[1, 5]         # p-value
  ))
}

results <- map(lineage_files, safely(analyze_lineage)) %>% 
  set_names(basename(lineage_files)) %>%
  transpose()
valid_results <- compact(results$result)
names(valid_results) <- map_chr(valid_results, ~ .x$lineage)

pdf("figures/Fig5C.pdf", height = 5, width = 4.5)
valid_results$Tumor$km_plot
dev.off()

## Figure 5D: ecotypes classification of TCGA patients-------
bulk_expr <- as.matrix(exprSet)

merged_obj_subset <- subset(merged_obj, Ki67_cluster!="Stromal-dominant")
pseudo_bulk <- AverageExpression(merged_obj_subset, group.by = "Ki67_cluster", assays = "RNA", slot = "data")
pseudo_bulk <- pseudo_bulk$RNA

common_genes <- intersect(rownames(pseudo_bulk), rownames(exprSet))
pb <- pseudo_bulk[common_genes, ]
bx <- exprSet[common_genes, ]

pb_z <- t(scale(t(pb)))
bx_z <- t(scale(t(bx)))

sample_class_pearson <- apply(bx_z, 2, function(x) {
  cors <- sapply(colnames(pb_z), function(cl) cor(x, pb_z[, cl], method = "pearson", use = "pairwise.complete.obs"))
  names(which.max(cors))
})

table(sample_class_pearson)

tcga_meta <- meta_combined
tcga_meta$sample_id <- rownames(tcga_meta)

pred_df <- data.frame(
  sample_id = names(sample_class_pearson),
  Ki67_like_cluster = sample_class_pearson
)

tcga_meta_pred <- tcga_meta %>%
  left_join(pred_df, by = "sample_id")

table(tcga_meta_pred$Ki67_like_cluster)

tcga_meta_pred <- tcga_meta_pred %>%
  mutate(
    OS_time = ifelse(!is.na(days_to_death), days_to_death, days_to_last_follow_up),
    OS_status = ifelse(vital_status == "Dead", 1, 0)
  )

tcga_meta_pred$Ki67_like_cluster <- factor(
  tcga_meta_pred$Ki67_like_cluster,
  levels = c("Tumor-dominant", "Immune-dominant", "Glial-dominant")
)

surv_obj <- Surv(time = tcga_meta_pred$OS_time, event = tcga_meta_pred$OS_status)
fit <- survfit(surv_obj ~ Ki67_like_cluster, data = tcga_meta_pred)

pdf("figures/Fig5D.pdf", height = 5.5, width = 5)
ggsurvplot(
  fit,
  data = tcga_meta_pred,
  pval = TRUE,
  risk.table = TRUE,
  legend.title = "Predicted Cluster",
  legend.labs = c("Tumor-dominant", "Immune-dominant", "Glial-dominant"),
  palette = cluster_colors,
  xlab = "Days",
  ylab = "Overall survival probability",
  title = "Survival by Predicted Ki67-like Cluster in TCGA Glioma"
)
dev.off()
