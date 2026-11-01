## Figure5
## Yijun Dong, 22/10/2025

## library-----
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
library(patchwork)
library(TCGAbiolinks)

## prepare data------
merged_obj <- readRDS("results/merged_obj_new.rds")

lgg_subtype <- TCGAquery_subtype(tumor = "lgg")
gbm_subtype <- TCGAquery_subtype(tumor = "gbm")
load("data/GBM-LGG-TCGA-exprSet_names_by_symbol.Rdata")

fpkm_data1 <- read_table("~/neurosurgery/NES/reference/CGGA/CGGA.mRNAseq_693.RSEM-genes.20200506.txt")
fpkm_data2 <- read_table("~/neurosurgery/NES/reference/CGGA/CGGA.mRNAseq_325.RSEM-genes.20200506.txt")

clinical_data1 <- read.table("~/neurosurgery/NES/reference/CGGA/CGGA.mRNAseq_693_clinical.20200506.txt", sep = "\t", header = TRUE)
clinical_data2 <- read.table("~/neurosurgery/NES/reference/CGGA/CGGA.mRNAseq_325_clinical.20200506.txt", sep = "\t", header = TRUE)

gene <- "MKI67"

## Figure 5A: Pan-glioma TCGA (GBM + LGG)-------------
barcode <- colnames(exprSet)
first_three <- unlist(lapply(strsplit(barcode, "-", fixed = TRUE),
                             function(x) paste(x[1:3], collapse = "-")))
colnames(exprSet) <- first_three

gbm_subtype$tumor_grade <- "GBM"
lgg_subtype$tumor_grade <- "LGG"

gene_expr_df <- data.frame(
  patient = colnames(exprSet),
  gene_expr = as.numeric(exprSet[gene, ])
)

meta_combined <- bind_rows(gbm_subtype, lgg_subtype) %>%
  left_join(gene_expr_df, by = "patient") %>%
  filter(!is.na(gene_expr)) %>%
  mutate(group = case_when(
    gene_expr < quantile(gene_expr, 0.5, na.rm = TRUE) ~ "Low",
    gene_expr >= quantile(gene_expr, 0.5, na.rm = TRUE) ~ "High",
    TRUE ~ NA_character_
  )) %>%
  mutate(
    OS_time_months = Survival..months.,
    OS_status = Vital.status..1.dead.
  )

head(meta_combined)

fit_tcga_pan <- survfit(Surv(OS_time_months, OS_status) ~ group, data = meta_combined)

pdf("figures/Fig5A.pdf", height = 5.5, width = 5)
ggsurvplot(
  fit_tcga_pan, data = meta_combined,
  pval = TRUE, risk.table = TRUE,
  censor.size = 3,
  surv.median.line = "hv",
  legend.labs = c("High", "Low"),
  xlab = "Time (Months)", ylab = "Overall survival probability",
  title = "TCGA Pan-Glioma (GBM + LGG)",
  palette = c("#8B235A", "#005787"),
  ggtheme = theme_classic()
)
dev.off()

## Figure 5B: Pan-glioma CGGA (GBM + LGG)-------------
identical_genes <- intersect(fpkm_data1$Gene_Name, fpkm_data2$Gene_Name)
fpkm_data1 <- fpkm_data1 %>% column_to_rownames("Gene_Name")
fpkm_data2 <- fpkm_data2 %>% column_to_rownames("Gene_Name")
fpkm_data <- cbind(fpkm_data1[identical_genes, ], fpkm_data2[identical_genes, ])

clinical_data <- rbind(clinical_data1, clinical_data2) %>%
  mutate(gene_expr = unlist(fpkm_data[gene, CGGA_ID])) %>%
  mutate(group = case_when(
    gene_expr < quantile(gene_expr, 0.5, na.rm = TRUE) ~ "Low",
    gene_expr >= quantile(gene_expr, 0.5, na.rm = TRUE) ~ "High",
    TRUE ~ NA_character_
  )) %>%
  mutate(OS_months = OS / 30.4375) # Convert to months

fit_cgga_pan <- survfit(Surv(OS_months, Censor..alive.0..dead.1.) ~ group, data = clinical_data)

pdf("figures/Fig5B.pdf", height = 5.5, width = 5)
ggsurvplot(
  fit_cgga_pan, data = clinical_data,
  pval = TRUE, risk.table = TRUE,
  censor.size = 3,
  surv.median.line = "hv",
  legend.labs = c("High", "Low"),
  xlab = "Time (Months)", ylab = "Overall survival probability",
  title = "CGGA Pan-Glioma (GBM + LGG)",
  palette = c("#FDDE00", "#1068B1"),
  ggtheme = theme_classic()
)
dev.off()

## Figure 5C: TCGA IDHmt Subgroup-----
tcga_lgg <- meta_combined %>% filter(c)

fit_tcga_lgg <- survfit(Surv(OS_time_months, OS_status) ~ group, data = tcga_lgg)

pdf("figures/Fig5C.pdf", height = 5.5, width = 5)
ggsurvplot(
  fit_tcga_lgg, data = tcga_lgg,
  pval = TRUE, risk.table = TRUE,
  censor.size = 3,
  surv.median.line = "hv",
  legend.labs = c("High", "Low"),
  xlab = "Time (Months)", ylab = "Overall survival probability",
  title = "TCGA Glioma (IDHmt Only)",
  palette = c("#8B235A", "#005787"),
  ggtheme = theme_classic()
)
dev.off()

## Figure 5D: TCGA GBM Subgroup-------
tcga_gbm <- meta_combined %>% filter(IDH.status == "WT")

fit_tcga_gbm <- survfit(Surv(OS_time_months, OS_status) ~ group, data = tcga_gbm)

pdf("figures/Fig5D.pdf", height = 5.5, width = 5)
ggsurvplot(
  fit_tcga_gbm, data = tcga_gbm,
  pval = TRUE, risk.table = TRUE,
  censor.size = 3,
  surv.median.line = "hv",
  legend.labs = c("High", "Low"),
  xlab = "Time (Months)", ylab = "Overall survival probability",
  title = "TCGA Glioma (IDHwt Only)",
  palette = c("#8B235A", "#005787"),
  ggtheme = theme_classic()
)
dev.off()

fit_cox <- coxph(Surv(OS_time_months, OS_status) ~ group, data = tcga_gbm)
summary(fit_cox)

## Figure 5E: CGGA LGG Subgroup----
cgga_lgg <- clinical_data %>% 
  filter(IDH_mutation_status == "Mutant")

fit_cgga_lgg <- survfit(Surv(OS_months, Censor..alive.0..dead.1.) ~ group, data = cgga_lgg)

pdf("figures/Fig5E.pdf", height = 5.5, width = 5)
ggsurvplot(
  fit_cgga_lgg, data = cgga_lgg, 
  pval = TRUE, risk.table = TRUE,
  censor.size = 3,
  surv.median.line = "hv",
  legend.labs = c("High", "Low"),
  xlab = "Time (Months)", ylab = "Overall survival probability",
  title = "CGGA Glioma (IDHmt Only)", 
  palette = c("#FDDE00", "#1068B1"),
  ggtheme = theme_classic()
)
dev.off()

fit_cox <- coxph(Surv(OS_months, Censor..alive.0..dead.1.) ~ group, data = cgga_lgg)
summary(fit_cox)

## Figure 5F: CGGA GBM Subgroup----
cgga_gbm <- clinical_data %>% 
  filter(IDH_mutation_status == "Wildtype")

fit_cgga_gbm <- survfit(Surv(OS_months, Censor..alive.0..dead.1.) ~ group, data = cgga_gbm)

pdf("figures/Fig5F.pdf", height = 5.5, width = 5)
ggsurvplot(
  fit_cgga_gbm, data = cgga_gbm, 
  pval = TRUE, risk.table = TRUE,
  censor.size = 3,
  surv.median.line = "hv",
  legend.labs = c("High", "Low"),
  xlab = "Time (Months)", ylab = "Overall survival probability",
  title = "CGGA Glioma (IDHwt Only)", 
  palette = c("#FDDE00", "#1068B1"),
  ggtheme = theme_classic()
)
dev.off()

print(fit_cgga_gbm)

fit_cox <- coxph(Surv(OS_months, Censor..alive.0..dead.1.) ~ group, data = cgga_gbm)
summary(fit_cox)

## Figure 5G: Multivariate Cox Regression Forest Plot----
cox_data <- meta_combined %>%
  select(OS_time_months, OS_status, group, Age..years.at.diagnosis., IDH.status, Grade, MGMT.promoter.status, Histology) %>% 
  na.omit() %>%
  mutate(
    Age = Age..years.at.diagnosis.,
    Grade = factor(Grade, levels = c("G2","G3","G4")),
    MGMT_promotor_status = factor(MGMT.promoter.status, levels = c("Methylated","Unmethylated")),
    MKI67_group = factor(group, levels = c("Low", "High")),
    IDH_status = factor(IDH.status, levels = c("Mutant", "WT"))
  )

res.cox <- coxph(Surv(OS_time_months, OS_status) ~ MKI67_group + Age + Grade + IDH_status + MGMT_promotor_status
                 , data = cox_data)

pdf("figures/Fig5G.pdf", height = 4, width = 10)
ggforest(res.cox, data = cox_data, 
         main = "Multivariate Cox Regression (TCGA)",
         cpositions = c(0.02, 0.22, 0.4), 
         fontsize = 1.0, 
         refLabel = "reference", noDigits = 4)
dev.off()

## Figure 5H: CIBERSORTx (Tumor Lineage)-----
fpkm_data_TCGA <- cbind(Gene = rownames(exprSet), exprSet)
fpkm_data_TCGA <- fpkm_data_TCGA[!duplicated(rownames(fpkm_data_TCGA)),] %>% na.omit()

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
    t() %>% as.data.frame() %>%
    dplyr::mutate(sample_id = gsub("\\.", "-", rownames(.))) %>%
    dplyr::select(sample_id, MKI67)
  
  merged_data <- merge(meta_combined, exp_data, by = "sample_id")
  mki67_values <- merged_data$MKI67[!is.na(merged_data$MKI67) & merged_data$MKI67 != 1]
  
  if (length(mki67_values) == 0) {
    message("Skipping lineage ", lineage_name, ": No valid MKI67 expression")
    return(NULL)
  }
  
  merged_data <- merged_data %>%
    mutate(group = case_when(
      MKI67 < quantile(MKI67, 0.5, na.rm = TRUE) ~ "Low",
      MKI67 >= quantile(MKI67, 0.5, na.rm = TRUE) ~ "High",
      TRUE ~ NA_character_ 
    ))
  
  fit <- survfit(Surv(OS_time_months, OS_status) ~ group, data = merged_data)
  km_plot <- ggsurvplot(
    fit, data = merged_data,
    pval = TRUE, risk.table = TRUE,
    surv.median.line = "hv",
    xlab = "Time (Months)", ylab = "Overall survival probability",
    title = paste("Survival by MKI67 in", lineage_name, "Cells"),
    legend.title = "MKI67",
    palette = c("#D55E00","#0072B2"),
    ggtheme = theme_classic()
  )
  
  cox_fit <- coxph(Surv(OS_time_months, OS_status) ~ group, data = merged_data)
  cox_summary <- summary(cox_fit)
  
  return(list(
    lineage = lineage_name,
    km_plot = km_plot,
    cox_hr = exp(cox_summary$coefficients[1, 1]), 
    cox_lcl = exp(cox_summary$conf.int[1, 3]),     
    cox_ucl = exp(cox_summary$conf.int[1, 4]),     
    cox_p = cox_summary$coefficients[1, 5]         
  ))
}

results <- map(lineage_files, safely(analyze_lineage)) %>% 
  set_names(basename(lineage_files)) %>%
  transpose()
valid_results <- compact(results$result)
names(valid_results) <- map_chr(valid_results, ~ .x$lineage)

pdf("figures/Fig5H.pdf", height = 5.5, width = 5)
print(valid_results$Tumor$km_plot)
dev.off()

## Figure 2D-2: Clinical Correlation Boxplots----
grade_color <- c("Mutant"="navy","WT"="firebrick3" )
grade_levels <- levels(factor(meta_combined$IDH.status))
grade_comparisons <- combn(grade_levels, 2, simplify = FALSE)

p_grade <- ggboxplot(meta_combined, x = "IDH.status", y = "gene_expr", fill = "IDH.status",
                     palette = grade_color, 
                     alpha = 0.6, 
                     size = 0.3,     
                     add = "jitter", 
                     add.params = list(size = 1.2, alpha = 0.5)) +
  stat_compare_means(comparisons = grade_comparisons,
                     method = "wilcox.test",
                     label = "p.signif", 
                     hide.ns = TRUE,
                     step.increase = 0.05) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_classic(base_size = 14) +
  labs(y = paste(gene, "Expression"), x = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

pdf("figures/Fig2D-2.pdf",height = 4, width = 5)
print(p_grade)
dev.off()
