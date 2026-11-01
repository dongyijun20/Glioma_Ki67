# Cellular Origins of Ki-67+ Cells in Gliomas

[![R-version](https://img.shields.io/badge/R-%3E%3D%204.3.0-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the analysis workflow and reproducible pipeline for the study:

> **Cellular Origins of Ki-67+ Cells in Gliomas**  
> *Yijun Dong, Kaishun Xia, Xiaopeng Guo, Wenbin Ma, Yu Wang\**  
> Department of Neurosurgery, Peking Union Medical College Hospital, Chinese Academy of Medical Sciences and Peking Union Medical College.

---

## Overview & Study Architecture

Ki-67 labeling index (encoded by *MKI67*) is the standard histological metric for tumor proliferation across human gliomas, yet bulk proliferation indices confound malignant neoplastic expansion with non-neoplastic stromal and immune cell proliferation.

In this study, we construct an integrated single-cell transcriptomic atlas spanning **40 glioma sample IDs representing 38 unique patient donors** across diverse histological and molecular subtypes (24 IDH-wildtype glioblastomas, 14 IDH-mutant lower-grade gliomas: 6 astrocytomas and 8 oligodendrogliomas; multi-region specimens from donor LGG-04 consolidated).

### Key Findings:
1. **Proliferative Heterogeneity**: Cycling cells identified through composite transcriptomic gating (*MKI67* raw transcript count >= 1 **OR** enhanced S/G2M cell-cycle phase labels; **23,029** union-positive cells, **8.97%** of **256,627** cells) reveal marked lineage shifts across molecular grades. While neoplastic cells drive proliferation in IDH-wildtype glioblastomas, microglial and myeloid cells contribute substantial fractions to the cycling compartment in IDH-mutant gliomas.
2. **Three Proliferative Ecotypes**: Ternary barycentric coordinate projection establishes three recurrent dominant-lineage ecotypes at the donor level:
   - **Tumor-dominant** (25/38 donors)
   - **Immune-dominant** (7/38 donors)
   - **Glial-dominant** (5/38 donors)  
   *(One specimen with prominent stromal proliferation, SF11488, is characterized as an exploratory outlier observation — not a fourth peer ecotype.)*
3. **Metaprogram Discovery**: Cross-dataset non-negative matrix factorization (GeneNMF, *k* = 10) uncovers recurring transcriptional programs underlying cycling glioma cells.
4. **Independent Prognostic Validation**: Bulk transcriptomic survival analyses across TCGA pan-glioma (*n* = 592 complete cases) and CGGA validation cohorts confirm *MKI67* as an independent predictor of overall survival after adjusting for age, grade, IDH status, and MGMT promoter methylation (median-split multivariate Cox; no continuous *MKI67* Cox in the revision).
5. **Multiplex Immunofluorescence (mIF)**: Confirmatory mIF assays independently validate the presence of both neoplastic (OLIG2+ Ki-67+) and microglial/macrophage (CD11b+ Ki-67+) proliferating subpopulations in glioblastoma tissue sections.

---

## Repository Structure

```text
Glioma_Ki67/
├── code_figure/                 # Main figure scripts (Figures 1–6)
├── data/                        # Metadata tables and public inputs (local / not always pushed)
├── scripts/                     # Pipeline helpers when present
├── revision/                    # Revision reproduction helpers (cluster-side)
├── figures/                     # Generated PDFs
├── results/                     # Intermediate RDS / tables
├── LICENSE
└── README.md
```

---


### Figure regeneration (revision pack)

Tracked plotters under `code_figure/` (cluster absolute inputs remain under `revision/`):

- `Figure2_review_assemble.R` — full Fig2 assembly; **Fig2B** left/right boxplots share log y-limits, identical 45° x labels, and cowplot `align="hv"` / `axis="tblr"` so baselines match. Vector `cairo_pdf` + 600 dpi PNG.
- `FigureS5_integration_diagnostics.R`, `FigureS6_Ki67_threshold_sensitivity.R`, `FigureS7_Ki67_definition_comparison.R` — supplementary diagnostics; **no under-plot caption strips** (captions stay in Word). Arial requested; cluster embeds Liberation Sans (Arial-metric).

Delivery copies (cluster): `revision/delivery_figures_YYYYMMDD/`.

## Installation & Prerequisites

The analyses are executed in **R (>= 4.3.0)**. Essential R dependencies include:

```r
install.packages(c("tidyverse", "cowplot", "ggtern", "patchwork", "ggpubr"))
BiocManager::install(c("SingleCellExperiment", "ComplexHeatmap"))
devtools::install_github("satijalab/seurat", "seurat5")
devtools::install_github("immunogenomics/harmony")
devtools::install_github("carmonalab/GeneNMF")
install.packages(c("survival", "survminer"))
```

---

## Step-by-Step Workflow (revision locks)

### Composite Ki-67+ calling (script_composite_union)
Union gating for revision statistics: raw *MKI67* count >= 1 **OR** `Phase_enhanced` in `{S, G2/M}`.  
Do **not** use `manuscript_text_union` (normalized MKI67>0 OR Seurat Phase S/G2M) for Fig3 clinical associations — that definition yields the incorrect age/sex/recurrence P values.

### Donor consolidation
Map `LGG-04-1`, `LGG-04-2`, `LGG-04-3` → donor `LGG-04` before patient-level tests (38 donors).

### Ecotype projection & bootstrap
Project donor proliferating fractions into ternary space (Tumor, Immune, Glial) with **10,000** within-donor bootstrap draws (seed 42). SF11488 is retained as a stromal-rich exploratory outlier.

### Survival
Median-split *MKI67* Kaplan–Meier and multivariate Cox (age, grade, IDH, MGMT). Continuous *MKI67* Cox / HR≈1.14 forest plots are **not** part of the revision figure code.

Example (when pipeline scripts are present):

```r
Rscript scripts/03_ecotype_classification.R --draws 10000 --seed 42
Rscript scripts/05_survival_and_cox.R
```

Fig3 clinical panels can be regenerated from precomputed donor metrics:

```bash
Rscript code_figure/Figure3.R
# or the revision review plotter (uses script_composite_union):
Rscript revision/20261023_reproduction/figure_review/plot_fig3_donor_review.R
```

---

## Data Availability

- **Single-cell RNA-seq datasets**: GEO accessions [GSE131928](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131928), [GSE182109](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182109), [GSE84465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84465), [GSE138794](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138794), and [GSE278456](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278456).
- **Bulk RNA-seq cohorts**: TCGA Pan-Glioma via [GDC](https://portal.gdc.cancer.gov/) and CGGA via [CGGA portal](http://www.cgga.org.cn/).

---

## Citation & Contact

If you use this codebase or methodology in your research, please cite:
> Dong Y, Xia K, Guo X, Ma W, Wang Y. Cellular Origins of Ki-67+ Cells in Gliomas. *BMC Cancer* (In Revision, 2026).

For inquiries or feedback, please contact:
- **Yu Wang, M.D., Ph.D.** (Corresponding author): `wangyu@pumch.cn`
- **Yijun Dong** (Lead author): `dongyijun20@mails.tsinghua.edu.cn`
