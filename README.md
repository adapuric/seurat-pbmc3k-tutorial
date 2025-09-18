# Seurat PBMC3k Tutorial

This repository contains my notes and code from working through the [Seurat PBMC3k tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html).

The goal is to practice single-cell RNA-seq analysis workflows using Seurat in R, understand the logic of each step, and build a reproducible pipeline.

---

## 🧬 Dataset
- 3k Peripheral Blood Mononuclear Cells (PBMCs) from 10x Genomics
- Download link: [10x Genomics dataset](https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_3k/pbmc_3k_filtered_gene_bc_matrices.tar.gz)

---

## 📦 Dependencies
- R (≥ 4.0)
- [Seurat](https://satijalab.org/seurat/)  
- dplyr, patchwork, ggplot2

Install with:
```r
install.packages(c("dplyr", "patchwork", "ggplot2"))
# Seurat
install.packages("Seurat")
