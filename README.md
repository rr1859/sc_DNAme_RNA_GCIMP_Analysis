# Dual capture of RNA and DNA methylation from longitudinal gliomas (single nuclei) (IN DEVELOPMENT)

This repository contains R scripts for analyzing single cell DNA methylation and RNA-seq data from IDH-mutant longitudinal gliomas
Note: final data will be added once finalized (manuscript in revision)

## RNA analysis

### Clustering and cell state assignment
We performed clustering analysis using Seurat and assigned cell state clusters for malignant cells based on Tirosh et al, Nature, 2016. 

### inferCNV
We use [inferCNV](https://github.com/broadinstitute/infercnv) to detect CNV events from snRNA-seq data (myeloid and oligodendrocyte cells were used reference non-malignant cells)

### Differential analysis
We performed differential gene expression analysis between G-CIMP high and low tumors using MAST (patient as latent variable) and without correction for patient information (FindMarkers function in Seurat)

---

## DNA methylation analysis

### CNV analysis
We detected large scale CNV events based on the log2 fold-change differences in 2MB bins in malignant cells compared to non-malignant cells

### G-CIMP score
We calculated the mean G-CIMP score for each cell based on the mean DNA methylation across 1kb regions around G-CIMP probes.

### PATH anaysis 
We used [PATH](https://github.com/landau-lab/PATH) (Schiffman et al, Nature Genetics, 2023) to infer G-CIMP and cell state heritability based on phylogenetic tress from DNA methylation data
 
