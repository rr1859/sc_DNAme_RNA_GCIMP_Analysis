library(infercnv)
library(tidyverse)
library(Seurat)

#gene information
genes=read.table("hg38/gencode.v29.basic.annotation_4col.bed.gz", stringsAsFactors = FALSE)
genes=genes[which(!duplicated(genes[,4])),]
rownames(genes) = genes[,4]
genes=genes[,1:3]

#load in RNA object and get raw counts for malignant, myeloid and oligo cells
obj = readRDS("SS2_RNA_obj.Rds")
cell_to_keep = colnames(subset(obj.2, cell_type_2 == "Malignant" | cell_type %in% c("Myeloid", "Oligo") & Patient %in% patient_order))
mat=obj@assays$RNA@counts[, cell_to_keep]


#annotation for cell type and patient - Myeloid and Oligo cells used as a refernce
annotation_df <- obj@meta.data %>% select(cell_type, Patient)
annotation_df$cell_type = ifelse(annotation_df$cell_type %in% c("Myeloid", "Oligo"), "normal", annotation_df$Patient)
#infercnv - re
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=mat,
                                    annotations_file=annotation_df %>% select(Patient),
                                    delim="\t",
                                    gene_order_file=genes,
                                    ref_group_names="normal") 

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir="all_samples_infercnv_w201_myeloid_oligo_ref",cluster_by_groups = TRUE,
                             denoise=TRUE,window_length = 201,
                             HMM=FALSE)
