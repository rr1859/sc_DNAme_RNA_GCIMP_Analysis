# ============================================================
# SS2 RNA Differential Expression and Heatmap Analysis
# Description: Perform Seurat-based differential gene expression 
#              analysis between GCIMP subtypes and visualize 
#              results using pheatmap.
# ============================================================

# ---- Load libraries ----
library(Seurat)        # For single-cell RNA-seq data analysis
library(dplyr)         # For data wrangling
library(MAST)          # For differential expression analysis
library(pheatmap)      # For creating heatmaps
library(viridisLite)   # For color palettes

# ---- Load Seurat object ----
obj = readRDS("SS2_RNA_seurat_obj.Rds")  # Load the single-cell Seurat object

# ---- Set cell identities ----
Idents(obj.2) = "GCIMP_int"

# ---- Subset object ----
# Remove normal samples and those with G-CIMP intermediate classification
obj.p = subset(obj.2, state_anno != "Normal" & GCIMP_int != "GCIMP_int")

# ---- Normalize and scale data ----
obj.p = ScaleData(obj.p, features = rownames(obj.p))

# ---- Differential expression analysis ----
# Compare G-CIMP_High vs GCIMP_Low, controlling for patient effects
diff_analysis_gcimp = FindMarkers(
  obj.p,
  ident.1 = "GCIMP_High",
  ident.2 = "GCIMP_Low",
  min.pct = 0.1,
  latent.vars = "Patient",
  test.use = "MAST"
)

#Compare G-CIMP_High vs GCIMP_Low uncorrected for patient specific effects
diff_analysis_gcimp_noc=FindMarkers(obj.p, ident.1 = "GCIMP_High", ident.2 = "GCIMP_Low", min.pct = 0.1)
diff_analysis_gcimp_noc$gene = rownames(diff_analysis_gcimp_noc)
diff_analysis_gcimp_noc_sig = subset(diff_analysis_gcimp_noc, abs(avg_log2FC) > 1 & p_val_adj < 0.01)

# ---- Compute average expression ----
# For genes with |log2FC| > 2 and adjusted p < 0.01
average_expr = AverageExpression(
  object = obj.p,
  group.by = c('GCIMP_int', 'sample_final'),
  features = unique(subset(diff_analysis_gcimp, abs(avg_log2FC) > 1 & p_val_adj < 0.01)$gene)
)$RNA

# ---- Prepare annotation dataframe ----
annotation_df = data.frame(
  annotation_filt_2 %>%
    group_by(Patient, gcimp_cat_int, Codel_status, Occurence, Sample_Final) %>%
    summarize(mean = mean(gcimp), count = n())
)

annotation_df = annotation_df[grep('GCIMP_High|GCIMP_Low', annotation_df$gcimp_cat_int), ]
rownames(annotation_df) = paste(annotation_df$gcimp_cat_int, annotation_df$Sample_Final, sep = "_")

# ---- Plot heatmap ----
pheatmap(
  average_expr[, rownames(annotation_df)], 
  scale = "row",
  clustering_method = "ward",
  color = magma(200),
  annotation_col = annotation_df %>% select(Patient, Codel_status, Occurence, gcimp_cat_int, mean),
  fontsize_row = 4
)

# ---- End of Script ----
