############################################################
# Load libraries
############################################################

library(data.table)     # Fast data loading and manipulation
library(dplyr)          # Data manipulation (filter, select, mutate, etc.)
library(GenomicRanges)  # For genomic coordinate handling
library(parallel)       # For parallel processing
library(stringr)        # For string operations (regex, pattern matching)
library(ggplot2)        # For visualization
library(parallelDist)   # Parallel distance calculations
library(fastcluster)    # Fast hierarchical clustering
library(hexbin)         # Hexagonal binning (large scatterplots)
library(umap)           # UMAP dimensionality reduction
library(impute)         # k-nearest neighbor imputation
library(matrixStats)    # Efficient row/column stats for matrices
library(harmony)        # Batch correction (Harmony integration)

############################################################
# Load methylation and annotation data
############################################################

load("DNAme_5kb_bins.RData")   # Load methylation matrix at 5kb resolution
site = "5kb_bins"              # Define bin size used

# Alternative: set to "100kb" for larger genomic bins
site = "100kb"

# Load sample annotation metadata (filtered for CNV and other criteria)
load("/gpfs/commons/groups/landau_lab/rraviram/Suva_lab_GBM/RRBS_scRNA/cov_files/DNAme_annotation_cnv_filt_Aug2024.RData")

############################################################
# Update file names in annotation to match binarized methylation data
############################################################

# Replace ".cov.gz" with ".binarize.<site>.cov" for matching column names
annotation_filt$meth_cell = gsub(".cov.gz", paste0(".binarize.", site, ".cov"), annotation_filt$og_file)

# For promoter/enhancer-level data, use ".1kb.cov"
if(site %in% c("promoters", "enhancers")){
  annotation_filt$meth_cell = gsub(".cov.gz", paste0(".binarize.", site, ".1kb.cov"), annotation_filt$og_file)
}

# Assign row names for easy reference
rownames(annotation_filt) = annotation_filt$meth_cell

############################################################
# Filter methylation data by available samples
############################################################

# Helper function: count NA values in each row
which_row_na = function(x){
  return(length(which(is.na(x))))
}

# Define threshold: remove rows with ≥ 20% missing values (adjustable)
threshold <- ncol(meth_3) * 0.2

# Remove those rows from the data
filtered_data_2 <- meth_3[-which(meth_na_rows >= threshold), ]

############################################################
# Identify variable genomic regions
############################################################

# Compute row variance (across cells)
row_vars = rowVars(as.matrix(filtered_data_2), na.rm = TRUE)
summary(row_vars)

############################################################
# Impute missing methylation values (KNN imputation)
############################################################

# Transpose matrix: impute.knn expects samples in rows
# Keep only regions with variance > 0.08
meth_an_impute = impute.knn(t(filtered_data_2[which(row_vars > 0.08), ]), k = 5)

# Convert imputed matrix back to data frame
meth_an_impute_data = as.data.frame(meth_an_impute$data)

############################################################
# Principal Component Analysis (PCA)
############################################################

pca = prcomp(as.data.frame(meth_an_impute_data), scale = TRUE)

# Extract PCA scores (PC coordinates per cell)
scores = as.data.frame(pca$x)

# Compute % variance explained by PC1 and PC2
var_pc1 = round(((pca$sdev)^2 / sum(pca$sdev^2))[1], 2)
var_pc2 = round(((pca$sdev)^2 / sum(pca$sdev^2))[2], 2)

############################################################
# Batch correction with Harmony
############################################################

# Integrate data by sample to remove batch effects
harmony_embeddings <- harmony::RunHarmony(
  scores[, 1:20],             # Use first 20 PCs
  annotation_filt,            # Metadata for batch variable
  'sample_final',             # Batch column to correct for
  verbose = FALSE
)

############################################################
# Dimensionality reduction (UMAP)
############################################################

umap_out = umap(harmony_embeddings)

# Create UMAP dataframe
umap_df = data.frame(umap_out$layout)
colnames(umap_df) = c("UMAP1", "UMAP2")

# Add metadata columns for visualization
umap_df$gcimp         = annotation_filt$gcimp
umap_df$Methylation   = annotation_filt$Methylation
umap_df$occ           = annotation_filt$Occurence
umap_df$sample_final  = annotation_filt$sample_final
umap_df$cell_type_2   = annotation_filt$cell_type_2
umap_df$cell_type     = annotation_filt$cell_type
umap_df$dname_amb     = annotation_filt$dname_amb
umap_df$Total_Reads   = annotation_filt$Total_Reads
umap_df$CpG_Sites     = annotation_filt$CpG_Sites
umap_df$stem_final    = annotation_filt$stem_final
umap_df$G1_S_genes12  = annotation_filt$G1_S_genes12
umap_df$state_anno    = annotation_filt$state_anno

############################################################
# Cluster cells in UMAP space
############################################################

kmeans_cluster = kmeans(umap_out$layout, centers = 3, nstart = 20)
umap_df$kmeans_cluster = kmeans_cluster$cluster

############################################################
# Visualization: UMAPs colored by various annotations
############################################################

# By G-CIMP status
ggplot(umap_df, aes(UMAP1, UMAP2)) +
  geom_point(size = 1, aes(colour = gcimp)) +
  theme_classic() + scale_colour_viridis_c() 

# By overall methylation level
ggplot(umap_df, aes(UMAP1, UMAP2)) +
  geom_point(size = 1, aes(colour = Methylation)) +
  theme_classic() + scale_colour_viridis_c() 

# By stemness, cell cycle, or state annotations
ggplot(umap_df, aes(UMAP1, UMAP2)) +
  geom_point(size = 1, aes(colour = stem_final)) +
  theme_classic() + scale_colour_viridis_c() 

ggplot(umap_df, aes(UMAP1, UMAP2)) +
  geom_point(size = 1, aes(colour = G1_S_genes12)) +
  theme_classic() + scale_colour_viridis_c() 

ggplot(umap_df, aes(UMAP1, UMAP2)) +
  geom_point(size = 1, aes(colour = state_anno)) +
  theme_classic() 

# By occurence (primary versus recurrent)
ggplot(umap_df, aes(UMAP1, UMAP2)) +
  geom_point(size = 1, aes(colour = occ)) +
  theme_classic() + scale_colour_manual(values = c("red", "darkred"))


