#load libraries
library(Seurat)
library(ggplot2)

# Read in counts and metadata file to create Seurat object
raw_counts = read.table("SS2_raw_counts.txt". header = T, row.names = 1, sep = "\t")
metadata = read.table("SS2_metadata.txt". header = T, row.names = 1, sep = "\t")
obj = CreateSeuratObject(counts = counts, meta.data = metadata)

#Normalization and clustering
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, features = VariableFeatures(object = obj))
obj <- FindNeighbors(obj, dims = 1:15)
obj <- FindClusters(obj, resolution = 0.3)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:10, n.neighbors = 15, seed.use = 24, spread = 0.5)


#score cells based on previous cell state signatures
gene_sets = readRDS("cell_state_gene_sets.Rds")
obj = AddModuleScore(obj, features = gene_sets, name = names(gene_sets))


#calculate stem_final and lineage score based
stem_score = function(x){
  return(x[1] - max(x[c(2,3)]))
}
module_scores_cells = obj@meta.data %>% select(stem_genes9,AC_genes10, OC_genes11)
module_scores_cells$stem_final = apply(module_scores_cells %>% select(stem_genes9,AC_genes10, OC_genes11), 1, stem_score)
module_scores_cells$lineage = rowMaxs(as.matrix(module_scores_cells %>% select(AC_genes10,OC_genes11)))
obj<- AddMetaData(
  object = obj,
  metadata = module_scores_cells$stem_final,
  col.name = 'stem_final'
)
obj<- AddMetaData(
  object = obj,
  metadata = module_scores_cells$lineage,
  col.name = 'lineage'
)

#Lineage plot
ggplot(subset(obj@meta.data, state_anno != "Normal"),aes(x = lineage, y = stem_final))+geom_point(size=2)+
  theme_classic()+
  xlab("Lineage score")+ylab("Stemness score")+
  theme(text = element_text(size=20),axis.text.x = element_text(angle=90, hjust=1), axis.text = element_text(color = "black"))
