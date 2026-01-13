

library(data.table)
library(dplyr)
library(GenomicRanges)
library(parallel)
library(stringr)
library(ggplot2)


sites = "gcimp684_500"
files = list.files( pattern = paste0("*binarize.", sites, ".cov$" ))
files = paste0(files)
f = as.list(files)
## Given GENE, add reads counts per GENE 
count_cpgs <- function(x,m) { k = dplyr::filter( m , V11 == x ) ; val = length(k$V5) ; return( c(x, val) ) } 
covs <- mclapply( f , function(x) { m <- fread(x) ; 
m$ID <- paste0( m$V1, m$V2, m$V3, m$V11) ; 
m <- m[ !duplicated(m$ID) , ]  ; 
new = m
new$V11 = paste(m$V8, m$V9, m$V10, sep = "_")
new$V12 = paste(m$V8, m$V9, m$V10, sep = "_")

out <- sapply( unique( new$V11 ) , count_cpgs , new  ) ;
out <- data.frame( t(out) );
#out <- data.frame( table(new$V12) );
#rownames(out) = out$Var1
colnames(out) <- c("GENE", basename(x) ) ; 
return(out) } 
, mc.cores = 50)
fun = function(x,y) merge(x,y, by = c("GENE") , all = T) 
## Reduce cell
results <- Reduce(fun, covs)
saveRDS(results, paste0("binarize_",sites , "_cpgs.rds"))

#Methylation 
add_reads <- function(x,m) { k = dplyr::filter( m , V11 == x ) ; val = ( sum(k$V5) / length(k$V5) )  ; return( c(x, val) ) } 

covs <- mclapply( f , function(x) { m <- fread(x) ; 
m$ID <- paste0( m$V1, m$V2, m$V3, m$V11) ; 
m <- m[ !duplicated(m$ID) , ]  ; 
new = m
new$V11 = paste(m$V8, m$V9, m$V10, sep = "_")
new$V12 = paste(m$V8, m$V9, m$V10, sep = "_")
out <- sapply( unique( new$V11 ) , add_reads , new  ) ;
out <- data.frame( t(out) );
colnames(out) <- c("GENE", basename(x) ) ; 
return(out) }, mc.cores = 500 
)
#save(covs, file = paste0("TSS_1kb_binarize_", sites, "_meth_covs.rds"))
fun = function(x,y) merge(x,y, by = c("GENE") , all = T) 

## Reduce cell
results <- Reduce(fun, covs)
saveRDS(results, paste0("binarize_", sites, "_meth.rds"))

meth = readRDS(paste0("TSS_1kb_binarize_", site, "_meth.rds"))
meth$GENE <- as.character(meth$GENE)
meth <- meth[ !is.na(meth$GENE) ,]
## A gene by cell matrix (inputs are # of meth)
rownames(meth) <- meth$GENE
meth <- meth[, -1]

annotation = read.table("annotation.txt", header = T, sep = "\t")

meth_2 = meth[,subset(annotation, cell_type_2 == "Malignant")$meth_cell]
meth_an = apply(meth_2, 2, function(x) as.numeric(x))
meth_an = as.data.frame(meth_an)
rownames(meth_an) = rownames(meth)

annotation$gcimp_score = rowMeans(meth_an, na.rm = T)          
                
