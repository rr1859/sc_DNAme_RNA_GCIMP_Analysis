
sites = "LTR_ERV"
library(data.table)
library(dplyr)
library(GenomicRanges)
library(parallel)
library(stringr)
library(ggplot2)
files = list.files( pattern = paste0("*binarize.", sites, ".cov$" ))
files = paste0(files)
f = as.list(files)
## Given GENE, add reads counts per GENE 
count_cpgs <- function(x,m) { k = dplyr::filter( m , V11 == x ) ; val = length(k$V5) ; return( c(x, val) ) } 
covs <- mclapply( f , function(x) { m <- fread(x) ; 
m$ID <- paste0( m$V1, m$V2, m$V3, m$V11) ; 
m <- m[ !duplicated(m$ID) , ]  ; 
new = m
if(sites %in% c("gene_body_noTSS1kb", "promoters.1kb", "promoters.enhancers.1kb")){
  new = m
  new$V12 = new$V11
} else if(sites != "promoters"){
  new = m
  new$V11 = paste(m$V8, m$V9, m$V10, sep = "_")
  new$V12 = paste(m$V8, m$V9, m$V10, sep = "_")
}
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
saveRDS(results, paste0("TSS_1kb_binarize_",sites , "_cpgs.rds"))
