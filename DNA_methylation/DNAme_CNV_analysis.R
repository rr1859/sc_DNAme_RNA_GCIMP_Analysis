library(GenomicRanges)
library(parallel)
library(stringr)
library(ggplot2)
library(dplyr)

load("DNAme_annotation_cnv_filt_Jan2024.RData")


bed_file=read.table("/gpfs/commons/groups/landau_lab/rraviram/hg38/hg38_20mb_5mb.bed", stringsAsFactors = FALSE)
bed_file.gr=GRanges(seqnames=bed_file[,1], IRanges(bed_file[,2], bed_file[,3]))
cpg_bin_counts=do.call(cbind, mclapply(1:nrow(annotation_filt), function(i){
#cpg_bin_counts=do.call(cbind, mclapply(1:10, function(i){
  print(i)
  cov_data=read.table(gsub(".cov.", ".cov.bin.", as.character(annotation_filt$og_file[i])), stringsAsFactors = FALSE)
  cov_data.gr=GRanges(seqnames=cov_data[,1], IRanges(cov_data[,2], cov_data[,3]), methylated=cov_data[,5], unmethylated=cov_data[,6], sample=cov_data[,7])
  cov_bed.idy=queryHits(findOverlaps(cov_data.gr, bed_file.gr))
  cov_bed.gr = cov_data.gr[unique(cov_bed.idy),]
  sample_total=sum(cov_bed.gr$methylated)+sum(cov_bed.gr$unmethylated)
  counts=unlist(lapply(1:length(bed_file.gr), function(i){
    cov_bed.sample.bin=cov_bed.gr[unique(queryHits(findOverlaps(cov_bed.gr, bed_file.gr[i,]))),]
    #print(cov_cpg.sample.bin)
    if(length(cov_bed.sample.bin) > 0){
      return(((sum(cov_bed.sample.bin$methylated)+sum(cov_bed.sample.bin$unmethylated))/sample_total)*10^6)
    } else {
      return(0)
    }
  }))
  return(counts)
}, mc.cores = 40))


cpg_bin_counts = t(cpg_bin_counts)
rownames(cpg_bin_counts) = annotation_filt$og_file


normal_cells=cpg_bin_counts[as.character(subset(annotation_filt, label == "Normal")$og_file),]
all_cells=rbind(cpg_bin_counts[as.character(subset(annotation_filt, label == "Normal")$og_file),],
                cpg_bin_counts[as.character(subset(annotation_filt, label != "Normal")$og_file),])

normal_cells_t=t(normal_cells)
median_value_normal=colMedians(normal_cells)
median_value_normal[which(median_value_normal == 0)] = 1
norm_malignant_fc=do.call(rbind, lapply(1:nrow(all_cells), function(k){
  if(k %% 100 == 0){
    print(k)
  }
  x_log=log(all_cells[k,]/median_value_normal,2)
  return(x_log)
}))

