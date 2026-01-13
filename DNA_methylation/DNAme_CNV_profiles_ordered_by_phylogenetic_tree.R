library(pheatmap)

load("phylo_CNV.RData")

#save a pdf of all heatmaps
pdf("All_CNV_20Mb_phylo_order_malignant.pdf", height = 14, width = 10)
for(patient_name in patient_order){
  time_points = unique(subset(annotation_filt_id, Patient == patient_name)$sample_final)
  names_both = NULL
  labels_both = NULL
  for(time_point in time_points){
    if(time_point %in% samples_order){
      edited_tree_path = files_edited_2[time_point]
      tree_path = files_2[time_point]
      tr <- read.newick(edited_tree_path)
      tr_unedit <- read.newick(tree_path)
      labels_orig = tr_unedit[["tip.label"]] #original unedited labels - RR
      labels <- tr_unedit[["tip.label"]]
      # get rid of index for some sample names
      labels <- gsub("_Index", "", labels)
      labels <- gsub("RRBS", "XRBS", labels)
      
      names <- annotation_filt[which(annotation_filt$DNA_id %in% labels),]
      names_select <- as.matrix(names[, c("DNA_id",features, "gcimp")])
      
      # get average/var gcimp
      avg_gcimp_score <- mean(as.numeric(names_select[,c("gcimp")]))
      var_gcimp_score <- var(as.numeric(names_select[,c("gcimp")]))
      
      # reorder matrix
      labels_all = labels
      labels = labels_all[which(labels %in% names_select[, 1])]
      labels_orig = labels_orig[which(labels_all %in% names_select[, 1])] #RR
      tree_unedited_labels = tr_unedit$tip.label #RR
      tree_labels_keep = tr$tip.label[which(tree_unedited_labels %in% labels_orig)] #RR
      tr <- keep.tip(tr, tree_labels_keep) #RR
      tr_unedit <- keep.tip(tr_unedit, labels_orig) #RR
      names_both = rbind(names_both,names_select)
      labels_both = append(labels_both, labels)
    }

  }
  names_both_malignant = names_both[which(names_both[,1] %in% as.character(annotation_filt_id$DNA_id)),]
  norm_malignant_fc_2_p= norm_malignant_fc_2[as.character(annotation_filt_id[as.character(labels_both[which(labels_both %in% names_both_malignant[, 1])]),]$og_file),1:587]
  rownames(norm_malignant_fc_2_p) = annotation_filt_id[labels_both[which(labels_both %in% names_both_malignant[, 1])],]$og_file
  print(pheatmap(norm_malignant_fc_2_p, annotation_col = col_colors, annotation_row=annotation_filt_2 %>% select(gcimp, stem_final, Occurence,Sample_Final),
              cluster_rows=F,cluster_cols=FALSE, labels_row=NULL, labels_col=c(""), fontsize_row=0.1, fontize=4, scale="none", 
              col= colorRampPalette(c("navy", "white", "firebrick3"))(100)))


}
dev.off()
