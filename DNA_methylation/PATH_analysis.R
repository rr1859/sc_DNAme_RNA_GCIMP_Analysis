# PATH analysis

library(ape)
library(phytools)
library(PATH)
library(castor)
library(qlcMatrix)
library(dplyr)
library(ggplot2)
library(patchwork)
library(reshape2)
library(ggpubr)


load("path_analysis.RData")

#' @param edited_tree_path Character string. File path to the edited tree object
#'   used for downstream analysis.
#' @param tree_path Character string. File path to the original (unedited) tree
#'   used as a reference.
#' @param features Character vector or data frame. Features to be analyzed or
#'   mapped onto the tree (e.g., stem_final, gcimp_score).
#' @param return_value Character or logical. Specifies what the function should
#'   return (e.g., "value" : z-score values, "plot" : auto-correlation plot, 
#' "auto_values": all values).
#'
#' @return Depends on `return_value`. May return a data frame, list, or
#'   visualization object.

path_analysis <- function(edited_tree_path, tree_path, features, return_value){
  tr <- read.newick(edited_tree_path)
  tr_unedit <- read.newick(tree_path)
  # remove normal cells
  normal <- tr[["tip.label"]][grep("Oligo|Myeloid|Neurons|Endothelial|Tcells", tr[["tip.label"]])]
  normal_unedit <- tr_unedit[["tip.label"]][grep("Oligo|Myeloid|Neurons|Endothelial|Tcells", tr[["tip.label"]])]
  # drop normal cells from both trees
  tr <- drop.tip(tr, normal)
  tr_unedit <- drop.tip(tr_unedit, normal_unedit)
  labels_orig = tr_unedit[["tip.label"]] #original unedited labels - RR
  labels <- tr_unedit[["tip.label"]]
  labels_noedit <- tr[["tip.label"]]
  
  # get rid of index for some sample names
  labels <- gsub("_Index", "", labels)
  labels <- gsub("RRBS", "XRBS", labels)
  names <- annotation_filt[which(annotation_filt$DNA_id %in% labels),]
  names_select <- as.matrix(names[, c("DNA_id",features, "gcimp")])
  
  # get average/var gcimp
  avg_gcimp_score <- mean(as.numeric(names_select[,c("gcimp")]))
  var_gcimp_score <- var(as.numeric(names_select[,c("gcimp")]))
  
  # reorder matrix
  labels_orig_filt = labels_orig[which(labels %in% names_select[, 1])] #RR
  tree_unedited_labels = tr_unedit$tip.label #RR
  tree_labels_keep = tr$tip.label[which(tree_unedited_labels %in% labels_orig_filt)] #RR
  tr <- keep.tip(tr, tree_labels_keep) #RR
  tr_unedit <- keep.tip(tr_unedit, labels_orig_filt) #RR

  matching_indices <- match(labels[which(labels %in% names_select[, 1])], names_select[, 1])
  genotype_matrix <- names_select[matching_indices, ]

  # remove names from genotype matrix
  genotype_mat_num <- genotype_matrix[,-c(1)]
  mat_num <- matrix(as.numeric(genotype_mat_num),
                    ncol = ncol(genotype_mat_num))
  colnames(mat_num) <- c(features, "gcimp")
  # Compute the phylogenetic node distances between cells in order to measure phylogenetic correlations.
  Winv <- inv_tree_dist(tr, node = TRUE, norm = TRUE)
  # Josh switched to exp.tree, switch back! alterations:
  #Winv <- exp.tree.dist(tr, node = TRUE, norm = TRUE)
  #Winv <- one_node.tree.dist(tr, norm = T)
  
  # Compute phylogenetic correlations between GBM modules. 
  modxcor <- xcor(mat_num, Winv)
  
  Idf <- reshape2::melt(modxcor$phy_cor, 
                        value.name = "I")
  Zdf <- reshape2::melt(modxcor$Z.score, 
                        value.name = "Z")
  
  df <- full_join(Idf, Zdf, by=c("Var1", "Var2"))
  
  df <- df %>% mutate(Var1 = as.factor(Var1), 
                      Var2 = as.factor(Var2))
  df$sample = sub(".*/([^/]+)_filter.phy.treefile", "\\1", tree_path)
  # Phylogenetic auto-correlation bar plot.
  maxz <- max(abs(df$Z))
  ## switch to non z-score value
  z_score_df <- df %>% dplyr::filter(Var1 == Var2) %>% dplyr::select("Z")
  #z_score_df <- df %>% dplyr::filter(Var1 == Var2) %>% dplyr::select("I")
  ##
  colnames(z_score_df) <- tree_path
  z_score_df[nrow(z_score_df)+1,1] <- avg_gcimp_score
  z_score_df[nrow(z_score_df)+1,1] <- var_gcimp_score
  auto_plot =ggplot(df,aes(x=Var1, y=Var2, fill=Z)) +
    geom_tile(col="white") +
    scale_fill_distiller(palette = 5, type = "div",
                         limits=c(-15,15)) +
    theme_classic() +
    scale_y_discrete(limits=rev) +
    labs(fill="Phylogenetic\ncorrelation\nz score") +
    xlab("Cell state") + ylab("Cell state") +
    theme(aspect.ratio = 1) +ggtitle(sub(".*/([^/]+)_filter.phy.treefile", "\\1", tree_path))
  if(return_value == "value"){
    return(z_score_df)
  } else if(return_value == "plot"){
    return(auto_plot)
  } else if(return_value == "auto_values"){
    return(df)
  }
}

#run path_analysis function for each sample
sample = "TKU3197"
path_analysis (edited_tree_path = subset(annotation, sample == "TKU3197")$edited_tree_path,
               tree_path = subset(annotation, sample == "TKU3197")$tree_path,
               features = c("stem_final", "gcimp"), 
               return_value = c("auto_values"))

