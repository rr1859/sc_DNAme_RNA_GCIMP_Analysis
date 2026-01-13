module load betools
for i in *binarize.cov; do bedtools intersect -a ${i} -b /gpfs/commons/groups/landau_lab/rraviram/Suva_lab_GBM/RRBS_scRNA/cov_files/probes_684_list_loci_500bprange_2.bed -wa -wb > ${i%.*}.gcimp684_500.cov; done

