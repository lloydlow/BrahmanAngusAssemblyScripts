#------------------------------------------------------
# Program name: brahman_angus_vcf_score_alt_allele_allChr_AngusRef.R
# Objective: load all chr for analysis from 
#   combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_<chr>.RData
#   and find overlap of selection interval with NCBI/EBI annotation
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(readr)
library(dplyr)
library(tidyr)
#library(reshape)
#library(refGenome)
library(ape)
library(GenomicRanges)
#library(IRanges)
library(stringr)

#ref is angus here
# path to filtered vcf where only wanted SNPs per chr are studied
dir1 <- 
  "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus/"

# load datasets
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus/combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_1.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus/combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_10.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus/combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_11.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus/combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_12.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus/combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_13.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus/combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_14.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus/combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_15.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus/combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_16.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus/combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_17.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus/combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_18.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus/combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_19.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus/combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_2.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus/combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_20.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus/combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_21.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus/combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_22.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus/combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_23.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus/combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_24.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus/combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_25.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus/combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_26.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus/combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_27.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus/combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_28.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus/combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_29.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus/combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_3.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus/combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_4.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus/combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_5.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus/combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_6.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus/combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_7.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus/combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_8.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus/combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_9.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus/combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_Y.RData")

#combine all 30 chr final SNPs file into one
merged_AngusRef <- do.call("rbind", list(combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_1,
                                combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_10,
                                combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_11,
                                combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_12,
                                combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_13,
                                combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_14,
                                combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_15,
                                combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_16,
                                combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_17,
                                combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_18,
                                combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_19,
                                combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_2,
                                combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_20,
                                combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_21,
                                combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_22,
                                combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_23,
                                combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_24,
                                combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_25,
                                combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_26,
                                combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_27,
                                combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_28,
                                combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_29,
                                combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_3,
                                combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_4,
                                combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_5,
                                combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_6,
                                combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_7,
                                combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_8,
                                combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_9,
                                combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_Y))

#use proportion per breed instead of proportion diff
proportion_diff_df_AngusRef <- merged_AngusRef %>% 
  group_by(CHROM,interval) %>% summarise(an_proportion_mean = mean(an_ALT_proportion),
                                         br_proportion_mean= mean(br_ALT_proportion),
                                         proportion_diff_mean = mean(proportion_diff),
                                         snp_counted = n())

proportion_diff_df_AngusRef_summ <- proportion_diff_df_AngusRef %>% filter(snp_counted >= 10)

#brahman percentile
quantile(proportion_diff_df_AngusRef_summ$br_proportion_mean, probs = seq(0,1,0.05))

#angus percentile
quantile(proportion_diff_df_AngusRef_summ$an_proportion_mean, probs = seq(0,1,0.05))

#proportion_diff_df_AngusRef$interval_no <- 1:nrow(proportion_diff_df_AngusRef)

proportion_diff_df_AngusRef_summ <- proportion_diff_df_AngusRef_summ %>% filter(an_proportion_mean <= 0.02) %>%
  filter(br_proportion_mean > 0.43) # %>% filter(snp_counted >= 10) #already done above

#split based on an1-an5 and br1-br6
proportion_diff_df_AngusRef_summ$interval <- gsub("\\(","",proportion_diff_df_AngusRef_summ$interval)
proportion_diff_df_AngusRef_summ$interval <- gsub("\\]","",proportion_diff_df_AngusRef_summ$interval)

proportion_diff_df_AngusRef_summ <- proportion_diff_df_AngusRef_summ %>% 
  separate(interval, c("start","end"),sep = ",",remove = FALSE,convert = TRUE) 

proportion_diff_df_AngusRef_summ$strand <- "*"
names(proportion_diff_df_AngusRef_summ)[1] <- "chr"

#write tsv to see what annotation features overlapped with
# options(scipen = 999)
# write.csv(proportion_diff_df_AngusRef_summ,file=paste0(dir1,"proportion_diff_df_AngusRef_summ.csv"),row.names = FALSE)

############################# EBI ###########################
#Angus EBI anno #chr only
ebi_anno_gtf <- read_tsv("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/EBI_annotation/ensembl_annotation_Paternal/Bos_indicus_x_bos_taurus_Pat.Bos_hybrid_PaternalHap_v2.0.94.chr.gtf",
                         comment = "#", col_names = FALSE,
                         col_types = list(col_character(),
                                          col_character(),
                                          col_character(),
                                          col_integer(),
                                          col_integer(),
                                          col_character(),
                                          col_character(),
                                          col_character(),
                                          col_character()))

colnames(ebi_anno_gtf) <- c("seqname","source","feature","start","end","score","strand","frame","attributes")

ebi_anno_gtf_gene <- ebi_anno_gtf %>% filter(feature == "gene") %>% select(seqname,start,end,strand,attributes)

names(ebi_anno_gtf_gene)[1] <- "chr"

#turn positively selected intervals and annotation into granges
#selection interval
selection_interval <- makeGRangesFromDataFrame(proportion_diff_df_AngusRef_summ, keep.extra.columns = TRUE, 
                                               seqnames.field="chr", start.field="start", 
                                               end.field="end", strand.field="strand")

#annotation - gene
gene <- makeGRangesFromDataFrame(ebi_anno_gtf_gene, keep.extra.columns = TRUE, 
                                 seqnames.field="chr", start.field="start", 
                                 end.field="end", strand.field="strand")

gene %>% as.data.frame() %>% dim()

subsetByOverlaps(selection_interval,gene)
subsetByOverlaps(gene,selection_interval)

ebi_overlap_selection_gene_df <- mergeByOverlaps(gene, selection_interval)

ebi_overlap_selection_gene_as_df <- as.data.frame(ebi_overlap_selection_gene_df)

#write out dataframe as csv
ebi_overlap_selection_gene_as_df <- ebi_overlap_selection_gene_as_df %>% 
  select(gene.seqnames:gene.attributes,selection_interval.seqnames:selection_interval.snp_counted) %>%
  arrange(gene.seqnames)

write.csv(ebi_overlap_selection_gene_as_df,file=paste0(dir1,"ebi_overlap_selection_gene_as_df_AngusRef.csv"),row.names = FALSE)

