#------------------------------------------------------
# Program name: brahman_angus_vcf_score_alt_allele_allChr_SevenBreeds_BrahmanRef.R
# Objective: load all chr for analysis from 
#   combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_<chr>.RData
#   and find overlap of selection interval with EBI annotation
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

#ref is brahman here
# path to filtered vcf where only wanted SNPs per chr are studied
dir1 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Brahman_SevenBreeds/"

# load datasets
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Brahman_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_1.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Brahman_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_10.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Brahman_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_11.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Brahman_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_12.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Brahman_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_13.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Brahman_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_14.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Brahman_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_15.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Brahman_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_16.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Brahman_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_17.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Brahman_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_18.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Brahman_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_19.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Brahman_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_2.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Brahman_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_20.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Brahman_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_21.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Brahman_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_22.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Brahman_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_23.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Brahman_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_24.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Brahman_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_25.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Brahman_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_26.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Brahman_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_27.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Brahman_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_28.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Brahman_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_29.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Brahman_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_3.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Brahman_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_4.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Brahman_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_5.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Brahman_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_6.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Brahman_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_7.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Brahman_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_8.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Brahman_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_9.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Brahman_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_X.RData")

#combine all 30 chr final SNPs file into one
merged <- do.call("rbind", list(combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_1,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_10,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_11,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_12,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_13,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_14,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_15,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_16,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_17,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_18,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_19,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_2,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_20,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_21,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_22,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_23,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_24,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_25,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_26,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_27,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_28,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_29,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_3,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_4,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_5,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_6,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_7,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_8,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_9,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_X))

#use proportion per breed instead of proportion diff
proportion_diff_df <- merged %>% 
  group_by(CHROM,interval) %>% summarise(an_proportion_mean = mean(an_ALT_proportion),
                                         br_proportion_mean= mean(br_ALT_proportion),
                                         ge_proportion_mean= mean(ge_ALT_proportion),
                                         he_proportion_mean= mean(he_ALT_proportion),
                                         re_proportion_mean= mean(re_ALT_proportion),
                                         sh_proportion_mean= mean(sh_ALT_proportion),
                                         si_proportion_mean= mean(si_ALT_proportion),
                                         snp_counted = n(),mean_diff_POS = mean(diff(POS)))

summary(proportion_diff_df$snp_counted) 

proportion_diff_df_summ <- proportion_diff_df %>% filter(snp_counted >= 10)

#all taurine combine mean proportion
proportion_diff_df_summ <- proportion_diff_df_summ %>% rowwise() %>% 
  mutate(tau_proportion_mean = sum(an_proportion_mean,ge_proportion_mean,he_proportion_mean,
                                    re_proportion_mean,sh_proportion_mean,si_proportion_mean)/6)
  
summary(proportion_diff_df_summ$mean_diff_POS)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 314    3947    5028    5321    6466   11020 

nrow(proportion_diff_df_summ)/nrow(proportion_diff_df)

#taurine percentile
quantile(proportion_diff_df_summ$tau_proportion_mean, probs = seq(0,1,0.05))
# 0%         5%        10%        15%        20%        25%        30%        35%        40% 
# 0.00000000 0.02996899 0.04941288 0.06843434 0.08578683 0.10347222 0.12157617 0.14027778 0.16025641 
# 45%        50%        55%        60%        65%        70%        75%        80%        85% 
# 0.17909635 0.19921875 0.21781863 0.23665751 0.25722222 0.27725524 0.29977778 0.32566919 0.35861111 
# 90%        95%       100% 
# 0.39759098 0.45494285 0.91792929 

#brahman percentile
quantile(proportion_diff_df_summ$br_proportion_mean, probs = seq(0,1,0.05))
# 0%        5%       10%       15%       20%       25%       30%       35%       40%       45% 
#   0.0000000 0.1000000 0.1272727 0.1463736 0.1611111 0.1750000 0.1875000 0.2000000 0.2100000 0.2210526 
# 50%       55%       60%       65%       70%       75%       80%       85%       90%       95% 
#   0.2314286 0.2434783 0.2555556 0.2681818 0.2818182 0.2958333 0.3142857 0.3352941 0.3650000 0.4117647 
# 100% 
# 0.8357143 

#histogram
tiff(filename = "hist_Brahman_vs_Angus_SevenBreeds_proportion_mean.tiff",width = 800)
par(mfrow=c(1,2))
hist(proportion_diff_df_summ$br_proportion_mean,xlab="Proportion mean of alternate allele in Brahman",
     main="",ylim=c(0,5000))
abline(v=0.1,lty="dashed",col="red")
arrows(0.1,4000,0,4000,lwd=0.75,col="red",length=0.1)

hist(proportion_diff_df_summ$tau_proportion_mean,xlab="Proportion mean of alternate allele in six taurine breeds",
     main="",ylim=c(0,5000))
abline(v=0.4,lty="dashed",col="red")
arrows(0.4,4000,0.95,4000,lwd=0.75,col="red",length=0.1)
par(mfrow=c(1,1))
dev.off()

proportion_diff_df_summ <- proportion_diff_df_summ %>% filter(br_proportion_mean <= 0.1) %>%
  filter(tau_proportion_mean >= 0.4) 

summary(proportion_diff_df_summ$tau_proportion_mean)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.4013  0.4359  0.4677  0.4839  0.5194  0.7298 
summary(proportion_diff_df_summ$br_proportion_mean)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.005556 0.052864 0.082576 0.071912 0.091827 0.100000 

#split based on all breeds
proportion_diff_df_summ$interval <- gsub("\\(","",proportion_diff_df_summ$interval)
proportion_diff_df_summ$interval <- gsub("\\]","",proportion_diff_df_summ$interval)

proportion_diff_df_summ <- proportion_diff_df_summ %>% 
  separate(interval, c("start","end"),sep = ",",remove = FALSE,convert = TRUE) 

proportion_diff_df_summ$strand <- "*"
names(proportion_diff_df_summ)[1] <- "chr"

# finding consecutive selection region
# vec <- c()
# for (i in 1:nrow(proportion_diff_df_summ)){
#   dummy <- proportion_diff_df_summ$start[i+1] - proportion_diff_df_summ$end[i]
#   vec <- c(vec,dummy)
# }

############################# EBI ###########################
# #Brahman EBI anno #chr only
# ebi_anno_gtf <- read_tsv("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/EBI_annotation/ensembl_annotation_Maternal_20190320/Bos_indicus_x_bos_taurus_Mom.UOA_Brahman_1.96.chr.gtf",
#                          comment = "#", col_names = FALSE,
#                          col_types = list(col_character(),
#                                           col_character(),
#                                           col_character(),
#                                           col_integer(),
#                                           col_integer(),
#                                           col_character(),
#                                           col_character(),
#                                           col_character(),
#                                           col_character()))
# 
# colnames(ebi_anno_gtf) <- c("seqname","source","feature","start","end","score","strand","frame","attributes")
# 
# ebi_anno_gtf_gene <- ebi_anno_gtf %>% filter(feature == "gene") %>% select(seqname,start,end,strand,attributes)
# 
# names(ebi_anno_gtf_gene)[1] <- "chr"
# 
# #turn positively selected intervals and annotation into granges
# #selection interval
# selection_interval <- makeGRangesFromDataFrame(proportion_diff_df_summ, keep.extra.columns = TRUE, 
#                                                seqnames.field="chr", start.field="start", 
#                                                end.field="end", strand.field="strand")
# 
# #annotation - gene
# gene <- makeGRangesFromDataFrame(ebi_anno_gtf_gene, keep.extra.columns = TRUE, 
#                                  seqnames.field="chr", start.field="start", 
#                                  end.field="end", strand.field="strand")
# 
# gene %>% as.data.frame() %>% dim()
# 
# # subsetByOverlaps(selection_interval,gene)
# # subsetByOverlaps(gene,selection_interval)
# 
# ebi_overlap_selection_gene_df <- mergeByOverlaps(gene, selection_interval)
# 
# ebi_overlap_selection_gene_as_df <- as.data.frame(ebi_overlap_selection_gene_df)
# 
# ebi_overlap_selection_gene_as_df <- ebi_overlap_selection_gene_as_df %>% select(gene.seqnames:gene.attributes,selection_interval.seqnames:selection_interval.tau_proportion_mean)
# 
# #separate the attributes column
# ebi_overlap_selection_gene_as_df <- ebi_overlap_selection_gene_as_df %>% 
#   separate(gene.attributes, c("gene_id","gene_version","gene_source","gene_biotype"),sep = ";",extra = "drop", fill = "right")
# 
# #tidy up the values in gene_id and gene_biotype
# ebi_overlap_selection_gene_as_df$gene_id <- gsub("gene_id \"","",ebi_overlap_selection_gene_as_df$gene_id)
# ebi_overlap_selection_gene_as_df$gene_id <- gsub("\\\"","",ebi_overlap_selection_gene_as_df$gene_id)
# 
# ebi_overlap_selection_gene_as_df$gene_biotype <- gsub("gene_biotype \"","",ebi_overlap_selection_gene_as_df$gene_biotype)
# ebi_overlap_selection_gene_as_df$gene_biotype <- gsub("\\\"","",ebi_overlap_selection_gene_as_df$gene_biotype)
# ebi_overlap_selection_gene_as_df$gene_biotype <- gsub("^ ","",ebi_overlap_selection_gene_as_df$gene_biotype)
# 
# #write out dataframe as csv
# ebi_overlap_selection_gene_as_df <- ebi_overlap_selection_gene_as_df %>% 
#   select(gene.seqnames:gene_id,gene_biotype,selection_interval.start:selection_interval.end,selection_interval.br_proportion_mean,selection_interval.snp_counted:selection_interval.tau_proportion_mean) %>%
#   arrange(gene.seqnames)
# 
# #ebi_overlap_100k_selection_gene_as_df_SevenBreeds_v2.csv is the 2nd version with tidied gene attributes
# write.csv(ebi_overlap_selection_gene_as_df,file=paste0(dir1,"ebi_overlap_100k_selection_gene_as_df_SevenBreeds_v2.csv"),
#           row.names = FALSE)

############################# EBI with gene names ###########################
#Brahman EBI anno #chr only
ebi_anno_gtf <- read_tsv("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/EBI_annotation/ensembl_annotation_Maternal_20190620/Bos_indicus_hybrid.UOA_Brahman_1.97.chr.gtf",
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
selection_interval <- makeGRangesFromDataFrame(proportion_diff_df_summ, keep.extra.columns = TRUE, 
                                               seqnames.field="chr", start.field="start", 
                                               end.field="end", strand.field="strand")

#annotation - gene
gene <- makeGRangesFromDataFrame(ebi_anno_gtf_gene, keep.extra.columns = TRUE, 
                                 seqnames.field="chr", start.field="start", 
                                 end.field="end", strand.field="strand")

gene %>% as.data.frame() %>% dim()

# subsetByOverlaps(selection_interval,gene)
# subsetByOverlaps(gene,selection_interval)

ebi_overlap_selection_gene_df <- mergeByOverlaps(gene, selection_interval)

ebi_overlap_selection_gene_as_df <- as.data.frame(ebi_overlap_selection_gene_df)

ebi_overlap_selection_gene_as_df <- ebi_overlap_selection_gene_as_df %>% select(gene.seqnames:gene.attributes,selection_interval.seqnames:selection_interval.tau_proportion_mean)

#separate the attributes column
ebi_overlap_selection_gene_as_df <- ebi_overlap_selection_gene_as_df %>% 
  separate(gene.attributes, c("gene_id","gene_version","gene_name","gene_source","gene_biotype"),sep = ";",extra = "drop", fill = "right")

#tidy up the values in gene_id, gene_name and gene_biotype
ebi_overlap_selection_gene_as_df$gene_id <- gsub("gene_id \"","",ebi_overlap_selection_gene_as_df$gene_id)
ebi_overlap_selection_gene_as_df$gene_id <- gsub("\\\"","",ebi_overlap_selection_gene_as_df$gene_id)

ebi_overlap_selection_gene_as_df$gene_biotype <- gsub("gene_biotype \"","",ebi_overlap_selection_gene_as_df$gene_biotype)
ebi_overlap_selection_gene_as_df$gene_biotype <- gsub("\\\"","",ebi_overlap_selection_gene_as_df$gene_biotype)
ebi_overlap_selection_gene_as_df$gene_biotype <- gsub("^ ","",ebi_overlap_selection_gene_as_df$gene_biotype)

ebi_overlap_selection_gene_as_df$gene_source <- gsub("gene_biotype \"","",ebi_overlap_selection_gene_as_df$gene_source)
ebi_overlap_selection_gene_as_df$gene_source <- gsub("\\\"","",ebi_overlap_selection_gene_as_df$gene_source)
ebi_overlap_selection_gene_as_df$gene_source <- gsub("^ ","",ebi_overlap_selection_gene_as_df$gene_source)

ebi_overlap_selection_gene_as_df$gene_name <- gsub("gene_name \"","",ebi_overlap_selection_gene_as_df$gene_name)
ebi_overlap_selection_gene_as_df$gene_name <- gsub("\\\"","",ebi_overlap_selection_gene_as_df$gene_name)

#if loop to check gene_source contains biotype and then modify biotype
for (i in 1:nrow(ebi_overlap_selection_gene_as_df)){
  if (ebi_overlap_selection_gene_as_df$gene_source[i] != "gene_source ensembl"){
    ebi_overlap_selection_gene_as_df$gene_biotype[i] <- ebi_overlap_selection_gene_as_df$gene_source[i]
  ebi_overlap_selection_gene_as_df$gene_name[i] <- "not available"
  }
}

#get rid of _ in gene_biotype
ebi_overlap_selection_gene_as_df$gene_biotype <- gsub("_"," ",ebi_overlap_selection_gene_as_df$gene_biotype)

#write out dataframe as csv
ebi_overlap_selection_gene_as_df <- ebi_overlap_selection_gene_as_df %>% 
  select(gene.seqnames:gene_id,gene_name,gene_biotype,selection_interval.start:selection_interval.end,selection_interval.br_proportion_mean,selection_interval.snp_counted:selection_interval.tau_proportion_mean) %>%
  arrange(gene.seqnames)

ebi_overlap_selection_gene_as_df$Brahman_Query <- paste0(ebi_overlap_selection_gene_as_df$gene.seqnames,":",
                                                         ebi_overlap_selection_gene_as_df$selection_interval.start,"-",
                                                         ebi_overlap_selection_gene_as_df$selection_interval.end)

#ebi_overlap_100k_selection_gene_as_df_SevenBreeds_v2.csv is the 2nd version with tidied gene attributes
write.csv(ebi_overlap_selection_gene_as_df,file=paste0(dir1,"ebi_overlap_100k_selection_gene_as_df_SevenBreeds_v2_withGeneName.csv"),
          row.names = FALSE)

############################# NCBI ###########################
# #Brahman NCBI anno #Initial analysis was done with this
# ncbi_anno_gff3 <- read.gff("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/NCBI/Brahman/GCF_003369695.1_UOA_Brahman_1_genomic.gff")
# 
# #rename chr name
# ncbi_anno_gff3$seqid <- as.character(ncbi_anno_gff3$seqid)
# 
# for (i in 40076:(40076+28)){
#   char <- paste0("NC_0",as.character(i),".1")
#   ncbi_anno_gff3$seqid[ncbi_anno_gff3$seqid == char] <- as.character(i-40075)
# }
# 
# ncbi_anno_gff3$seqid[ncbi_anno_gff3$seqid == "NC_040105.1"] <- "X"
# 
# #remove the unplaced bcos not studying them
# ncbi_anno_gff3 <- ncbi_anno_gff3 %>% filter(!str_detect(seqid, "NW"))
# 
# ncbi_anno_gff3_gene <- ncbi_anno_gff3 %>% filter(type == "gene") %>% select(seqid,start,end,strand,
#                                                                             attributes)
# names(ncbi_anno_gff3_gene)[1] <- "chr"
# 
# #turn positively selected intervals and annotation into granges
# #selection interval
# selection_interval <- makeGRangesFromDataFrame(proportion_diff_df_summ, keep.extra.columns = TRUE, 
#                                                seqnames.field="chr", start.field="start", 
#                                                end.field="end", strand.field="strand")
# 
# #annotation - gene
# gene <- makeGRangesFromDataFrame(ncbi_anno_gff3_gene, keep.extra.columns = TRUE, 
#                                  seqnames.field="chr", start.field="start", 
#                                  end.field="end", strand.field="strand")
# 
# gene %>% as.data.frame() %>% dim()
# 
# subsetByOverlaps(selection_interval,gene)
# subsetByOverlaps(gene,selection_interval)
# 
# overlap_selection_gene_df <- mergeByOverlaps(gene, selection_interval)
# 
# overlap_selection_gene_as_df <- as.data.frame(overlap_selection_gene_df)
# 
# overlap_selection_gene_as_df <- overlap_selection_gene_as_df %>% select(gene.seqnames:gene.attributes,selection_interval.seqnames:selection_interval.tau_proportion_mean)
# 
# #separate the attributes column
# overlap_selection_gene_as_df <- overlap_selection_gene_as_df %>% 
#   separate(gene.attributes, c("ID","Dbxref","Name","gbkey","gene","gene_biotype"),
#            sep = ";",extra = "drop", fill = "right")
# 
# #tidy up the values in gene_id and gene_biotype
# overlap_selection_gene_as_df$ID <- gsub("ID=","",overlap_selection_gene_as_df$ID)
# overlap_selection_gene_as_df$Dbxref <- gsub("Dbxref=GeneID:","",overlap_selection_gene_as_df$Dbxref)
# overlap_selection_gene_as_df$Name <- gsub("Name=","",overlap_selection_gene_as_df$Name)
# overlap_selection_gene_as_df$gbkey <- gsub("gbkey=","",overlap_selection_gene_as_df$gbkey)
# overlap_selection_gene_as_df$gene <- gsub("gene=","",overlap_selection_gene_as_df$gene)
# overlap_selection_gene_as_df$gene_biotype <- gsub("gene_biotype=","",overlap_selection_gene_as_df$gene_biotype)
# 
# #write out dataframe as csv
# overlap_selection_gene_as_df <- overlap_selection_gene_as_df %>% 
#   select(gene.seqnames:Name,gene_biotype,selection_interval.start:selection_interval.end,selection_interval.br_proportion_mean,selection_interval.snp_counted:selection_interval.tau_proportion_mean) %>%
#   arrange(gene.seqnames)
# 
# #ncbi_overlap_100k_selection_gene_as_df_SevenBreeds_v2.csv is the 2rd version. This is because ensembl anno missed a lot
# write.csv(overlap_selection_gene_as_df,file=paste0(dir1,"ncbi_overlap_100k_selection_gene_as_df_SevenBreeds_v2.csv"),
#           row.names = FALSE)

##### Overlap with UMD3.1 QTL #####
#read in qtl umd3.1
qtl_gff3 <- read.gff("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/QTLs/QTL_UMD_3.1.gff")

#rename chr name
qtl_gff3$seqid <- as.character(qtl_gff3$seqid)

qtl_gff3$seqid <- gsub("Chr.","",qtl_gff3$seqid)

names(qtl_gff3)[1] <- "chr"

#work with complete cases in start and end
compl_cases <- complete.cases(qtl_gff3$start,qtl_gff3$end)
qtl_gff3 <- qtl_gff3[compl_cases,]

#need a strand
qtl_gff3$strand <- rep("*",length(qtl_gff3$strand))
qtl_gff3 <- as.tbl(qtl_gff3)

#read is mashmap brahman equivalent umd3.1 region
brahman_mashmap_umd3_1_1 <- read_delim("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/QTLs/sweep/script/all_mashmap",
             " ", col_names = FALSE)
names(brahman_mashmap_umd3_1_1) <- c("Query_name","Query_length","Query_start","Query_end","Orientation",
                                                     "Ref_name","Ref_length","Ref_start","Ref_end","Percentid")

#remove dup
brahman_mashmap_umd3_1_1_modi <- brahman_mashmap_umd3_1_1 %>% distinct(Query_name,.keep_all = TRUE)

#write it out to try my awk script
# write_tsv(brahman_mashmap_umd3_1_1_modi,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/QTLs/sweep/script/brahman_mashmap_umd3_1_1_modi.tsv")

brahman_mashmap_umd3_1_1_modi <- brahman_mashmap_umd3_1_1_modi %>%
  arrange(Ref_name,Ref_start) %>% mutate(alignedL = Ref_end - Ref_start +1) %>%
  mutate(alignProportion = alignedL/Query_length)

brahman_mashmap_umd3_1_1_modi$Ref_name <- gsub("chr","",brahman_mashmap_umd3_1_1_modi$Ref_name)

#need a strand
brahman_mashmap_umd3_1_1_modi$strand <- rep("*",length(brahman_mashmap_umd3_1_1_modi$Ref_name))

#turn brahman equivalent umd3.1 region and QTL
#brahman equivalent umd3.1 region
brahman_equivalent_umd3.1_GR <- makeGRangesFromDataFrame(brahman_mashmap_umd3_1_1_modi, keep.extra.columns = TRUE,
                                               seqnames.field="Ref_name", start.field="Ref_start",
                                               end.field="Ref_end", strand.field="strand")

#QTL
QTL_GR <- makeGRangesFromDataFrame(qtl_gff3, keep.extra.columns = TRUE,
                                 seqnames.field="chr", start.field="start",
                                 end.field="end",strand.field="strand")

#overlap selection interval transformed to umd3.1.1 with QTL
overlap_selection_QTL_df <- mergeByOverlaps(brahman_equivalent_umd3.1_GR,QTL_GR)

overlap_selection_QTL_as_df <- as.data.frame(overlap_selection_QTL_df) #this will only work with tible was used for merge overlap

overlap_selection_QTL_as_df <- overlap_selection_QTL_as_df %>% 
  select(brahman_equivalent_umd3.1_GR.seqnames:brahman_equivalent_umd3.1_GR.end,Query_name,QTL_GR.seqnames:QTL_GR.end,QTL_GR.type,QTL_GR.attributes)

overlap_selection_QTL_as_df <- overlap_selection_QTL_as_df %>% 
  separate(QTL_GR.attributes, c("QTL_ID","Name","Abbrev","PUBMED_ID","trait_ID","trait","breed","FlankMarkers"),
           sep = ";",extra = "drop", fill = "right",remove = FALSE)

overlap_selection_QTL_as_df$Name <- gsub("Name=","",overlap_selection_QTL_as_df$Name)
overlap_selection_QTL_as_df$Abbrev <- gsub("Abbrev=","",overlap_selection_QTL_as_df$Abbrev)
overlap_selection_QTL_as_df$PUBMED_ID <- gsub("PUBMED_ID=","",overlap_selection_QTL_as_df$PUBMED_ID)
overlap_selection_QTL_as_df$trait_ID <- gsub("trait_ID=","",overlap_selection_QTL_as_df$trait_ID)
overlap_selection_QTL_as_df$trait <- gsub("trait=","",overlap_selection_QTL_as_df$trait)
overlap_selection_QTL_as_df$breed <- gsub("breed=","",overlap_selection_QTL_as_df$breed)
overlap_selection_QTL_as_df$QTL_ID <- gsub("QTL_ID=","",overlap_selection_QTL_as_df$QTL_ID)

#get rid of bad breed names 208,230, get rid of flank markers, arrange by chr
#rename brahman_equivalent_umd3.1_GR.seqnames,start,end
overlap_selection_QTL_as_df$breed[c(208,230)] <- c("NA","NA")
overlap_selection_QTL_as_df <- overlap_selection_QTL_as_df %>% 
  select(-FlankMarkers) %>% arrange(brahman_equivalent_umd3.1_GR.seqnames)
names(overlap_selection_QTL_as_df)[c(1,2,3,4)] <- c("UMD3.1_equivalent_chr","UMD3.1_equivalent_start","UMD3.1_equivalent_end","Brahman_Query")

#write out dataframe as csv
write.csv(overlap_selection_QTL_as_df,file=paste0(dir1,"overlap_selection_QTL_as_df.csv"),
          row.names = FALSE)

#merging selective sweep interval with QTL interval
sweep_QTL_combined <- merge(ebi_overlap_selection_gene_as_df,overlap_selection_QTL_as_df,by="Brahman_Query")

sweep_QTL_combined <- sweep_QTL_combined %>% arrange(gene.seqnames)

#write out dataframe as csv
write.csv(sweep_QTL_combined,file=paste0(dir1,"sweep_QTL_combined.csv"),
          row.names = FALSE)

#breakdown of biotype genes in sweep intervals
ebi_overlap_selection_gene_as_df %>% distinct(gene_id, .keep_all = TRUE) %>% 
  group_by(gene_biotype) %>% summarise(count = n()) %>% mutate(proportion = count/sum(count))
# A tibble: 8 x 3
# gene_biotype   count proportion
# <chr>          <int>      <dbl>
#   1 lncRNA            11    0.0859 
# 2 miRNA              2    0.0156 
# 3 misc RNA           2    0.0156 
# 4 protein coding   102    0.797  
# 5 pseudogene         1    0.00781
# 6 rRNA               2    0.0156 
# 7 snoRNA             4    0.0312 
# 8 snRNA              4    0.0312 

length(unique(ebi_overlap_selection_gene_as_df$Brahman_Query))
#60

length(unique(overlap_selection_QTL_as_df$Brahman_Query))
#53

base::setdiff(unique(ebi_overlap_selection_gene_as_df$Brahman_Query),unique(overlap_selection_QTL_as_df$Brahman_Query))
base::intersect(unique(ebi_overlap_selection_gene_as_df$Brahman_Query),unique(overlap_selection_QTL_as_df$Brahman_Query))
base::union(unique(ebi_overlap_selection_gene_as_df$Brahman_Query),unique(overlap_selection_QTL_as_df$Brahman_Query))
base::setdiff(unique(overlap_selection_QTL_as_df$Brahman_Query),unique(ebi_overlap_selection_gene_as_df$Brahman_Query))

#understanding the overlap of selection interval and qtl interval
library(VennDiagram)
venn.diagram(
  x = list(unique(overlap_selection_QTL_as_df$Brahman_Query), unique(ebi_overlap_selection_gene_as_df$Brahman_Query)),
  cat.pos = 0,category.names = c("QTL" , "Selection"),fill = c('turquoise', 'red'),cex=2,cat.cex=1.5,
  filename = 'QTLVsSelection_venn.png')
