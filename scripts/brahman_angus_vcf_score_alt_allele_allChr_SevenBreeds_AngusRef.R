#------------------------------------------------------
# Program name: brahman_angus_vcf_score_alt_allele_allChr_SevenBreeds_AngusRef.R
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

#ref is angus here
# path to filtered vcf where only wanted SNPs per chr are studied
dir1 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus_SevenBreeds/"

# load datasets
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_1.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_10.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_11.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_12.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_13.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_14.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_15.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_16.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_17.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_18.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_19.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_2.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_20.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_21.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_22.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_23.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_24.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_25.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_26.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_27.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_28.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_29.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_3.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_4.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_5.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_6.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_7.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_8.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_9.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_Y.RData")

#combine all 30 chr final SNPs file into one
merged <- do.call("rbind", list(combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_1,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_10,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_11,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_12,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_13,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_14,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_15,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_16,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_17,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_18,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_19,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_2,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_20,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_21,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_22,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_23,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_24,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_25,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_26,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_27,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_28,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_29,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_3,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_4,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_5,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_6,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_7,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_8,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_9,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3_Y))

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

#all non-angus combine mean proportion
proportion_diff_df_summ <- proportion_diff_df_summ %>% rowwise() %>% 
  mutate(nonangus_proportion_mean = sum(ge_proportion_mean,he_proportion_mean,
                                   re_proportion_mean,sh_proportion_mean,si_proportion_mean)/5)

summary(proportion_diff_df_summ$mean_diff_POS)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 119.8  3613.4  4743.2  5034.9  6201.8 11060.1 

nrow(proportion_diff_df_summ)/nrow(proportion_diff_df)

#nonangus percentile
quantile(proportion_diff_df_summ$nonangus_proportion_mean, probs = seq(0,1,0.05))
# 0%         5%        10%        15%        20%        25%        30%        35%        40%        45%        50% 
# 0.01623932 0.05098915 0.06051587 0.06785714 0.07458333 0.08055556 0.08680556 0.09249557 0.09823232 0.10416667 0.11023392 
# 55%        60%        65%        70%        75%        80%        85%        90%        95%       100% 
#   0.11710526 0.12395833 0.13163580 0.13976559 0.15000000 0.16208709 0.17722222 0.19745532 0.23379430 0.57350427 

#angus percentile
quantile(proportion_diff_df_summ$an_proportion_mean, probs = seq(0,1,0.05))
# 0%          5%         10%         15%         20%         25%         30%         35%         40%         45% 
#   0.000000000 0.005555556 0.013888889 0.020833333 0.027777778 0.034552846 0.041666667 0.047619048 0.053921569 0.061111111 
# 50%         55%         60%         65%         70%         75%         80%         85%         90%         95% 
#   0.068627451 0.076388889 0.083333333 0.094002016 0.104166667 0.116666667 0.131944444 0.151041667 0.177083333 0.219696970 
# 100% 
# 0.600000000 

#so choose from  0.006 (in angus) to  0.2 (in nonangus) 2nd window
#so choose from  0.01 (in angus) to  0.17 (in nonangus) 3rd window

#histogram
tiff(filename = "hist_Angus_vs_nonAngus_SevenBreeds_proportion_mean.tiff",width = 800)
par(mfrow=c(1,2))
hist(proportion_diff_df_summ$an_proportion_mean,xlab="Proportion mean of alternate allele in Angus",
     main="",ylim=c(0,8000))
abline(v=0.006,lty="dashed",col="red")
arrows(0.006,4000,0,4000,lwd=0.75,col="red",length=0.05)

hist(proportion_diff_df_summ$nonangus_proportion_mean,xlab="Proportion mean of alternate allele in other five taurine breeds",
     main="",ylim=c(0,8000))
abline(v=0.2,lty="dashed",col="red")
arrows(0.2,4000,0.55,4000,lwd=0.75,col="red",length=0.05)
par(mfrow=c(1,1))
dev.off()

proportion_diff_df_summ <- proportion_diff_df_summ %>% filter(an_proportion_mean <= 0.01) %>%
  filter(nonangus_proportion_mean >= 0.17) 

summary(proportion_diff_df_summ$nonangus_proportion_mean)

summary(proportion_diff_df_summ$an_proportion_mean)

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
