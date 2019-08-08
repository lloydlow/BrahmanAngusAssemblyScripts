#!/usr/bin/env Rscript
#------------------------------------------------------
# Program name: brahman_angus_vcf_score_alt_allele_AngusRef_Rscript.R
# Objective: from vcf annotated with annovar, count
#           alt alleles
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(readr)
library(dplyr)
library(tidyr)

args <- commandArgs(TRUE)
vcfFile <- args[1]

#get chr no
chr <- gsub("combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_","",vcfFile)
chr <- gsub(".vcf","",chr)

# path to annotated vcf
dir1 <- 
  "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus/"

# reading combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1.vcf
# path1 <- paste0(dir1,"combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1.vcf")

path1 <- paste0(dir1,vcfFile)

combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1 <- 
  read_tsv(path1,col_names = FALSE)
names(combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1) <- 
  c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","an1","an2","an3","an4","an5",
    "an6","br1","br2","br3","br4","br5")

#split based on an1-an5 and br1-br6
combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split <- combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1 %>% 
  separate(an1, c("an1_GT","an1_AD","an1_DP","an1_GQ","an1_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(an2, c("an2_GT","an2_AD","an2_DP","an2_GQ","an2_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(an3, c("an3_GT","an3_AD","an3_DP","an3_GQ","an3_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(an4, c("an4_GT","an4_AD","an4_DP","an4_GQ","an4_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(an5, c("an5_GT","an5_AD","an5_DP","an5_GQ","an5_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(an6, c("an6_GT","an6_AD","an6_DP","an6_GQ","an6_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(br1, c("br1_GT","br1_AD","br1_DP","br1_GQ","br1_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(br2, c("br2_GT","br2_AD","br2_DP","br2_GQ","br2_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(br3, c("br3_GT","br3_AD","br3_DP","br3_GQ","br3_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(br4, c("br4_GT","br4_AD","br4_DP","br4_GQ","br4_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(br5, c("br5_GT","br5_AD","br5_DP","br5_GQ","br5_PL"),sep = ":",extra = "drop", fill = "right")

#complete calls for all genotypes
combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split <- combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split %>% 
  filter(an1_GT != "./.") %>%
  filter(an2_GT != "./.") %>%
  filter(an3_GT != "./.") %>%
  filter(an4_GT != "./.") %>%
  filter(an5_GT != "./.") %>%
  filter(an6_GT != "./.") %>%
  filter(br1_GT != "./.") %>%
  filter(br2_GT != "./.") %>%
  filter(br3_GT != "./.") %>%
  filter(br4_GT != "./.") %>%
  filter(br5_GT != "./.")

#filter for at least 5 reads mapped
combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split <- combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split %>% 
  filter(an1_DP > 4) %>%
  filter(an2_DP > 4) %>%
  filter(an3_DP > 4) %>%
  filter(an4_DP > 4) %>%
  filter(an5_DP > 4) %>%
  filter(an6_DP > 4) 
# filter(br1_DP > 4) %>%
# filter(br2_DP > 4) %>%
# filter(br3_DP > 4) %>%
# filter(br4_DP > 4) %>%
# filter(br5_DP > 4)

#split GT column into 2 (e.g. "an1_GT_A","an1_GT_B")
combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT <- combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split %>%
  separate(an1_GT, c("an1_GT_A","an1_GT_B"),sep = "/",convert = TRUE) %>%
  separate(an2_GT, c("an2_GT_A","an2_GT_B"),sep = "/",convert = TRUE) %>%
  separate(an3_GT, c("an3_GT_A","an3_GT_B"),sep = "/",convert = TRUE) %>%
  separate(an4_GT, c("an4_GT_A","an4_GT_B"),sep = "/",convert = TRUE) %>%
  separate(an5_GT, c("an5_GT_A","an5_GT_B"),sep = "/",convert = TRUE) %>%
  separate(an6_GT, c("an6_GT_A","an6_GT_B"),sep = "/",convert = TRUE) %>%
  separate(br1_GT, c("br1_GT_A","br1_GT_B"),sep = "/",convert = TRUE) %>%
  separate(br2_GT, c("br2_GT_A","br2_GT_B"),sep = "/",convert = TRUE) %>%
  separate(br3_GT, c("br3_GT_A","br3_GT_B"),sep = "/",convert = TRUE) %>%
  separate(br4_GT, c("br4_GT_A","br4_GT_B"),sep = "/",convert = TRUE) %>%
  separate(br5_GT, c("br5_GT_A","br5_GT_B"),sep = "/",convert = TRUE)

#convert non-zero to 1 for allele column
combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$an1_GT_A[combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$an1_GT_A != 0] <- 1
combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$an1_GT_B[combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$an1_GT_B != 0] <- 1
combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$an2_GT_A[combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$an2_GT_A != 0] <- 1
combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$an2_GT_B[combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$an2_GT_B != 0] <- 1
combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$an3_GT_A[combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$an3_GT_A != 0] <- 1
combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$an3_GT_B[combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$an3_GT_B != 0] <- 1
combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$an4_GT_A[combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$an4_GT_A != 0] <- 1
combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$an4_GT_B[combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$an4_GT_B != 0] <- 1
combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$an5_GT_A[combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$an5_GT_A != 0] <- 1
combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$an5_GT_B[combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$an5_GT_B != 0] <- 1
combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$an6_GT_A[combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$an6_GT_A != 0] <- 1
combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$an6_GT_B[combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$an6_GT_B != 0] <- 1
combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$br1_GT_A[combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$br1_GT_A != 0] <- 1
combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$br1_GT_B[combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$br1_GT_B != 0] <- 1
combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$br2_GT_A[combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$br2_GT_A != 0] <- 1
combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$br2_GT_B[combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$br2_GT_B != 0] <- 1
combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$br3_GT_A[combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$br3_GT_A != 0] <- 1
combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$br3_GT_B[combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$br3_GT_B != 0] <- 1
combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$br4_GT_A[combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$br4_GT_A != 0] <- 1
combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$br4_GT_B[combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$br4_GT_B != 0] <- 1
combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$br5_GT_A[combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$br5_GT_A != 0] <- 1
combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$br5_GT_B[combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT$br5_GT_B != 0] <- 1

#create new var for scoring number of ALT allele
combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT <- combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT %>%
  mutate(an1_GT_ALT_count = an1_GT_A + an1_GT_B) %>%
  mutate(an2_GT_ALT_count = an2_GT_A + an2_GT_B) %>%
  mutate(an3_GT_ALT_count = an3_GT_A + an3_GT_B) %>%
  mutate(an4_GT_ALT_count = an4_GT_A + an4_GT_B) %>%
  mutate(an5_GT_ALT_count = an5_GT_A + an5_GT_B) %>%
  mutate(an6_GT_ALT_count = an6_GT_A + an6_GT_B) %>%
  mutate(br1_GT_ALT_count = br1_GT_A + br1_GT_B) %>%
  mutate(br2_GT_ALT_count = br2_GT_A + br2_GT_B) %>%
  mutate(br3_GT_ALT_count = br3_GT_A + br3_GT_B) %>%
  mutate(br4_GT_ALT_count = br4_GT_A + br4_GT_B) %>%
  mutate(br5_GT_ALT_count = br5_GT_A + br5_GT_B)

#create breed proportion variable
combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro <- combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT %>%
  rowwise() %>%
  mutate(an_ALT_proportion = sum(an1_GT_ALT_count,an2_GT_ALT_count,an3_GT_ALT_count,an4_GT_ALT_count,an5_GT_ALT_count,an6_GT_ALT_count)/12,
         br_ALT_proportion = sum(br1_GT_ALT_count,br2_GT_ALT_count,br3_GT_ALT_count,br4_GT_ALT_count,br5_GT_ALT_count)/10)

#create proportion difference var, angus minus brahman
combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro <- combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro %>%
  mutate(proportion_diff = an_ALT_proportion - br_ALT_proportion)

#create interval to summarise proportion
combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro$interval <- 
  with(combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro, cut(POS, seq(1, max(POS), by = 100e3)))

#use proportion per breed instead of proportion diff
proportion_diff_df <- combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro %>% 
  group_by(interval) %>% summarise(an_proportion_mean = mean(an_ALT_proportion),
                                   br_proportion_mean= mean(br_ALT_proportion),
                                   proportion_diff_mean = mean(proportion_diff),
                                   snp_counted = n())

proportion_diff_df$interval_no <- 1:nrow(proportion_diff_df)

#simple plot to detect any "big" difference
#blue is brahman, black is angus!
plot(proportion_diff_df$interval_no,proportion_diff_df$an_proportion_mean,type="l",col="black",
     xlab = "Window of 100kb interval",ylab = "Proportion of ALT allele")
lines(proportion_diff_df$interval_no,proportion_diff_df$br_proportion_mean,col="blue")

#any consecutive low Brahman ALT allele ref?
proportion_diff_df_summ <- proportion_diff_df %>% filter(an_proportion_mean <= 0.1) %>%
  filter(br_proportion_mean > 0.5) %>% filter(snp_counted >= 10)

proportion_diff_df_summ$chr <- chr

#write out tables per chr
#final filtered vcf file with new columns from my split
assign(paste('combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_', chr, sep=''), 
       combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro)

path2 <- paste0(dir1,'combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_',chr,".RData")

save(list = paste0('combinedBreedRefAngus_filteredSNP.Angus_multianno_subset_1_split_GT_pro_', chr),file = 
       path2)

#proportion_diff_df_summ of monomorphic brahman #not used so commented out below
# assign(paste('proportion_diff_df_summ_', chr, sep=''), 
#        proportion_diff_df_summ)
# 
# path3 <- paste0(dir1,'proportion_diff_df_summ_',chr,".RData")
# 
# save(list = paste0('proportion_diff_df_summ_',chr),file = 
#        path3)
