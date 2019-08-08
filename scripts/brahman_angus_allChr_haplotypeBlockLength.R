#------------------------------------------------------
# Program name: brahman_angus_allChr_haplotypeBlockLength.R
# Objective: load all chr for analysis and 
#   created a function haploLength to look at consecutive 
#   homozygous run
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(readr)
library(dplyr)
library(tidyr)

#brahman
# path to filtered vcf where only wanted SNPs per chr are studied
dir1 <- 
  "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/"

# load datasets
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_1.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_10.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_11.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_12.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_13.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_14.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_15.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_16.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_17.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_18.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_19.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_20.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_21.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_22.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_23.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_24.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_25.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_26.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_27.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_28.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_29.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_4.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_5.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_6.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_7.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_8.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_9.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_X.RData")


df.list <- list(combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_1,
                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_10,
                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_11,
                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_12,
                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_13,
                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_14,
                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_15,
                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_16,
                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_17,
                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_18,
                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_19,
                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2,
                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_20,
                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_21,
                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_22,
                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_23,
                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_24,
                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_25,
                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_26,
                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_27,
                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_28,
                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_29,
                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3,
                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_4,
                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_5,
                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_6,
                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_7,
                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_8,
                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_9,
                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_X)

#loop thro df.list
haplo_list <- lapply(df.list,haploLength)

haplo_vec <- unlist(haplo_list, recursive = TRUE)

summary(haplo_vec)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1    1627    6946   14344   17730 1896803 

#functions
#turn storing haplotype block length into a function
haploLength <- function(df){
  #implement searching for consecutive homozygous snps
  df <- df %>%
    mutate(all_an_alt = an1_GT_ALT_count + an2_GT_ALT_count + an3_GT_ALT_count + 
             an4_GT_ALT_count + an5_GT_ALT_count + an6_GT_ALT_count)
  
  #loop thro df and save POS of consecutive 0 (i.e. homozygous)
  #make rle object
  runs = rle(df$all_an_alt == 0)
  
  #indices of runs with length more than 1
  myruns = which(runs$values == TRUE & runs$lengths > 1)
  
  #extract end positions
  runs.lengths.cumsum = cumsum(runs$lengths)
  ends = runs.lengths.cumsum[myruns]
  
  #find start position
  newindex = ifelse(myruns>1, myruns-1, 0)
  starts = runs.lengths.cumsum[newindex] + 1
  if (0 %in% newindex) starts = c(1,starts)
  
  #store length of shared haplotype blocks
  vec <- c()
  
  for (i in 1:length(ends)){
    block_length <- df$POS[ends[i]] -
      df$POS[starts[i]]
    vec <- c(vec,block_length)
  }
  
  return(vec)
}
