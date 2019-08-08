#------------------------------------------------------
# Program name: brahman_angus_Fst.R
# Objective: Analyse Fst tables from vcftools
#   
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(readr)
library(dplyr)
library(tidyr)

#Angus ref
# path to taurine_vs_brahman_Wind_RefAngus.weir.fst
dir1 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/FST/SevenBreedsAngus/"

# reading taurine_vs_brahman_Wind_RefAngus.windowed.weir.fst
path1 <- paste0(dir1,"taurine_vs_brahman_Wind_RefAngus.windowed.weir.fst")

taurine_vs_brahman_Wind_RefAngus_fst <- read_tsv(path1,col_names = TRUE, col_types = "ciiidd")

taurine_vs_brahman_Wind_RefAngus_fst_chr7 <- taurine_vs_brahman_Wind_RefAngus_fst %>% filter(CHROM == "7") %>% 
  filter(BIN_START > 66100000-1e6) %>% filter(BIN_START < 66200000+1e6)

plot(taurine_vs_brahman_Wind_RefAngus_fst_chr7$BIN_START,
     taurine_vs_brahman_Wind_RefAngus_fst_chr7$MEAN_FST,type="l")
abline(v=c(66133220,66187098))

plot(taurine_vs_brahman_Wind_RefAngus_fst_chr7$BIN_START,
     taurine_vs_brahman_Wind_RefAngus_fst_chr7$WEIGHTED_FST,type="l")
abline(v=c(66133220,66187098))

#Angus ref
#path to taurine_vs_brahman_RefAngus.weir.fst
dir1 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/FST/SevenBreedsAngus/"

# reading taurine_vs_brahman_RefAngus.weir.fst
path1 <- paste0(dir1,"taurine_vs_brahman_RefAngus.weir.fst")

taurine_vs_brahman_RefAngus_fst <- read_tsv(path1,col_names = TRUE, col_types = "cid")

taurine_vs_brahman_RefAngus_fst_chr7 <- taurine_vs_brahman_RefAngus_fst %>% filter(CHROM == "7") %>% 
  filter(POS > 66100000) %>% filter(POS < 66200000) %>%
  filter(WEIR_AND_COCKERHAM_FST != "NaN") 

plot(taurine_vs_brahman_RefAngus_fst_chr7$POS,taurine_vs_brahman_RefAngus_fst_chr7$WEIR_AND_COCKERHAM_FST,
     type="l")
abline(v=c(66133220,66187098))

#############################################################
#Brahman ref
# path to angus_vs_brahman.windowed.weir.fst
dir1 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/FST/SevenBreedsBrahman/"

# reading bostaurus_angus_bionano_NCBI_full_corrected_gapfill_arrow_fil_chronly.coor
path1 <- paste0(dir1,"angus_vs_brahman.windowed.weir.fst")

angus_vs_brahman_fst <- read_tsv(path1,col_names = TRUE, col_types = "ciiidd")

angus_vs_brahman_fst_chr7 <- angus_vs_brahman_fst %>% filter(CHROM == "7") %>% 
  filter(BIN_START > 66100000-1e6) %>% filter(BIN_START < 66200000+1e6)

plot(angus_vs_brahman_fst_chr7$BIN_START,angus_vs_brahman_fst_chr7$MEAN_FST,type="l")
abline(v=c(66133220,66187098))

plot(angus_vs_brahman_fst_chr7$BIN_START,angus_vs_brahman_fst_chr7$WEIGHTED_FST,type="l")
abline(v=c(66133220,66187098))

#Brahman ref
# path to taurine_vs_brahman_Wind.windowed.weir_RefBrahman.fst
dir1 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/FST/SevenBreedsBrahman/"

# reading bostaurus_angus_bionano_NCBI_full_corrected_gapfill_arrow_fil_chronly.coor
path1 <- paste0(dir1,"taurine_vs_brahman_Wind.windowed.weir_RefBrahman.fst")

taurine_vs_brahman_Wind_fst <- read_tsv(path1,col_names = TRUE, col_types = "ciiidd")

taurine_vs_brahman_Wind_fst_chr7 <- taurine_vs_brahman_Wind_fst %>% filter(CHROM == "7") %>% 
  filter(BIN_START > 66100000-1e6) %>% filter(BIN_START < 66200000+1e6)

plot(taurine_vs_brahman_Wind_fst_chr7$BIN_START,taurine_vs_brahman_Wind_fst_chr7$MEAN_FST,type="l")
abline(v=c(66133220,66187098))

plot(taurine_vs_brahman_Wind_fst_chr7$BIN_START,taurine_vs_brahman_Wind_fst_chr7$WEIGHTED_FST,type="l")
abline(v=c(66133220,66187098))

#Brahman ref
#path to taurine_vs_brahman.weir_RefBrahman.fst
dir1 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/FST/SevenBreedsBrahman/"

# reading taurine_vs_brahman.weir_RefBrahman.fst
path1 <- paste0(dir1,"taurine_vs_brahman.weir_RefBrahman.fst")

taurine_vs_brahman_RefBrahman_fst <- read_tsv(path1,col_names = TRUE, col_types = "cid")

taurine_vs_brahman_RefBrahman_fst_chr7 <- taurine_vs_brahman_RefBrahman_fst %>% filter(CHROM == "7") %>% 
  filter(POS > 66100000) %>% filter(POS < 66200000) %>%
  filter(WEIR_AND_COCKERHAM_FST != "NaN") 

plot(taurine_vs_brahman_RefBrahman_fst_chr7$POS,taurine_vs_brahman_RefBrahman_fst_chr7$WEIR_AND_COCKERHAM_FST,
     type="l")
abline(v=c(66133220,66187098))
