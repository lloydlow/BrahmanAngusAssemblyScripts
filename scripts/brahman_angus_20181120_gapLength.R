#------------------------------------------------------
# Program name: brahman_angus_20181120_gapLength.R
# Objective: summarise assembly statistics from the
#           sequential improvement using gapLength output
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

########## checking sequential improvement in assembly versions #############

#angus

library(readr)
library(dplyr)
library(ggplot2)

#calculate N50 function
N50 <- function(assembly_df){
  vect <- c()
  for (i in 1:nrow(assembly_df)){
    N50length <- 0.5*sum(as.numeric(assembly_df$length))
    vect <- c(vect,assembly_df$length[i])
    cumsumlength <- sum(as.numeric(vect))
    if (cumsumlength >= N50length){
      print(assembly_df$length[i])
      break
    } 
  }
}

# path to all gapLength results
dir1 <- 
  "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/fasta_file_order_checkBases_N/"

# reading bostaurus_angus_No_Ns.rls
path1 <- paste0(dir1,"bostaurus_angus_No_Ns.rls")

bostaurus_angus_No_Ns <- read_tsv(path1,col_names = FALSE)
names(bostaurus_angus_No_Ns) <- c("scaffold","nothing","gap","length","perc_gap")

bostaurus_angus_No_Ns <- bostaurus_angus_No_Ns %>%
  dplyr::select(scaffold,gap,length) %>% arrange(desc(length))

sum(as.numeric(bostaurus_angus_No_Ns$length)) - 
  sum(as.numeric(bostaurus_angus_No_Ns$gap))

N50(bostaurus_angus_No_Ns)
nrow(bostaurus_angus_No_Ns)
sum(as.numeric(bostaurus_angus_No_Ns$length))

# reading f1_sire_salsa_No_Ns.rls
path1 <- paste0(dir1,"f1_sire_salsa_No_Ns.rls")

f1_sire_salsa_No_Ns <- read_tsv(path1,col_names = FALSE)
names(f1_sire_salsa_No_Ns) <- c("scaffold","nothing","gap","length","perc_gap")

f1_sire_salsa_No_Ns <- f1_sire_salsa_No_Ns %>%
  dplyr::select(scaffold,gap,length) %>% arrange(desc(length))

sum(as.numeric(f1_sire_salsa_No_Ns$length)) - 
  sum(as.numeric(f1_sire_salsa_No_Ns$gap))

N50(f1_sire_salsa_No_Ns)
nrow(f1_sire_salsa_No_Ns)
sum(as.numeric(f1_sire_salsa_No_Ns$length))

# reading bostaurus_angus_bionano_NCBI_full_No_Ns.rls
path1 <- paste0(dir1,"bostaurus_angus_bionano_NCBI_full_No_Ns.rls")

bostaurus_angus_bionano_NCBI_full_No_Ns <- read_tsv(path1,col_names = FALSE)
names(bostaurus_angus_bionano_NCBI_full_No_Ns) <- c("scaffold","nothing","gap","length","perc_gap")

bostaurus_angus_bionano_NCBI_full_No_Ns <- bostaurus_angus_bionano_NCBI_full_No_Ns %>%
  dplyr::select(scaffold,gap,length) %>% arrange(desc(length))

sum(as.numeric(bostaurus_angus_bionano_NCBI_full_No_Ns$length)) - 
  sum(as.numeric(bostaurus_angus_bionano_NCBI_full_No_Ns$gap))

N50(bostaurus_angus_bionano_NCBI_full_No_Ns)
nrow(bostaurus_angus_bionano_NCBI_full_No_Ns)
sum(as.numeric(bostaurus_angus_bionano_NCBI_full_No_Ns$length))

# reading bostaurus_angus_bionano_NCBI_full_corrected_gapfill_arrow_fil_chronly_No_Ns.rls
path1 <- paste0(dir1,"bostaurus_angus_bionano_NCBI_full_corrected_gapfill_arrow_fil_chronly_No_Ns.rls")

bostaurus_angus_bionano_NCBI_full_corrected_gapfill_arrow_fil_chronly_No_Ns <- read_tsv(path1,col_names = FALSE)
names(bostaurus_angus_bionano_NCBI_full_corrected_gapfill_arrow_fil_chronly_No_Ns) <- c("scaffold","nothing","gap","length","perc_gap")

bostaurus_angus_bionano_NCBI_full_corrected_gapfill_arrow_fil_chronly_No_Ns <- bostaurus_angus_bionano_NCBI_full_corrected_gapfill_arrow_fil_chronly_No_Ns %>%
  dplyr::select(scaffold,gap,length) %>% arrange(desc(length))

sum(as.numeric(bostaurus_angus_bionano_NCBI_full_corrected_gapfill_arrow_fil_chronly_No_Ns$length)) - 
  sum(as.numeric(bostaurus_angus_bionano_NCBI_full_corrected_gapfill_arrow_fil_chronly_No_Ns$gap))

N50(bostaurus_angus_bionano_NCBI_full_corrected_gapfill_arrow_fil_chronly_No_Ns)
nrow(bostaurus_angus_bionano_NCBI_full_corrected_gapfill_arrow_fil_chronly_No_Ns)
sum(as.numeric(bostaurus_angus_bionano_NCBI_full_corrected_gapfill_arrow_fil_chronly_No_Ns$length))

#brahman
# reading bostaurus_brahma_No_Ns.rls
path1 <- paste0(dir1,"bostaurus_brahma_No_Ns.rls")

bostaurus_brahma_No_Ns <- read_tsv(path1,col_names = FALSE)
names(bostaurus_brahma_No_Ns) <- c("scaffold","nothing","gap","length","perc_gap")

bostaurus_brahma_No_Ns <- bostaurus_brahma_No_Ns %>%
  dplyr::select(scaffold,gap,length) %>% arrange(desc(length))

sum(as.numeric(bostaurus_brahma_No_Ns$length)) - 
  sum(as.numeric(bostaurus_brahma_No_Ns$gap))

N50(bostaurus_brahma_No_Ns)
nrow(bostaurus_brahma_No_Ns)
sum(as.numeric(bostaurus_brahma_No_Ns$length))

# reading f1_dam_salsa_No_Ns.rls
path1 <- paste0(dir1,"f1_dam_salsa_No_Ns.rls")

f1_dam_salsa_No_Ns <- read_tsv(path1,col_names = FALSE)
names(f1_dam_salsa_No_Ns) <- c("scaffold","nothing","gap","length","perc_gap")

f1_dam_salsa_No_Ns <- f1_dam_salsa_No_Ns %>%
  dplyr::select(scaffold,gap,length) %>% arrange(desc(length))

sum(as.numeric(f1_dam_salsa_No_Ns$length)) - 
  sum(as.numeric(f1_dam_salsa_No_Ns$gap))

N50(f1_dam_salsa_No_Ns)
nrow(f1_dam_salsa_No_Ns)
sum(as.numeric(f1_dam_salsa_No_Ns$length))

# reading bostaurus_brahma_bionano_NCBI_full_No_Ns.rls
path1 <- paste0(dir1,"bostaurus_brahma_bionano_NCBI_full_No_Ns.rls")

bostaurus_brahma_bionano_NCBI_full_No_Ns <- read_tsv(path1,col_names = FALSE)
names(bostaurus_brahma_bionano_NCBI_full_No_Ns) <- c("scaffold","nothing","gap","length","perc_gap")

bostaurus_brahma_bionano_NCBI_full_No_Ns <- bostaurus_brahma_bionano_NCBI_full_No_Ns %>%
  dplyr::select(scaffold,gap,length) %>% arrange(desc(length))

sum(as.numeric(bostaurus_brahma_bionano_NCBI_full_No_Ns$length)) - 
  sum(as.numeric(bostaurus_brahma_bionano_NCBI_full_No_Ns$gap))

N50(bostaurus_brahma_bionano_NCBI_full_No_Ns)
nrow(bostaurus_brahma_bionano_NCBI_full_No_Ns)
sum(as.numeric(bostaurus_brahma_bionano_NCBI_full_No_Ns$length))

# reading bostaurus_brahma_bionano_NCBI_full_corrected_gapfill_arrow_fil_withM_chronly_No_Ns.rls
path1 <- paste0(dir1,"bostaurus_brahma_bionano_NCBI_full_corrected_gapfill_arrow_fil_withM_chronly_No_Ns.rls")

bostaurus_brahma_bionano_NCBI_full_corrected_gapfill_arrow_fil_withM_chronly_No_Ns <- read_tsv(path1,col_names = FALSE)
names(bostaurus_brahma_bionano_NCBI_full_corrected_gapfill_arrow_fil_withM_chronly_No_Ns) <- c("scaffold","nothing","gap","length","perc_gap")

bostaurus_brahma_bionano_NCBI_full_corrected_gapfill_arrow_fil_withM_chronly_No_Ns <- bostaurus_brahma_bionano_NCBI_full_corrected_gapfill_arrow_fil_withM_chronly_No_Ns %>%
  dplyr::select(scaffold,gap,length) %>% arrange(desc(length))

sum(as.numeric(bostaurus_brahma_bionano_NCBI_full_corrected_gapfill_arrow_fil_withM_chronly_No_Ns$length)) - 
  sum(as.numeric(bostaurus_brahma_bionano_NCBI_full_corrected_gapfill_arrow_fil_withM_chronly_No_Ns$gap))

N50(bostaurus_brahma_bionano_NCBI_full_corrected_gapfill_arrow_fil_withM_chronly_No_Ns)
nrow(bostaurus_brahma_bionano_NCBI_full_corrected_gapfill_arrow_fil_withM_chronly_No_Ns)
sum(as.numeric(bostaurus_brahma_bionano_NCBI_full_corrected_gapfill_arrow_fil_withM_chronly_No_Ns$length))

