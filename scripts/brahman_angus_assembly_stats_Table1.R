#------------------------------------------------------
# Program name: brahman_angus_assembly_stats_Table1.R
# Objective: All calculations relevant to assembly stats
#           for Brahman and Angus paper
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(readr)
library(dplyr)

#calculate N50 function #revised with arranging the df
N50 <- function(assembly_df){
  assembly_df <- assembly_df %>% arrange(desc(length))
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

#angus
# path to all gapLength results
dir1 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/fasta_file_order_checkBases_N/"

# reading bostaurus_angus_No_Ns.rls
path1 <- paste0(dir1,"bostaurus_angus_No_Ns.rls")

df <- read_tsv(path1,col_names = FALSE)
names(df) <- c("scaffold","nothing","gap","length","perc_gap")

df <- df %>% dplyr::select(scaffold,gap,length) %>% arrange(desc(length))

#ungapped length
sum(as.numeric(df$length)) - sum(as.numeric(df$gap))
#2573806573

N50(df)
#29439230

nrow(df)
#1747

#no of gap
#0

# reading f1_sire_salsa_No_Ns.rls
path1 <- paste0(dir1,"f1_sire_salsa_No_Ns.rls")

df <- read_tsv(path1,col_names = FALSE)
names(df) <- c("scaffold","nothing","gap","length","perc_gap")

df <- df %>% dplyr::select(scaffold,gap,length) %>% arrange(desc(length))

#ungapped length
sum(as.numeric(df$length)) - sum(as.numeric(df$gap))
#2573806573

N50(df)
#104553360

nrow(df)
#1515

#no of gap
#wc -l $PWD/f1_sire_salsa.coor
#235 /fast/users/a1223107/gap_genome_analysis/species/BrahmanAngus/f1_sire_salsa.coor

# reading bostaurus_angus_bionano_NCBI_full_No_Ns.rls
path1 <- paste0(dir1,"bostaurus_angus_bionano_NCBI_full_No_Ns.rls")

df <- read_tsv(path1,col_names = FALSE)
names(df) <- c("scaffold","nothing","gap","length","perc_gap")

df <- df %>% dplyr::select(scaffold,gap,length) %>% arrange(desc(length))

#ungapped length
sum(as.numeric(df$length)) - sum(as.numeric(df$gap))
#2573806573

N50(df)
#35241510

nrow(df)
#1595

#no of gap
#wc -l $PWD/bostaurus_angus_bionano_NCBI_full.coor 
#181 /fast/users/a1223107/gap_genome_analysis/species/BrahmanAngus/bostaurus_angus_bionano_NCBI_full.coor

# reading bostaurus_angus_bionano_NCBI_full_corrected_gapfill_arrow_fil_No_Ns.rls
path1 <- paste0(dir1,"bostaurus_angus_bionano_NCBI_full_corrected_gapfill_arrow_fil_No_Ns.rls")

df <- read_tsv(path1,col_names = FALSE)
names(df) <- c("scaffold","nothing","gap","length","perc_gap")

df <- df %>% dplyr::select(scaffold,gap,length) %>% arrange(desc(length))

#ungapped length
sum(as.numeric(df$length)) - sum(as.numeric(df$gap))
#2574284011

N50(df)
#102763607

nrow(df)
#1435

#no of gap
#wc -l $PWD/bostaurus_angus_bionano_NCBI_full_corrected_gapfill_arrow_fil.coor 
#277 /fast/users/a1223107/gap_genome_analysis/species/BrahmanAngus/bostaurus_angus_bionano_NCBI_full_corrected_gapfill_arrow_fil.coor

#brahman
# reading bostaurus_brahma_No_Ns.rls
path1 <- paste0(dir1,"bostaurus_brahma_No_Ns.rls")

df <- read_tsv(path1,col_names = FALSE)
names(df) <- c("scaffold","nothing","gap","length","perc_gap")

df <- df %>% dplyr::select(scaffold,gap,length) %>% arrange(desc(length))

#ungapped length
sum(as.numeric(df$length)) - sum(as.numeric(df$gap))
#2678769087

N50(df)
#23445275

nrow(df)
#1585

#no of gap
#0

# reading f1_dam_salsa_No_Ns.rls
path1 <- paste0(dir1,"f1_dam_salsa_No_Ns.rls")

df <- read_tsv(path1,col_names = FALSE)
names(df) <- c("scaffold","nothing","gap","length","perc_gap")

df <- df %>% dplyr::select(scaffold,gap,length) %>% arrange(desc(length))

#ungapped length
sum(as.numeric(df$length)) - sum(as.numeric(df$gap))
#2678769087

N50(df)
#72609562

nrow(df)
#1370

#no of gap
#wc -l $PWD/f1_dam_salsa.coor
#216 /fast/users/a1223107/gap_genome_analysis/species/BrahmanAngus/f1_dam_salsa.coor

# reading bostaurus_brahma_bionano_NCBI_full_No_Ns.rls
path1 <- paste0(dir1,"bostaurus_brahma_bionano_NCBI_full_No_Ns.rls")

df <- read_tsv(path1,col_names = FALSE)
names(df) <- c("scaffold","nothing","gap","length","perc_gap")

df <- df %>% dplyr::select(scaffold,gap,length) %>% arrange(desc(length))

#ungapped length
sum(as.numeric(df$length)) - sum(as.numeric(df$gap))
#2678769087

N50(df)
#31743416

nrow(df)
#1353

#no of gap
#wc -l $PWD/bostaurus_brahma_bionano_NCBI_full.coor
#268 /fast/users/a1223107/gap_genome_analysis/species/BrahmanAngus/bostaurus_brahma_bionano_NCBI_full.coor

# reading bostaurus_brahma_bionano_NCBI_full_corrected_gapfill_arrow_fil_withM_No_Ns.rls
path1 <- paste0(dir1,"bostaurus_brahma_bionano_NCBI_full_corrected_gapfill_arrow_fil_withM_No_Ns.rls")

df <- read_tsv(path1,col_names = FALSE)
names(df) <- c("scaffold","nothing","gap","length","perc_gap")

df <- df %>% dplyr::select(scaffold,gap,length) %>% arrange(desc(length))

#ungapped length
sum(as.numeric(df$length)) - sum(as.numeric(df$gap))
#2679332898

N50(df)
#104466507

nrow(df)
#1251

#no of gap
#wc -l $PWD/bostaurus_brahma_bionano_NCBI_full_corrected_gapfill_arrow_fil_withM.coor
#302 /fast/users/a1223107/gap_genome_analysis/species/BrahmanAngus/bostaurus_brahma_bionano_NCBI_full_corrected_gapfill_arrow_fil_withM.coor


