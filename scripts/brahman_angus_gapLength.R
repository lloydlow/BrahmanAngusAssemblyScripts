#------------------------------------------------------
# Program name: brahman_angus_gapLength.R
# Objective: check assembly gap and length of the
#           various versions
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(readr)
library(dplyr)

# reading dam_scaffold_salsa_only.ref_No_Ns.rls
dir2 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/fasta_file_order_checkBases_N/"
path1 <- paste0(dir2,"dam_scaffold_salsa_only.ref_No_Ns.rls")

dam_scaffold_salsa_only.ref_No_Ns <- read_tsv(path1,col_names = FALSE)
names(dam_scaffold_salsa_only.ref_No_Ns) <- c("scaffold","nothing","gap","length","perc_gap")

dam_scaffold_salsa_only.ref_No_Ns <- dam_scaffold_salsa_only.ref_No_Ns %>%
  select(scaffold,gap,length)

sum(as.numeric(dam_scaffold_salsa_only.ref_No_Ns$length)) - 
  sum(as.numeric(dam_scaffold_salsa_only.ref_No_Ns$gap))

# reading sire_scaffold_salsa_only.complete_No_Ns.rls
dir2 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/fasta_file_order_checkBases_N/"
path2 <- paste0(dir2,"sire_scaffold_salsa_only.complete_No_Ns.rls")

sire_scaffold_salsa_only.complete_No_Ns <- read_tsv(path2,col_names = FALSE)
names(sire_scaffold_salsa_only.complete_No_Ns) <- c("scaffold","nothing","gap","length","perc_gap")

sire_scaffold_salsa_only.complete_No_Ns <- sire_scaffold_salsa_only.complete_No_Ns %>%
  select(scaffold,gap,length)

sum(as.numeric(sire_scaffold_salsa_only.complete_No_Ns$length)) - 
  sum(as.numeric(sire_scaffold_salsa_only.complete_No_Ns$gap))

# reading angus_3ddna2_v1.complete_No_Ns.rls
dir2 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/fasta_file_order_checkBases_N/"
path3 <- paste0(dir2,"angus_3ddna2_v1.complete_No_Ns.rls")

angus_3ddna2_v1.complete_No_Ns <- read_tsv(path3,col_names = FALSE)
names(angus_3ddna2_v1.complete_No_Ns) <- c("scaffold","nothing","gap","length","perc_gap")

angus_3ddna2_v1.complete_No_Ns <- angus_3ddna2_v1.complete_No_Ns %>%
  select(scaffold,gap,length)

sum(as.numeric(angus_3ddna2_v1.complete_No_Ns$length)) - 
  sum(as.numeric(angus_3ddna2_v1.complete_No_Ns$gap))

# reading dam_assembly.cleaned_No_Ns.rls
dir2 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/fasta_file_order_checkBases_N/"
path4 <- paste0(dir2,"dam_assembly.cleaned_No_Ns.rls")

dam_assembly.cleaned_No_Ns <- read_tsv(path4,col_names = FALSE)
names(dam_assembly.cleaned_No_Ns) <- c("scaffold","nothing","gap","length","perc_gap")

dam_assembly.cleaned_No_Ns <- dam_assembly.cleaned_No_Ns %>%
  select(scaffold,gap,length)

sum(as.numeric(dam_assembly.cleaned_No_Ns$length)) - 
  sum(as.numeric(dam_assembly.cleaned_No_Ns$gap))

# reading sire_assembly.cleaned_No_Ns.rls
dir2 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/fasta_file_order_checkBases_N/"
path5 <- paste0(dir2,"sire_assembly.cleaned_No_Ns.rls")

sire_assembly.cleaned_No_Ns <- read_tsv(path5,col_names = FALSE)
names(sire_assembly.cleaned_No_Ns) <- c("scaffold","nothing","gap","length","perc_gap")

sire_assembly.cleaned_No_Ns <- sire_assembly.cleaned_No_Ns %>%
  select(scaffold,gap,length)

sum(as.numeric(sire_assembly.cleaned_No_Ns$length)) - 
  sum(as.numeric(sire_assembly.cleaned_No_Ns$gap))


# reading dam_scaffold_salsa_only.complete_No_Ns.rls
dir2 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/fasta_file_order_checkBases_N/"
path6 <- paste0(dir2,"dam_scaffold_salsa_only.complete_No_Ns.rls")

dam_scaffold_salsa_only.complete_No_Ns <- read_tsv(path6,col_names = FALSE)
names(dam_scaffold_salsa_only.complete_No_Ns) <- c("scaffold","nothing","gap","length","perc_gap")

dam_scaffold_salsa_only.complete_No_Ns <- dam_scaffold_salsa_only.complete_No_Ns %>%
  select(scaffold,gap,length)

sum(as.numeric(dam_scaffold_salsa_only.complete_No_Ns$length)) - 
  sum(as.numeric(dam_scaffold_salsa_only.complete_No_Ns$gap))

# f1_dam_salsa_No_Ns.rls
dir2 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/fasta_file_order_checkBases_N/"
path7 <- paste0(dir2,"f1_dam_salsa_No_Ns.rls")

f1_dam_salsa_No_Ns <- read_tsv(path7,col_names = FALSE)
names(f1_dam_salsa_No_Ns) <- c("scaffold","nothing","gap","length","perc_gap")

f1_dam_salsa_No_Ns <- f1_dam_salsa_No_Ns %>%
  select(scaffold,gap,length)

sum(as.numeric(f1_dam_salsa_No_Ns$length)) - 
  sum(as.numeric(f1_dam_salsa_No_Ns$gap))
