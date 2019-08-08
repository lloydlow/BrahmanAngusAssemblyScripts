# brahman_angus_scaffold_gaplength_20180315.R

library(readr)
library(dplyr)

# reading *No_Ns.rls file to ensure sequence manipulation gives 
# the same number of bases

#1. bostaurus_angus_No_Ns.rls
dir1 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/brahman_Angus_scaffolding/"
path1 <- paste0(dir1,"bostaurus_angus_No_Ns.rls")

bostaurus_angus_No_Ns <- read_tsv(path1,col_names = FALSE)
names(bostaurus_angus_No_Ns) <- c("scaffold","nothing","gap","length","perc_gap")

bostaurus_angus_No_Ns <- bostaurus_angus_No_Ns %>%
  select(scaffold,gap,length)

sum(as.numeric(bostaurus_angus_No_Ns$length)) - 
  sum(as.numeric(bostaurus_angus_No_Ns$gap))
nrow(bostaurus_angus_No_Ns)

#2. bostaurus_brahma_No_Ns.rls
dir2 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/brahman_Angus_scaffolding/"
path2 <- paste0(dir2,"bostaurus_brahma_No_Ns.rls")

bostaurus_brahma_No_Ns <- read_tsv(path2,col_names = FALSE)
names(bostaurus_brahma_No_Ns) <- c("scaffold","nothing","gap","length","perc_gap")

bostaurus_brahma_No_Ns <- bostaurus_brahma_No_Ns %>%
  select(scaffold,gap,length)

sum(as.numeric(bostaurus_brahma_No_Ns$length)) - 
  sum(as.numeric(bostaurus_brahma_No_Ns$gap))
nrow(bostaurus_brahma_No_Ns)

#3. f1_dam_3ddna_No_Ns.rls
dir3 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/brahman_Angus_scaffolding/"
path3 <- paste0(dir3,"f1_dam_3ddna_No_Ns.rls")

f1_dam_3ddna_No_Ns <- read_tsv(path3,col_names = FALSE)
names(f1_dam_3ddna_No_Ns) <- c("scaffold","nothing","gap","length","perc_gap")

f1_dam_3ddna_No_Ns <- f1_dam_3ddna_No_Ns %>%
  select(scaffold,gap,length)

sum(as.numeric(f1_dam_3ddna_No_Ns$length)) - 
  sum(as.numeric(f1_dam_3ddna_No_Ns$gap))
nrow(f1_dam_3ddna_No_Ns)

#4. f1_dam_3ddna_v2_No_Ns.rls
dir4 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/brahman_Angus_scaffolding/"
path4 <- paste0(dir4,"f1_dam_3ddna_v2_No_Ns.rls")

f1_dam_3ddna_v2_No_Ns <- read_tsv(path4,col_names = FALSE)
names(f1_dam_3ddna_v2_No_Ns) <- c("scaffold","nothing","gap","length","perc_gap")

f1_dam_3ddna_v2_No_Ns <- f1_dam_3ddna_v2_No_Ns %>%
  select(scaffold,gap,length)

sum(as.numeric(f1_dam_3ddna_v2_No_Ns$length)) - 
  sum(as.numeric(f1_dam_3ddna_v2_No_Ns$gap))
nrow(f1_dam_3ddna_v2_No_Ns)

#5. f1_dam_phase_No_Ns.rls
dir5 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/brahman_Angus_scaffolding/"
path5 <- paste0(dir5,"f1_dam_phase_No_Ns.rls")

f1_dam_phase_No_Ns <- read_tsv(path5,col_names = FALSE)
names(f1_dam_phase_No_Ns) <- c("scaffold","nothing","gap","length","perc_gap")

f1_dam_phase_No_Ns <- f1_dam_phase_No_Ns %>%
  select(scaffold,gap,length)

sum(as.numeric(f1_dam_phase_No_Ns$length)) - 
  sum(as.numeric(f1_dam_phase_No_Ns$gap))
nrow(f1_dam_phase_No_Ns)

#6. f1_dam_salsa_No_Ns.rls
dir6 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/brahman_Angus_scaffolding/"
path6 <- paste0(dir6,"f1_dam_salsa_No_Ns.rls")

f1_dam_salsa_No_Ns <- read_tsv(path6,col_names = FALSE)
names(f1_dam_salsa_No_Ns) <- c("scaffold","nothing","gap","length","perc_gap")

f1_dam_salsa_No_Ns <- f1_dam_salsa_No_Ns %>%
  select(scaffold,gap,length)

sum(as.numeric(f1_dam_salsa_No_Ns$length)) - 
  sum(as.numeric(f1_dam_salsa_No_Ns$gap))
nrow(f1_dam_salsa_No_Ns)

#7. f1_sire_3ddna_No_Ns.rls
dir7 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/brahman_Angus_scaffolding/"
path7 <- paste0(dir7,"f1_sire_3ddna_No_Ns.rls")

f1_sire_3ddna_No_Ns <- read_tsv(path7,col_names = FALSE)
names(f1_sire_3ddna_No_Ns) <- c("scaffold","nothing","gap","length","perc_gap")

f1_sire_3ddna_No_Ns <- f1_sire_3ddna_No_Ns %>%
  select(scaffold,gap,length)

sum(as.numeric(f1_sire_3ddna_No_Ns$length)) - 
  sum(as.numeric(f1_sire_3ddna_No_Ns$gap))
nrow(f1_sire_3ddna_No_Ns)

#8. f1_sire_3ddna_v2_No_Ns.rls
dir8 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/brahman_Angus_scaffolding/"
path8 <- paste0(dir8,"f1_sire_3ddna_v2_No_Ns.rls")

f1_sire_3ddna_v2_No_Ns <- read_tsv(path8,col_names = FALSE)
names(f1_sire_3ddna_v2_No_Ns) <- c("scaffold","nothing","gap","length","perc_gap")

f1_sire_3ddna_v2_No_Ns <- f1_sire_3ddna_v2_No_Ns %>%
  select(scaffold,gap,length)

sum(as.numeric(f1_sire_3ddna_v2_No_Ns$length)) - 
  sum(as.numeric(f1_sire_3ddna_v2_No_Ns$gap))
nrow(f1_sire_3ddna_v2_No_Ns)

#9. f1_sire_phase_No_Ns.rls
dir9 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/brahman_Angus_scaffolding/"
path9 <- paste0(dir9,"f1_sire_phase_No_Ns.rls")

f1_sire_phase_No_Ns <- read_tsv(path9,col_names = FALSE)
names(f1_sire_phase_No_Ns) <- c("scaffold","nothing","gap","length","perc_gap")

f1_sire_phase_No_Ns <- f1_sire_phase_No_Ns %>%
  select(scaffold,gap,length)

sum(as.numeric(f1_sire_phase_No_Ns$length)) - 
  sum(as.numeric(f1_sire_phase_No_Ns$gap))
nrow(f1_sire_phase_No_Ns)

#10. f1_sire_salsa_No_Ns.rls
dir10 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/brahman_Angus_scaffolding/"
path10 <- paste0(dir10,"f1_sire_salsa_No_Ns.rls")

f1_sire_salsa_No_Ns <- read_tsv(path10,col_names = FALSE)
names(f1_sire_salsa_No_Ns) <- c("scaffold","nothing","gap","length","perc_gap")

f1_sire_salsa_No_Ns <- f1_sire_salsa_No_Ns %>%
  select(scaffold,gap,length)

sum(as.numeric(f1_sire_salsa_No_Ns$length)) - 
  sum(as.numeric(f1_sire_salsa_No_Ns$gap))
nrow(f1_sire_salsa_No_Ns)

#11. dam_best_scaffold_reference_No_Ns.rls
dir11 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/brahman_Angus_scaffolding/"
path11 <- paste0(dir11,"dam_best_scaffold_reference_No_Ns.rls")

dam_best_scaffold_reference_No_Ns <- read_tsv(path11,col_names = FALSE)
names(dam_best_scaffold_reference_No_Ns) <- c("scaffold","nothing","gap","length","perc_gap")

dam_best_scaffold_reference_No_Ns <- dam_best_scaffold_reference_No_Ns %>%
  select(scaffold,gap,length)

sum(as.numeric(dam_best_scaffold_reference_No_Ns$length)) - 
  sum(as.numeric(dam_best_scaffold_reference_No_Ns$gap))
nrow(dam_best_scaffold_reference_No_Ns)

#12. sire_best_scaffold_reference_No_Ns.rls
dir12 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/brahman_Angus_scaffolding/"
path12 <- paste0(dir12,"sire_best_scaffold_reference_No_Ns.rls")

sire_best_scaffold_reference_No_Ns <- read_tsv(path12,col_names = FALSE)
names(sire_best_scaffold_reference_No_Ns) <- c("scaffold","nothing","gap","length","perc_gap")

sire_best_scaffold_reference_No_Ns <- sire_best_scaffold_reference_No_Ns %>%
  select(scaffold,gap,length)

sum(as.numeric(sire_best_scaffold_reference_No_Ns$length)) - 
  sum(as.numeric(sire_best_scaffold_reference_No_Ns$gap))
nrow(sire_best_scaffold_reference_No_Ns)
