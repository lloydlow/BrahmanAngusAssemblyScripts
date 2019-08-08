#------------------------------------------------------
# Program name: brahman_angus_bionano_align_orderSNP.R
# Objective: check tab file from derek alignAndOrderSnpProbes.pl
#           using the bionano scaffolds as input
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(readr)
library(dplyr)

##### Dam #####
# reading HD probes from bionano_NCBI_fasta, not scaffolded fasta is not here
dir1 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/alignAndOrderSnp_result/"
path1 <- paste0(dir1,"bostaurus_brahma_bionano_NCBI_full.HDProbes.tab")

f1_dam_bionano_HDProbes_tab <- read_tsv(path1,col_names = FALSE)
names(f1_dam_bionano_HDProbes_tab) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

#subs the _arrow to |arrow #only for bionano stuff
f1_dam_bionano_HDProbes_tab$scaffold <- gsub("_arrow","|arrow",f1_dam_bionano_HDProbes_tab$scaffold)

#check each chr is made up of how many scaffolds
f1_dam_bionano_HDProbes_tab_by_chrScaff <- f1_dam_bionano_HDProbes_tab %>% group_by(chromosome,scaffold) %>% 
  summarise(n=n()) %>% arrange(chromosome,desc(n))

#scaffold with the most probe hits
f1_dam_bionano_HDProbes_tab_by_chrScaff_mostprobes <- f1_dam_bionano_HDProbes_tab_by_chrScaff %>% 
  filter(scaffold != "*") %>% group_by(scaffold) %>%
  slice(which.max(n)) %>% arrange(desc(n))

# reading rc probes
dir2 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/alignAndOrderSnp_result/"
path2 <- paste0(dir1,"bostaurus_brahma_bionano_NCBI_full.rcmap.tab")

f1_dam_bionano_fasta_rcmap_tab <- read_tsv(path2,col_names = FALSE)
names(f1_dam_bionano_fasta_rcmap_tab) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

#subs the _arrow to |arrow #only for bionano stuff
f1_dam_bionano_fasta_rcmap_tab$scaffold <- gsub("_arrow","|arrow",f1_dam_bionano_fasta_rcmap_tab$scaffold)

#check each chr is made up of how many scaffolds
f1_dam_bionano_fasta_rcmap_tab_by_chrScaff <- f1_dam_bionano_fasta_rcmap_tab %>% group_by(chromosome,scaffold) %>% 
  summarise(n=n()) %>% arrange(chromosome,desc(n))

#scaffold with the most probe hits
f1_dam_bionano_fasta_rcmap_tab_by_chrScaff_mostprobes <- f1_dam_bionano_fasta_rcmap_tab_by_chrScaff %>% 
  filter(scaffold != "*") %>% group_by(scaffold) %>%
  slice(which.max(n)) %>% arrange(desc(n))

save(f1_dam_bionano_HDProbes_tab_by_chrScaff_mostprobes,f1_dam_bionano_fasta_rcmap_tab_by_chrScaff_mostprobes,
     file = "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/order_bionano_scaffolds_to_chr_level/f1_dam_bionano_fasta_rcmap_tab_by_chrScaff_mostprobes.RData")

save(f1_dam_bionano_HDProbes_tab,f1_dam_bionano_fasta_rcmap_tab,
     file="/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/order_bionano_scaffolds_to_chr_level/f1_dam_bionano_HD_rc_Probes_tab.RData")

##### Sire #####
# reading HD probes from bionano_NCBI_fasta, not scaffolded fasta is not here
dir1 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/alignAndOrderSnp_result/"
path1 <- paste0(dir1,"bostaurus_angus_bionano_NCBI_full.HDProbes.tab")

f1_sire_bionano_HDProbes_tab <- read_tsv(path1,col_names = FALSE)
names(f1_sire_bionano_HDProbes_tab) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

#subs the _arrow to |arrow #only for bionano stuff
f1_sire_bionano_HDProbes_tab$scaffold <- gsub("_arrow","|arrow",f1_sire_bionano_HDProbes_tab$scaffold)

#check each chr is made up of how many scaffolds
f1_sire_bionano_HDProbes_tab_by_chrScaff <- f1_sire_bionano_HDProbes_tab %>% group_by(chromosome,scaffold) %>% 
  summarise(n=n()) %>% arrange(chromosome,desc(n))

#scaffold with the most probe hits
f1_sire_bionano_HDProbes_tab_by_chrScaff_mostprobes <- f1_sire_bionano_HDProbes_tab_by_chrScaff %>% 
  filter(scaffold != "*") %>% group_by(scaffold) %>%
  slice(which.max(n)) %>% arrange(desc(n))

# reading rc probes
dir2 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/alignAndOrderSnp_result/"
path2 <- paste0(dir1,"bostaurus_angus_bionano_NCBI_full.rcmap.tab")

f1_sire_bionano_fasta_rcmap_tab <- read_tsv(path2,col_names = FALSE)
names(f1_sire_bionano_fasta_rcmap_tab) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

#subs the _arrow to |arrow #only for bionano stuff
f1_sire_bionano_fasta_rcmap_tab$scaffold <- gsub("_arrow","|arrow",f1_sire_bionano_fasta_rcmap_tab$scaffold)

#check each chr is made up of how many scaffolds
f1_sire_bionano_fasta_rcmap_tab_by_chrScaff <- f1_sire_bionano_fasta_rcmap_tab %>% group_by(chromosome,scaffold) %>% 
  summarise(n=n()) %>% arrange(chromosome,desc(n))

#scaffold with the most probe hits
f1_sire_bionano_fasta_rcmap_tab_by_chrScaff_mostprobes <- f1_sire_bionano_fasta_rcmap_tab_by_chrScaff %>% 
  filter(scaffold != "*") %>% group_by(scaffold) %>%
  slice(which.max(n)) %>% arrange(desc(n))

save(f1_sire_bionano_HDProbes_tab_by_chrScaff_mostprobes,f1_sire_bionano_fasta_rcmap_tab_by_chrScaff_mostprobes,
     file = "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/order_bionano_scaffolds_to_chr_level/f1_sire_bionano_fasta_rcmap_tab_by_chrScaff_mostprobes.RData")

save(f1_sire_bionano_HDProbes_tab,f1_sire_bionano_fasta_rcmap_tab,
     file="/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/order_bionano_scaffolds_to_chr_level/f1_sire_bionano_HD_rc_Probes_tab.RData")
