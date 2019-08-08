#------------------------------------------------------
# Program name: brahman_angus_bionano_corrected_cut.R
# Objective: read in before and after bionano cut
#           correction and make change to scaffold order
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(readr)
library(dplyr)

##### Dam #####

dir1 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/conflict_resolution/"

#uncorrected
path1 <- paste0(dir1,"Brahma_bionano_cut_uncorrected.txt")

Brahma_bionano_cut_uncorrected <- read_tsv(path1,col_names = TRUE)

#corrected
path2 <- paste0(dir1,"Brahma_bionano_cut_corrected.txt")

Brahma_bionano_cut_corrected <- read_tsv(path2,col_names = TRUE)

#scaffold order book pre making agp
dir2 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/"
path3 <- paste0(dir2,"scaffold_order_bionano_brahma.txt")

scaffold_order_bionano_brahma <- read_tsv(path3,col_names = FALSE)
colnames(scaffold_order_bionano_brahma) <- c("chr","canu","start","end","orientation")

#loop thro bionano cut contig and modify scaffold order book
#1. changed cut and corrected while in the scaffold order bionano book
for (i in 1:nrow(Brahma_bionano_cut_uncorrected)){
  if (Brahma_bionano_cut_uncorrected$change[i] == "yes"){
    for (j in 1:nrow(scaffold_order_bionano_brahma)){
      if (scaffold_order_bionano_brahma$canu[j] == Brahma_bionano_cut_uncorrected$canu[i] && 
          scaffold_order_bionano_brahma$start[j] == Brahma_bionano_cut_uncorrected$start[i] && 
          scaffold_order_bionano_brahma$end[j] == Brahma_bionano_cut_uncorrected$end[i]){
        
        #perform the correction here
        scaffold_order_bionano_brahma$start[j] <- Brahma_bionano_cut_corrected$start[i]
        scaffold_order_bionano_brahma$end[j] <- Brahma_bionano_cut_corrected$end[i]
        #print(j)
      } 
    }
  }
}

#write out tab delimited cut corrected file
write_tsv(scaffold_order_bionano_brahma,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_bionano_brahma_corrected.txt")

#2. changed cut with only partial cut in the scaffold order bionano book
# for (i in 1:nrow(Brahma_bionano_cut_corrected)){
#   for (j in 1:nrow(scaffold_order_bionano_brahma)){
#     Brahma_bionano_cut_corrected$canu[i] == scaffold_order_bionano_brahma$canu[j]
#   }
# }
#above is unfinished, decided to do in Excel as only 11 in this category

#3. cut but totally not in scaffold order bionano book
#again not that many, marked as "blanks" in cut_check.xlsx

##### Sire #####
dir1 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/conflict_resolution/"

#uncorrected
path1 <- paste0(dir1,"Angus_bionano_cut_uncorrected.txt")

Angus_bionano_cut_uncorrected <- read_tsv(path1,col_names = TRUE)

#corrected
path2 <- paste0(dir1,"Angus_bionano_cut_corrected.txt")

Angus_bionano_cut_corrected <- read_tsv(path2,col_names = TRUE)

#scaffold order book pre making agp
dir2 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/"
path3 <- paste0(dir2,"scaffold_order_bionano_angus.txt")

scaffold_order_bionano_angus <- read_tsv(path3,col_names = FALSE)
colnames(scaffold_order_bionano_angus) <- c("chr","canu","start","end","orientation")

#loop thro bionano cut contig and modify scaffold order book
#1. changed cut and corrected while in the scaffold order bionano book
for (i in 1:nrow(Angus_bionano_cut_uncorrected)){
  if (Angus_bionano_cut_uncorrected$change[i] == "yes"){
    for (j in 1:nrow(scaffold_order_bionano_angus)){
      if (scaffold_order_bionano_angus$canu[j] == Angus_bionano_cut_uncorrected$canu[i] && 
          scaffold_order_bionano_angus$start[j] == Angus_bionano_cut_uncorrected$start[i] && 
          scaffold_order_bionano_angus$end[j] == Angus_bionano_cut_uncorrected$end[i]){
        
        #perform the correction here
        scaffold_order_bionano_angus$start[j] <- Angus_bionano_cut_corrected$start[i]
        scaffold_order_bionano_angus$end[j] <- Angus_bionano_cut_corrected$end[i]
        #print(j)
      } 
    }
  }
}

#write out tab delimited cut corrected file
write_tsv(scaffold_order_bionano_angus,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_bionano_angus_corrected.txt")
