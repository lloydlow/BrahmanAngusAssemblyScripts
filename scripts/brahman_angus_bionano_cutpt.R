#------------------------------------------------------
# Program name: brahman_angus_bionano_cutpt.R
# Objective: make cutpt file to be looped over to get
#           cov
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(readr)
library(dplyr)

#Dam
#read in cutpt_brahma, testing with cut point not corrected
cutpt_brahma <- read_tsv("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/conflict_resolution/cutpt_brahma.txt",
                         col_names = TRUE)
names(cutpt_brahma) <- c("canu","startseq","endseq")

cutpt_brahma <- cutpt_brahma %>% mutate(start = endseq - 1e5, end = endseq + 1e5)

selc <- c()
for (i in 1:nrow(cutpt_brahma)){
  if (nrow(cutpt_brahma) != i && cutpt_brahma$canu[i] == cutpt_brahma$canu[i+1]){
    holder <- TRUE
    selc <- c(selc,holder)
  } else {
    holder <- FALSE
    selc <- c(selc,holder)
  }
}

cutpt_brahma <- cutpt_brahma[selc,]

for (i in 1:nrow(cutpt_brahma)){
  if (cutpt_brahma$start[i] < 0){
    cutpt_brahma$start[i] <- cutpt_brahma$startseq[i]
  }
}

write_tsv(cutpt_brahma,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/conflict_resolution/cutpt_for_looping/cutpt_brahma.tsv",
          col_names = FALSE)

#Sire
#read in cutpt_angus, testing with cut point not corrected
cutpt_angus <- read_tsv("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/conflict_resolution/cutpt_angus.txt",
                        col_names = TRUE)
names(cutpt_angus) <- c("canu","startseq","endseq")

cutpt_angus <- cutpt_angus %>% mutate(start = endseq - 1e5, end = endseq + 1e5)

selc <- c()
for (i in 1:nrow(cutpt_angus)){
  if (nrow(cutpt_angus) != i && cutpt_angus$canu[i] == cutpt_angus$canu[i+1]){
    holder <- TRUE
    selc <- c(selc,holder)
  } else {
    holder <- FALSE
    selc <- c(selc,holder)
  }
}

cutpt_angus <- cutpt_angus[selc,]

for (i in 1:nrow(cutpt_angus)){
  if (cutpt_angus$start[i] < 0){
    cutpt_angus$start[i] <- cutpt_angus$startseq[i]
  }
}

write_tsv(cutpt_angus,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/conflict_resolution/cutpt_for_looping/cutpt_angus.tsv",
          col_names = FALSE)
