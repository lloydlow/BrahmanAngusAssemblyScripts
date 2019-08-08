#------------------------------------------------------
# Program name: brahman_angus_polishing_2rounds.R
# Objective: script to check what has actually changed
#           in each round of polishing
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

#gff3 format
#c("seqid","source","type","start","end","score","strand","phase","attributes")

#angus
#changed Y to 30
# path to all arrow polish results
dir1 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/arrow_polish/"

# reading bostaurus_angus_arrow1_clean.gff
path1 <- paste0(dir1,"bostaurus_angus_arrow1_clean.gff")

bostaurus_angus_No_Ns <- read_tsv(path1,col_names = FALSE)
names(bostaurus_angus_No_Ns) <- c("seqid","source","type","start","end","score","strand","phase","attributes")



