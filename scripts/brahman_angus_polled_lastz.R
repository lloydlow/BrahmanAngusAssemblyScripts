#------------------------------------------------------
# Program name: brahman_angus_polled_lastz.R
# Objective: format=rdotplot option can be used to
#           show dotplot 
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(readr)
library(dplyr)

dir1 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/genomeAlignmentToCattle/lastzCattle/"

#read each line of the lastz rdotplot output file
path1 <- paste0(dir1,"Angus_vs_Brahman_polled_status1")
#files <- readLines(path1)

dotplot <- read_tsv(path1,comment = "NA",skip=1)
names(dotplot) <- c("brahman_chr1","angus_chr1")

# plot(dotplot$brahman_chr1,dotplot$angus_chr1,xlab = "brahman_chr1",ylab="angus_chr1")

brahman_1Mb_1 <- dotplot %>% filter(brahman_chr1 <= 1e6)

png(filename = "dotplot.png")
plot(dotplot$brahman_chr1,dotplot$brahman_chr1,xlab = "brahman_chr1",ylab="angus_chr1")
dev.off()
