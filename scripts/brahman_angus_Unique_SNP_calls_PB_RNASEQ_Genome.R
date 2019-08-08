#------------------------------------------------------
# Program name: brahman_angus_Unique_SNP_calls_PB_RNASEQ_Genome.R
# Objective: replot liz unique SNP calls that she did in 
#          Excel
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(readr)
library(xlsx)
library(dplyr)
library(ggplot2)

dir1 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/isoseq/IsoPhase/April2019_IsoPhase_with_Lloyd/"
path1 <- paste0(dir1,"Unique_SNP_calls_PB_RNASEQ_Genome.xlsx")

Unique_SNP_calls <- read.xlsx2(path1,1,stringsAsFactors = FALSE)
Unique_SNP_calls$count <- as.numeric(Unique_SNP_calls$count)
Unique_SNP_calls$totalByType <- as.numeric(Unique_SNP_calls$totalByType)

Unique_SNP_calls <- Unique_SNP_calls %>% mutate(percentage = round(count/totalByType*100))

tiff(filename = "FigFinal_Unique_SNP_calls_PB_RNASEQ_Genome.tiff",width = 600, height = 600)
g <- ggplot(Unique_SNP_calls, aes(SNP, percentage,fill=type)) 
g <- g + geom_bar(position="dodge",stat="identity") +labs(x="SNP",y="Percentage of unique SNP calls")
g <- g + theme_bw()
g
dev.off()
