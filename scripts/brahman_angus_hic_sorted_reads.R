#------------------------------------------------------
# Program name: brahman_angus_hic_sorted_reads.R
# Objective: investigate the kmer used to separate
#           haplotypes
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(readr)
library(dplyr)
#library(data.table)

angus_kmer <- read_delim("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/hic_stuff/trio_hic/bostaurus_angus.kmers", 
                         " ",col_names = FALSE)

colnames(angus_kmer) <- c("kmer","freq")

brahman_kmer <- read_delim("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/hic_stuff/trio_hic/bostaurus_brahma.kmers",
                         " ",col_names = FALSE)

#angus_kmer <- fread("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/hic_stuff/trio_hic/bostaurus_angus.kmers")

bostaurus_angus_kmer_sorted_uniqcount <- 
  read_delim("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/hic_stuff/trio_hic/bostaurus_angus_kmer_sorted.uniqcount",
             " ",col_names = FALSE)

colnames(bostaurus_angus_kmer_sorted_uniqcount) <- c("count","kmer_freq")

par(mai=c(1.02,2,0.82,0.42))
plot(bostaurus_angus_kmer_sorted_uniqcount$kmer_freq,bostaurus_angus_kmer_sorted_uniqcount$count/1e6,
     type="l",xlab = expression(paste("Count, ",i)), ylab = "sum of particular count in million",
     main = expression(paste(sum(k[i],i==11,N==99)," versus Count i")))
#abline(v=20)

bostaurus_brahma_kmer_sorted_uniqcount <- 
  read_delim("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/hic_stuff/trio_hic/bostaurus_brahma_kmer_sorted.uniqcount",
             " ",col_names = FALSE)

colnames(bostaurus_brahma_kmer_sorted_uniqcount) <- c("count","kmer_freq")

plot(bostaurus_brahma_kmer_sorted_uniqcount$kmer_freq,bostaurus_brahma_kmer_sorted_uniqcount$count/1e6,
     type="l",xlab = expression(paste("Count, ",i)), 
     ylab = expression(paste(sum(k[i],i==11,N==99),"  ,where ",k[i],
                             " refers to kmer with Count i")))
#abline(v=20)

##### this is for the MAPQ based bwa alignment sorting #####
mapq_angus_read1 <- read_tsv("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/hic_stuff/trio_hic/mapq/mapq_angus_read1",
                             col_names = FALSE)

colnames(mapq_angus_read1) <- c("count","mapq")

mapq_angus_read1_no0_no60 <- mapq_angus_read1[2:60,]

plot(mapq_angus_read1_no0_no60$mapq,mapq_angus_read1_no0_no60$count, type="l")
