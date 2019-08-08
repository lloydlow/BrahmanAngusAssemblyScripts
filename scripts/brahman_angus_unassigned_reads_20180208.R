# This script here is to look at the 1% unassigned reads in 
# the brahman and angus assembly. These reads are never used in 
# the assembly so potentially some valuable info may be lost here.

library(xlsx)
library(readr)
library(dplyr)
library(ggplot2)

# reading *No_Ns.rls file generated from neither.fa
dir1 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/onePercentUnassignedReads/"
path1 <- paste0(dir1,"neither_No_Ns.rls")

neither_No_Ns <- read_tsv(path1,col_names = FALSE)
names(neither_No_Ns) <- c("scaffold","nothing","gap","length","perc_gap")

neither_No_Ns <- neither_No_Ns %>%
  select(scaffold,gap,length)

#neither_No_Ns %>% filter(length >= 5000)

#reading *No_Ns.rls from brahman only reads
path2 <- paste0(dir1,"brahman_reads_No_Ns.rls")

brahman_reads_No_Ns <- read_tsv(path2,col_names = FALSE)

names(brahman_reads_No_Ns) <- c("scaffold","nothing","gap","length","perc_gap")

brahman_reads_No_Ns <- brahman_reads_No_Ns %>%
  select(scaffold,gap,length)

#reading *No_Ns.rls from angus only reads
path3 <- paste0(dir1,"angus_reads_No_Ns.rls")

angus_reads_No_Ns <- read_tsv(path3,col_names = FALSE)

names(angus_reads_No_Ns) <- c("scaffold","nothing","gap","length","perc_gap")

angus_reads_No_Ns <- angus_reads_No_Ns %>%
  select(scaffold,gap,length)

summary(neither_No_Ns$length)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#35     367     831    1266    1677   40545 

summary(brahman_reads_No_Ns$length)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#35    4339    8973   10860   15657  137782 

summary(angus_reads_No_Ns$length)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#35    5065    9996   11662   16567  172752 

summary(assembly_reads$length)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#35    4672    9466   11246   16109  172752 

assembly_reads <- rbind(brahman_reads_No_Ns,angus_reads_No_Ns)

all_reads <- rbind(assembly_reads,neither_No_Ns)

label_used_reads <- rep("Brahman + Angus",nrow(assembly_reads))

label_unused_reads <- rep("Unused",nrow(neither_No_Ns))

all_labels <- c(label_used_reads,label_unused_reads)

all_reads$label <- all_labels

#Length comparisons
#boxplot(neither_No_Ns$length,brahman_reads_No_Ns$length,angus_reads_No_Ns$length)

#boxplot(neither_No_Ns$length,assembly_reads$length)

png(filename = "unassigned_length_comparison.png",height = 600,width = 600)
g <- ggplot(all_reads, aes(x=as.factor(label),y=length))
g <- g + geom_boxplot(fill="slateblue", alpha=0.2)
g <- g + xlab("Reads") + ylab("Read length") + theme(text = element_text(size=20))
g
dev.off()

#wilcox.test(neither_No_Ns$length,brahman_reads_No_Ns$length)
#data:  neither_No_Ns$length and brahman_reads_No_Ns$length
#W = 2.9358e+12, p-value < 2.2e-16
#wilcox.test(neither_No_Ns$length,angus_reads_No_Ns$length)
#data:  neither_No_Ns$length and angus_reads_No_Ns$length
#W = 2.1608e+12, p-value < 2.2e-16
#wilcox.test(brahman_reads_No_Ns$length,angus_reads_No_Ns$length)
#data:  brahman_reads_No_Ns$length and angus_reads_No_Ns$length
#W = 1.213e+14, p-value < 2.2e-16

#Only this Mann Whitney test is used
wilcox.test(length ~ label,data = all_reads)
#wilcox.test(assembly_reads$length,neither_No_Ns$length)
#W = 9.5919e+13, p-value < 2.2e-16

#filter the fasta id of those >= 2kb
neither_atleast_2kb <- neither_No_Ns %>% filter(length >= 2000) %>%
  select(scaffold)
save(neither_atleast_2kb,file = "neither_atleast_2kb.RData")
outpath <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/onePercentUnassignedReads/neither_subset_fasta"
write_lines(neither_atleast_2kb$scaffold,outpath)

# neither_atleast_4kb <- neither_No_Ns %>% filter(length >= 4000) %>%
#   select(scaffold)
# save(neither_atleast_4kb,file = "neither_atleast_4kb.RData")
save(neither_No_Ns,file="neither_No_Ns.RData")

