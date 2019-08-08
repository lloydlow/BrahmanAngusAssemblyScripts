# This script is to study how many hits the same genomic region
# (e.g. chr 1 same coordinates overlap) in the brahman and angus
# nonHits that hits umd3.1 at >= 80% identity and >= 80% query overlap

library(readr)
library(dplyr)
library(ggplot2)
library(easyGgplot2)
library(stringr)

dir_hit_umd31 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/onePercentUnassignedReads/"

#file_path <- paste0(dir_hit_umd31,"umd31_nonHits4kb.blnm6qcov")
file_path <- paste0(dir_hit_umd31,"ARS24_nonHits2kb.blnm6qcov")

#Bcos it is not functionised, still var name is hit_umd31
hit_umd31 <- read_tsv(file_path, col_names = FALSE)

#used chromosome instead of sseqid
names(hit_umd31) <- c("qseqid","chromosome","pident","length",
                                     "mismatch", "gapopen", "qstart", "qend", 
                                     "sstart", "send", "evalue", "bitscore", 
                                     "qcovs")

#Option 2: those that hits at qcov >= 95% and percent id >= 85% 
hit_umd31 <- hit_umd31 %>% arrange(chromosome,sstart) #%>% 
#  filter(qcovs >= 95) %>%
#  filter(pident >= 85)

summary(hit_umd31$qcovs)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#80.00   95.00   99.00   96.24   99.00  100.00 
summary(hit_umd31$pident)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#80.00   84.25   86.32   86.08   88.14   91.99

#remove "chromosome_" for easy plot
hit_umd31$chromosome <- gsub("chromosome_","",hit_umd31$chromosome)

hit_umd31[hit_umd31$chromosome == "X",2] <- "30"

#optional code to get rid of leftover & keep the ori
hit_umd31_ori <- hit_umd31
hit_umd31 <- hit_umd31 %>% filter(!str_detect(chromosome,"Left"))

hit_umd31$chromosome <- as.integer(hit_umd31$chromosome)

hit_umd31 <- hit_umd31 %>% arrange(chromosome)

#hit_umd31[hit_umd31$chromosome == 30,2] <- "X"

hit_umd31_X <- hit_umd31 %>% filter(chromosome == 30)
ggplot2.histogram(data=hit_umd31_X,xName= 'sstart',binwidth = 0.1e6)
hit_umd31_1 <- hit_umd31 %>% filter(chromosome == 1)
ggplot2.histogram(data=hit_umd31_1,xName= 'sstart',binwidth = 0.1e6)
hit_umd31_7 <- hit_umd31 %>% filter(chromosome == 7)
ggplot2.histogram(data=hit_umd31_7,xName= 'sstart',binwidth = 0.11e6)

ggplot2.histogram(data=hit_umd31, xName= 'sstart', xtitle="Position",
                  groupName='chromosome', legendPosition="right",
                  faceting=TRUE, facetingVarNames="chromosome",
                  binwidth = 0.1e6,yShowTitle=FALSE,yShowTickLabel=FALSE,
                  hideAxisTicks=TRUE) 

#How long are non-Hits that mapped? this is for >= 4kb
#load("neither_No_Ns.RData")
neither_atleast_4kb_length <- neither_No_Ns %>% filter(length >= 4000) %>%
     select(scaffold,length)

neither_atleast_4kb_length$scaffold_change <- gsub("[[:space:]]RQ.*$","",
                                                   neither_atleast_4kb_length$scaffold)

hit_selc <- neither_atleast_4kb_length$scaffold_change %in% hit_umd31$qseqid
hit_only_DF <- neither_atleast_4kb_length[hit_selc,]
summary(hit_only_DF$length)
sum(hit_only_DF$length)/1e6

#How long are non-Hits that mapped? this is for >= 2kb, ARS24
load("neither_No_Ns.RData")
neither_atleast_2kb_length <- neither_No_Ns %>% filter(length >= 2000) %>%
  select(scaffold,length)

neither_atleast_2kb_length$scaffold_change <- gsub("[[:space:]]RQ.*$","",
                                                   neither_atleast_2kb_length$scaffold)

hit_selc <- neither_atleast_2kb_length$scaffold_change %in% hit_umd31$qseqid
hit_only_DF <- neither_atleast_2kb_length[hit_selc,]
summary(hit_only_DF$length)
sum(hit_only_DF$length)/1e6

#finding universal non-hits and blast nr
#need to use the intersect var from brahman_angus_unassigned_reads_20180213.R
the_intersect_noRQ <- gsub("[[:space:]]RQ.*$","",
                           the_intersect)
universal_non_hit_selc <- !the_intersect_noRQ %in% hit_umd31_ori$qseqid
the_intersect_nonHit <- the_intersect[universal_non_hit_selc]
write_lines(the_intersect_nonHit,
            "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/onePercentUnassignedReads/universal_non_hit_2kb_ars24")


  