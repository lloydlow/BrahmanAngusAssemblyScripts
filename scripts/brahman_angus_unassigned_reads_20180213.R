# This script is to look at those reads that dont meet the
# 80% query cov and 80% identity
# subset these reads out with seqtk subseq and then blast to
# cattle UMD3.1

library(readr)
library(dplyr)
library(VennDiagram)

#load("neither_atleast_4kb.RData")
load("neither_atleast_2kb.RData")

dir_angus <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/brahman_Angus_BlastDB/bostaurus_angus/output/"

dir_brahman <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/brahman_Angus_BlastDB/bostaurus_brahman/output/"

#file_path_angus <- paste0(dir_angus,"bostaurus_angus_neither4kb.blnm6qcov")
file_path_angus <- paste0(dir_angus,"bostaurus_angus_neither2kb.blnm6qcov")

#file_path_brahman <- paste0(dir_brahman,"bostaurus_brahma_neither4kb.blnm6qcov")
file_path_brahman <- paste0(dir_brahman,"bostaurus_brahma_neither2kb.blnm6qcov")

#read in angus
unassigned_hitting_angus <- read_tsv(file_path_angus,col_names = FALSE)
names(unassigned_hitting_angus) <- c("qseqid","sseqid","pident","length",
                                     "mismatch", "gapopen", "qstart", "qend", 
                                     "sstart", "send", "evalue", "bitscore", 
                                     "qcovs")

#read in brahman
unassigned_hitting_brahman <- read_tsv(file_path_brahman,col_names = FALSE)
names(unassigned_hitting_brahman) <- c("qseqid","sseqid","pident","length",
                                       "mismatch", "gapopen", "qstart", "qend", 
                                       "sstart", "send", "evalue", "bitscore", 
                                       "qcovs")

#one DF name substitution to make it easier for remaining analysis
neither_atleast_4kb <- neither_atleast_2kb

#change scaffold for matching
neither_atleast_4kb$scaffold_change <- gsub("[[:space:]]RQ.*$","",neither_atleast_4kb$scaffold)

#Look for non-hits in angus
angus_selc <- neither_atleast_4kb$scaffold_change %in% unassigned_hitting_angus$qseqid 

angus_non_hits <- neither_atleast_4kb[!angus_selc,]

#Look for non-hits in brahman
brahman_selc <- neither_atleast_4kb$scaffold_change %in% unassigned_hitting_brahman$qseqid 

brahman_non_hits <- neither_atleast_4kb[!brahman_selc,]

#look for overlap of non-hits in both angus and brahman
different_in_angus <- setdiff(angus_non_hits$scaffold,brahman_non_hits$scaffold)
different_in_brahman <- setdiff(brahman_non_hits$scaffold,angus_non_hits$scaffold)
the_intersect <- intersect(angus_non_hits$scaffold,brahman_non_hits$scaffold)

#write_lines(the_intersect,
#            "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/onePercentUnassignedReads/nonHits4kb")
write_lines(the_intersect,
            "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/onePercentUnassignedReads/nonHits2kb")

grid.newpage()
draw.pairwise.venn(nrow(brahman_non_hits), nrow(angus_non_hits), length(the_intersect),
                   category = c("Brahman", "Angus"),
                   fill = c("light blue", "pink"), alpha = rep(0.5, 2))
#                   lty = rep("blank", 2), fill = c("light blue", "pink"), 
#                   alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2))
