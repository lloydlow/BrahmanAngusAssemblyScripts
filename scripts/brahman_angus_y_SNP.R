#------------------------------------------------------
# Program name: brahman_angus_y_SNP.R
# Objective: which contigs might be chrY using Ben  
#           probes
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(readr)
library(dplyr)
library(VennDiagram)

load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/brahman_Angus_mashmap/contig_scaff_UCDv25.RData")

# reading blastn m6 file dam with y_snp_probes.fa 
dir1 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/brahman_Angus_BlastDB/bostaurus_brahman/output/"
path1 <- paste0(dir1,"dam_assembly_cleaned_ySNP.blnm6")

dam_assembly_cleaned_ySNP <- read_tsv(path1,col_names = FALSE)
names(dam_assembly_cleaned_ySNP) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                      "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs")

#loop through blast snp DF to get the mapping to dominette chr
ref <- c()
for (i in 1:nrow(dam_assembly_cleaned_ySNP)){
  logic_selc <- dam_mashmap_modi$Query_name %in% dam_assembly_cleaned_ySNP$sseqid[i]
  ref <- c(ref,dam_mashmap_modi$Ref_name[logic_selc])
}
dam_assembly_cleaned_ySNP$Ref_name <- ref

dam_assembly_cleaned_ySNP <- dam_assembly_cleaned_ySNP %>% filter(pident == 100 & qcovs == 100)
dam_probes <- dam_assembly_cleaned_ySNP$qseqid

# reading blastn m6 file sire with y_snp_probes.fa 
dir2 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/brahman_Angus_BlastDB/bostaurus_angus/output/"
path2 <- paste0(dir2,"sire_assembly_cleaned_ySNP.blnm6")

sire_assembly_cleaned_ySNP <- read_tsv(path2,col_names = FALSE)
names(sire_assembly_cleaned_ySNP) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                      "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs")

#loop through blast snp DF to get the mapping to dominette chr
ref2 <- c()
for (i in 1:nrow(sire_assembly_cleaned_ySNP)){
  logic_selc2 <- sire_mashmap_modi$Query_name %in% sire_assembly_cleaned_ySNP$sseqid[i]
  if (sum(logic_selc2) == 0) {ref2 <- c(ref2,"missing")}
  else ref2 <- c(ref2,sire_mashmap_modi$Ref_name[logic_selc2])
}
sire_assembly_cleaned_ySNP$Ref_name <- ref2

sire_assembly_cleaned_ySNP <- sire_assembly_cleaned_ySNP %>% filter(pident == 100 & qcovs == 100)
sire_probes <- sire_assembly_cleaned_ySNP$qseqid
  
# reading blastn m6 file sire with y_snp_probes.fa 
dir3 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/brahman_Angus_BlastDB/latest_hereford/output/"
path3 <- paste0(dir3,"UCDv25_ySNP.blnm6")

UCDv25_ySNP <- read_tsv(path3,col_names = FALSE)
names(UCDv25_ySNP) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                       "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs")
UCDv25_ySNP <- UCDv25_ySNP %>% filter(pident == 100 & qcovs == 100)
dominette_probe <- UCDv25_ySNP$qseqid

# reading blastn m6 file baylorY with y_snp_probes.fa 
dir4 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/brahman_Angus_BlastDB/CM001061/output/"
path4 <- paste0(dir4,"BtauY_ySNP.blnm6")

BaylorY_ySNP <- read_tsv(path4,col_names = FALSE)
names(BaylorY_ySNP) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                        "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs")
BaylorY_ySNP <- BaylorY_ySNP %>% filter(pident == 100 & qcovs == 100)
BaylorY_probe <- BaylorY_ySNP$qseqid

#probes that mapped to sire only
sire_probes_only <- setdiff(setdiff(sire_probes,dam_probes),dominette_probe)

#probes that mapped overlapped in dam, sire, dominette
overlap_all_probes <- Reduce(intersect,list(dam_probes,sire_probes,dominette_probe))

venn.diagram(
  x = list(dam_probes, sire_probes, dominette_probe),
  category.names = c("Dam" , "Sire" , "Dominette"),fill = c('yellow', 'purple', 'green'),
  filename = 'venn_dam_sire_dom.png')

#find which contigs have overlap_all_probes
#sire
sire_selc <- sire_assembly_cleaned_ySNP$qseqid %in% overlap_all_probes
sire_assembly_cleaned_ySNP_overlap_all <- sire_assembly_cleaned_ySNP[sire_selc,]

#dam
dam_selc <- dam_assembly_cleaned_ySNP$qseqid %in% overlap_all_probes
dam_assembly_cleaned_ySNP_overlap_all <- dam_assembly_cleaned_ySNP[dam_selc,]

#dominette
dom_selc <- UCDv25_ySNP$qseqid %in% overlap_all_probes
UCDv25_ySNP_overlap_all <- UCDv25_ySNP[dom_selc,]

#overlap summary 
sire_assembly_cleaned_ySNP_overlap_all %>% 
  group_by(Ref_name) %>% summarise(count = n())

dam_assembly_cleaned_ySNP_overlap_all %>% 
  group_by(Ref_name) %>% summarise(count = n())

#read in Ben's reads_y.found matched Baylor Y
#*********** need to map the probes on my own ***************
reads_y <- read_csv("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/y_SNPs/reads_y.found",
                    col_names = TRUE)

#sire_probe_baylorY
sire_probe_baylorY <- intersect(reads_y$name,sire_probes_only)
venn.diagram(
  x = list(reads_y$name, sire_probes_only),
  category.names = c("Baylor Y" , "Sire only"),fill = c('yellow', 'purple'),
  filename = 'venn_baylorY_sire.png')

#use sire_probe_baylorY to filter sire_assembly_cleaned_ySNP
sire_assembly_cleaned_ySNP_baylor_selc <- sire_assembly_cleaned_ySNP$qseqid %in% sire_probe_baylorY
sire_assembly_cleaned_ySNP_baylor <- sire_assembly_cleaned_ySNP[sire_assembly_cleaned_ySNP_baylor_selc,]

#how many contigs in the missing
missing_sire_contig_probe <- sire_assembly_cleaned_ySNP_baylor %>% 
  filter(Ref_name == "missing") %>% group_by(sseqid) %>%
  summarise(count = n()) %>% arrange(desc(count))

#read in gapLength file for this sire contigs
sire_assembly.cleaned_No_Ns <- read_tsv("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/fasta_file_order_checkBases_N/sire_assembly.cleaned_No_Ns.rls",
         col_names = FALSE)
names(sire_assembly.cleaned_No_Ns) <- c("name","nothing","gap","length","nothing2")

length_contig <- c()
for (i in 1:nrow(missing_sire_contig_probe)){
  logic_selc <- sire_assembly.cleaned_No_Ns$name %in% missing_sire_contig_probe$sseqid[i]
  length_contig <- c(length_contig,sire_assembly.cleaned_No_Ns$length[logic_selc])
}
missing_sire_contig_probe$length <- length_contig

venn.diagram(
  x = list(dam_probes, sire_probes, dominette_probe, reads_y$name),
  category.names = c("Dam" , "Sire" , "Dominette", "Baylor Y v1"),fill = c('yellow', 'purple', 'green', 'blue'),
  imagetype = "png", filename = 'venn_dam_sire_dom_baylorY.png')

venn.diagram(
  x = list(dam_probes, sire_probes, dominette_probe, BaylorY_probe),
  category.names = c("Dam" , "Sire" , "Dominette", "Baylor Y v2"),fill = c('yellow', 'purple', 'green', 'red'),
  imagetype = "png",filename = 'venn_dam_sire_dom_baylorYv2.png')


#probes in Baylor Y but not our sire
BaylorY_setdiff <- setdiff(BaylorY_probe,sire_probes)
setdiff_logic1 <- BaylorY_ySNP$qseqid %in% BaylorY_setdiff
BaylorY_ySNP_setdiffBaylorY <- BaylorY_ySNP[setdiff_logic1,]
BaylorY_ySNP_setdiffBaylorY %>% group_by(qseqid) %>% summarise(count = n()) %>% arrange(desc(count))
