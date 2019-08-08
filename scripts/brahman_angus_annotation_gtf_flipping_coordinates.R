#------------------------------------------------------
# Program name: brahman_angus_annotation_gtf_flipping_coordinates.R
# Objective: create reverse comp gtf coordinates for ribbon plot
#           
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(readr)
library(dplyr)
library(ggplot2)
library(refGenome)

# names(brahmanGTF) <- c("seqname","source","feature","start","end","score","strand","frame","attribute")

#using refGenome to read in genome gtf

#Brahman
# create ensemblGenome object for storing Ensembl genomic annotation data
ens <- ensemblGenome()

setwd("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/EBI_annotation/ensembl_annotation_Maternal_20190320/")

#Brahman EBI first draft v1.96
read.gtf(ens,"Bos_indicus_x_bos_taurus_Mom.UOA_Brahman_1.96.gtf")

# counts all annotations on each seqname
tableSeqids(ens)

# # create table of genes
Brahman_gene <- getGenePositions(ens)
dim(Brahman_gene)
# 
# # gene IDs are unique
# length(Brahman_gene$gene_id)
# length(unique(Brahman_gene$gene_id))
# 
# # use dplyr to create more summaries
# # number of genes on each seqname
# Brahman_gene_groupByChrGeneBiotype <- Brahman_gene %>% group_by(seqid,gene_biotype) %>% summarise(count=n())
# 
# #group all Un as a single term
# Brahman_gene_groupByChrGeneBiotype$seqid[grepl("^P.*",Brahman_gene_groupByChrGeneBiotype$seqid)] <- "Un"
# 
# Brahman_gene_groupByChrGeneBiotype$seqid[Brahman_gene_groupByChrGeneBiotype$seqid == "X"] <- "XY"
# 
# Brahman_gene_groupByChrGeneBiotype <- Brahman_gene_groupByChrGeneBiotype %>% group_by(seqid)
# 
# #order_seqid
# order_seqid <- c("1","2","3","4","5","6","7","8","9","10","11","12","13",
#                  "14","15","16","17","18","19","20","21","22","23","24",
#                  "25","26","27","28","29","XY","Un")
# 
# Brahman_gene_groupByChrGeneBiotype$seqid <- factor(Brahman_gene_groupByChrGeneBiotype$seqid, levels = order_seqid)
# 
# Brahman_gene_groupByChrGeneBiotype$breed <- "Brahman"
# 
# #protein_coding only
# Brahman_gene_groupByChrGeneBiotype_protein_coding <- Brahman_gene_groupByChrGeneBiotype %>%
#   filter(gene_biotype == "protein_coding")

#get chr15 and rev comp annotation
#length of chr15 84270218
#region of interest is 3748952 to 5140465
left <- 3748952
right <- 5140465
length <- 5140465 - 3748952 + 1
selectedRegionName <- "UOA_Brahman_1_chr15_3748952_5140465_rc"

Brahman_gene_chr15 <- Brahman_gene %>% filter(seqid == 15) %>% 
  filter(start >= left) %>% filter(end <= right)

#new zoom in coordinates
Brahman_gene_chr15 <- Brahman_gene_chr15 %>% mutate(newstart = start -left + 1) %>%
  mutate(newend = end -left + 1)

#rev comp
Brahman_gene_chr15 <- Brahman_gene_chr15 %>% mutate(newstart = length - newstart +1) %>% 
  mutate(newend = length - newend +1)
Brahman_gene_chr15$newstrand <- ifelse(Brahman_gene_chr15$strand == "+","-","+")

Brahman_gene_chr15$selectedRegionName <- selectedRegionName

newcoorDF <- Brahman_gene_chr15 %>% 
  select(selectedRegionName,newend,newstart,gene_id,score,newstrand,gene_biotype)

#write out the results as bed for overlap in ribbon
write_tsv(newcoorDF,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/misassembly_issues/interesting_region/brahman_chr15_FADS2P1/synteny_draw/UOA_Brahman_1_chr15_3748952_5140465_rc.bed",
          col_names = FALSE)

#Angus
# create ensemblGenome object for storing Ensembl genomic annotation data
ens <- ensemblGenome()

setwd("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/EBI_annotation/ensembl_annotation_Paternal_20190320/")

#Brahman EBI first draft v1.96
read.gtf(ens,"Bos_indicus_x_bos_taurus_Pat.UOA_Angus_1.96.gtf")

# counts all annotations on each seqname
tableSeqids(ens)

# # create table of genes
Angus_gene <- getGenePositions(ens)
dim(Angus_gene)

#get chr15 annotation
#length of chr15 84012196
#region of interest is 78799177 to 80168904
left <- 78799177
right <- 80168904
length <- 80168904 - 78799177 + 1
selectedRegionName <- "UOA_Angus_1_chr15_78799177_80168904"

Angus_gene_chr15 <- Angus_gene %>% filter(seqid == 15) %>% 
  filter(start >= left) %>% filter(end <= right)

#new zoom in coordinates
Angus_gene_chr15 <- Angus_gene_chr15 %>% mutate(newstart = start -left + 1) %>%
  mutate(newend = end -left + 1)

Angus_gene_chr15$selectedRegionName <- selectedRegionName

newcoorAngusDF <- Angus_gene_chr15 %>% 
  select(selectedRegionName,newstart,newend,gene_id,score,strand,gene_biotype)

#write out the results as bed for overlap in ribbon
write_tsv(newcoorAngusDF,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/misassembly_issues/interesting_region/brahman_chr15_FADS2P1/synteny_draw/UOA_Angus_1_chr15_78799177_80168904.bed",
          col_names = FALSE)

#Herefore
# create ensemblGenome object for storing Ensembl genomic annotation data
ens <- ensemblGenome()

setwd("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/EBI_annotation/arsucd1.2/")

#Brahman EBI first draft v1.96
read.gtf(ens,"Bos_taurus.ARS-UCD1.2.95.gtf")

# counts all annotations on each seqname
tableSeqids(ens)

# # create table of genes
Hereford_gene <- getGenePositions(ens)
dim(Hereford_gene)

#get chr15 annotation
#length of chr15 85007780
#region of interest is 78791037 to 80120961
left <- 78791037
right <- 80120961
length <- 80120961 - 78791037 + 1
selectedRegionName <- "ARSUCD1_2_chr15_78791037_80120961"

Hereford_gene_chr15 <- Hereford_gene %>% filter(seqid == 15) %>% 
  filter(start >= left) %>% filter(end <= right)

#new zoom in coordinates
Hereford_gene_chr15 <- Hereford_gene_chr15 %>% mutate(newstart = start -left + 1) %>%
  mutate(newend = end -left + 1)

Hereford_gene_chr15$selectedRegionName <- selectedRegionName

newcoorHerefordDF <- Hereford_gene_chr15 %>% 
  select(selectedRegionName,newstart,newend,gene_id,score,strand,gene_biotype)

#write out the results as bed for overlap in ribbon
write_tsv(newcoorHerefordDF,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/misassembly_issues/interesting_region/brahman_chr15_FADS2P1/synteny_draw/ARSUCD1_2_chr15_78791037_80120961.bed",
          col_names = FALSE)

