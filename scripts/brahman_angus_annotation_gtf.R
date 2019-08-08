#------------------------------------------------------
# Program name: brahman_angus_annotation_gtf.R
# Objective: check annotation features
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

setwd("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/EBI_annotation/ensembl_annotation_Maternal/")

#Brahman EBI first draft v1.94
read.gtf(ens,"Bos_indicus_x_bos_taurus_Mom.UOA_Brahman_1.94.gtf")

# counts all annotations on each seqname
tableSeqids(ens)

# create table of genes
Brahman_gene <- getGenePositions(ens)
dim(Brahman_gene)

# gene IDs are unique
length(Brahman_gene$gene_id)
length(unique(Brahman_gene$gene_id))

# use dplyr to create more summaries
# number of genes on each seqname
Brahman_gene_groupByChrGeneBiotype <- Brahman_gene %>% group_by(seqid,gene_biotype) %>% summarise(count=n())

#group all Un as a single term
Brahman_gene_groupByChrGeneBiotype$seqid[grepl("^P.*",Brahman_gene_groupByChrGeneBiotype$seqid)] <- "Un"

Brahman_gene_groupByChrGeneBiotype$seqid[Brahman_gene_groupByChrGeneBiotype$seqid == "X"] <- "XY"

Brahman_gene_groupByChrGeneBiotype <- Brahman_gene_groupByChrGeneBiotype %>% group_by(seqid)

#order_seqid
order_seqid <- c("1","2","3","4","5","6","7","8","9","10","11","12","13",
                 "14","15","16","17","18","19","20","21","22","23","24",
                 "25","26","27","28","29","XY","Un")

Brahman_gene_groupByChrGeneBiotype$seqid <- factor(Brahman_gene_groupByChrGeneBiotype$seqid, levels = order_seqid)

Brahman_gene_groupByChrGeneBiotype$breed <- "Brahman"

#protein_coding only
Brahman_gene_groupByChrGeneBiotype_protein_coding <- Brahman_gene_groupByChrGeneBiotype %>%
  filter(gene_biotype == "protein_coding")

#plot protein_coding only
tiff(filename = "Annotation_Brahman_protein_coding_byChr.tiff",width = 800)
g <- ggplot(Brahman_gene_groupByChrGeneBiotype_protein_coding, aes(x = seqid, y = count))
g <- g + geom_bar(aes(fill = gene_biotype), stat = 'identity', position = 'dodge')
g <- g + xlab("chromosome") + ylab("count")
g
dev.off()

#miRNA and lncRNA 
Brahman_gene_groupByChrGeneBiotype_miRNA_lncRNA <- Brahman_gene_groupByChrGeneBiotype %>%
  filter(gene_biotype == "miRNA" | gene_biotype == "lncRNA")

#plot miRNA and lncRNA 
tiff(filename = "Annotation_Brahman_miRNA_lncRNA_byChr.tiff",width = 800)
g <- ggplot(Brahman_gene_groupByChrGeneBiotype_miRNA_lncRNA, aes(x = seqid, y = count))
g <- g + geom_bar(aes(fill = gene_biotype), stat = 'identity', position = 'dodge')
g <- g + xlab("chromosome") + ylab("count")
g
dev.off()

#Angus
# create ensemblGenome object for storing Ensembl genomic annotation data
ens <- ensemblGenome()

setwd("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/EBI_annotation/ensembl_annotation_Paternal/")

#Angus EBI first draft v0.94.
read.gtf(ens,"Bos_indicus_x_bos_taurus_Pat.Bos_hybrid_PaternalHap_v2.0.94.gtf")

# counts all annotations on each seqname
tableSeqids(ens)

# create table of genes
Angus_gene <- getGenePositions(ens)
dim(Angus_gene)

# gene IDs are unique
length(Angus_gene$gene_id)
length(unique(Angus_gene$gene_id))

# use dplyr to create more summaries
# number of genes on each seqname
Angus_gene_groupByChrGeneBiotype <- Angus_gene %>% group_by(seqid,gene_biotype) %>% summarise(count=n())

#group all Un as a single term
Angus_gene_groupByChrGeneBiotype$seqid[grepl("^P.*",Angus_gene_groupByChrGeneBiotype$seqid)] <- "Un"

Angus_gene_groupByChrGeneBiotype$seqid[Angus_gene_groupByChrGeneBiotype$seqid == "Y"] <- "XY"

Angus_gene_groupByChrGeneBiotype <- Angus_gene_groupByChrGeneBiotype %>% group_by(seqid)

#order_seqid
order_seqid <- c("1","2","3","4","5","6","7","8","9","10","11","12","13",
                 "14","15","16","17","18","19","20","21","22","23","24",
                 "25","26","27","28","29","XY","Un")

Angus_gene_groupByChrGeneBiotype$seqid <- factor(Angus_gene_groupByChrGeneBiotype$seqid, levels = order_seqid)

Angus_gene_groupByChrGeneBiotype$breed <- "Angus"

#protein_coding only
Angus_gene_groupByChrGeneBiotype_protein_coding <- Angus_gene_groupByChrGeneBiotype %>%
  filter(gene_biotype == "protein_coding")

#plot protein_coding only
tiff(filename = "Annotation_Angus_protein_coding_byChr.tiff",width = 800)
g <- ggplot(Angus_gene_groupByChrGeneBiotype_protein_coding, aes(x = seqid, y = count))
g <- g + geom_bar(aes(fill = gene_biotype), stat = 'identity', position = 'dodge')
g <- g + xlab("chromosome") + ylab("count")
g
dev.off()

#miRNA and lncRNA 
Angus_gene_groupByChrGeneBiotype_miRNA_lncRNA <- Angus_gene_groupByChrGeneBiotype %>%
  filter(gene_biotype == "miRNA" | gene_biotype == "lncRNA")

#plot miRNA and lncRNA 
tiff(filename = "Annotation_Angus_miRNA_lncRNA_byChr.tiff",width = 800)
g <- ggplot(Angus_gene_groupByChrGeneBiotype_miRNA_lncRNA, aes(x = seqid, y = count))
g <- g + geom_bar(aes(fill = gene_biotype), stat = 'identity', position = 'dodge')
g <- g + xlab("chromosome") + ylab("count")
g
dev.off()

#Hereford
# create ensemblGenome object for storing Ensembl genomic annotation data
ens <- ensemblGenome()

setwd("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/EBI_annotation/arsucd1.2/")

#Hereford EBI arsucd1.2
read.gtf(ens,"Bos_taurus.ARS-UCD1.2.95.gtf")

# counts all annotations on each seqname
tableSeqids(ens)

# create table of genes
Hereford_gene <- getGenePositions(ens)
dim(Hereford_gene)

# gene IDs are unique
length(Hereford_gene$gene_id)
length(unique(Hereford_gene$gene_id))

# use dplyr to create more summaries
# number of genes on each seqname
Hereford_gene_groupByChrGeneBiotype <- Hereford_gene %>% group_by(seqid,gene_biotype) %>% summarise(count=n())

#group all Un as a single term
Hereford_gene_groupByChrGeneBiotype$seqid[grepl("^N.*",Hereford_gene_groupByChrGeneBiotype$seqid)] <- "Un"

Hereford_gene_groupByChrGeneBiotype$seqid[Hereford_gene_groupByChrGeneBiotype$seqid == "X"] <- "XY"

Hereford_gene_groupByChrGeneBiotype <- Hereford_gene_groupByChrGeneBiotype %>% group_by(seqid)

#order_seqid
order_seqid <- c("1","2","3","4","5","6","7","8","9","10","11","12","13",
                 "14","15","16","17","18","19","20","21","22","23","24",
                 "25","26","27","28","29","XY","Un")

Hereford_gene_groupByChrGeneBiotype$seqid <- factor(Hereford_gene_groupByChrGeneBiotype$seqid, levels = order_seqid)

Hereford_gene_groupByChrGeneBiotype$breed <- "Hereford"

#protein_coding only
Hereford_gene_groupByChrGeneBiotype_protein_coding <- Angus_gene_groupByChrGeneBiotype %>%
  filter(gene_biotype == "protein_coding")

#plot protein_coding only
tiff(filename = "Annotation_Angus_protein_coding_byChr.tiff",width = 800)
g <- ggplot(Hereford_gene_groupByChrGeneBiotype_protein_coding, aes(x = seqid, y = count))
g <- g + geom_bar(aes(fill = gene_biotype), stat = 'identity', position = 'dodge')
g <- g + xlab("chromosome") + ylab("count")
g
dev.off()

#miRNA and lncRNA 
Hereford_gene_groupByChrGeneBiotype_miRNA_lncRNA <- Hereford_gene_groupByChrGeneBiotype %>%
  filter(gene_biotype == "miRNA" | gene_biotype == "lncRNA")

#plot miRNA and lncRNA 
tiff(filename = "Annotation_Angus_miRNA_lncRNA_byChr.tiff",width = 800)
g <- ggplot(Hereford_gene_groupByChrGeneBiotype_miRNA_lncRNA, aes(x = seqid, y = count))
g <- g + geom_bar(aes(fill = gene_biotype), stat = 'identity', position = 'dodge')
g <- g + xlab("chromosome") + ylab("count")
g
dev.off()

##### Combine 3 breeds together #####
merged_breed_groupByChrGeneBiotype <- rbind(Angus_gene_groupByChrGeneBiotype,Brahman_gene_groupByChrGeneBiotype)
merged_breed_groupByChrGeneBiotype <- rbind(merged_breed_groupByChrGeneBiotype,Hereford_gene_groupByChrGeneBiotype)

#protein_coding df
merged_breed_groupByChrGeneBiotype_protein_coding <- merged_breed_groupByChrGeneBiotype %>% 
  filter(gene_biotype == "protein_coding")

#Plot protein_coding per breed
tiff(filename = "Annotation_merged_breed_protein_coding.tiff",width = 600, height = 600)
g <- ggplot(merged_breed_groupByChrGeneBiotype_protein_coding, aes(seqid, count, fill=breed)) 
g <- g + geom_bar(position="dodge",stat="identity") +labs(x="chromosome",y="count")
g <- g + theme_bw()
g
dev.off()

#miRNA df
merged_breed_groupByChrGeneBiotype_miRNA <- merged_breed_groupByChrGeneBiotype %>% 
  filter(gene_biotype == "miRNA")

#Plot protein_coding per breed
tiff(filename = "Annotation_merged_breed_miRNA.tiff",width = 600, height = 600)
g <- ggplot(merged_breed_groupByChrGeneBiotype_miRNA, aes(seqid, count, fill=breed)) 
g <- g + geom_bar(position="dodge",stat="identity") +labs(x="chromosome",y="count")
g <- g + theme_bw()
g
dev.off()

#lncRNA df
merged_breed_groupByChrGeneBiotype_lncRNA <- merged_breed_groupByChrGeneBiotype %>% 
  filter(gene_biotype == "lncRNA") %>% filter(breed != "Hereford")

#Plot lncRNA per breed
tiff(filename = "Annotation_merged_breed_lncRNA.tiff",width = 600, height = 600)
g <- ggplot(merged_breed_groupByChrGeneBiotype_lncRNA, aes(seqid, count, fill=breed)) 
g <- g + geom_bar(position="dodge",stat="identity") +labs(x="chromosome",y="count")
g <- g + theme_bw()
g
dev.off()

