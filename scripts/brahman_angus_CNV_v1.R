#------------------------------------------------------
# Program name: brahman_angus_CNV_v1.R
# Objective: analyse Derek CN per window, overlapped with  
#     gene anno, Vst     
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)

#read in *bed and *melt files of CNVRs
# path to CNV results
dir1 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/CopyNumberVariation/20190530/vst_gene_lists/"

# angus
# reading data files
path1 <- paste0(dir1,"angus.cnvrs_windows_vst_genes.tvi.bed")
path2 <- paste0(dir1,"angus.cnvrs_windows_vst_genes.tvi.melt")

angus.cnvrs_windows_vst_genes.tvi.bed <- read_tsv(path1,col_names = FALSE, col_types = "cdddc")
colnames(angus.cnvrs_windows_vst_genes.tvi.bed) <- c("chr","start","end","Vst","gene")

angus.cnvrs_windows_vst_genes.tvi.melt <- read_tsv(path2,col_names = TRUE)

angus.cnvrs_windows_vst_genes.tvi.melt$Gene <- gsub("\\(","",angus.cnvrs_windows_vst_genes.tvi.melt$Gene)
angus.cnvrs_windows_vst_genes.tvi.melt$Gene <- gsub("\\)","",angus.cnvrs_windows_vst_genes.tvi.melt$Gene)

angus.cnvrs_windows_vst_genes.tvi.melt <- angus.cnvrs_windows_vst_genes.tvi.melt %>% 
  separate(Gene, c("gene","Vst"),sep = " ",extra = "drop", fill = "right", convert = TRUE)

# brahman
# reading data files
path1 <- paste0(dir1,"brahman.cnvrs_windows_vst_genes.tvi.bed")
path2 <- paste0(dir1,"brahman.cnvrs_windows_vst_genes.tvi.melt")

brahman.cnvrs_windows_vst_genes.tvi.bed <- read_tsv(path1,col_names = FALSE, col_types = "cdddc")
colnames(brahman.cnvrs_windows_vst_genes.tvi.bed) <- c("chr","start","end","Vst","gene")

brahman.cnvrs_windows_vst_genes.tvi.melt <- read_tsv(path2,col_names = TRUE)

brahman.cnvrs_windows_vst_genes.tvi.melt$Gene <- gsub("\\(","",brahman.cnvrs_windows_vst_genes.tvi.melt$Gene)
brahman.cnvrs_windows_vst_genes.tvi.melt$Gene <- gsub("\\)","",brahman.cnvrs_windows_vst_genes.tvi.melt$Gene)

brahman.cnvrs_windows_vst_genes.tvi.melt <- brahman.cnvrs_windows_vst_genes.tvi.melt %>% 
  separate(Gene, c("gene","Vst"),sep = " ",extra = "drop", fill = "right", convert = TRUE)

# arsucd
# reading data files #this one with ensembl r96 anno
path1 <- paste0(dir1,"arsucd.cnvrs_windows_vst_ensembl.bed")
path2 <- paste0(dir1,"arsucd.cnvrs_windows_vst_ensembl.melt")

arsucd.cnvrs_windows_vst_genes.tvi.bed <- read_tsv(path1,col_names = FALSE, col_types = "cdddc")
colnames(arsucd.cnvrs_windows_vst_genes.tvi.bed) <- c("chr","start","end","Vst","gene")

arsucd.cnvrs_windows_vst_genes.tvi.melt <- read_tsv(path2,col_names = TRUE)

arsucd.cnvrs_windows_vst_genes.tvi.melt$Gene <- gsub("\\(","",arsucd.cnvrs_windows_vst_genes.tvi.melt$Gene)
arsucd.cnvrs_windows_vst_genes.tvi.melt$Gene <- gsub("\\)","",arsucd.cnvrs_windows_vst_genes.tvi.melt$Gene)

arsucd.cnvrs_windows_vst_genes.tvi.melt <- arsucd.cnvrs_windows_vst_genes.tvi.melt %>% 
  separate(Gene, c("gene","Vst"),sep = " ",extra = "drop", fill = "right", convert = TRUE)

############################# Annotation ###########################
# path to Angus ensembl r96 anno
dir2 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/EBI_annotation/ensembl_annotation_Paternal_20190620/"

#Angus EBI anno #chr only
ebi_angus_anno_gtf <- read_tsv(paste0(dir2,"Bos_taurus_hybrid.UOA_Angus_1.97.chr.gtf"),
                               comment = "#", col_names = FALSE,
                               col_types = list(col_character(),
                                                col_character(),
                                                col_character(),
                                                col_integer(),
                                                col_integer(),
                                                col_character(),
                                                col_character(),
                                                col_character(),
                                                col_character()))

colnames(ebi_angus_anno_gtf) <- c("seqname","source","feature","start","end","score","strand","frame","attributes")

ebi_angus_anno_gtf_gene <- ebi_angus_anno_gtf %>% as.tbl() %>% filter(feature == "gene") %>% dplyr::select(seqname,start,end,strand,attributes)

names(ebi_angus_anno_gtf_gene)[1] <- "chr"

ebi_angus_anno_gtf_gene <- ebi_angus_anno_gtf_gene %>% 
  separate(attributes, c("gene_id","gene_version","gene_name"),sep = ";",extra = "drop", fill = "right")

ebi_angus_anno_gtf_gene$gene_id <- gsub("gene_id \"","",ebi_angus_anno_gtf_gene$gene_id)
ebi_angus_anno_gtf_gene$gene_id <- gsub("\"","",ebi_angus_anno_gtf_gene$gene_id)

#"gene_sourceensembl" means no gene name
ebi_angus_anno_gtf_gene$gene_name <- gsub("gene_name \"","",ebi_angus_anno_gtf_gene$gene_name, fixed = TRUE)
ebi_angus_anno_gtf_gene$gene_name <- gsub("\"","",ebi_angus_anno_gtf_gene$gene_name, fixed = TRUE)
ebi_angus_anno_gtf_gene$gene_name <- gsub(" ","",ebi_angus_anno_gtf_gene$gene_name, fixed = TRUE)

ebi_angus_anno_gtf_gene <- ebi_angus_anno_gtf_gene %>% dplyr::select(chr:gene_id,gene_name)

# path to Brahman ensembl r96 anno
dir3 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/EBI_annotation/ensembl_annotation_Maternal_20190620/"

#Brahman EBI anno #chr only
ebi_brahman_anno_gtf <- read_tsv(paste0(dir3,"Bos_indicus_hybrid.UOA_Brahman_1.97.chr.gtf"),
                                 comment = "#", col_names = FALSE,
                                 col_types = list(col_character(),
                                                  col_character(),
                                                  col_character(),
                                                  col_integer(),
                                                  col_integer(),
                                                  col_character(),
                                                  col_character(),
                                                  col_character(),
                                                  col_character()))

colnames(ebi_brahman_anno_gtf) <- c("seqname","source","feature","start","end","score","strand","frame","attributes")

ebi_brahman_anno_gtf_gene <- ebi_brahman_anno_gtf %>% as.tbl() %>% filter(feature == "gene") %>% dplyr::select(seqname,start,end,strand,attributes)

names(ebi_brahman_anno_gtf_gene)[1] <- "chr"

ebi_brahman_anno_gtf_gene <- ebi_brahman_anno_gtf_gene %>% 
  separate(attributes, c("gene_id","gene_version","gene_name"),sep = ";",extra = "drop", fill = "right")

ebi_brahman_anno_gtf_gene$gene_id <- gsub("gene_id \"","",ebi_brahman_anno_gtf_gene$gene_id)
ebi_brahman_anno_gtf_gene$gene_id <- gsub("\"","",ebi_brahman_anno_gtf_gene$gene_id)

#"gene_sourceensembl" means no gene name
ebi_brahman_anno_gtf_gene$gene_name <- gsub("gene_name \"","",ebi_brahman_anno_gtf_gene$gene_name, fixed = TRUE)
ebi_brahman_anno_gtf_gene$gene_name <- gsub("\"","",ebi_brahman_anno_gtf_gene$gene_name, fixed = TRUE)
ebi_brahman_anno_gtf_gene$gene_name <- gsub(" ","",ebi_brahman_anno_gtf_gene$gene_name, fixed = TRUE)

ebi_brahman_anno_gtf_gene <- ebi_brahman_anno_gtf_gene %>% dplyr::select(chr:gene_id,gene_name)

# path to Arsucd ensembl r96 anno
dir4 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/EBI_annotation/arsucd1.2/"

#Arsucd EBI anno #chr only
ebi_arsucd_anno_gtf <- read_tsv(paste0(dir4,"Bos_taurus.ARS-UCD1.2.96.chr.gtf"),
                                comment = "#", col_names = FALSE,
                                col_types = list(col_character(),
                                                 col_character(),
                                                 col_character(),
                                                 col_integer(),
                                                 col_integer(),
                                                 col_character(),
                                                 col_character(),
                                                 col_character(),
                                                 col_character()))

colnames(ebi_arsucd_anno_gtf) <- c("seqname","source","feature","start","end","score","strand","frame","attributes")

ebi_arsucd_anno_gtf_gene <- ebi_arsucd_anno_gtf %>% as.tbl() %>% filter(feature == "gene") %>% dplyr::select(seqname,start,end,strand,attributes)

names(ebi_arsucd_anno_gtf_gene)[1] <- "chr"

ebi_arsucd_anno_gtf_gene <- ebi_arsucd_anno_gtf_gene %>% 
  separate(attributes, c("gene_id","gene_version","gene_name"),sep = ";",extra = "drop", fill = "right")

ebi_arsucd_anno_gtf_gene$gene_id <- gsub("gene_id \"","",ebi_arsucd_anno_gtf_gene$gene_id)
ebi_arsucd_anno_gtf_gene$gene_id <- gsub("\"","",ebi_arsucd_anno_gtf_gene$gene_id)

#"gene_sourceensembl" means no gene name
ebi_arsucd_anno_gtf_gene$gene_name <- gsub("gene_name \"","",ebi_arsucd_anno_gtf_gene$gene_name, fixed = TRUE)
ebi_arsucd_anno_gtf_gene$gene_name <- gsub("\"","",ebi_arsucd_anno_gtf_gene$gene_name, fixed = TRUE)
ebi_arsucd_anno_gtf_gene$gene_name <- gsub(" ","",ebi_arsucd_anno_gtf_gene$gene_name, fixed = TRUE)

ebi_arsucd_anno_gtf_gene <- ebi_arsucd_anno_gtf_gene %>% dplyr::select(chr:gene_id,gene_name)

############################# Overlap bed with Annotation- Manhattan ###########################
#angus
#loop thro bed to get annotation
df <- data.frame(strand=character(),gene_name=character(),stringsAsFactors=FALSE)

for(i in 1:nrow(angus.cnvrs_windows_vst_genes.tvi.bed)){
  row_i <- which(angus.cnvrs_windows_vst_genes.tvi.bed$gene[i] == ebi_angus_anno_gtf_gene$gene_id)
  df_add <- ebi_angus_anno_gtf_gene[row_i,c(4,6)]
  df <- rbind(df,df_add)
}

#combine bed with annotation
angus.cnvrs_windows_vst_genes.tvi.bed_comb <- cbind(angus.cnvrs_windows_vst_genes.tvi.bed,df)

#order chr
angus.cnvrs_windows_vst_genes.tvi.bed_comb$chr <- as.character(angus.cnvrs_windows_vst_genes.tvi.bed_comb$chr)

order_chr <- c(as.character(1:29),"Y")

angus.cnvrs_windows_vst_genes.tvi.bed_comb$chr <- factor(angus.cnvrs_windows_vst_genes.tvi.bed_comb$chr,
                                     levels = order_chr)

#plot
tiff(filename = "FigFinal_angus_CNV_Manhattan_like.tiff",width = 1000,height = 300)
g <- ggplot(data = angus.cnvrs_windows_vst_genes.tvi.bed_comb, aes(x = chr,y = Vst))
g <- g + geom_jitter(aes(colour = Vst),width = 0.2) + labs(x = "Chromosome") + scale_y_continuous(expand = c(0, 0.1)) + theme_bw()
g 
dev.off()

#brahman
#loop thro bed to get annotation
df <- data.frame(strand=character(),gene_name=character(),stringsAsFactors=FALSE)

for(i in 1:nrow(brahman.cnvrs_windows_vst_genes.tvi.bed)){
  row_i <- which(brahman.cnvrs_windows_vst_genes.tvi.bed$gene[i] == ebi_brahman_anno_gtf_gene$gene_id)
  df_add <- ebi_brahman_anno_gtf_gene[row_i,c(4,6)]
  df <- rbind(df,df_add)
}

#combine bed with annotation
brahman.cnvrs_windows_vst_genes.tvi.bed_comb <- cbind(brahman.cnvrs_windows_vst_genes.tvi.bed,df)

#order chr
brahman.cnvrs_windows_vst_genes.tvi.bed_comb$chr <- as.character(brahman.cnvrs_windows_vst_genes.tvi.bed_comb$chr)

order_chr <- c(as.character(1:29),"X")

brahman.cnvrs_windows_vst_genes.tvi.bed_comb$chr <- factor(brahman.cnvrs_windows_vst_genes.tvi.bed_comb$chr,
                                                           levels = order_chr)

#plot
tiff(filename = "FigFinal_brahman_CNV_Manhattan_like.tiff",width = 1000,height = 300)
g <- ggplot(data = brahman.cnvrs_windows_vst_genes.tvi.bed_comb, aes(x = chr,y = Vst))
g <- g + geom_jitter(aes(colour = Vst),width = 0.2) + labs(x = "Chromosome") + scale_y_continuous(expand = c(0, 0.1)) + theme_bw()
g
dev.off()

#arsucd
#loop thro bed to get annotation
df <- data.frame(strand=character(),gene_name=character(),stringsAsFactors=FALSE)

for(i in 1:nrow(arsucd.cnvrs_windows_vst_genes.tvi.bed)){
  row_i <- which(arsucd.cnvrs_windows_vst_genes.tvi.bed$gene[i] == ebi_arsucd_anno_gtf_gene$gene_id)
  df_add <- ebi_arsucd_anno_gtf_gene[row_i,c(4,6)]
  df <- rbind(df,df_add)
}

#combine bed with annotation
arsucd.cnvrs_windows_vst_genes.tvi.bed_comb <- cbind(arsucd.cnvrs_windows_vst_genes.tvi.bed,df)

#order chr
arsucd.cnvrs_windows_vst_genes.tvi.bed_comb$chr <- as.character(arsucd.cnvrs_windows_vst_genes.tvi.bed_comb$chr)

order_chr <- c(as.character(1:29),"X")

arsucd.cnvrs_windows_vst_genes.tvi.bed_comb$chr <- factor(arsucd.cnvrs_windows_vst_genes.tvi.bed_comb$chr,
                                                          levels = order_chr)

#plot
tiff(filename = "FigFinal_arsucd_CNV_Manhattan_like.tiff",width = 1000,height = 300)
g <- ggplot(data = arsucd.cnvrs_windows_vst_genes.tvi.bed_comb, aes(x = chr,y = Vst))
g <- g + geom_jitter(aes(colour = Vst),width = 0.2) + labs(x = "Chromosome") + scale_y_continuous(expand = c(0, 0.1)) + theme_bw()
g
dev.off()

############################# Overlap bed with Annotation - Boxplot ###########################
#modify melt table for plotting
angus.cnvrs_windows_vst_genes.tvi.melt <- angus.cnvrs_windows_vst_genes.tvi.melt %>% 
  mutate(Population = ifelse(Pop > 1, "Indicine","Taurine"))

angus.cnvrs_windows_vst_genes.tvi.melt$Population <- as.factor(angus.cnvrs_windows_vst_genes.tvi.melt$Population)

angus.cnvrs_windows_vst_genes.tvi.melt_merge <- merge(angus.cnvrs_windows_vst_genes.tvi.melt,
                                                      angus.cnvrs_windows_vst_genes.tvi.bed_comb,by="gene")

colnames(angus.cnvrs_windows_vst_genes.tvi.melt_merge)[2] <- "Vst"

#looking at significant genes in Angus ASM
#filter for Vst more than X among those with > 1.5 CN
#if X is 0.3, we have 19 genes, we can filter out Y
angus.cnvrs_windows_vst_genes.tvi.melt_sig_genes <- 
  angus.cnvrs_windows_vst_genes.tvi.melt_merge %>% filter(Vst > 0.3) %>% filter(chr != "Y")

#change name for better plot
# ENSBIXG00000001701 is WC1-a
# ENSBIXG00000002096 is TMPRSS11D
# ENSBIXG00000004482 is T-cell receptor alpha chain
# ENSBIXG00000014359 is WC1-b
# ENSBIXG00000030034 is HSFY-a
# ENSBIXG00000030039 is HSFY-b
# ENSBIXG00000030042 is Ubiquitin-conjugating enzyme E2D-a 
# ENSBIXG00000030047 is Ubiquitin-conjugating enzyme E2D-b
# ENSBIXG00000030050 is HSFY-c
# ENSBIXG00000030054 is TSPY
# ENSBIXG00000030057 is ZNF280BY

angus.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene_name <- rep(NA,nrow(angus.cnvrs_windows_vst_genes.tvi.melt_sig_genes))

angus.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene_name[angus.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene == "ENSBIXG00000001701"] <- "WC1-a"
angus.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene_name[angus.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene == "ENSBIXG00000002096"] <- "TMPRSS11D"
angus.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene_name[angus.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene == "ENSBIXG00000004482"] <- "T-cell receptor alpha chain"
angus.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene_name[angus.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene == "ENSBIXG00000014359"] <- "WC1-b"
# angus.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene_name[angus.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene == "ENSBIXG00000030034"] <- "HSFY-a"
# angus.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene_name[angus.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene == "ENSBIXG00000030039"] <- "HSFY-b"
# angus.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene_name[angus.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene == "ENSBIXG00000030042"] <- "Ubiquitin -conjugating enzyme E2D-a"
# angus.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene_name[angus.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene == "ENSBIXG00000030047"] <- "Ubiquitin -conjugating enzyme E2D-b"
# angus.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene_name[angus.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene == "ENSBIXG00000030050"] <- "HSFY-c"
# angus.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene_name[angus.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene == "ENSBIXG00000030054"] <- "TSPY"
# angus.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene_name[angus.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene == "ENSBIXG00000030057"] <- "ZNF280BY"

angus.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene_name <- as.factor(angus.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene_name)

order_gene <- c("WC1-a",
                "WC1-b",
                "TMPRSS11D",
                "T-cell receptor alpha chain")

angus.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene_name <- factor(angus.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene_name,
                                                                      levels = order_gene)

tiff(filename = "FigFinal_angus_boxplot_Vst0.3CN1.5.tiff",width = 350,height = 400)
bp <- ggplot(angus.cnvrs_windows_vst_genes.tvi.melt_sig_genes, aes(x=Population, y=CN, group=Population)) + geom_boxplot(aes(fill=Population))
bp <- bp + facet_grid(. ~ gene_name, labeller=label_wrap_gen(width = 10, multi_line = TRUE)) + ylab("Normalized copy number") 
bp <- bp + theme_bw() + theme(axis.text.x = element_blank(),axis.text.y = element_text(color = "black"), axis.ticks.x=element_blank()) + xlab("")
bp
dev.off()

#brahman
brahman.cnvrs_windows_vst_genes.tvi.melt <- brahman.cnvrs_windows_vst_genes.tvi.melt %>% 
  mutate(Population = ifelse(Pop > 1, "Indicine","Taurine"))

brahman.cnvrs_windows_vst_genes.tvi.melt$Population <- as.factor(brahman.cnvrs_windows_vst_genes.tvi.melt$Population)

brahman.cnvrs_windows_vst_genes.tvi.melt_merge <- merge(brahman.cnvrs_windows_vst_genes.tvi.melt,
                                                        brahman.cnvrs_windows_vst_genes.tvi.bed_comb,by="gene")

colnames(brahman.cnvrs_windows_vst_genes.tvi.melt_merge)[2] <- "Vst"

#looking at significant genes in brahman ASM
#filter for Vst more than X among those with > 1.5 CN
#if X is 0.3, we have 19 genes, we can filter out Y
brahman.cnvrs_windows_vst_genes.tvi.melt_sig_genes <- 
  brahman.cnvrs_windows_vst_genes.tvi.melt_merge %>% filter(Vst > 0.3) %>% filter(chr != "X")

#change name for better plot
# ENSBIXG00005002807 is Ubiquitin-conjugating enzyme E2D3
# ENSBIXG00005017008 is lncRNA on chr27
# ENSBIXG00005017017 is enteric beta-defensin precursor
# ENSBIXG00005017240 is beta-defensin-like precursor
# ENSBIXG00005031373 is keratin-associated protein 9-1
# ENSBIXG00005031375 is keratin-associated protein 9-2

brahman.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene_name <- rep(NA,nrow(brahman.cnvrs_windows_vst_genes.tvi.melt_sig_genes))

brahman.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene_name[brahman.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene == "ENSBIXG00005002807"] <- "Ubiquitin- conjugating enzyme E2D3"
brahman.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene_name[brahman.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene == "ENSBIXG00005017008"] <- "lncRNA on chr27"
brahman.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene_name[brahman.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene == "ENSBIXG00005017017"] <- "enteric beta-defensin precursor"
brahman.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene_name[brahman.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene == "ENSBIXG00005017240"] <- "beta-defensin- like precursor"
brahman.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene_name[brahman.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene == "ENSBIXG00005031373"] <- "keratin- associated protein 9-1"
brahman.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene_name[brahman.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene == "ENSBIXG00005031375"] <- "keratin- associated protein 9-2"

brahman.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene_name <- as.factor(brahman.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene_name)

order_gene <- c("Ubiquitin- conjugating enzyme E2D3",
                "enteric beta-defensin precursor",
                "beta-defensin- like precursor",
                "keratin- associated protein 9-1",
                "keratin- associated protein 9-2",
                "lncRNA on chr27")

brahman.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene_name <- factor(brahman.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene_name,
                                                                       levels = order_gene)

tiff(filename = "FigFinal_brahman_boxplot_Vst0.3CN1.5.tiff",width = 500,height = 300)
bp <- ggplot(brahman.cnvrs_windows_vst_genes.tvi.melt_sig_genes, aes(x=Population, y=CN, group=Population)) + geom_boxplot(aes(fill=Population))
bp <- bp + facet_grid(. ~ gene_name, labeller=label_wrap_gen(width = 10, multi_line = TRUE)) + ylab("Normalized copy number") 
bp <- bp + theme_bw() + theme(axis.text.x = element_blank(),axis.text.y = element_text(color = "black"), axis.ticks.x=element_blank()) + xlab("")
bp
dev.off()

#arsucd
arsucd.cnvrs_windows_vst_genes.tvi.melt <- arsucd.cnvrs_windows_vst_genes.tvi.melt %>% 
  mutate(Population = ifelse(Pop > 1, "Indicine","Taurine"))

arsucd.cnvrs_windows_vst_genes.tvi.melt$Population <- as.factor(arsucd.cnvrs_windows_vst_genes.tvi.melt$Population)

arsucd.cnvrs_windows_vst_genes.tvi.melt_merge <- merge(arsucd.cnvrs_windows_vst_genes.tvi.melt,
                                                       arsucd.cnvrs_windows_vst_genes.tvi.bed_comb,by="gene")

colnames(arsucd.cnvrs_windows_vst_genes.tvi.melt_merge)[2] <- "Vst"

#looking at significant genes in arsucd ASM
#filter for Vst more than X among those with > 1.5 CN
#if X is 0.3, we have 19 genes, we can filter out Y
arsucd.cnvrs_windows_vst_genes.tvi.melt_sig_genes <- 
  arsucd.cnvrs_windows_vst_genes.tvi.melt_merge %>% filter(Vst > 0.3) %>% filter(chr != "X")

#change name for better plot  
# ENSBTAG00000001925 is TMPRSS11D
# ENSBTAG00000015575 is tachykinin receptor 1
# ENSBTAG00000023157 is putative protein FAM90A12P
# ENSBTAG00000037937 is interferon-induced very large GTPase 1-a
# ENSBTAG00000050853 is lncRNA on chr21
# ENSBTAG00000052878 is interferon-induced very large GTPase 1-b
# ENSBTAG00000052940 is olfactory receptor 145-like
# ENSBTAG00000053555 is beta-defensin-like precursor

arsucd.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene_name <- rep(NA,nrow(arsucd.cnvrs_windows_vst_genes.tvi.melt_sig_genes))

arsucd.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene_name[arsucd.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene == "ENSBTAG00000001925"] <- "TMPRSS11D"
arsucd.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene_name[arsucd.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene == "ENSBTAG00000015575"] <- "tachykinin receptor 1"
arsucd.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene_name[arsucd.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene == "ENSBTAG00000023157"] <- "putative protein FAM90A12P"
arsucd.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene_name[arsucd.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene == "ENSBTAG00000037937"] <- "interferon -induced very large GTPase 1-a"
arsucd.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene_name[arsucd.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene == "ENSBTAG00000050853"] <- "lncRNA on chr21"
arsucd.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene_name[arsucd.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene == "ENSBTAG00000052878"] <- "interferon -induced very large GTPase 1-b"
arsucd.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene_name[arsucd.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene == "ENSBTAG00000052940"] <- "olfactory receptor 145-like"
arsucd.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene_name[arsucd.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene == "ENSBTAG00000053555"] <- "beta-defensin- like precursor"

arsucd.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene_name <- as.factor(arsucd.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene_name)

order_gene <- c("TMPRSS11D",
                "tachykinin receptor 1",
                "putative protein FAM90A12P",
                "interferon -induced very large GTPase 1-a",
                "interferon -induced very large GTPase 1-b",
                "olfactory receptor 145-like",
                "beta-defensin- like precursor",
                "lncRNA on chr21")

arsucd.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene_name <- factor(arsucd.cnvrs_windows_vst_genes.tvi.melt_sig_genes$gene_name,
                                                                      levels = order_gene)

tiff(filename = "FigFinal_arsucd_boxplot_Vst0.3CN1.5.tiff",width = 650,height = 400)
bp <- ggplot(arsucd.cnvrs_windows_vst_genes.tvi.melt_sig_genes, aes(x=Population, y=CN, group=Population)) + geom_boxplot(aes(fill=Population))
bp <- bp + facet_grid(. ~ gene_name, labeller=label_wrap_gen(width = 10, multi_line = TRUE)) + ylab("Normalized copy number") 
bp <- bp + theme_bw() + theme(axis.text.x = element_blank(),axis.text.y = element_text(color = "black"), axis.ticks.x=element_blank()) + xlab("")
bp
dev.off()

############################# Extra Unused ###########################

# angus.cnvrs_windows_vst_genes.tvi.melt$chr <- c()
# angus.cnvrs_windows_vst_genes.tvi.melt$gene_name <- c()
# 
# angus.cnvrs_windows_vst_genes.tvi.melt$gene_name <- angus.cnvrs_windows_vst_genes.tvi.melt$gene
# 
# vec <- c()
# vec2 <- c()
# for (i in 1:nrow(angus.cnvrs_windows_vst_genes.tvi.melt)){
#   row_i <- which(angus.cnvrs_windows_vst_genes.tvi.melt$gene_name[i] == angus.cnvrs_windows_vst_genes.tvi.bed_comb$gene)
#   vec <- angus.cnvrs_windows_vst_genes.tvi.bed_comb$gene_name[row_i]
#   
#   vec2 <- as.character(angus.cnvrs_windows_vst_genes.tvi.bed_comb$chr[row_i])
#   angus.cnvrs_windows_vst_genes.tvi.melt$chr <- vec2
#   
#   if (vec != "gene_sourceensembl"){
#     angus.cnvrs_windows_vst_genes.tvi.melt$gene_name[i] <- vec
#   }
# }
# 
# angus.cnvrs_windows_vst_genes.tvi.melt$gene_name <- as.factor(angus.cnvrs_windows_vst_genes.tvi.melt$gene_name)

#looking at significant genes in Angus ASM
#filter for Vst more than X among those with > 1.5 CN
#if X is 0.75, we have 6 genes, of which five on the Y chr, 1 is TMPRSS11D on chr 6
# angus.cnvrs_windows_vst_genes.tvi.melt_sig_genes <- 
#   angus.cnvrs_windows_vst_genes.tvi.melt %>% filter(Vst > 0.75)
# 
# tiff(filename = "FigFinal_boxplot_Vst0.75CN1.5.tiff",width = 800,height = 300)
# bp <- ggplot(angus.cnvrs_windows_vst_genes.tvi.melt_sig_genes, aes(x=Population, y=CN, group=Population)) + geom_boxplot(aes(fill=Population))
# bp <- bp + facet_grid(. ~ gene_name) + ylab("Normalized copy number")
# bp <- bp + theme_bw()
# bp
# dev.off()