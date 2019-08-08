#------------------------------------------------------
# Program name: brahman_angus_IsoPhase_allelic_imbalance_v2.R
# Objective: analyse any allelic imbalance specific to selected  
#         gene
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(epade)
library(tidyr)
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(ggbiplot)
library(preprocessCore)

#analysing one gene PB.13819

# path to isophase results
dir1 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/isoseq/IsoPhase/April2019_IsoPhase_with_Lloyd/"

#Is the transcript (at least one of the isoforms) annotated by SQANTI2
# path to annotated vcf is dir1

# reading all.hq.5merge.collapsed.rep_classification.txt
path2 <- paste0(dir1,"all.hq.5merge.collapsed.rep_classification.txt")

classified_iso_to_gene_trans <- read_tsv(path2,col_names = TRUE)

################################################################################################
#####Plot trancript categories
#manipulate table to become category, coding and count
classified_iso_to_gene_trans_plot_struc <- 
  as.data.frame.matrix(table(classified_iso_to_gene_trans$structural_category,classified_iso_to_gene_trans$coding))

category <- rownames(classified_iso_to_gene_trans_plot_struc)

classified_iso_to_gene_trans_plot_struc <- cbind(classified_iso_to_gene_trans_plot_struc,category)
classified_iso_to_gene_trans_plot_struc <- as.tbl(classified_iso_to_gene_trans_plot_struc)

classified_iso_to_gene_trans_plot_struc <- classified_iso_to_gene_trans_plot_struc %>% gather("coding","count",1:2)

classified_iso_to_gene_trans_plot_struc$coding <- gsub("_","-",classified_iso_to_gene_trans_plot_struc$coding)
classified_iso_to_gene_trans_plot_struc$category <- gsub("_"," ",classified_iso_to_gene_trans_plot_struc$category)

#convert to percentage
classified_iso_to_gene_trans_plot_struc <- classified_iso_to_gene_trans_plot_struc %>% mutate(count = count/sum(count)*100)

#order cat
order_cat <- c("full-splice match","incomplete-splice match","novel in catalog","novel not in catalog","antisense","intergenic","genic")
classified_iso_to_gene_trans_plot_struc$category <- factor(classified_iso_to_gene_trans_plot_struc$category,levels = order_cat)

#plot #+ theme(legend.position="top")
tiff(filename = "FigFinal_transcript_structural_category.tiff",width = 300,height = 250)
g <- ggplot(data = classified_iso_to_gene_trans_plot_struc, aes(x = category,y = count,fill = coding)) +
  geom_bar(stat = "identity",position="stack") + ylab("% transcripts") + xlab("Structural category")
g <- g + coord_flip()
g <- g + theme_bw()
g <- g + theme(plot.title = element_text(hjust = 0.5),
               axis.text.x = element_text(angle = 0, hjust = 0.5,color = "black"),
               axis.text.y = element_text(color = "black")) 
g
dev.off()

# classified_iso_to_gene_trans_plot_struc %>% filter(coding == "coding") %>% summarise(total = sum(count))
# total
# 83.29004
# classified_iso_to_gene_trans_plot_struc %>% filter(category != "full-splice match") %>% filter(category != "incomplete-splice match") %>% summarise(total = sum(count))
# total
# 67.7167

#hist of transcript length
tiff(filename = "FigFinal_transcript_length.tiff",width = 400,height = 250)
g <- ggplot(classified_iso_to_gene_trans, aes(x=length)) + geom_histogram(binwidth = 100,color="#F8766D", fill="white")
g <- g + ylab("Count") + xlab("Transcript length (bp)") + scale_x_continuous(breaks = seq(0, 12500, 1000))
g <- g + theme(plot.title = element_text(hjust = 0.5),
               axis.text.x = element_text(angle = 0, hjust = 0.5,color = "black"),
               axis.text.y = element_text(color = "black")) 
g <- g + theme_bw()
g
dev.off()

# getmode <- function(v) {
#      uniqv <- unique(v)
#      uniqv[which.max(tabulate(match(v, uniqv)))]
# }
# getmode(classified_iso_to_gene_trans$length)
# [1] 4125

################################################################################################
#Isoform separated 
path3 <- paste0(dir1,"evaled_isophase.demux_iso_count.txt")

evaled_isophase.demux_iso_count <- read_tsv(path3)

#####Normalization
#by lib size, TPM
evaled_isophase.demux_iso_count_TPM <- evaled_isophase.demux_iso_count %>%
  mutate(heart_p0_new = (heart_p0/(sum(heart_p0)+sum(heart_p1)))*1e6) %>%
  mutate(heart_p1_new = (heart_p1/(sum(heart_p0)+sum(heart_p1)))*1e6) %>%
  mutate(liver_p0_new = (liver_p0/(sum(liver_p0)+sum(liver_p1)))*1e6) %>%
  mutate(liver_p1_new = (liver_p1/(sum(liver_p0)+sum(liver_p1)))*1e6) %>%
  mutate(kidney_p0_new = (kidney_p0/(sum(kidney_p0)+sum(kidney_p1)))*1e6) %>%
  mutate(kidney_p1_new = (kidney_p1/(sum(kidney_p0)+sum(kidney_p1)))*1e6) %>%
  mutate(brain_p0_new = (brain_p0/(sum(brain_p0)+sum(brain_p1)))*1e6) %>%
  mutate(brain_p1_new = (brain_p1/(sum(brain_p0)+sum(brain_p1)))*1e6) %>%
  mutate(lung_p0_new = (lung_p0/(sum(lung_p0)+sum(lung_p1)))*1e6) %>%
  mutate(lung_p1_new = (lung_p1/(sum(lung_p0)+sum(lung_p1)))*1e6) %>%
  mutate(muscle_p0_new = (muscle_p0/(sum(muscle_p0)+sum(muscle_p1)))*1e6) %>%
  mutate(muscle_p1_new = (muscle_p1/(sum(muscle_p0)+sum(muscle_p1)))*1e6) %>%
  mutate(placenta_p0_new = (placenta_p0/(sum(placenta_p0)+sum(placenta_p1)))*1e6) %>%
  mutate(placenta_p1_new = (placenta_p1/(sum(placenta_p0)+sum(placenta_p1)))*1e6)

#remove original count columns, replace with and rename the lib size TPM normalized columns
ori_col_names <- colnames(evaled_isophase.demux_iso_count_TPM)[5:18]
evaled_isophase.demux_iso_count_TPM <- evaled_isophase.demux_iso_count_TPM %>% select(isoform:phase1,heart_p0_new:placenta_p1_new)
colnames(evaled_isophase.demux_iso_count_TPM)[5:18] <- ori_col_names

################################################################################################
#####different transcript utilization
#get total transcript of both alelles per transcript
evaled_isophase.demux_iso_count_TPM_transcript_use <- evaled_isophase.demux_iso_count_TPM %>% 
  mutate(heart_total_transcript = heart_p0+heart_p1) %>%
  mutate(liver_total_transcript = liver_p0+liver_p1) %>%
  mutate(kidney_total_transcript = kidney_p0+kidney_p1) %>%
  mutate(brain_total_transcript = brain_p0+brain_p1) %>%
  mutate(lung_total_transcript = lung_p0+lung_p1) %>%
  mutate(muscle_total_transcript = muscle_p0+muscle_p1) %>%
  mutate(placenta_total_transcript = placenta_p0+placenta_p1)

#get abs difference for each transcript in the two alleles
evaled_isophase.demux_iso_count_TPM_transcript_use <- 
  evaled_isophase.demux_iso_count_TPM_transcript_use %>% 
  mutate(diff_heart = heart_p0 - heart_p1) %>%
  mutate(diff_liver = liver_p0 - liver_p1) %>%
  mutate(diff_kidney = kidney_p0 - kidney_p1) %>%
  mutate(diff_brain = brain_p0 - brain_p1) %>%
  mutate(diff_lung = lung_p0 - lung_p1) %>%
  mutate(diff_muscle = muscle_p0 - muscle_p1) %>%
  mutate(diff_placenta = placenta_p0 - placenta_p1)

#5806 different genes, 52270 transcripts
evaled_isophase.demux_iso_count_TPM_transcript_use_summ <-
  evaled_isophase.demux_iso_count_TPM_transcript_use %>% group_by(locus) %>%
  dplyr::summarise(heart_usage = sum(heart_total_transcript),liver_usage = sum(liver_total_transcript),
                   kidney_usage = sum(kidney_total_transcript),brain_usage = sum(brain_total_transcript),
                   lung_usage = sum(lung_total_transcript),muscle_usage = sum(muscle_total_transcript),
                   placenta_usage = sum(placenta_total_transcript))

#loop thro to get total transcript column per transcript
totaltranscript_heart_vec <- c()
totaltranscript_liver_vec <- c()
totaltranscript_kidney_vec <- c()
totaltranscript_brain_vec <- c()
totaltranscript_lung_vec <- c()
totaltranscript_muscle_vec <- c()
totaltranscript_placenta_vec <- c()

for (i in 1:nrow(evaled_isophase.demux_iso_count_TPM_transcript_use)){
  locus <- evaled_isophase.demux_iso_count_TPM_transcript_use$locus[i]
  
  logic1 <- evaled_isophase.demux_iso_count_TPM_transcript_use_summ$locus == locus
  
  heart_value <- evaled_isophase.demux_iso_count_TPM_transcript_use_summ$heart_usage[logic1]
  totaltranscript_heart_vec <- c(totaltranscript_heart_vec,heart_value)
  
  liver_value <- evaled_isophase.demux_iso_count_TPM_transcript_use_summ$liver_usage[logic1]
  totaltranscript_liver_vec <- c(totaltranscript_liver_vec,liver_value)
  
  kidney_value <- evaled_isophase.demux_iso_count_TPM_transcript_use_summ$kidney_usage[logic1]
  totaltranscript_kidney_vec <- c(totaltranscript_kidney_vec,kidney_value)
  
  brain_value <- evaled_isophase.demux_iso_count_TPM_transcript_use_summ$brain_usage[logic1]
  totaltranscript_brain_vec <- c(totaltranscript_brain_vec,brain_value)
  
  lung_value <- evaled_isophase.demux_iso_count_TPM_transcript_use_summ$lung_usage[logic1]
  totaltranscript_lung_vec <- c(totaltranscript_lung_vec,lung_value)
  
  muscle_value <- evaled_isophase.demux_iso_count_TPM_transcript_use_summ$muscle_usage[logic1]
  totaltranscript_muscle_vec <- c(totaltranscript_muscle_vec,muscle_value)
  
  placenta_value <- evaled_isophase.demux_iso_count_TPM_transcript_use_summ$placenta_usage[logic1]
  totaltranscript_placenta_vec <- c(totaltranscript_placenta_vec,placenta_value)
}

#Add the total transcript per gene to the df
evaled_isophase.demux_iso_count_TPM_transcript_use <- cbind(evaled_isophase.demux_iso_count_TPM_transcript_use,
                                                            totaltranscript_heart_vec,
                                                              totaltranscript_liver_vec,
                                                              totaltranscript_kidney_vec,
                                                              totaltranscript_brain_vec,
                                                              totaltranscript_lung_vec, 
                                                              totaltranscript_muscle_vec,
                                                              totaltranscript_placenta_vec)

#get the proportion diff over the total transcript
evaled_isophase.demux_iso_count_TPM_transcript_use <- evaled_isophase.demux_iso_count_TPM_transcript_use %>%
  mutate(diff_total_ratio_heart = diff_heart/totaltranscript_heart_vec) %>%
  mutate(diff_total_ratio_liver = diff_liver/totaltranscript_liver_vec) %>%
  mutate(diff_total_ratio_kidney = diff_kidney/totaltranscript_kidney_vec) %>%
  mutate(diff_total_ratio_brain = diff_brain/totaltranscript_brain_vec) %>%
  mutate(diff_total_ratio_lung = diff_lung/totaltranscript_lung_vec) %>%
  mutate(diff_total_ratio_muscle = diff_muscle/totaltranscript_muscle_vec) %>%
  mutate(diff_total_ratio_placenta = diff_placenta/totaltranscript_placenta_vec)

#test on the highly expressed brain genes
evaled_isophase.demux_iso_count_TPM_transcript_use_brain_5000 <- 
  evaled_isophase.demux_iso_count_TPM_transcript_use %>% as.tbl() %>% filter(totaltranscript_brain_vec >= 5000) %>%
  select(isoform:phase1,brain_p0,brain_p1,brain_total_transcript,diff_brain,totaltranscript_brain_vec,diff_total_ratio_brain) %>% arrange(locus)

evaled_isophase.demux_iso_count_TPM_transcript_use_brain_5000 %>% filter(diff_total_ratio_brain > 0.2)

################################################################################################
#pregnancy-associated glycoprotein 1-like
#annotation of ENSBIXT00005011437 differs from the best matched NCBI annotated brahman protein
PB.13819 <- classified_iso_to_gene_trans %>% filter(str_detect(isoform, pattern="PB.13819"))
table(PB.13819$associated_transcript)

#confirming iso table that expression only seen in brain and placenta
evaled_isophase.demux_iso_count_TPM %>% filter(locus == "PB.13819_size455") %>% 
  select(heart_p0:placenta_p1) %>% summarise_all(list(sum))

evaled_isophase.demux_iso_count_TPM_PB.13819_size455 <- evaled_isophase.demux_iso_count_TPM %>% filter(locus == "PB.13819_size455")

#reshaping the iso table for plot
iso_reshape_PB.13819_size455 <- evaled_isophase.demux_iso_count_TPM_PB.13819_size455 %>% gather("tissue","count",5:18)
iso_reshape_PB.13819_size455 <- iso_reshape_PB.13819_size455 %>% 
  separate(tissue, c("tissue","breed"),sep = "_",extra = "drop", fill = "right")
iso_reshape_PB.13819_size455$breed <- gsub("p0","brahman",iso_reshape_PB.13819_size455$breed) 
iso_reshape_PB.13819_size455$breed <- gsub("p1","angus",iso_reshape_PB.13819_size455$breed)

#checking gather, separate worked!
iso_reshape_PB.13819_size455_summ <- iso_reshape_PB.13819_size455 %>% group_by(tissue,breed) %>% dplyr::summarise(total = sum(count))
iso_reshape_PB.13819_size455_summ_brain_placenta <- iso_reshape_PB.13819_size455_summ %>% filter(tissue == "brain" | tissue == "placenta")

#barplot just 2 D (count per tissue per breed)
tiff(filename = "FigFinal_PB.13819_size455_iso_count_per_tissue_per_breed_2D.tiff",width = 400,height = 400)
g <- ggplot(data = iso_reshape_PB.13819_size455_summ_brain_placenta, aes(x = tissue, y = total,fill = breed)) +
  geom_bar(stat = "identity",position="dodge") + ylab("Normalized transcript count per million") + xlab("Tissue")
g <- g + theme_bw(base_size = 12)
g <- g + theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.5,color = "black"),axis.text.y = element_text(color = "black"))
g
dev.off()

#filter for just brain and placenta in iso table for pregnancy associated glycoprotein
iso_reshape_PB.13819_size455_brain_placenta <- iso_reshape_PB.13819_size455 %>% filter(tissue == "brain" | tissue == "placenta")

#make unique tissue breed label
iso_reshape_PB.13819_size455_brain_placenta$tissue_breed <- paste0(iso_reshape_PB.13819_size455_brain_placenta$tissue,iso_reshape_PB.13819_size455_brain_placenta$breed)

#remove PB.13819. label and sort by transcript number
iso_reshape_PB.13819_size455_brain_placenta$isoform <- gsub("PB.13819.","",iso_reshape_PB.13819_size455_brain_placenta$isoform)
iso_reshape_PB.13819_size455_brain_placenta$isoform <- as.integer(iso_reshape_PB.13819_size455_brain_placenta$isoform)
iso_reshape_PB.13819_size455_brain_placenta <- iso_reshape_PB.13819_size455_brain_placenta %>% arrange(isoform)

#manipulate to matrix with proper col and row names
iso_reshape_PB.13819_size455_brain_placenta_mx <- spread(iso_reshape_PB.13819_size455_brain_placenta, isoform, count)
iso_reshape_PB.13819_size455_brain_placenta_mx <- as.matrix(t(iso_reshape_PB.13819_size455_brain_placenta_mx[7:25]))
# colnames(iso_reshape_PB.13819_size455_brain_placenta_mx) <- c("brainangus","brainbrahman","placentaangus","placentabrahman")
colnames(iso_reshape_PB.13819_size455_brain_placenta_mx) <- c("brain","brain","placenta","placenta")
rownames(iso_reshape_PB.13819_size455_brain_placenta_mx) <- paste("transcript",rownames(iso_reshape_PB.13819_size455_brain_placenta_mx),sep = " ")

#barplot 3 D (count per tissue per breed)
tiff(filename = "FigFinal_PB.13819_size455_iso_count_per_tissue_per_breed_3D.tiff",width = 650,height = 650)
bar3d.ade(iso_reshape_PB.13819_size455_brain_placenta_mx,wall = 2,alpha=0.5,
          col=c("coral1","#48d1cc","coral1","#48d1cc"),ylab="Normalized transcript count per million")
dev.off()

#check all walls
for(i in 0:6){
  bar3d.ade(iso_reshape_PB.13819_size455_brain_placenta_mx,wall = i,alpha=0.5,
            col=c("coral1","#48d1cc","coral1","#48d1cc"),ylab="Full length transcript count")
}

################################################################################################
#PB.13819.5 getting SNPs from Liz

dir2 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/isoseq/IsoPhase/April2019_IsoPhase_with_Lloyd/PB.13819.5/"

path4 <- paste0(dir2,"phased.partial.cleaned.vcf")

PB.13819.5_phased.partial.cleaned.vcf <- read_tsv(path4,comment = "##")

PB.13819.5_phased.partial.cleaned.vcf %>% select(`#CHROM`:FORMAT,PB.13819.5)

################################################################################################
#####PB.14210_size2844
evaled_isophase.demux_iso_count_TPM %>% filter(locus == "PB.14210_size2844") %>% 
  select(heart_p0:placenta_p1) %>% summarise_all(list(sum))

evaled_isophase.demux_iso_count_TPM_PB.14210_size2844 <- evaled_isophase.demux_iso_count_TPM %>% filter(locus == "PB.14210_size2844")

#reshaping the iso table for plot
iso_reshape_PB.14210_size2844 <- evaled_isophase.demux_iso_count_TPM_PB.14210_size2844 %>% gather("tissue","count",5:18)
iso_reshape_PB.14210_size2844 <- iso_reshape_PB.14210_size2844 %>% 
  separate(tissue, c("tissue","breed"),sep = "_",extra = "drop", fill = "right")
iso_reshape_PB.14210_size2844$breed <- gsub("p0","brahman",iso_reshape_PB.14210_size2844$breed) 
iso_reshape_PB.14210_size2844$breed <- gsub("p1","angus",iso_reshape_PB.14210_size2844$breed)

#checking gather, separate worked!
iso_reshape_PB.14210_size2844_summ <- iso_reshape_PB.14210_size2844 %>% group_by(tissue,breed) %>% dplyr::summarise(total = sum(count))
# iso_reshape_PB.14210_size2844_summ_brain_placenta <- iso_reshape_PB.14210_size2844_summ %>% filter(tissue == "brain" | tissue == "placenta")

#barplot just 2 D (count per tissue per breed)
tiff(filename = "PB.14210_size2844_iso_count_per_tissue_per_breed_2D.tiff",width = 400,height = 400)
g <- ggplot(data = iso_reshape_PB.14210_size2844_summ, aes(x = tissue, y = total,fill = breed)) +
  geom_bar(stat = "identity",position="dodge") + ylab("Normalized transcript count per million") + xlab("Tissue")
g <- g + theme_bw(base_size = 12)
g <- g + theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.5,color = "black"),axis.text.y = element_text(color = "black"))
g
dev.off()

#filter for just brain and placenta in iso table for pregnancy associated glycoprotein
# iso_reshape_PB.14210_size2844_brain_placenta <- iso_reshape_PB.14210_size2844 %>% filter(tissue == "brain" | tissue == "placenta")

#make unique tissue breed label
iso_reshape_PB.14210_size2844$tissue_breed <- paste0(iso_reshape_PB.14210_size2844$tissue,iso_reshape_PB.14210_size2844$breed)

#remove PB.14210. label and sort by transcript number
iso_reshape_PB.14210_size2844$isoform <- gsub("PB.14210.","",iso_reshape_PB.14210_size2844$isoform)
iso_reshape_PB.14210_size2844$isoform <- as.integer(iso_reshape_PB.14210_size2844$isoform)
iso_reshape_PB.14210_size2844 <- iso_reshape_PB.14210_size2844 %>% arrange(isoform)

#manipulate to matrix with proper col and row names
iso_reshape_PB.14210_size2844_mx <- spread(iso_reshape_PB.14210_size2844, isoform, count)

name_for_col <- iso_reshape_PB.14210_size2844_mx$tissue_breed

iso_reshape_PB.14210_size2844_mx <- as.matrix(t(iso_reshape_PB.14210_size2844_mx[7:39]))

colnames(iso_reshape_PB.14210_size2844_mx) <- name_for_col

rownames(iso_reshape_PB.14210_size2844_mx) <- paste("transcript",rownames(iso_reshape_PB.14210_size2844_mx),sep = " ")

#barplot 3 D (count per tissue per breed)
tiff(filename = "iso_count_per_tissue_per_breed_3D.tiff",width = 650,height = 650)
bar3d.ade(iso_reshape_PB.14210_size2844_mx,wall = 2,alpha=0.5,ylab="Normalized transcript count per million")
dev.off()

################################################################################################
#####PB.6234_size1635
evaled_isophase.demux_iso_count_TPM %>% filter(locus == "PB.6234_size1635") %>% 
  select(heart_p0:placenta_p1) %>% summarise_all(list(sum))

evaled_isophase.demux_iso_count_TPM_PB.6234_size1635 <- evaled_isophase.demux_iso_count_TPM %>% filter(locus == "PB.6234_size1635")

#reshaping the iso table for plot
iso_reshape_PB.6234_size1635 <- evaled_isophase.demux_iso_count_TPM_PB.6234_size1635 %>% gather("tissue","count",5:18)
iso_reshape_PB.6234_size1635 <- iso_reshape_PB.6234_size1635 %>% 
  separate(tissue, c("tissue","breed"),sep = "_",extra = "drop", fill = "right")
iso_reshape_PB.6234_size1635$breed <- gsub("p0","brahman",iso_reshape_PB.6234_size1635$breed) 
iso_reshape_PB.6234_size1635$breed <- gsub("p1","angus",iso_reshape_PB.6234_size1635$breed)

#checking gather, separate worked!
iso_reshape_PB.6234_size1635_summ <- iso_reshape_PB.6234_size1635 %>% group_by(tissue,breed) %>% dplyr::summarise(total = sum(count))
# iso_reshape_PB.6234_size1635_summ_brain_placenta <- iso_reshape_PB.6234_size1635_summ %>% filter(tissue == "brain" | tissue == "placenta")

#barplot just 2 D (count per tissue per breed)
tiff(filename = "PB.6234_size1635_iso_count_per_tissue_per_breed_2D.tiff",width = 400,height = 400)
g <- ggplot(data = iso_reshape_PB.6234_size1635_summ, aes(x = tissue, y = total,fill = breed)) +
  geom_bar(stat = "identity",position="dodge") + ylab("Normalized transcript count per million") + xlab("Tissue")
g <- g + theme_bw(base_size = 12)
g <- g + theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.5,color = "black"),axis.text.y = element_text(color = "black"))
g
dev.off()

#filter for just brain and placenta in iso table for pregnancy associated glycoprotein
# iso_reshape_PB.6234_size1635_brain_placenta <- iso_reshape_PB.6234_size1635 %>% filter(tissue == "brain" | tissue == "placenta")

#make unique tissue breed label
iso_reshape_PB.6234_size1635$tissue_breed <- paste0(iso_reshape_PB.6234_size1635$tissue,iso_reshape_PB.6234_size1635$breed)

#remove PB.6234. label and sort by transcript number
iso_reshape_PB.6234_size1635$isoform <- gsub("PB.6234.","",iso_reshape_PB.6234_size1635$isoform)
iso_reshape_PB.6234_size1635$isoform <- as.integer(iso_reshape_PB.6234_size1635$isoform)
iso_reshape_PB.6234_size1635 <- iso_reshape_PB.6234_size1635 %>% arrange(isoform)

#manipulate to matrix with proper col and row names
iso_reshape_PB.6234_size1635_mx <- spread(iso_reshape_PB.6234_size1635, isoform, count)

name_for_col <- iso_reshape_PB.6234_size1635_mx$tissue_breed

iso_reshape_PB.6234_size1635_mx <- as.matrix(t(iso_reshape_PB.6234_size1635_mx[7:39]))

colnames(iso_reshape_PB.6234_size1635_mx) <- name_for_col

rownames(iso_reshape_PB.6234_size1635_mx) <- paste("transcript",rownames(iso_reshape_PB.6234_size1635_mx),sep = " ")

#barplot 3 D (count per tissue per breed)
tiff(filename = "iso_count_per_tissue_per_breed_3D.tiff",width = 650,height = 650)
bar3d.ade(iso_reshape_PB.6234_size1635_mx,wall = 2,alpha=0.5,ylab="Normalized transcript count per million")
dev.off()

################################################################################################
#####PB.1922_size10529
evaled_isophase.demux_iso_count_TPM %>% filter(locus == "PB.1922_size10529") %>% 
  select(heart_p0:placenta_p1) %>% summarise_all(list(sum))

evaled_isophase.demux_iso_count_TPM_PB.1922_size10529 <- evaled_isophase.demux_iso_count_TPM %>% filter(locus == "PB.1922_size10529")

#reshaping the iso table for plot
iso_reshape_PB.1922_size10529 <- evaled_isophase.demux_iso_count_TPM_PB.1922_size10529 %>% gather("tissue","count",5:18)
iso_reshape_PB.1922_size10529 <- iso_reshape_PB.1922_size10529 %>% 
  separate(tissue, c("tissue","breed"),sep = "_",extra = "drop", fill = "right")
iso_reshape_PB.1922_size10529$breed <- gsub("p0","brahman",iso_reshape_PB.1922_size10529$breed) 
iso_reshape_PB.1922_size10529$breed <- gsub("p1","angus",iso_reshape_PB.1922_size10529$breed)

#checking gather, separate worked!
iso_reshape_PB.1922_size10529_summ <- iso_reshape_PB.1922_size10529 %>% group_by(tissue,breed) %>% dplyr::summarise(total = sum(count))
# iso_reshape_PB.1922_size10529_summ_brain_placenta <- iso_reshape_PB.1922_size10529_summ %>% filter(tissue == "brain" | tissue == "placenta")

#barplot just 2 D (count per tissue per breed)
tiff(filename = "FigFinal_PB.1922_size10529_iso_count_per_tissue_per_breed_2D.tiff",width = 400,height = 400)
g <- ggplot(data = iso_reshape_PB.1922_size10529_summ, aes(x = tissue, y = total,fill = breed)) +
  geom_bar(stat = "identity",position="dodge") + ylab("Normalized transcript count per million") + xlab("Tissue")
g <- g + theme_bw(base_size = 12)
g <- g + theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.5,color = "black"),axis.text.y = element_text(color = "black"))
g
dev.off()

#filter for just brain and placenta in iso table for pregnancy associated glycoprotein
# iso_reshape_PB.1922_size10529_brain_placenta <- iso_reshape_PB.1922_size10529 %>% filter(tissue == "brain" | tissue == "placenta")

#make unique tissue breed label
iso_reshape_PB.1922_size10529$tissue_breed <- paste0(iso_reshape_PB.1922_size10529$tissue,iso_reshape_PB.1922_size10529$breed)

#remove PB.1922. label and sort by transcript number
iso_reshape_PB.1922_size10529$isoform <- gsub("PB.1922.","",iso_reshape_PB.1922_size10529$isoform)
iso_reshape_PB.1922_size10529$isoform <- as.integer(iso_reshape_PB.1922_size10529$isoform)
iso_reshape_PB.1922_size10529 <- iso_reshape_PB.1922_size10529 %>% arrange(isoform)

#manipulate to matrix with proper col and row names
iso_reshape_PB.1922_size10529_mx <- spread(iso_reshape_PB.1922_size10529, isoform, count)

name_for_col <- iso_reshape_PB.1922_size10529_mx$tissue_breed

#check here for total columns and adjust accordingly
iso_reshape_PB.1922_size10529_mx <- as.matrix(t(iso_reshape_PB.1922_size10529_mx[7:46]))

colnames(iso_reshape_PB.1922_size10529_mx) <- name_for_col

rownames(iso_reshape_PB.1922_size10529_mx) <- paste("transcript",rownames(iso_reshape_PB.1922_size10529_mx),sep = " ")

#barplot 3 D (count per tissue per breed)
tiff(filename = "PB.1922_size10529_iso_count_per_tissue_per_breed_3D.tiff",width = 650,height = 650)
bar3d.ade(iso_reshape_PB.1922_size10529_mx,wall = 2,alpha=0.5,ylab="Normalized transcript count per million")
dev.off()

################################################################################################
#####PB.10404_size272
evaled_isophase.demux_iso_count_TPM %>% filter(locus == "PB.10404_size272") %>% 
  select(heart_p0:placenta_p1) %>% summarise_all(list(sum))

evaled_isophase.demux_iso_count_TPM_PB.10404_size272 <- evaled_isophase.demux_iso_count_TPM %>% filter(locus == "PB.10404_size272")

#reshaping the iso table for plot
iso_reshape_PB.10404_size272 <- evaled_isophase.demux_iso_count_TPM_PB.10404_size272 %>% gather("tissue","count",5:18)
iso_reshape_PB.10404_size272 <- iso_reshape_PB.10404_size272 %>% 
  separate(tissue, c("tissue","breed"),sep = "_",extra = "drop", fill = "right")
iso_reshape_PB.10404_size272$breed <- gsub("p0","brahman",iso_reshape_PB.10404_size272$breed) 
iso_reshape_PB.10404_size272$breed <- gsub("p1","angus",iso_reshape_PB.10404_size272$breed)

#checking gather, separate worked!
iso_reshape_PB.10404_size272_summ <- iso_reshape_PB.10404_size272 %>% group_by(tissue,breed) %>% dplyr::summarise(total = sum(count))
# iso_reshape_PB.10404_size272_summ_brain_placenta <- iso_reshape_PB.10404_size272_summ %>% filter(tissue == "brain" | tissue == "placenta")

#barplot just 2 D (count per tissue per breed)
tiff(filename = "FigFinal_PB.10404_size272_iso_count_per_tissue_per_breed_2D.tiff",width = 400,height = 500)
g <- ggplot(data = iso_reshape_PB.10404_size272_summ, aes(x = tissue, y = total,fill = breed)) +
  geom_bar(stat = "identity",position="dodge") + ylab("Normalized transcript count per million") + xlab("Tissue")
g <- g + theme_bw(base_size = 12)
g <- g + theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.5,color = "black"),axis.text.y = element_text(color = "black"))
g
dev.off()

#filter for just brain and placenta in iso table for pregnancy associated glycoprotein
# iso_reshape_PB.10404_size272_brain_placenta <- iso_reshape_PB.10404_size272 %>% filter(tissue == "brain" | tissue == "placenta")

#make unique tissue breed label
iso_reshape_PB.10404_size272$tissue_breed <- paste0(iso_reshape_PB.10404_size272$tissue,iso_reshape_PB.10404_size272$breed)

#remove PB.10404. label and sort by transcript number
iso_reshape_PB.10404_size272$isoform <- gsub("PB.10404.","",iso_reshape_PB.10404_size272$isoform)
iso_reshape_PB.10404_size272$isoform <- as.integer(iso_reshape_PB.10404_size272$isoform)
iso_reshape_PB.10404_size272 <- iso_reshape_PB.10404_size272 %>% arrange(isoform)

#manipulate to matrix with proper col and row names
iso_reshape_PB.10404_size272_mx <- spread(iso_reshape_PB.10404_size272, isoform, count)

name_for_col <- iso_reshape_PB.10404_size272_mx$tissue_breed

#check here for total columns and adjust accordingly
iso_reshape_PB.10404_size272_mx <- as.matrix(t(iso_reshape_PB.10404_size272_mx[7:29]))

colnames(iso_reshape_PB.10404_size272_mx) <- name_for_col

rownames(iso_reshape_PB.10404_size272_mx) <- paste("transcript",rownames(iso_reshape_PB.10404_size272_mx),sep = " ")

#barplot 3 D (count per tissue per breed)
tiff(filename = "FigFinal_PB.10404_size272_iso_count_per_tissue_per_breed_3D.tiff",width = 600,height = 600)
bar3d.ade(iso_reshape_PB.10404_size272_mx,wall = 2,alpha=0.5,ylab="Normalized transcript count per million")
dev.off()

#condense to major transcripts
iso_reshape_PB.10404_size272_df <- as.tbl(as.data.frame(iso_reshape_PB.10404_size272_mx))

iso_reshape_PB.10404_size272_df <- iso_reshape_PB.10404_size272_df %>% 
  mutate(heart_combined_transcript = heartangus+heartbrahman) %>%
  mutate(liver_combined_transcript = liverangus+liverbrahman) %>%
  mutate(kidney_combined_transcript = kidneyangus+kidneybrahman) %>%
  mutate(brain_combined_transcript = brainangus+brainbrahman) %>%
  mutate(lung_combined_transcript = lungangus+lungbrahman) %>%
  mutate(muscle_combined_transcript = muscleangus+musclebrahman) %>%
  mutate(placenta_combined_transcript = placentaangus+placentabrahman)

iso_reshape_PB.10404_size272_df <- iso_reshape_PB.10404_size272_df %>% 
  mutate(heart_trans_prop = heart_combined_transcript/sum(heart_combined_transcript)*100) %>%
  mutate(liver_trans_prop = liver_combined_transcript/sum(liver_combined_transcript)*100) %>%
  mutate(kidney_trans_prop = kidney_combined_transcript/sum(kidney_combined_transcript)*100) %>%
  mutate(brain_trans_prop = brain_combined_transcript/sum(brain_combined_transcript)*100) %>%
  mutate(lung_trans_prop = lung_combined_transcript/sum(lung_combined_transcript)*100) %>%
  mutate(muscle_trans_prop = muscle_combined_transcript/sum(muscle_combined_transcript)*100) %>%
  mutate(placenta_trans_prop = placenta_combined_transcript/sum(placenta_combined_transcript)*100)

iso_reshape_PB.10404_size272_df$sum <- iso_reshape_PB.10404_size272_df$heart_trans_prop + 
  iso_reshape_PB.10404_size272_df$liver_trans_prop + iso_reshape_PB.10404_size272_df$kidney_trans_prop +
  iso_reshape_PB.10404_size272_df$brain_trans_prop + iso_reshape_PB.10404_size272_df$lung_trans_prop +
  iso_reshape_PB.10404_size272_df$muscle_trans_prop + iso_reshape_PB.10404_size272_df$placenta_trans_prop

iso_reshape_PB.10404_size272_df <- iso_reshape_PB.10404_size272_df %>% mutate(sum_prop = sum/sum(sum)*100)

iso_reshape_PB.10404_size272_df$trans_no <- 1:nrow(iso_reshape_PB.10404_size272_df)

iso_reshape_PB.10404_size272_df <- iso_reshape_PB.10404_size272_df %>% arrange(desc(sum_prop))

