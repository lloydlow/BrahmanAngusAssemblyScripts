#------------------------------------------------------
# Program name: brahman_angus_SV_Assemblytics_2.R
# Objective: take 12 DF with SV (6 types) specific to
#         Angus and Brahman and overlap with annotation,
#         (and maybe repeatmasker out and cnv derek stuff)
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(readr)
library(dplyr)
library(ggplot2)
library(VennDiagram)
library(tidyr)
library(easyGgplot2)
library(stringr)
library(ape)
library(GenomicRanges)

load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/SV6typesBrahmanAngus.RData")

# path where SV results will be sent
dir1 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/SV/"

#first test the overlap with gene anno and see if it works
#then try to loop over the 6 SV types
#ARS-UCD1.2 r96 EBI anno #chr only
arsucd_anno_gtf <- read_tsv("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/EBI_annotation/arsucd1.2/Bos_taurus.ARS-UCD1.2.96.chr.gtf",
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

colnames(arsucd_anno_gtf) <- c("seqname","source","feature","start","end","score","strand","frame","attributes")

arsucd_anno_gtf_gene <- arsucd_anno_gtf %>% filter(feature == "gene") %>% select(seqname,start,end,strand,attributes)

names(arsucd_anno_gtf_gene)[1] <- "chr"

arsucd_anno_gtf_exon <- arsucd_anno_gtf %>% filter(feature == "exon") %>% select(seqname,start,end,strand,attributes)

names(arsucd_anno_gtf_exon)[1] <- "chr"

#combine all Brahman 6 SV types
ARSUCD_as_REF_Brahman_6SVtypes <- rbind(ARSUCD_VS_ANGUS_Deletion_BrahmanOnly_DF,
                                        ARSUCD_VS_ANGUS_Insertion_BrahmanOnly_DF,
                                        ARSUCD_VS_ANGUS_Repeat_contraction_BrahmanOnly_DF,
                                        ARSUCD_VS_ANGUS_Repeat_expansion_BrahmanOnly_DF,
                                        ARSUCD_VS_ANGUS_Tandem_contraction_BrahmanOnly_DF,
                                        ARSUCD_VS_ANGUS_Tandem_expansion_BrahmanOnly_DF)

#combine all Angus 6 SV types
ARSUCD_as_REF_Angus_6SVtypes <- rbind(ARSUCD_VS_ANGUS_Deletion_AngusOnly_DF,
                                      ARSUCD_VS_ANGUS_Insertion_AngusOnly_DF,
                                      ARSUCD_VS_ANGUS_Repeat_contraction_AngusOnly_DF,
                                      ARSUCD_VS_ANGUS_Repeat_expansion_AngusOnly_DF,
                                      ARSUCD_VS_ANGUS_Tandem_contraction_AngusOnly_DF,
                                      ARSUCD_VS_ANGUS_Tandem_expansion_AngusOnly_DF)

#turn sv intervals and annotation into granges
#sv
sv_interval_Brahman <- makeGRangesFromDataFrame(ARSUCD_as_REF_Brahman_6SVtypes, keep.extra.columns = TRUE, 
                                               seqnames.field="chr", start.field="start", 
                                               end.field="end", strand.field="strand")

sv_interval_Angus <- makeGRangesFromDataFrame(ARSUCD_as_REF_Angus_6SVtypes, keep.extra.columns = TRUE, 
                                                seqnames.field="chr", start.field="start", 
                                                end.field="end", strand.field="strand")

#annotation - gene
gene <- makeGRangesFromDataFrame(arsucd_anno_gtf_gene, keep.extra.columns = TRUE, 
                                 seqnames.field="chr", start.field="start", 
                                 end.field="end", strand.field="strand")

gene %>% as.data.frame() %>% dim()

# subsetByOverlaps(sv_interval_Brahman,gene)
# subsetByOverlaps(gene,sv_interval_Brahman)

#Brahman
ebi_overlap_selection_gene_Brahman_df <- mergeByOverlaps(gene, sv_interval_Brahman)

ebi_overlap_selection_gene_Brahman_as_df <- as.data.frame(ebi_overlap_selection_gene_Brahman_df)

#Angus
ebi_overlap_selection_gene_Angus_df <- mergeByOverlaps(gene, sv_interval_Angus)

ebi_overlap_selection_gene_Angus_as_df <- as.data.frame(ebi_overlap_selection_gene_Angus_df)

################################################################################################
##### Brahman #####
#some tidy up to the merge overlap table
#go thro and remove unwanted columns, arrange by sv_interval_Brahman.type, 
#sv_interval_Brahman.seqnames:sv_interval_Brahman.end is the UCD1.2 equivalent Brahman coor
#sv_interval_Brahman.query_chr:sv_interval_Brahman.query_end is the actual Brahman coor in frozen fasta
ebi_overlap_selection_gene_Brahman_as_df <- ebi_overlap_selection_gene_Brahman_as_df %>% 
  select(gene.seqnames:gene.end,gene.attributes,sv_interval_Brahman.seqnames:sv_interval_Brahman.end,
         sv_interval_Brahman.size:sv_interval_Brahman.type,sv_interval_Brahman.query_chr:sv_interval_Brahman.query_end,
         sv_interval_Brahman.uniqueName)

ebi_overlap_selection_gene_Brahman_as_df <- ebi_overlap_selection_gene_Brahman_as_df %>%
  arrange(sv_interval_Brahman.type,gene.seqnames)

ebi_overlap_selection_gene_Brahman_as_df <- ebi_overlap_selection_gene_Brahman_as_df %>%
  separate(gene.attributes, c("gene_id","gene_version","gene_name","gene_source","gene_biotype"),
           sep = ";",extra = "drop", fill = "right")

ebi_overlap_selection_gene_Brahman_as_df$gene_id <- gsub("gene_id \"","",ebi_overlap_selection_gene_Brahman_as_df$gene_id)
ebi_overlap_selection_gene_Brahman_as_df$gene_id <- gsub("\"","",ebi_overlap_selection_gene_Brahman_as_df$gene_id)

#Deletion
ebi_overlap_selection_gene_Brahman_as_df_Deletion <- ebi_overlap_selection_gene_Brahman_as_df %>% 
  filter(sv_interval_Brahman.type == "Deletion")

#need to unique the ensembl id before panther pathway
ebi_overlap_selection_gene_Brahman_as_df_Deletion_uniq <- ebi_overlap_selection_gene_Brahman_as_df_Deletion %>% 
  distinct(gene_id,.keep_all = TRUE)

#write out
write.csv(ebi_overlap_selection_gene_Brahman_as_df_Deletion, 
          file=paste0(dir1,"ebi_overlap_selection_gene_Brahman_as_df_Deletion.csv"),row.names = FALSE)
write.csv(ebi_overlap_selection_gene_Brahman_as_df_Deletion_uniq, 
          file=paste0(dir1,"ebi_overlap_selection_gene_Brahman_as_df_Deletion_uniq.csv"),row.names = FALSE)

#Insertion
ebi_overlap_selection_gene_Brahman_as_df_Insertion <- ebi_overlap_selection_gene_Brahman_as_df %>% 
  filter(sv_interval_Brahman.type == "Insertion")

#need to unique the ensembl id before panther pathway
ebi_overlap_selection_gene_Brahman_as_df_Insertion_uniq <- ebi_overlap_selection_gene_Brahman_as_df_Insertion %>% 
  distinct(gene_id,.keep_all = TRUE)

#write out
write.csv(ebi_overlap_selection_gene_Brahman_as_df_Insertion, 
          file=paste0(dir1,"ebi_overlap_selection_gene_Brahman_as_df_Insertion.csv"),row.names = FALSE)
write.csv(ebi_overlap_selection_gene_Brahman_as_df_Insertion_uniq, 
          file=paste0(dir1,"ebi_overlap_selection_gene_Brahman_as_df_Insertion_uniq.csv"),row.names = FALSE)

#Repeat_contraction
ebi_overlap_selection_gene_Brahman_as_df_Repeat_contraction <- ebi_overlap_selection_gene_Brahman_as_df %>% 
  filter(sv_interval_Brahman.type == "Repeat_contraction")

#need to unique the ensembl id before panther pathway
ebi_overlap_selection_gene_Brahman_as_df_Repeat_contraction_uniq <- ebi_overlap_selection_gene_Brahman_as_df_Repeat_contraction %>% 
  distinct(gene_id,.keep_all = TRUE)

#write out
write.csv(ebi_overlap_selection_gene_Brahman_as_df_Repeat_contraction, 
          file=paste0(dir1,"ebi_overlap_selection_gene_Brahman_as_df_Repeat_contraction.csv"),row.names = FALSE)
write.csv(ebi_overlap_selection_gene_Brahman_as_df_Repeat_contraction_uniq, 
          file=paste0(dir1,"ebi_overlap_selection_gene_Brahman_as_df_Repeat_contraction_uniq.csv"),row.names = FALSE)

#Repeat_expansion
ebi_overlap_selection_gene_Brahman_as_df_Repeat_expansion <- ebi_overlap_selection_gene_Brahman_as_df %>% 
  filter(sv_interval_Brahman.type == "Repeat_expansion")

#need to unique the ensembl id before panther pathway
ebi_overlap_selection_gene_Brahman_as_df_Repeat_expansion_uniq <- ebi_overlap_selection_gene_Brahman_as_df_Repeat_expansion %>% 
  distinct(gene_id,.keep_all = TRUE)

#write out
write.csv(ebi_overlap_selection_gene_Brahman_as_df_Repeat_expansion, 
          file=paste0(dir1,"ebi_overlap_selection_gene_Brahman_as_df_Repeat_expansion.csv"),row.names = FALSE)
write.csv(ebi_overlap_selection_gene_Brahman_as_df_Repeat_expansion_uniq, 
          file=paste0(dir1,"ebi_overlap_selection_gene_Brahman_as_df_Repeat_expansion_uniq.csv"),row.names = FALSE)

#Tandem_contraction
ebi_overlap_selection_gene_Brahman_as_df_Tandem_contraction <- ebi_overlap_selection_gene_Brahman_as_df %>% 
  filter(sv_interval_Brahman.type == "Tandem_contraction")

#need to unique the ensembl id before panther pathway
ebi_overlap_selection_gene_Brahman_as_df_Tandem_contraction_uniq <- ebi_overlap_selection_gene_Brahman_as_df_Tandem_contraction %>% 
  distinct(gene_id,.keep_all = TRUE)

#Tandem_expansion
ebi_overlap_selection_gene_Brahman_as_df_Tandem_expansion <- ebi_overlap_selection_gene_Brahman_as_df %>% 
  filter(sv_interval_Brahman.type == "Tandem_expansion")

#need to unique the ensembl id before panther pathway
ebi_overlap_selection_gene_Brahman_as_df_Tandem_expansion_uniq <- ebi_overlap_selection_gene_Brahman_as_df_Tandem_expansion %>% 
  distinct(gene_id,.keep_all = TRUE)

#write out
write.csv(ebi_overlap_selection_gene_Brahman_as_df_Tandem_expansion, 
          file=paste0(dir1,"ebi_overlap_selection_gene_Brahman_as_df_Tandem_expansion.csv"),row.names = FALSE)
write.csv(ebi_overlap_selection_gene_Brahman_as_df_Tandem_expansion_uniq, 
          file=paste0(dir1,"ebi_overlap_selection_gene_Brahman_as_df_Tandem_expansion_uniq.csv"),row.names = FALSE)

#write out
write.csv(ebi_overlap_selection_gene_Brahman_as_df_Tandem_contraction, 
          file=paste0(dir1,"ebi_overlap_selection_gene_Brahman_as_df_Tandem_contraction.csv"),row.names = FALSE)
write.csv(ebi_overlap_selection_gene_Brahman_as_df_Tandem_contraction_uniq, 
          file=paste0(dir1,"ebi_overlap_selection_gene_Brahman_as_df_Tandem_contraction_uniq.csv"),row.names = FALSE)

##### Angus #####
ebi_overlap_selection_gene_Angus_as_df <- ebi_overlap_selection_gene_Angus_as_df %>% 
  select(gene.seqnames:gene.end,gene.attributes,sv_interval_Angus.seqnames:sv_interval_Angus.end,
         sv_interval_Angus.size:sv_interval_Angus.type,sv_interval_Angus.query_chr:sv_interval_Angus.query_end,
         sv_interval_Angus.uniqueName)

ebi_overlap_selection_gene_Angus_as_df <- ebi_overlap_selection_gene_Angus_as_df %>% 
  arrange(sv_interval_Angus.type,gene.seqnames)

ebi_overlap_selection_gene_Angus_as_df <- ebi_overlap_selection_gene_Angus_as_df %>%
  separate(gene.attributes, c("gene_id","gene_version","gene_name","gene_source","gene_biotype"),
           sep = ";",extra = "drop", fill = "right")

ebi_overlap_selection_gene_Angus_as_df$gene_id <- gsub("gene_id \"","",ebi_overlap_selection_gene_Angus_as_df$gene_id)
ebi_overlap_selection_gene_Angus_as_df$gene_id <- gsub("\"","",ebi_overlap_selection_gene_Angus_as_df$gene_id)

#Deletion
ebi_overlap_selection_gene_Angus_as_df_Deletion <- ebi_overlap_selection_gene_Angus_as_df %>% 
  filter(sv_interval_Angus.type == "Deletion")

#need to unique the ensembl id before panther pathway
ebi_overlap_selection_gene_Angus_as_df_Deletion_uniq <- ebi_overlap_selection_gene_Angus_as_df_Deletion %>% 
  distinct(gene_id,.keep_all = TRUE)

#write out
write.csv(ebi_overlap_selection_gene_Angus_as_df_Deletion, 
          file=paste0(dir1,"ebi_overlap_selection_gene_Angus_as_df_Deletion.csv"),row.names = FALSE)
write.csv(ebi_overlap_selection_gene_Angus_as_df_Deletion_uniq, 
          file=paste0(dir1,"ebi_overlap_selection_gene_Angus_as_df_Deletion_uniq.csv"),row.names = FALSE)

#Insertion
ebi_overlap_selection_gene_Angus_as_df_Insertion <- ebi_overlap_selection_gene_Angus_as_df %>% 
  filter(sv_interval_Angus.type == "Insertion")

#need to unique the ensembl id before panther pathway
ebi_overlap_selection_gene_Angus_as_df_Insertion_uniq <- ebi_overlap_selection_gene_Angus_as_df_Insertion %>% 
  distinct(gene_id,.keep_all = TRUE)

#write out
write.csv(ebi_overlap_selection_gene_Angus_as_df_Insertion, 
          file=paste0(dir1,"ebi_overlap_selection_gene_Angus_as_df_Insertion.csv"),row.names = FALSE)
write.csv(ebi_overlap_selection_gene_Angus_as_df_Insertion_uniq, 
          file=paste0(dir1,"ebi_overlap_selection_gene_Angus_as_df_Insertion_uniq.csv"),row.names = FALSE)

#Repeat_contraction
ebi_overlap_selection_gene_Angus_as_df_Repeat_contraction <- ebi_overlap_selection_gene_Angus_as_df %>% 
  filter(sv_interval_Angus.type == "Repeat_contraction")

#need to unique the ensembl id before panther pathway
ebi_overlap_selection_gene_Angus_as_df_Repeat_contraction_uniq <- ebi_overlap_selection_gene_Angus_as_df_Repeat_contraction %>% 
  distinct(gene_id,.keep_all = TRUE)

#write out
write.csv(ebi_overlap_selection_gene_Angus_as_df_Repeat_contraction, 
          file=paste0(dir1,"ebi_overlap_selection_gene_Angus_as_df_Repeat_contraction.csv"),row.names = FALSE)
write.csv(ebi_overlap_selection_gene_Angus_as_df_Repeat_contraction_uniq, 
          file=paste0(dir1,"ebi_overlap_selection_gene_Angus_as_df_Repeat_contraction_uniq.csv"),row.names = FALSE)

#Repeat_expansion
ebi_overlap_selection_gene_Angus_as_df_Repeat_expansion <- ebi_overlap_selection_gene_Angus_as_df %>% 
  filter(sv_interval_Angus.type == "Repeat_expansion")

#need to unique the ensembl id before panther pathway
ebi_overlap_selection_gene_Angus_as_df_Repeat_expansion_uniq <- ebi_overlap_selection_gene_Angus_as_df_Repeat_expansion %>% 
  distinct(gene_id,.keep_all = TRUE)

#write out
write.csv(ebi_overlap_selection_gene_Angus_as_df_Repeat_expansion, 
          file=paste0(dir1,"ebi_overlap_selection_gene_Angus_as_df_Repeat_expansion.csv"),row.names = FALSE)
write.csv(ebi_overlap_selection_gene_Angus_as_df_Repeat_expansion_uniq, 
          file=paste0(dir1,"ebi_overlap_selection_gene_Angus_as_df_Repeat_expansion_uniq.csv"),row.names = FALSE)

#Tandem_contraction
ebi_overlap_selection_gene_Angus_as_df_Tandem_contraction <- ebi_overlap_selection_gene_Angus_as_df %>% 
  filter(sv_interval_Angus.type == "Tandem_contraction")

#need to unique the ensembl id before panther pathway
ebi_overlap_selection_gene_Angus_as_df_Tandem_contraction_uniq <- ebi_overlap_selection_gene_Angus_as_df_Tandem_contraction %>% 
  distinct(gene_id,.keep_all = TRUE)

#write out
write.csv(ebi_overlap_selection_gene_Angus_as_df_Tandem_contraction, 
          file=paste0(dir1,"ebi_overlap_selection_gene_Angus_as_df_Tandem_contraction.csv"),row.names = FALSE)
write.csv(ebi_overlap_selection_gene_Angus_as_df_Tandem_contraction_uniq, 
          file=paste0(dir1,"ebi_overlap_selection_gene_Angus_as_df_Tandem_contraction_uniq.csv"),row.names = FALSE)

#Tandem_expansion
ebi_overlap_selection_gene_Angus_as_df_Tandem_expansion <- ebi_overlap_selection_gene_Angus_as_df %>% 
  filter(sv_interval_Angus.type == "Tandem_expansion")

#need to unique the ensembl id before panther pathway
ebi_overlap_selection_gene_Angus_as_df_Tandem_expansion_uniq <- ebi_overlap_selection_gene_Angus_as_df_Tandem_expansion %>% 
  distinct(gene_id,.keep_all = TRUE)

#write out
write.csv(ebi_overlap_selection_gene_Angus_as_df_Tandem_expansion, 
          file=paste0(dir1,"ebi_overlap_selection_gene_Angus_as_df_Tandem_expansion.csv"),row.names = FALSE)
write.csv(ebi_overlap_selection_gene_Angus_as_df_Tandem_expansion_uniq, 
          file=paste0(dir1,"ebi_overlap_selection_gene_Angus_as_df_Tandem_expansion_uniq.csv"),row.names = FALSE)

############################################################################################
#Any SV hitting FASD2 using Brahman specific SV? Ans: No
ebi_overlap_selection_gene_Brahman_as_df %>% filter(gene.seqnames == "15") %>% 
  filter(sv_interval_Brahman.start >= 78791037) %>% filter(sv_interval_Brahman.end <= 80120961)

ebi_overlap_selection_gene_Angus_as_df %>% filter(gene.seqnames == "15") %>% 
  filter(sv_interval_Angus.start >= 78791037) %>% filter(sv_interval_Angus.end <= 80120961)

#Any SV in IGF1R? Ans: Yes
ebi_overlap_selection_gene_Brahman_as_df %>% filter(gene_id == "ENSBTAG00000021527") %>%
  select(sv_interval_Brahman.size:sv_interval_Brahman.query_end)
# sv_interval_Brahman.size sv_interval_Brahman.type sv_interval_Brahman.query_chr sv_interval_Brahman.query_start
# 1                     1758                 Deletion                            21                         7925379
# 2                       90                 Deletion                            21                         7991613
# 3                      108                Insertion                            21                         8013450
# 4                      122         Tandem_expansion                            21                         7973777
# sv_interval_Brahman.query_end
# 1                       7925417
# 2                       7991613
# 3                       8013560
# 4                       7973936

ebi_overlap_selection_gene_Angus_as_df %>% filter(gene_id == "ENSBTAG00000021527") %>%
  select(sv_interval_Angus.size:sv_interval_Angus.query_end)
# sv_interval_Angus.size sv_interval_Angus.type sv_interval_Angus.query_chr sv_interval_Angus.query_start
# 1                    639     Tandem_contraction                          21                       7259618
# 2                    131     Tandem_contraction                          21                       7338285
# sv_interval_Angus.query_end
# 1                     7260732
# 2                     7338587

#how many SVs hit genes per SV type?
#how many unique genes got hit by SV?
# Deletion Insertion Repeat_contraction Repeat_expansion Tandem_contraction Tandem_expansion

##########################################################################################
#Brahman
#Deletion
num_sv_Brahman_Deletion <- ebi_overlap_selection_gene_Brahman_as_df_Deletion %>% distinct(sv_interval_Brahman.uniqueName) %>% nrow()
denum_sv_Brahman_Deletion <- nrow(ARSUCD_VS_ANGUS_Deletion_BrahmanOnly_DF)
uniqgene_Brahman_Deletion <- nrow(ebi_overlap_selection_gene_Brahman_as_df_Deletion_uniq)
totalgenes_Brahman_Deletion <- gene %>% as.data.frame() %>% nrow()

ebi_overlap_selection_gene_exon_Brahman_as_df_Deletion <- ebi_overlap_selection_gene_exon_Brahman_as_df %>% 
  filter(sv_interval_Brahman.type == "Deletion")

print(paste0("For Brahman, ",num_sv_Brahman_Deletion," out of ",denum_sv_Brahman_Deletion," Deletion SVs or ",round(num_sv_Brahman_Deletion/denum_sv_Brahman_Deletion * 100,1),"% of the Brahman specific SV overlapped with ",uniqgene_Brahman_Deletion," genes or ",round(uniqgene_Brahman_Deletion/totalgenes_Brahman_Deletion * 100,1),"% of all genes."))

#Brahman
#Insertion
num_sv_Brahman_Insertion <- ebi_overlap_selection_gene_Brahman_as_df_Insertion %>% distinct(sv_interval_Brahman.uniqueName) %>% nrow()
denum_sv_Brahman_Insertion <- nrow(ARSUCD_VS_ANGUS_Insertion_BrahmanOnly_DF)
uniqgene_Brahman_Insertion <- nrow(ebi_overlap_selection_gene_Brahman_as_df_Insertion_uniq)
totalgenes_Brahman_Insertion <- gene %>% as.data.frame() %>% nrow()

ebi_overlap_selection_gene_exon_Brahman_as_df_Insertion <- ebi_overlap_selection_gene_exon_Brahman_as_df %>% 
  filter(sv_interval_Brahman.type == "Insertion")

print(paste0("For Brahman, ",num_sv_Brahman_Insertion," out of ",denum_sv_Brahman_Insertion," Insertion SVs or ",round(num_sv_Brahman_Insertion/denum_sv_Brahman_Insertion * 100,1),"% of the Brahman specific SV overlapped with ",uniqgene_Brahman_Insertion," genes or ",round(uniqgene_Brahman_Insertion/totalgenes_Brahman_Insertion * 100,1),"% of all genes."))

#Brahman
#Repeat_contraction
num_sv_Brahman_Repeat_contraction <- ebi_overlap_selection_gene_Brahman_as_df_Repeat_contraction %>% distinct(sv_interval_Brahman.uniqueName) %>% nrow()
denum_sv_Brahman_Repeat_contraction <- nrow(ARSUCD_VS_ANGUS_Repeat_contraction_BrahmanOnly_DF)
uniqgene_Brahman_Repeat_contraction <- nrow(ebi_overlap_selection_gene_Brahman_as_df_Repeat_contraction_uniq)
totalgenes_Brahman_Repeat_contraction <- gene %>% as.data.frame() %>% nrow()

ebi_overlap_selection_gene_exon_Brahman_as_df_Repeat_contraction <- ebi_overlap_selection_gene_exon_Brahman_as_df %>% 
  filter(sv_interval_Brahman.type == "Repeat_contraction")

print(paste0("For Brahman, ",num_sv_Brahman_Repeat_contraction," out of ",denum_sv_Brahman_Repeat_contraction," Repeat_contraction SVs or ",round(num_sv_Brahman_Repeat_contraction/denum_sv_Brahman_Repeat_contraction * 100,1),"% of the Brahman specific SV overlapped with ",uniqgene_Brahman_Repeat_contraction," genes or ",round(uniqgene_Brahman_Repeat_contraction/totalgenes_Brahman_Repeat_contraction * 100,1),"% of all genes."))

#Brahman
#Repeat_expansion
num_sv_Brahman_Repeat_expansion <- ebi_overlap_selection_gene_Brahman_as_df_Repeat_expansion %>% distinct(sv_interval_Brahman.uniqueName) %>% nrow()
denum_sv_Brahman_Repeat_expansion <- nrow(ARSUCD_VS_ANGUS_Repeat_expansion_BrahmanOnly_DF)
uniqgene_Brahman_Repeat_expansion <- nrow(ebi_overlap_selection_gene_Brahman_as_df_Repeat_expansion_uniq)
totalgenes_Brahman_Repeat_expansion <- gene %>% as.data.frame() %>% nrow()

ebi_overlap_selection_gene_exon_Brahman_as_df_Repeat_expansion <- ebi_overlap_selection_gene_exon_Brahman_as_df %>% 
  filter(sv_interval_Brahman.type == "Repeat_expansion")

print(paste0("For Brahman, ",num_sv_Brahman_Repeat_expansion," out of ",denum_sv_Brahman_Repeat_expansion," Repeat_expansion SVs or ",round(num_sv_Brahman_Repeat_expansion/denum_sv_Brahman_Repeat_expansion * 100,1),"% of the Brahman specific SV overlapped with ",uniqgene_Brahman_Repeat_expansion," genes or ",round(uniqgene_Brahman_Repeat_expansion/totalgenes_Brahman_Repeat_expansion * 100,1),"% of all genes."))

#Brahman
#Tandem_contraction
num_sv_Brahman_Tandem_contraction <- ebi_overlap_selection_gene_Brahman_as_df_Tandem_contraction %>% distinct(sv_interval_Brahman.uniqueName) %>% nrow()
denum_sv_Brahman_Tandem_contraction <- nrow(ARSUCD_VS_ANGUS_Tandem_contraction_BrahmanOnly_DF)
uniqgene_Brahman_Tandem_contraction <- nrow(ebi_overlap_selection_gene_Brahman_as_df_Tandem_contraction_uniq)
totalgenes_Brahman_Tandem_contraction <- gene %>% as.data.frame() %>% nrow()

ebi_overlap_selection_gene_exon_Brahman_as_df_Tandem_contraction <- ebi_overlap_selection_gene_exon_Brahman_as_df %>% 
  filter(sv_interval_Brahman.type == "Tandem_contraction")

print(paste0("For Brahman, ",num_sv_Brahman_Tandem_contraction," out of ",denum_sv_Brahman_Tandem_contraction," Tandem_contraction SVs or ",round(num_sv_Brahman_Tandem_contraction/denum_sv_Brahman_Tandem_contraction * 100,1),"% of the Brahman specific SV overlapped with ",uniqgene_Brahman_Tandem_contraction," genes or ",round(uniqgene_Brahman_Tandem_contraction/totalgenes_Brahman_Tandem_contraction * 100,1),"% of all genes."))

#Brahman
#Tandem_expansion
num_sv_Brahman_Tandem_expansion <- ebi_overlap_selection_gene_Brahman_as_df_Tandem_expansion %>% distinct(sv_interval_Brahman.uniqueName) %>% nrow()
denum_sv_Brahman_Tandem_expansion <- nrow(ARSUCD_VS_ANGUS_Tandem_expansion_BrahmanOnly_DF)
uniqgene_Brahman_Tandem_expansion <- nrow(ebi_overlap_selection_gene_Brahman_as_df_Tandem_expansion_uniq)
totalgenes_Brahman_Tandem_expansion <- gene %>% as.data.frame() %>% nrow()

ebi_overlap_selection_gene_exon_Brahman_as_df_Tandem_expansion <- ebi_overlap_selection_gene_exon_Brahman_as_df %>% 
  filter(sv_interval_Brahman.type == "Tandem_expansion")

print(paste0("For Brahman, ",num_sv_Brahman_Tandem_expansion," out of ",denum_sv_Brahman_Tandem_expansion," Tandem_expansion SVs or ",round(num_sv_Brahman_Tandem_expansion/denum_sv_Brahman_Tandem_expansion * 100,1),"% of the Brahman specific SV overlapped with ",uniqgene_Brahman_Tandem_expansion," genes or ",round(uniqgene_Brahman_Tandem_expansion/totalgenes_Brahman_Tandem_expansion * 100,1),"% of all genes."))

#Angus
#Deletion
num_sv_Angus_Deletion <- ebi_overlap_selection_gene_Angus_as_df_Deletion %>% distinct(sv_interval_Angus.uniqueName) %>% nrow()
denum_sv_Angus_Deletion <- nrow(ARSUCD_VS_ANGUS_Deletion_AngusOnly_DF)
uniqgene_Angus_Deletion <- nrow(ebi_overlap_selection_gene_Angus_as_df_Deletion_uniq)
totalgenes_Angus_Deletion <- gene %>% as.data.frame() %>% nrow()

ebi_overlap_selection_gene_exon_Angus_as_df_Deletion <- ebi_overlap_selection_gene_exon_Angus_as_df %>% 
  filter(sv_interval_Angus.type == "Deletion")

print(paste0("For Angus, ",num_sv_Angus_Deletion," out of ",denum_sv_Angus_Deletion," Deletion SVs or ",round(num_sv_Angus_Deletion/denum_sv_Angus_Deletion * 100,1),"% of the Angus specific SV overlapped with ",uniqgene_Angus_Deletion," genes or ",round(uniqgene_Angus_Deletion/totalgenes_Angus_Deletion * 100,1),"% of all genes."))

#Angus
#Insertion
num_sv_Angus_Insertion <- ebi_overlap_selection_gene_Angus_as_df_Insertion %>% distinct(sv_interval_Angus.uniqueName) %>% nrow()
denum_sv_Angus_Insertion <- nrow(ARSUCD_VS_ANGUS_Insertion_AngusOnly_DF)
uniqgene_Angus_Insertion <- nrow(ebi_overlap_selection_gene_Angus_as_df_Insertion_uniq)
totalgenes_Angus_Insertion <- gene %>% as.data.frame() %>% nrow()

ebi_overlap_selection_gene_exon_Angus_as_df_Insertion <- ebi_overlap_selection_gene_exon_Angus_as_df %>% 
  filter(sv_interval_Angus.type == "Insertion")

print(paste0("For Angus, ",num_sv_Angus_Insertion," out of ",denum_sv_Angus_Insertion," Insertion SVs or ",round(num_sv_Angus_Insertion/denum_sv_Angus_Insertion * 100,1),"% of the Angus specific SV overlapped with ",uniqgene_Angus_Insertion," genes or ",round(uniqgene_Angus_Insertion/totalgenes_Angus_Insertion * 100,1),"% of all genes."))

#Angus
#Repeat_contraction
num_sv_Angus_Repeat_contraction <- ebi_overlap_selection_gene_Angus_as_df_Repeat_contraction %>% distinct(sv_interval_Angus.uniqueName) %>% nrow()
denum_sv_Angus_Repeat_contraction <- nrow(ARSUCD_VS_ANGUS_Repeat_contraction_AngusOnly_DF)
uniqgene_Angus_Repeat_contraction <- nrow(ebi_overlap_selection_gene_Angus_as_df_Repeat_contraction_uniq)
totalgenes_Angus_Repeat_contraction <- gene %>% as.data.frame() %>% nrow()

ebi_overlap_selection_gene_exon_Angus_as_df_Repeat_contraction <- ebi_overlap_selection_gene_exon_Angus_as_df %>% 
  filter(sv_interval_Angus.type == "Repeat_contraction")

print(paste0("For Angus, ",num_sv_Angus_Repeat_contraction," out of ",denum_sv_Angus_Repeat_contraction," Repeat_contraction SVs or ",round(num_sv_Angus_Repeat_contraction/denum_sv_Angus_Repeat_contraction * 100,1),"% of the Angus specific SV overlapped with ",uniqgene_Angus_Repeat_contraction," genes or ",round(uniqgene_Angus_Repeat_contraction/totalgenes_Angus_Repeat_contraction * 100,1),"% of all genes."))

#Angus
#Repeat_expansion
num_sv_Angus_Repeat_expansion <- ebi_overlap_selection_gene_Angus_as_df_Repeat_expansion %>% distinct(sv_interval_Angus.uniqueName) %>% nrow()
denum_sv_Angus_Repeat_expansion <- nrow(ARSUCD_VS_ANGUS_Repeat_expansion_AngusOnly_DF)
uniqgene_Angus_Repeat_expansion <- nrow(ebi_overlap_selection_gene_Angus_as_df_Repeat_expansion_uniq)
totalgenes_Angus_Repeat_expansion <- gene %>% as.data.frame() %>% nrow()

ebi_overlap_selection_gene_exon_Angus_as_df_Repeat_expansion <- ebi_overlap_selection_gene_exon_Angus_as_df %>% 
  filter(sv_interval_Angus.type == "Repeat_expansion")

print(paste0("For Angus, ",num_sv_Angus_Repeat_expansion," out of ",denum_sv_Angus_Repeat_expansion," Repeat_expansion SVs or ",round(num_sv_Angus_Repeat_expansion/denum_sv_Angus_Repeat_expansion * 100,1),"% of the Angus specific SV overlapped with ",uniqgene_Angus_Repeat_expansion," genes or ",round(uniqgene_Angus_Repeat_expansion/totalgenes_Angus_Repeat_expansion * 100,1),"% of all genes."))

#Angus
#Tandem_contraction
num_sv_Angus_Tandem_contraction <- ebi_overlap_selection_gene_Angus_as_df_Tandem_contraction %>% distinct(sv_interval_Angus.uniqueName) %>% nrow()
denum_sv_Angus_Tandem_contraction <- nrow(ARSUCD_VS_ANGUS_Tandem_contraction_AngusOnly_DF)
uniqgene_Angus_Tandem_contraction <- nrow(ebi_overlap_selection_gene_Angus_as_df_Tandem_contraction_uniq)
totalgenes_Angus_Tandem_contraction <- gene %>% as.data.frame() %>% nrow()

ebi_overlap_selection_gene_exon_Angus_as_df_Tandem_contraction <- ebi_overlap_selection_gene_exon_Angus_as_df %>% 
  filter(sv_interval_Angus.type == "Tandem_contraction")

print(paste0("For Angus, ",num_sv_Angus_Tandem_contraction," out of ",denum_sv_Angus_Tandem_contraction," Tandem_contraction SVs or ",round(num_sv_Angus_Tandem_contraction/denum_sv_Angus_Tandem_contraction * 100,1),"% of the Angus specific SV overlapped with ",uniqgene_Angus_Tandem_contraction," genes or ",round(uniqgene_Angus_Tandem_contraction/totalgenes_Angus_Tandem_contraction * 100,1),"% of all genes."))

#Angus
#Tandem_expansion
num_sv_Angus_Tandem_expansion <- ebi_overlap_selection_gene_Angus_as_df_Tandem_expansion %>% distinct(sv_interval_Angus.uniqueName) %>% nrow()
denum_sv_Angus_Tandem_expansion <- nrow(ARSUCD_VS_ANGUS_Tandem_expansion_AngusOnly_DF)
uniqgene_Angus_Tandem_expansion <- nrow(ebi_overlap_selection_gene_Angus_as_df_Tandem_expansion_uniq)
totalgenes_Angus_Tandem_expansion <- gene %>% as.data.frame() %>% nrow()

ebi_overlap_selection_gene_exon_Angus_as_df_Tandem_expansion <- ebi_overlap_selection_gene_exon_Angus_as_df %>% 
  filter(sv_interval_Angus.type == "Tandem_expansion")

print(paste0("For Angus, ",num_sv_Angus_Tandem_expansion," out of ",denum_sv_Angus_Tandem_expansion," Tandem_expansion SVs or ",round(num_sv_Angus_Tandem_expansion/denum_sv_Angus_Tandem_expansion * 100,1),"% of the Angus specific SV overlapped with ",uniqgene_Angus_Tandem_expansion," genes or ",round(uniqgene_Angus_Tandem_expansion/totalgenes_Angus_Tandem_expansion * 100,1),"% of all genes."))

##### Overlap with QTL will light up so many hits like Xmas tree! Not worth doing

##### Testing repeat expansion and tandem expansion overlap in Brahman
#turn sv intervals and annotation into granges
#sv
Repeat_expansion_Brahman <- makeGRangesFromDataFrame(ebi_overlap_selection_gene_Brahman_as_df_Repeat_expansion, keep.extra.columns = TRUE, 
                                                seqnames.field="sv_interval_Brahman.query_chr", start.field="sv_interval_Brahman.query_start", 
                                                end.field="sv_interval_Brahman.query_end")

Tandem_expansion_Brahman <- makeGRangesFromDataFrame(ebi_overlap_selection_gene_Brahman_as_df_Tandem_expansion, keep.extra.columns = TRUE, 
                                                     seqnames.field="sv_interval_Brahman.query_chr", start.field="sv_interval_Brahman.query_start", 
                                                     end.field="sv_interval_Brahman.query_end")

#Brahman Repeat_expansion overlaps Tandem_expansion
Repeat_expansion_overlaps_Tandem_expansion_Brahman_df <- mergeByOverlaps(Repeat_expansion_Brahman, Tandem_expansion_Brahman)

Repeat_expansion_overlaps_Tandem_expansion_Brahman_as_df <- as.data.frame(Repeat_expansion_overlaps_Tandem_expansion_Brahman_df)
#note that there is no overlap

################################################################################################


#extra
# #Brahman
# #Deletion
# ebi_overlap_selection_gene_Brahman_as_df_Deletion %>% distinct(sv_interval_Brahman.uniqueName) %>% nrow()
# #[1] 1311
# nrow(ARSUCD_VS_ANGUS_Deletion_BrahmanOnly_DF)
# #[1] 7181
# 1311/7181
# #[1] 0.1825651
# 
# nrow(ebi_overlap_selection_gene_Brahman_as_df_Deletion_uniq)
# gene %>% as.data.frame() %>% dim()
# 972/27270
# #[1] 0.03564356
# 
# #Insertion
# ebi_overlap_selection_gene_Brahman_as_df_Insertion %>% distinct(sv_interval_Brahman.uniqueName) %>% nrow()
# #[1] 1440
# nrow(ARSUCD_VS_ANGUS_Insertion_BrahmanOnly_DF)
# #[1] 8529
# 1440/8529
# #[1] 0.1688357
# 
# nrow(ebi_overlap_selection_gene_Brahman_as_df_Insertion_uniq)
# 1017/27270
# #[1] 0.03729373
# 
# #forget about repeat for now
# #Angus
# #Deletion
# ebi_overlap_selection_gene_Angus_as_df_Deletion %>% distinct(sv_interval_Angus.uniqueName) %>% nrow()
# #[1] 414
# nrow(ARSUCD_VS_ANGUS_Deletion_AngusOnly_DF)
# #[1] 2394
# 414/2394
# #[1] 0.1729323
# 
# nrow(ebi_overlap_selection_gene_Angus_as_df_Deletion_uniq)
# gene %>% as.data.frame() %>% dim()
# 361/27270
# #[1] 0.01323799
# 
# #Insertion
# ebi_overlap_selection_gene_Angus_as_df_Insertion %>% distinct(sv_interval_Angus.uniqueName) %>% nrow()
# #[1] 562
# nrow(ARSUCD_VS_ANGUS_Insertion_AngusOnly_DF)
# #[1] 3326
# 
# 562/3326
# #[1] 0.1689717
# 
# nrow(ebi_overlap_selection_gene_Angus_as_df_Insertion_uniq)
# 458/27270
# #[1] 0.01679501
