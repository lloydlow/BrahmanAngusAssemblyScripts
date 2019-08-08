#------------------------------------------------------
# Program name: brahman_angus_SV_Assemblytics_4.R
# Objective: overlap brahman and angus SV to get count
#   of exonic, intronic, intergenic, promoter.
#   main diff to v3 is that this one use reduce_ranges=FALSE
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(GenomicFeatures)
library(systemPipeR)
library(tidyr)
library(dplyr)
library(ggplot2)

# arsucd12r96_txdb <- makeTxDbFromGFF(file = "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/EBI_annotation/arsucd1.2/Bos_taurus.ARS-UCD1.2.96.chr.gtf",
#                         format="gtf",dataSource = "ensembl",organism = "Bos taurus")
# 
# saveDb(arsucd12r96_txdb, file="/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/EBI_annotation/arsucd1.2/arsucd12r96.sqlite")

arsucd12r96_txdb <- loadDb("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/EBI_annotation/arsucd1.2/arsucd12r96.sqlite")

#note that arsucd1.2 has 27270 genes
#systemPipeR way of extracting features
#https://mbio.asm.org/content/9/1/e01650-17#DC4 FIGS2 has some examples
# myfeatures <- c("tx_type", "promoter", "intron", "exon", "cds", "intergenic")

feat <- genFeatures(arsucd12r96_txdb, featuretype="all", reduce_ranges=FALSE, upstream=1000, downstream=0)

feat$exon

names(feat)

# ## Obtain feature lists by genes, here for promoter
# split(feat$promoter, unlist(mcols(feat$promoter)$feature_by))
# 
# ## Return all features in single GRanges object
# unlist(feat)

#GenomicFeatures way of extracting features #found out the exons number ranges not consistent. why?
# exonic_parts2 <- exonicParts(arsucd12r96_txdb, linked.to.single.gene.only=TRUE)
# exonic_parts2

#load the SV intervals
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/SV6typesBrahmanAngus.RData")

# path where SV results will be sent
dir1 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/SV/"

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

#Brahman
#exon
ebi_overlap_selection_exon_Brahman_df <- mergeByOverlaps(feat$exon, sv_interval_Brahman)

exon_Brahman <- table(ebi_overlap_selection_exon_Brahman_df$type)

exon_Brahman_df <- data.frame(rbind(exon_Brahman))

#Brahman
#intron
ebi_overlap_selection_intron_Brahman_df <- mergeByOverlaps(feat$intron, sv_interval_Brahman)

intron_Brahman <- table(ebi_overlap_selection_intron_Brahman_df$type)

intron_Brahman_df <- data.frame(rbind(intron_Brahman))

#Brahman
#promoter
ebi_overlap_selection_promoter_Brahman_df <- mergeByOverlaps(feat$promoter, sv_interval_Brahman)

promoter_Brahman <- table(ebi_overlap_selection_promoter_Brahman_df$type)

promoter_Brahman_df <- data.frame(rbind(promoter_Brahman))

#Brahman
#cds
ebi_overlap_selection_cds_Brahman_df <- mergeByOverlaps(feat$cds, sv_interval_Brahman)

cds_Brahman <- table(ebi_overlap_selection_cds_Brahman_df$type)

cds_Brahman_df <- data.frame(rbind(cds_Brahman))

#Brahman
#intergenic
ebi_overlap_selection_intergenic_Brahman_df <- mergeByOverlaps(feat$intergenic, sv_interval_Brahman)

intergenic_Brahman <- table(ebi_overlap_selection_intergenic_Brahman_df$type)

intergenic_Brahman_df <- data.frame(rbind(intergenic_Brahman))

# #Brahman
# #fiveUTR
# ebi_overlap_selection_fiveUTR_Brahman_df <- mergeByOverlaps(feat$fiveUTR, sv_interval_Brahman)
# 
# fiveUTR_Brahman <- table(ebi_overlap_selection_fiveUTR_Brahman_df$type)
# 
# fiveUTR_Brahman_df <- data.frame(rbind(fiveUTR_Brahman))
# 
# #Brahman
# #threeUTR
# ebi_overlap_selection_threeUTR_Brahman_df <- mergeByOverlaps(feat$threeUTR, sv_interval_Brahman)
# 
# threeUTR_Brahman <- table(ebi_overlap_selection_threeUTR_Brahman_df$type)
# 
# threeUTR_Brahman_df <- data.frame(rbind(threeUTR_Brahman))

################################################################################################
#Angus
#exon
ebi_overlap_selection_exon_Angus_df <- mergeByOverlaps(feat$exon, sv_interval_Angus)

exon_Angus <- table(ebi_overlap_selection_exon_Angus_df$type)

exon_Angus_df <- data.frame(rbind(exon_Angus))

#Angus
#intron
ebi_overlap_selection_intron_Angus_df <- mergeByOverlaps(feat$intron, sv_interval_Angus)

intron_Angus <- table(ebi_overlap_selection_intron_Angus_df$type)

intron_Angus_df <- data.frame(rbind(intron_Angus))

#Angus
#promoter
ebi_overlap_selection_promoter_Angus_df <- mergeByOverlaps(feat$promoter, sv_interval_Angus)

promoter_Angus <- table(ebi_overlap_selection_promoter_Angus_df$type)

promoter_Angus_df <- data.frame(rbind(promoter_Angus))

#Angus
#cds
ebi_overlap_selection_cds_Angus_df <- mergeByOverlaps(feat$cds, sv_interval_Angus)

cds_Angus <- table(ebi_overlap_selection_cds_Angus_df$type)

cds_Angus_df <- data.frame(rbind(cds_Angus))

#Angus
#intergenic
ebi_overlap_selection_intergenic_Angus_df <- mergeByOverlaps(feat$intergenic, sv_interval_Angus)

intergenic_Angus <- table(ebi_overlap_selection_intergenic_Angus_df$type)

intergenic_Angus_df <- data.frame(rbind(intergenic_Angus))

# #Angus
# #fiveUTR
# ebi_overlap_selection_fiveUTR_Angus_df <- mergeByOverlaps(feat$fiveUTR, sv_interval_Angus)
# 
# fiveUTR_Angus <- table(ebi_overlap_selection_fiveUTR_Angus_df$type)
# 
# fiveUTR_Angus_df <- data.frame(rbind(fiveUTR_Angus))
# 
# #Angus
# #threeUTR
# ebi_overlap_selection_threeUTR_Angus_df <- mergeByOverlaps(feat$threeUTR, sv_interval_Angus)
# 
# threeUTR_Angus <- table(ebi_overlap_selection_threeUTR_Angus_df$type)
# 
# threeUTR_Angus_df <- data.frame(rbind(threeUTR_Angus))

################################################################################################
#combine all counts to do barchart 
#fiveUTR_red_Brahman_df,threeUTR_red_Brahman_df,fiveUTR_red_Angus_df,threeUTR_red_Angus_df
all_feat_6SVtypes <- rbind(exon_Brahman_df,
                           intron_Brahman_df,
                           promoter_Brahman_df,
                           cds_Brahman_df,
                           intergenic_Brahman_df,
                           exon_Angus_df,
                           intron_Angus_df,
                           promoter_Angus_df,
                           cds_Angus_df,
                           intergenic_Angus_df)

all_feat_6SVtypes$type <- rep(c("exon","intron","promoter","CDS","intergenic"),2)
all_feat_6SVtypes$breed <- c(rep("brahman",5),rep("angus",5))

#reshape wide to long
all_feat_6SVtypes_long <- gather(all_feat_6SVtypes, SV, count, Deletion:Tandem_expansion, factor_key=TRUE)
all_feat_6SVtypes_long$type <- as.factor(all_feat_6SVtypes_long$type)
all_feat_6SVtypes_long

all_feat_6SVtypes_long$type <- factor(all_feat_6SVtypes_long$type, 
                                      levels = c("exon","intron","promoter","CDS","intergenic"))

all_feat_6SVtypes_long$SV <- gsub("_"," ",all_feat_6SVtypes_long$SV)

all_feat_6SVtypes_long %>% group_by(breed,SV) %>% summarise(total = sum(count))

#plot faceted barplot #scales="free_y" in facet_grid # + guides(fill=FALSE)
tiff(filename = "FigFinal_SV_vs_annotation_type_reduce_ranges_false.tiff", width = 1000, height = 400)
g <- ggplot(all_feat_6SVtypes_long, aes(x=type,y=count,color=type)) + geom_bar(stat="identity")
g <- g + facet_grid(breed ~ SV) + ylim(0,6000) + scale_y_continuous(trans='log10')
g <- g + theme_bw()
g <- g + theme(strip.text.x = element_text(size = 10),strip.text.y = element_text(size = 10))
g <- g + theme(axis.text.x = element_text(size = 10, angle = 45, color = "black"),axis.text.y = element_text(color = "black"))
g <- g + scale_x_discrete(name="Annotation type")
g <- g + labs(y=expression("SV count in log"[10]*" scale"))
g
dev.off()

################################################################################################

