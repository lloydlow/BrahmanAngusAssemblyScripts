#------------------------------------------------------
# Program name: brahman_angus_vcf_score_alt_allele_allChr_BrahmanRef.R
# Objective: load all chr for analysis from 
#   combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_<chr>.RData
#   and find overlap of selection interval with NCBI/EBI annotation
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(readr)
library(dplyr)
library(tidyr)
#library(reshape)
#library(refGenome)
library(ape)
library(GenomicRanges)
#library(IRanges)
library(stringr)

#ref is brahman here
# path to filtered vcf where only wanted SNPs per chr are studied
dir1 <- 
  "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/"

# load datasets
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_1.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_10.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_11.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_12.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_13.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_14.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_15.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_16.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_17.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_18.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_19.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_20.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_21.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_22.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_23.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_24.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_25.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_26.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_27.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_28.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_29.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_4.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_5.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_6.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_7.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_8.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_9.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_X.RData")

#combine all 30 chr final SNPs file into one

merged <- do.call("rbind", list(combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_1,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_10,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_11,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_12,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_13,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_14,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_15,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_16,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_17,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_18,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_19,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_20,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_21,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_22,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_23,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_24,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_25,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_26,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_27,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_28,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_29,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_3,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_4,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_5,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_6,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_7,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_8,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_9,
                                combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_X))

#use proportion per breed instead of proportion diff
proportion_diff_df <- merged %>% 
  group_by(CHROM,interval) %>% summarise(an_proportion_mean = mean(an_ALT_proportion),
                                   br_proportion_mean= mean(br_ALT_proportion),
                                   proportion_diff_mean = mean(proportion_diff),
                                   snp_counted = n(),mean_diff_POS = mean(diff(POS)))

proportion_diff_df_summ <- proportion_diff_df %>% filter(snp_counted >= 10)

summary(proportion_diff_df_summ$mean_diff_POS)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 567.7  3582.7  4632.6  4953.7  6020.0 11022.6 

#brahman percentile
quantile(proportion_diff_df_summ$br_proportion_mean, probs = seq(0,1,0.05))

#angus percentile
quantile(proportion_diff_df_summ$an_proportion_mean, probs = seq(0,1,0.05))

#histogram
tiff(filename = "hist_Brahman_vs_Angus_proportion_mean.tiff",width = 800)
par(mfrow=c(1,2))
hist(proportion_diff_df_summ$br_proportion_mean,xlab="Proportion mean of alternate allele in Brahman",
     main="",ylim=c(0,5000))
abline(v=0.1,lty="dashed",col="red")

hist(proportion_diff_df_summ$an_proportion_mean,xlab="Proportion mean of alternate allele in Angus",
     main="",ylim=c(0,5000))
abline(v=0.5,lty="dashed",col="red")
par(mfrow=c(1,1))
dev.off()

#proportion_diff_df$interval_no <- 1:nrow(proportion_diff_df)

proportion_diff_df_summ <- proportion_diff_df %>% filter(br_proportion_mean <= 0.1) %>%
  filter(an_proportion_mean > 0.5) %>% filter(snp_counted >= 10)

#split based on an1-an5 and br1-br6
proportion_diff_df_summ$interval <- gsub("\\(","",proportion_diff_df_summ$interval)
proportion_diff_df_summ$interval <- gsub("\\]","",proportion_diff_df_summ$interval)

proportion_diff_df_summ <- proportion_diff_df_summ %>% 
  separate(interval, c("start","end"),sep = ",",remove = FALSE,convert = TRUE) 

proportion_diff_df_summ$strand <- "*"
names(proportion_diff_df_summ)[1] <- "chr"

#write tsv to see what annotation features overlapped with
# options(scipen = 999)
# write.csv(proportion_diff_df_summ,file=paste0(dir1,"proportion_diff_df_summ.csv"),row.names = FALSE)

############################# NCBI ###########################
#Brahman NCBI anno #Initial analysis was done with this
ncbi_anno_gff3 <- read.gff("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/NCBI/Brahman/GCF_003369695.1_UOA_Brahman_1_genomic.gff")

#rename chr name
ncbi_anno_gff3$seqid <- as.character(ncbi_anno_gff3$seqid)

for (i in 40076:(40076+28)){
  char <- paste0("NC_0",as.character(i),".1")
  ncbi_anno_gff3$seqid[ncbi_anno_gff3$seqid == char] <- as.character(i-40075)
}

ncbi_anno_gff3$seqid[ncbi_anno_gff3$seqid == "NC_040105.1"] <- "X"

#remove the unplaced bcos not studying them
ncbi_anno_gff3 <- ncbi_anno_gff3 %>% filter(!str_detect(seqid, "NW"))

ncbi_anno_gff3_gene <- ncbi_anno_gff3 %>% filter(type == "gene") %>% select(seqid,start,end,strand,
                                                                            attributes)
names(ncbi_anno_gff3_gene)[1] <- "chr"

#turn positively selected intervals and annotation into granges
#selection interval
selection_interval <- makeGRangesFromDataFrame(proportion_diff_df_summ, keep.extra.columns = TRUE, 
                                               seqnames.field="chr", start.field="start", 
                                               end.field="end", strand.field="strand")

#annotation - gene
gene <- makeGRangesFromDataFrame(ncbi_anno_gff3_gene, keep.extra.columns = TRUE, 
                                 seqnames.field="chr", start.field="start", 
                                 end.field="end", strand.field="strand")

gene %>% as.data.frame() %>% dim()

subsetByOverlaps(selection_interval,gene)
subsetByOverlaps(gene,selection_interval)

overlap_selection_gene_df <- mergeByOverlaps(gene, selection_interval)

overlap_selection_gene_as_df <- as.data.frame(overlap_selection_gene_df)

#write out dataframe as csv
overlap_selection_gene_as_df <- overlap_selection_gene_as_df %>% 
  select(gene.seqnames:gene.attributes,selection_interval.seqnames:selection_interval.snp_counted)

write.csv(overlap_selection_gene_as_df,file=paste0(dir1,"overlap_selection_gene_as_df.csv"),row.names = FALSE)

############################# EBI ###########################
#Brahman EBI anno #chr only
ebi_anno_gtf <- read_tsv("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/EBI_annotation/ensembl_annotation_Maternal/Bos_indicus_x_bos_taurus_Mom.UOA_Brahman_1.94.chr.gtf",
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

colnames(ebi_anno_gtf) <- c("seqname","source","feature","start","end","score","strand","frame","attributes")

ebi_anno_gtf_gene <- ebi_anno_gtf %>% filter(feature == "gene") %>% select(seqname,start,end,strand,attributes)

names(ebi_anno_gtf_gene)[1] <- "chr"

#turn positively selected intervals and annotation into granges
#selection interval
selection_interval <- makeGRangesFromDataFrame(proportion_diff_df_summ, keep.extra.columns = TRUE, 
                                               seqnames.field="chr", start.field="start", 
                                               end.field="end", strand.field="strand")

#annotation - gene
gene <- makeGRangesFromDataFrame(ebi_anno_gtf_gene, keep.extra.columns = TRUE, 
                                 seqnames.field="chr", start.field="start", 
                                 end.field="end", strand.field="strand")

gene %>% as.data.frame() %>% dim()

subsetByOverlaps(selection_interval,gene)
subsetByOverlaps(gene,selection_interval)

ebi_overlap_selection_gene_df <- mergeByOverlaps(gene, selection_interval)

ebi_overlap_selection_gene_as_df <- as.data.frame(ebi_overlap_selection_gene_df)

#write out dataframe as csv
ebi_overlap_selection_gene_as_df <- ebi_overlap_selection_gene_as_df %>% 
  select(gene.seqnames:gene.attributes,selection_interval.seqnames:selection_interval.snp_counted) %>%
  arrange(gene.seqnames)

write.csv(ebi_overlap_selection_gene_as_df,file=paste0(dir1,"ebi_overlap_selection_gene_as_df.csv"),row.names = FALSE)
