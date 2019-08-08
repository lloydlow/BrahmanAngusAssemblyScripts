#------------------------------------------------------
# Program name: brahman_angus_IsoPhase_correlate_RNASeq.R
# Objective: check how well IsoPhase results correlate with 
#         Cynthia RNASeq CPM
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------
library(corrplot)
library(dplyr)
library(readr)

load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/isoseq/IsoPhase/April2019_IsoPhase_with_Lloyd/evaled_isophase.demux_hap_count_TPM.RData")

evaled_isophase.demux_hap_count_TPM_corr <- evaled_isophase.demux_hap_count_TPM

evaled_isophase.demux_hap_count_TPM_corr$locus <- gsub("_.*","",evaled_isophase.demux_hap_count_TPM_corr$locus, perl = TRUE)

evaled_isophase.demux_hap_count_TPM_corr <- evaled_isophase.demux_hap_count_TPM_corr %>% 
  arrange(locus) %>% select(locus, liver_total_transcript,brain_total_transcript:placenta_total_transcript)

#read in Cynthia CPM sample97-pacbio-5tissues.csv 
# path to isophase results
dir1 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/isoseq/IsoPhase/April2019_IsoPhase_with_Lloyd/"

# reading all.hq.5merge.collapsed.rep_classification.txt
path1 <- paste0(dir1,"sample97-pacbio-5tissues.csv")

RNASeq_CPM_F97 <- read_csv(path1)
colnames(RNASeq_CPM_F97)[1] <- "locus" 

#loop thro Isophase to get RNASeq CPM
df <- data.frame(brain=numeric(),liver=numeric(),lung=numeric(),muscle=numeric(),placenta=numeric(),
                 stringsAsFactors=FALSE)

for(i in 1:nrow(evaled_isophase.demux_hap_count_TPM_corr)){
  row_i <- which(evaled_isophase.demux_hap_count_TPM_corr$locus[i] == RNASeq_CPM_F97$locus)
  df_add <- RNASeq_CPM_F97[row_i,2:6]
  df <- rbind(df,df_add)
}

#combine IsoPhase with RNASeq CPM
evaled_isophase.demux_hap_count_TPM_corr_comb <- cbind(evaled_isophase.demux_hap_count_TPM_corr,
                                                       df)

#brain pearson correlation
cor(evaled_isophase.demux_hap_count_TPM_corr_comb$brain_total_transcript,
    evaled_isophase.demux_hap_count_TPM_corr_comb$brain)  

#liver pearson correlation
cor(evaled_isophase.demux_hap_count_TPM_corr_comb$liver_total_transcript,
    evaled_isophase.demux_hap_count_TPM_corr_comb$liver)  

#lung pearson correlation
cor(evaled_isophase.demux_hap_count_TPM_corr_comb$lung_total_transcript,
    evaled_isophase.demux_hap_count_TPM_corr_comb$lung)  

#muscle pearson correlation
cor(evaled_isophase.demux_hap_count_TPM_corr_comb$muscle_total_transcript,
    evaled_isophase.demux_hap_count_TPM_corr_comb$muscle)  

#placenta pearson correlation
cor(evaled_isophase.demux_hap_count_TPM_corr_comb$placenta_total_transcript,
    evaled_isophase.demux_hap_count_TPM_corr_comb$placenta)  

#corrplot
evaled_isophase.demux_hap_count_TPM_corr_comb_mx <- evaled_isophase.demux_hap_count_TPM_corr_comb[,2:11]
evaled_isophase.demux_hap_count_TPM_corr_comb_mx <- as.matrix(evaled_isophase.demux_hap_count_TPM_corr_comb_mx)

rownames(evaled_isophase.demux_hap_count_TPM_corr_comb_mx) <- evaled_isophase.demux_hap_count_TPM_corr_comb[,1]

res3 <- cor(evaled_isophase.demux_hap_count_TPM_corr_comb_mx)
round(res3, 2)

corrplot(res3, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45)
