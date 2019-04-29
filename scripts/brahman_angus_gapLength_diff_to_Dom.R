#------------------------------------------------------
# Program name: brahman_angus_gapLength_diff_to_Dom.R
# Objective: compare final length, gaps of brahman and
#           angus and compare with Dominette
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(readr)
library(dplyr)
library(reshape2)
library(ggplot2)

# reading ARS_UCD1v25_No_Ns.rls
dir1 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/fasta_file_order_checkBases_N/"
path1 <- paste0(dir1,"ARS_UCD1v25_No_Ns.rls")

ARS_UCD1v25_No_Ns <- read_tsv(path1,col_names = FALSE)
names(ARS_UCD1v25_No_Ns) <- c("scaffold","nothing","gap","length","perc_gap")

ARS_UCD1v25_No_Ns <- ARS_UCD1v25_No_Ns %>% select(scaffold,gap,length)
names(ARS_UCD1v25_No_Ns) <- c("chromosome","dominette_gap","dominette_length")

# reading dam_hd_rc_dom_v2_No_Ns.rls
# path2 <- paste0(dir1,"dam_hd_rc_dom_v2_No_Ns.rls")
# reading bostaurus_brahma_bionano_NCBI_full_corrected_gapfill_arrow_No_Ns.rls
path2 <- paste0(dir1,"bostaurus_brahma_bionano_NCBI_full_corrected_gapfill_arrow_No_Ns.rls")

dam_hd_rc_dom_No_Ns <- read_tsv(path2,col_names = FALSE)
names(dam_hd_rc_dom_No_Ns) <- c("scaffold","nothing","gap","length","perc_gap")

dam_hd_rc_dom_No_Ns <- dam_hd_rc_dom_No_Ns %>% select(scaffold,gap,length)
names(dam_hd_rc_dom_No_Ns) <- c("chromosome","brahman_gap","brahman_length")

# reading sire_hd_rc_dom_v2_No_Ns.rls
# path3 <- paste0(dir1,"sire_hd_rc_dom_v2_No_Ns.rls")
# reading bostaurus_angus_bionano_NCBI_full_corrected_No_Ns.rls
path3 <- paste0(dir1,"bostaurus_angus_bionano_NCBI_full_corrected_No_Ns.rls")

sire_hd_rc_dom_No_Ns <- read_tsv(path3,col_names = FALSE)
names(sire_hd_rc_dom_No_Ns) <- c("scaffold","nothing","gap","length","perc_gap")

sire_hd_rc_dom_No_Ns <- sire_hd_rc_dom_No_Ns %>% select(scaffold,gap,length)
names(sire_hd_rc_dom_No_Ns) <- c("chromosome","angus_gap","angus_length")

#change Y to X first
sire_hd_rc_dom_No_Ns$chromosome <- gsub("Y","X",sire_hd_rc_dom_No_Ns$chromosome)

#merge brahman with dominette
dom_brahman_DF <- merge(ARS_UCD1v25_No_Ns,dam_hd_rc_dom_No_Ns)

#merge angus to both brahman and dominette
dom_brahman_angus_DF <- merge(dom_brahman_DF,sire_hd_rc_dom_No_Ns)

dom_brahman_angus_DF <- dom_brahman_angus_DF %>% arrange(as.integer(chromosome)) 
#%>% mutate(base_diff_Mb = (dominette_length - brahman_length)/1e6) %>% mutate(percent_diff = ((dominette_length - brahman_length)/dominette_length)*100)

#create DF suitable for plotting
df = melt(data.frame(ARSUCD1.2=as.numeric(dom_brahman_angus_DF$dominette_length), 
                     UOA_Brahman_1=as.numeric(dom_brahman_angus_DF$brahman_length),
                     UOA_Angus_1=as.numeric(dom_brahman_angus_DF$angus_length),
                     chromosome=c(1:29,"X")),
          variable.name="genome")

#plot barplot
df$chromosome <- factor(df$chromosome, levels = df$chromosome[1:30])

tiff(filename = "Brahman_angus_vs_dominette_length_comparison.tiff",width = 600, height = 600)
options(scipen=10000)
g <- ggplot(df, aes(chromosome, value/1e6,fill=genome)) 
g <- g + geom_bar(position="dodge",stat="identity") +labs(x="Chromosome",y="Chromosome length (Mbp)")
g <- g + theme_bw()
g
dev.off()

