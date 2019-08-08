#------------------------------------------------------
# Program name: brahman_angus_repeatmasker.R
# Objective: check frozen brahman/angus repeat masker output
#           to make plots
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------
library(readr)
library(dplyr)
library(ggplot2)

#target assemblies: frozenbuff (maybe remove later?), umd3.1, arsucd1.2, brahman, angus

####### frozen buffalo #########
## Example processing of *rm.out below
#need to process *rm.out file to replace space with tab
# remove first 3 lines with tail -n
# sed -e 's/  */	/g' < GCF_000471725.1_UMD_CASPUR_WB_2.0_rm.out > tmp
# sed 's/^	*//' < tmp > tmp1
# mv tmp1 GCF_000471725.1_UMD_CASPUR_WB_2.0_rm.out; rm tmp

# path to repeatmasker.out folder
dir1 <- "/Users/lloyd/Documents/lloyd_2017/Research/RiverBuffalo/buffalo_NextGenAssembly/buffalo_NextGenAssembly/repeatmasker_output/"

#read in frozen buff rm.out
path1 <- paste0(dir1,"GCF_003121395.1_ASM312139v1_rm.out")
frozenbuff_rm <- read_tsv(path1, col_names = FALSE)

names(frozenbuff_rm) <- c("score", "divergence", "deletion", "insertion", "query_sequence", "query_begin", 
                          "query_end", "query_left", "strand","repeat", "family", "repeat_begin", 
                          "repeat_end","repeat_left","ID")

#convert query_sequence to proper name e.g. 1, 2, 3, ... X
chr1_24 <- c("NC_037545.1","NC_037546.1","NC_037547.1","NC_037548.1","NC_037549.1","NC_037550.1",
             "NC_037551.1","NC_037552.1","NC_037553.1","NC_037554.1","NC_037555.1","NC_037556.1",
             "NC_037557.1","NC_037558.1","NC_037559.1","NC_037560.1","NC_037561.1","NC_037562.1",
             "NC_037563.1","NC_037564.1","NC_037565.1","NC_037566.1","NC_037567.1","NC_037568.1")
chrX <- "NC_037569.1"
chr1_X <- c(chr1_24,chrX)

frozenbuff_rm$query_sequence[!(frozenbuff_rm$query_sequence %in% chr1_X)] <- "Unplaced"

#loop thro character vector representing chromosomes to sub in "1","2",... for easier chr identification
vect1 <- c()
for (i in 1:length(chr1_24)){
  vect1 <- frozenbuff_rm$query_sequence %in% chr1_24[i]
  frozenbuff_rm$query_sequence[vect1] <- as.character(i)
}

frozenbuff_rm$query_sequence[frozenbuff_rm$query_sequence %in% chrX] <- "X"

#repeat_begin and repeat_left have silly ()
new_repeat_begin <- gsub("\\(","",frozenbuff_rm$repeat_begin)
new_repeat_begin <- gsub("\\)","",new_repeat_begin)
new_repeat_begin <- as.integer(new_repeat_begin)

new_repeat_left <- gsub("\\(","",frozenbuff_rm$repeat_left)
new_repeat_left <- gsub("\\)","",new_repeat_left)
new_repeat_left <- as.integer(new_repeat_left)

frozenbuff_rm$new_repeat_begin <- new_repeat_begin
frozenbuff_rm$new_repeat_end <- frozenbuff_rm$repeat_end
frozenbuff_rm$new_repeat_left <- new_repeat_left

#loop thro and add query_seq_align_len and perc_id
frozenbuff_rm_subset <- frozenbuff_rm %>% 
  mutate(query_align_len = query_end - query_begin + 1, perc_id = 1 - (divergence/100)) %>%
  filter(perc_id > 0.6) %>%
  select("query_sequence","query_begin","query_end","query_left","strand","repeat","family","query_align_len")

frozenbuff_rm_subset$assembly <- rep("UOA_WB_1",nrow(frozenbuff_rm_subset))

#top align length of repeat family with highest presence
frozenbuff_rm_subset_topRepeatLen <- frozenbuff_rm_subset %>% 
  group_by(family) %>% summarise(total_len = sum(query_align_len)) %>% 
  arrange(desc(total_len))
#############################################################################################

#read in umd3.1 rm.out
path1 <- paste0(dir1,"GCA_000003055.5_Bos_taurus_UMD_3.1.1_rm.out")
umd3_1 <- read_tsv(path1, col_names = FALSE)

names(umd3_1) <- c("score", "divergence", "deletion", "insertion", "query_sequence", "query_begin", 
                          "query_end", "query_left", "strand","repeat", "family", "repeat_begin", 
                          "repeat_end","repeat_left","ID")

#convert query_sequence to proper name e.g. 1, 2, 3, ... X
chr1_29 <- c("GK000001.2","GK000002.2","GK000003.2","GK000004.2","GK000005.2","GK000006.2",
             "GK000007.2","GK000008.2","GK000009.2","GK000010.2","GK000011.2","GK000012.2",
             "GK000013.2","GK000014.2","GK000015.2","GK000016.2","GK000017.2","GK000018.2",
             "GK000019.2","GK000020.2","GK000021.2","GK000022.2","GK000023.2","GK000024.2",
             "GK000025.2","GK000026.2","GK000027.2","GK000028.2","GK000029.2")
chrX <- "GK000030.2"
chr1_X <- c(chr1_29,chrX)

umd3_1$query_sequence[!(umd3_1$query_sequence %in% chr1_X)] <- "Unplaced"

#loop thro character vector representing chromosomes to sub in "1","2",... for easier chr identification
vect1 <- c()
for (i in 1:length(chr1_29)){
  vect1 <- umd3_1$query_sequence %in% chr1_29[i]
  umd3_1$query_sequence[vect1] <- as.character(i)
}

umd3_1$query_sequence[umd3_1$query_sequence %in% chrX] <- "X"

#repeat_begin and repeat_left have silly ()
new_repeat_begin <- gsub("\\(","",umd3_1$repeat_begin)
new_repeat_begin <- gsub("\\)","",new_repeat_begin)
new_repeat_begin <- as.integer(new_repeat_begin)

new_repeat_left <- gsub("\\(","",umd3_1$repeat_left)
new_repeat_left <- gsub("\\)","",new_repeat_left)
new_repeat_left <- as.integer(new_repeat_left)

umd3_1$new_repeat_begin <- new_repeat_begin
umd3_1$new_repeat_end <- umd3_1$repeat_end
umd3_1$new_repeat_left <- new_repeat_left

#loop thro and add query_seq_align_len and perc_id
umd3_1_subset <- umd3_1 %>% 
  mutate(query_align_len = query_end - query_begin + 1, perc_id = 1 - (divergence/100)) %>%
  filter(perc_id > 0.6) %>%
  select("query_sequence","query_begin","query_end","query_left","strand","repeat","family","query_align_len")

umd3_1_subset$assembly <- rep("UMD3.1.1",nrow(umd3_1_subset))

#top align length of repeat family with highest presence
umd3_1_subset_topRepeatLen <- umd3_1_subset %>% 
  group_by(family) %>% summarise(total_len = sum(query_align_len)) %>% 
  arrange(desc(total_len))

#############################################################################################

#read in arsucd1.2 rm.out
path1 <- paste0(dir1,"GCF_002263795.1_ARS-UCD1.2_rm.out")
arsucd1_2 <- read_tsv(path1, col_names = FALSE)

#has column 16 containing *, need to be removed
names(arsucd1_2) <- c("score", "divergence", "deletion", "insertion", "query_sequence", "query_begin", 
                   "query_end", "query_left", "strand","repeat", "family", "repeat_begin", 
                   "repeat_end","repeat_left","ID","removed")

arsucd1_2 <- arsucd1_2 %>% select(score:ID)

#convert query_sequence to proper name e.g. 1, 2, 3, ... X
chr1_29 <- c("NC_037328.1","NC_037329.1","NC_037330.1","NC_037331.1","NC_037332.1","NC_037333.1",
             "NC_037334.1","NC_037335.1","NC_037336.1","NC_037337.1","NC_037338.1","NC_037339.1",
             "NC_037340.1","NC_037341.1","NC_037342.1","NC_037343.1","NC_037344.1","NC_037345.1",
             "NC_037346.1","NC_037347.1","NC_037348.1","NC_037349.1","NC_037350.1","NC_037351.1",
             "NC_037352.1","NC_037353.1","NC_037354.1","NC_037355.1","NC_037356.1")
chrX <- "NC_037357.1"
chr1_X <- c(chr1_29,chrX)

arsucd1_2$query_sequence[!(arsucd1_2$query_sequence %in% chr1_X)] <- "Unplaced"

#loop thro character vector representing chromosomes to sub in "1","2",... for easier chr identification
vect1 <- c()
for (i in 1:length(chr1_29)){
  vect1 <- arsucd1_2$query_sequence %in% chr1_29[i]
  arsucd1_2$query_sequence[vect1] <- as.character(i)
}

arsucd1_2$query_sequence[arsucd1_2$query_sequence %in% chrX] <- "X"

#repeat_begin and repeat_left have silly ()
new_repeat_begin <- gsub("\\(","",arsucd1_2$repeat_begin)
new_repeat_begin <- gsub("\\)","",new_repeat_begin)
new_repeat_begin <- as.integer(new_repeat_begin)

new_repeat_left <- gsub("\\(","",arsucd1_2$repeat_left)
new_repeat_left <- gsub("\\)","",new_repeat_left)
new_repeat_left <- as.integer(new_repeat_left)

arsucd1_2$new_repeat_begin <- new_repeat_begin
arsucd1_2$new_repeat_end <- arsucd1_2$repeat_end
arsucd1_2$new_repeat_left <- new_repeat_left

#loop thro and add query_seq_align_len and perc_id
arsucd1_2_subset <- arsucd1_2 %>% 
  mutate(query_align_len = query_end - query_begin + 1, perc_id = 1 - (divergence/100)) %>%
  filter(perc_id > 0.6) %>%
  select("query_sequence","query_begin","query_end","query_left","strand","repeat","family","query_align_len")

arsucd1_2_subset$assembly <- rep("ARS-UCD1.2",nrow(arsucd1_2_subset))

#top align length of repeat family with highest presence
arsucd1_2_subset_topRepeatLen <- arsucd1_2_subset %>% 
  group_by(family) %>% summarise(total_len = sum(query_align_len)) %>% 
  arrange(desc(total_len))

#############################################################################################

#read in angus rm.out
path1 <- paste0(dir1,"angus_repeatmasking.out")
angus <- read_tsv(path1, col_names = FALSE, col_types = cols(X5 = col_character()))

#has column 16 containing *, need to be removed
names(angus) <- c("score", "divergence", "deletion", "insertion", "query_sequence", "query_begin", 
                      "query_end", "query_left", "strand","repeat", "family", "repeat_begin", 
                      "repeat_end","repeat_left","ID","removed")

angus <- angus %>% select(score:ID)

#convert query_sequence to proper name e.g. 1, 2, 3, ... X
chr1_29 <- c("1","2","3","4","5","6",
             "7","8","9","10","11","12",
             "13","14","15","16","17","18",
             "19","20","21","22","23","24",
             "25","26","27","28","29")
chrX <- "Y"
chr1_X <- c(chr1_29,chrX)

angus$query_sequence[!(angus$query_sequence %in% chr1_X)] <- "Unplaced"

#loop thro character vector representing chromosomes to sub in "1","2",... for easier chr identification
vect1 <- c()
for (i in 1:length(chr1_29)){
  vect1 <- angus$query_sequence %in% chr1_29[i]
  angus$query_sequence[vect1] <- as.character(i)
}

angus$query_sequence[angus$query_sequence %in% chrX] <- "Y"

#repeat_begin and repeat_left have silly ()
new_repeat_begin <- gsub("\\(","",angus$repeat_begin)
new_repeat_begin <- gsub("\\)","",new_repeat_begin)
new_repeat_begin <- as.integer(new_repeat_begin)

new_repeat_left <- gsub("\\(","",angus$repeat_left)
new_repeat_left <- gsub("\\)","",new_repeat_left)
new_repeat_left <- as.integer(new_repeat_left)

angus$new_repeat_begin <- new_repeat_begin
angus$new_repeat_end <- angus$repeat_end
angus$new_repeat_left <- new_repeat_left

#loop thro and add query_seq_align_len and perc_id
angus_subset <- angus %>% 
  mutate(query_align_len = query_end - query_begin + 1, perc_id = 1 - (divergence/100)) %>%
  filter(perc_id > 0.6) %>%
  select("query_sequence","query_begin","query_end","query_left","strand","repeat","family","query_align_len")

angus_subset$assembly <- rep("Angus",nrow(angus_subset))

#top align length of repeat family with highest presence
angus_subset_topRepeatLen <- angus_subset %>% 
  group_by(family) %>% summarise(total_len = sum(query_align_len)) %>% 
  arrange(desc(total_len))

#############################################################################################

#read in brahman rm.out
path1 <- paste0(dir1,"brahman_repeatmasking.out")
brahman <- read_tsv(path1, col_names = FALSE, col_types = cols(X5 = col_character()))

#has column 16 containing *, need to be removed
names(brahman) <- c("score", "divergence", "deletion", "insertion", "query_sequence", "query_begin", 
                  "query_end", "query_left", "strand","repeat", "family", "repeat_begin", 
                  "repeat_end","repeat_left","ID","removed")

brahman <- brahman %>% select(score:ID)

#convert query_sequence to proper name e.g. 1, 2, 3, ... X
chr1_29 <- c("1","2","3","4","5","6",
             "7","8","9","10","11","12",
             "13","14","15","16","17","18",
             "19","20","21","22","23","24",
             "25","26","27","28","29")
chrX <- "X"
chr1_X <- c(chr1_29,chrX)

brahman$query_sequence[!(brahman$query_sequence %in% chr1_X)] <- "Unplaced"

#loop thro character vector representing chromosomes to sub in "1","2",... for easier chr identification
vect1 <- c()
for (i in 1:length(chr1_29)){
  vect1 <- brahman$query_sequence %in% chr1_29[i]
  brahman$query_sequence[vect1] <- as.character(i)
}

brahman$query_sequence[brahman$query_sequence %in% chrX] <- "X"

#repeat_begin and repeat_left have silly ()
new_repeat_begin <- gsub("\\(","",brahman$repeat_begin)
new_repeat_begin <- gsub("\\)","",new_repeat_begin)
new_repeat_begin <- as.integer(new_repeat_begin)

new_repeat_left <- gsub("\\(","",brahman$repeat_left)
new_repeat_left <- gsub("\\)","",new_repeat_left)
new_repeat_left <- as.integer(new_repeat_left)

brahman$new_repeat_begin <- new_repeat_begin
brahman$new_repeat_end <- brahman$repeat_end
brahman$new_repeat_left <- new_repeat_left

#loop thro and add query_seq_align_len and perc_id
brahman_subset <- brahman %>% 
  mutate(query_align_len = query_end - query_begin + 1, perc_id = 1 - (divergence/100)) %>%
  filter(perc_id > 0.6) %>%
  select("query_sequence","query_begin","query_end","query_left","strand","repeat","family","query_align_len")

brahman_subset$assembly <- rep("Brahman",nrow(brahman_subset))

#top align length of repeat family with highest presence
brahman_subset_topRepeatLen <- brahman_subset %>% 
  group_by(family) %>% summarise(total_len = sum(query_align_len)) %>% 
  arrange(desc(total_len))

#############################################################################################

#Combine assembly datasets
umd3_1_subset_noUnplaced <- umd3_1_subset %>% filter(query_sequence != "Unplaced")
arsucd1_2_subset_noUnplaced <- arsucd1_2_subset %>% filter(query_sequence != "Unplaced")
angus_subset_noUnplaced <- angus_subset %>% filter(query_sequence != "Unplaced")
brahman_subset_noUnplaced <- brahman_subset %>% filter(query_sequence != "Unplaced")

#all_spp_rm_subset <- rbind(umd3_1_subset,arsucd1_2_subset,angus_subset,brahman_subset)
all_spp_rm_subset <- rbind(umd3_1_subset_noUnplaced,arsucd1_2_subset_noUnplaced,angus_subset_noUnplaced,brahman_subset_noUnplaced)

#order spp
order_spp <- c("Angus","Brahman","ARS-UCD1.2","UMD3.1.1")
all_spp_rm_subset$assembly <- factor(all_spp_rm_subset$assembly, 
                                         levels = order_spp)

all_spp_rm_subset_nofil_summ <- all_spp_rm_subset %>% group_by(family,assembly) %>% summarise(count = n()) %>%
  arrange(desc(count)) 

all_spp_rm_subset_summ <- all_spp_rm_subset %>% group_by(family,assembly) %>% summarise(count = n()) %>%
  arrange(desc(count)) %>% 
  filter(family == "LINE/L1" | family == "LINE/RTE-BovB" | family == "Satellite/centr")

#barplot of different repeat family count (FILTERED for 3 families)
# tiff(filename = "RepeatsFamiliyBrahmanAngus.tiff",width = 600)
# g <- ggplot(data = all_spp_rm_subset_summ, aes(x = family, y = count,fill = assembly)) + 
#   geom_bar(stat = "identity",position="dodge") + ylab("Count") + xlab("Repeat family")
# g <- g + theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.5))
# #g <- g + scale_fill_manual("Assembly", values = c("ARS1" = "maroon", "CASPUR_WB_2" = "turquoise", "UOA_WB_1" = "dodgerblue3"))
# g <- g + theme_bw(base_size = 12)
# g
# dev.off()

#### trying violin plot to show query_align_len distribution 
all_spp_rm_subset_violin <- all_spp_rm_subset %>% group_by(assembly) %>% 
  filter(family == "LINE/L1" | family == "LINE/RTE-BovB" | family == "Satellite/centr") %>%
  filter(query_align_len > 2.5e3)

#The shape of violin starts and ends with the IQR ends
tiff(filename = "RepeatsFamiliyByLengthBrahmanAngus.tiff",width = 900)
bp <- ggplot(all_spp_rm_subset_violin, aes(x=assembly, y=query_align_len, group=assembly))
bp <- bp + geom_violin(aes(fill=assembly)) + guides(fill=FALSE)
bp <- bp + scale_y_log10(labels = scales::comma,breaks=c(0,5000,10000,15000,20000,25000,30000,40000,50000,60000,80000,100000,120000,150000,250000,500000))
bp <- bp + facet_grid(. ~ family) + theme_bw(base_size = 13.5) 
bp <- bp + theme(strip.text.x = element_text(size = 20),element_line(colour = "black"))
bp <- bp + ylab("Length of reference matched to repeat (bp)") + xlab("Assembly")
bp
dev.off()

#for supplementary
#Combine umd3_1_subset,arsucd1_2_subset,angus_subset,brahman_subset to look at unplaced vs in chr
all_spp_rm_subset_supp <- rbind(umd3_1_subset,arsucd1_2_subset,angus_subset,brahman_subset)

all_spp_rm_subset_supp_summ <- all_spp_rm_subset_supp %>% group_by(query_sequence,family,assembly) %>% summarise(count = n()) %>%
  filter(family == "LINE/L1" | family == "LINE/RTE-BovB" | family == "Satellite/centr") %>%
  mutate(placement = ifelse(query_sequence == "Unplaced","Unplaced","Chromosome"))

#barplot of different repeat family count by sequence placement status
tiff(filename = "RepeatsFamiliyBrahmanAngus.tiff",width = 600)
h <- ggplot(data = all_spp_rm_subset_supp_summ, aes(x = family, y = log10(count),fill = assembly)) + 
  geom_bar(stat = "identity",position="dodge") + ylab(expression(log[10]*" count")) + xlab("Repeat family")
h <- h + facet_grid(~placement) 
h <- h + theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.5))
#h <- h + scale_fill_manual("Assembly", values = c("ARS1" = "maroon", "UOA_WB_1" = "dodgerblue3"))
h <- h + theme_bw(base_size = 12)
h
dev.off()


#angus
#Aligned unplaced length in angus to satellite/centromeric
angus_subset_unplaced_cent_alignLength <- angus_subset %>% filter(query_sequence == "Unplaced" & family == "Satellite/centr")

#angus unplaced length is 96948465
sum(angus_subset_unplaced_cent_alignLength$query_align_len)/96948465
#20 Mb of sequence or 21% is alignable to satellite and centromeric repeats

#brahman
#Aligned unplaced length in brahman to satellite/centromeric
brahman_subset_unplaced_cent_alignLength <- brahman_subset %>% filter(query_sequence == "Unplaced" & family == "Satellite/centr")

#brahman unplaced length is 56786952
sum(brahman_subset_unplaced_cent_alignLength$query_align_len)/56786952
#7.9 Mb of sequence or 14% is alignable to satellite and centromeric repeats

#Seems like brahman might just have less sat/centr repeats, why? Bcos it has no Y seq
#Aligned in chr length in angus to satellite/centromeric
angus_subset_chr_cent_alignLength <- angus_subset %>% filter(query_sequence != "Unplaced" & family == "Satellite/centr")
sum(angus_subset_chr_cent_alignLength$query_align_len)/1e6
#so total 20 Mb (unplaced) + 1.6 Mb (chr) = 21.6 Mb

#Aligned in chr length in brahman to satellite/centromeric
brahman_subset_chr_cent_alignLength <- brahman_subset %>% filter(query_sequence != "Unplaced" & family == "Satellite/centr")
sum(brahman_subset_chr_cent_alignLength$query_align_len)/1e6
#so total 7.9 Mb (unplaced) + 1.4 Mb (chr) = 9.3 Mb

#what about dominette, which also has chr X like brahman
arsucd1_2_subset_chr_cent_alignLength <- arsucd1_2_subset %>% filter(family == "Satellite/centr")
sum(arsucd1_2_subset_chr_cent_alignLength$query_align_len)/1e6
#Dominette has 48.75535 Mb, wow!

##### <START> PERCENTAGE OF ANGUS and BRAHMAN REPEATS in the GENOME #####
#loop angus and calc align query length but dont filter percent id to get %repeat
angus_subset_nofilpercid <- angus %>% 
  mutate(query_align_len = query_end - query_begin + 1, perc_id = 1 - (divergence/100)) %>%
  select("query_sequence","query_begin","query_end","query_left","strand","repeat","family","query_align_len")

#angus repeat of all classes straight from repeat masker/genome size
sum(angus_subset_nofilpercid$query_align_len)/2580764822
#0.4913026

#loop brahman and calc align query length but dont filter percent id to get %repeat
brahman_subset_nofilpercid <- brahman %>% 
  mutate(query_align_len = query_end - query_begin + 1, perc_id = 1 - (divergence/100)) %>%
  select("query_sequence","query_begin","query_end","query_left","strand","repeat","family","query_align_len")

#brahman repeat of all classes straight from repeat masker/genome size
sum(brahman_subset_nofilpercid$query_align_len)/2680953056
#0.4927449
##### <END> PERCENTAGE OF ANGUS and BRAHMAN REPEATS in the GENOME #####

##### <START> PERCENTAGE OF TOP TWO repeat families #####
sum(angus_subset_topRepeatLen$total_len[1:2])/2580764822
#0.2466591
sum(brahman_subset_topRepeatLen$total_len[1:2])/2680953056
#0.2529637
##### <END> PERCENTAGE OF TOP TWO repeat families #####

##### <START> Highest repeat family in the unplaced #####
angus_subset_Unplaced <- angus_subset %>% filter(query_sequence == "Unplaced")
angus_subset_Unplaced_topRepeatLen <- angus_subset_Unplaced %>% group_by(family) %>% 
  summarise(total_len = sum(query_align_len)) %>% 
  arrange(desc(total_len))

brahman_subset_Unplaced <- brahman_subset %>% filter(query_sequence == "Unplaced")
brahman_subset_Unplaced_topRepeatLen <- brahman_subset_Unplaced %>% group_by(family) %>% 
  summarise(total_len = sum(query_align_len)) %>% 
  arrange(desc(total_len))
##### <END> Highest repeat family in the unplaced #####

##### <START> centromeric repeats #####
#toggle on-off %>% filter(query_align_len > 10e3) to get the results for paper
brahman_subset_centromere <- brahman_subset %>% filter(family == "Satellite/centr") %>% 
  filter(query_sequence != "Unplaced") %>% filter(query_align_len > 10e3) 
brahman_subset_centromere$query_left <- gsub("\\(","",brahman_subset_centromere$query_left)
brahman_subset_centromere$query_left <- gsub("\\)","",brahman_subset_centromere$query_left)
brahman_subset_centromere$query_left <- as.integer(brahman_subset_centromere$query_left)

brahman_subset_centromere_fil <- brahman_subset_centromere %>% 
  filter(query_begin < 100000 | query_left < 100000) %>% 
  group_by(query_sequence) %>% summarise(count = n()) %>% arrange(as.integer(query_sequence))
nrow(brahman_subset_centromere_fil)

angus_subset_centromere <- angus_subset %>% filter(family == "Satellite/centr") %>% 
  filter(query_sequence != "Unplaced") %>% filter(query_align_len > 10e3)
angus_subset_centromere$query_left <- gsub("\\(","",angus_subset_centromere$query_left)
angus_subset_centromere$query_left <- gsub("\\)","",angus_subset_centromere$query_left)
angus_subset_centromere$query_left <- as.integer(angus_subset_centromere$query_left)

angus_subset_centromere_fil <- angus_subset_centromere %>% 
  filter(query_begin < 100000 | query_left < 100000) %>% 
  group_by(query_sequence) %>% summarise(count = n()) %>% arrange(as.integer(query_sequence))
nrow(angus_subset_centromere_fil)
##### <END> centromeric repeats #####

##### <START> telomeric repeats #####
arsucd1_2_subset_telomere <- arsucd1_2_subset %>% filter(`repeat` == "(TTAGGG)n") %>% 
  filter(query_sequence != "Unplaced") %>% filter(query_align_len >= 60)

frozenbuff_rm_subset_telomere <- frozenbuff_rm_subset %>% filter(`repeat` == "(TTAGGG)n") %>% 
  filter(query_sequence != "Unplaced") %>% filter(query_align_len >= 60)

brahman_subset_telomere <- brahman_subset %>% filter(`repeat` == "(TTAGGG)n") %>% 
  filter(query_sequence != "Unplaced") %>% filter(query_align_len >= 60)
brahman_subset_telomere
brahman_subset %>% filter(`repeat` == "(TTAGGG)n") %>% 
  filter(query_sequence == "Unplaced") %>% filter(query_align_len >= 60)

angus_subset_telomere <- angus_subset %>% filter(`repeat` == "(TTAGGG)n") %>% 
  filter(query_sequence != "Unplaced") %>% filter(query_align_len >= 60)
angus_subset_telomere
angus_subset %>% filter(`repeat` == "(TTAGGG)n") %>% 
  filter(query_sequence == "Unplaced") %>% filter(query_align_len >= 60)
##### <END> telomeric repeats #####

