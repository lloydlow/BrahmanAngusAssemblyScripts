#------------------------------------------------------
# Program name: brahman_angus_bionano_order_scaffold3.R
# Objective: read the final agp scaffold order book
#           to check that contigs use to make chr
#           makes sense
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(readr)
library(dplyr)

load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/order_bionano_scaffolds_to_chr_level/scaff_order_DF.RData")

##### Dam #####
dam_scaffold_order_book <- read_tsv("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_bionano_brahma_corrected.txt",
                                    col_names = TRUE)
#colnames is true for above bcos it was read into Excel multiple times. For agp build, no header
#chr	canu	start	end	orientation

bostaurus_brahma_No_Ns <- read_tsv("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/fasta_file_order_checkBases_N/bostaurus_brahma_No_Ns.rls",col_names = FALSE)
names(bostaurus_brahma_No_Ns) <- c("canu","nothing","gap","length","perc_gap")

bostaurus_brahma_No_Ns <- bostaurus_brahma_No_Ns %>%
  select(canu,gap,length)

bostaurus_brahma_No_Ns$canu <- gsub("\\|","_",bostaurus_brahma_No_Ns$canu)

dam_scaffold_order_book_groupContig <- dam_scaffold_order_book %>% mutate(total = end - start + 1) %>% 
  group_by(canu) %>% summarise(sumtotal = sum(total))

merge_dam_df <- merge(bostaurus_brahma_No_Ns,dam_scaffold_order_book_groupContig)

#which contig is inconsistent?
merge_dam_df %>% filter(length != sumtotal) %>% mutate(difference = length - sumtotal) %>% arrange(canu)

#which scaffold not in merge_*_df, i.e. not already in scaffold order book
merge_dam_df$shortid <- gsub("_.*","",merge_dam_df$canu)

#make the table that has corresponding bionano scaff having same shortid
dam_scaff_order_DF$shortid <- gsub("\\|.*","",dam_scaff_order_DF$component_id)

#logic to get scaffolds not in our existing contigs
logic_dam <- !dam_scaff_order_DF$shortid %in% merge_dam_df$shortid
dam_scaff_order_DF_filtered <- dam_scaff_order_DF[logic_dam,]
unplaced <- unique(dam_scaff_order_DF_filtered$object)
unplaced <- gsub("-","",unplaced)
unplaced <- gsub(":","_",unplaced)
unplaced <- gsub("\\|","_",unplaced)
unplaced <- paste0(unplaced,".fa")
write(unplaced,"unplaced_file")

#this is the one to make agp
write_tsv(dam_scaffold_order_book,
          path = "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_bionano_brahma_corrected.tsv",
          col_names = FALSE)

##### Sire #####
sire_scaffold_order_book <- read_tsv("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_bionano_angus_corrected.txt",
                                    col_names = TRUE)
#colnames is true for above bcos it was read into Excel multiple times. For agp build, no header
#chr	canu	start	end	orientation

bostaurus_angus_No_Ns <- read_tsv("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/fasta_file_order_checkBases_N/bostaurus_angus_No_Ns.rls",col_names = FALSE)
names(bostaurus_angus_No_Ns) <- c("canu","nothing","gap","length","perc_gap")

bostaurus_angus_No_Ns <- bostaurus_angus_No_Ns %>%
  select(canu,gap,length)

bostaurus_angus_No_Ns$canu <- gsub("\\|","_",bostaurus_angus_No_Ns$canu)

sire_scaffold_order_book_groupContig <- sire_scaffold_order_book %>% mutate(total = end - start + 1) %>% 
  group_by(canu) %>% summarise(sumtotal = sum(total))

merge_sire_df <- merge(bostaurus_angus_No_Ns,sire_scaffold_order_book_groupContig)

#which contig is inconsistent?
merge_sire_df %>% filter(length != sumtotal) %>% mutate(difference = length - sumtotal) %>% arrange(canu)

#which scaffold not in merge_*_df, i.e. not already in scaffold order book
merge_sire_df$shortid <- gsub("_.*","",merge_sire_df$canu)

#make the table that has corresponding bionano scaff having same shortid
sire_scaff_order_DF$shortid <- gsub("\\|.*","",sire_scaff_order_DF$component_id)

#logic to get scaffolds not in our existing contigs
logic_sire <- !sire_scaff_order_DF$shortid %in% merge_sire_df$shortid
sire_scaff_order_DF_filtered <- sire_scaff_order_DF[logic_sire,]
unplaced <- unique(sire_scaff_order_DF_filtered$object)
unplaced <- gsub("-","",unplaced)
unplaced <- gsub(":","_",unplaced)
unplaced <- gsub("\\|","_",unplaced)
unplaced <- paste0(unplaced,".fa")
write(unplaced,"unplaced_file")

#this is the one to make agp
write_tsv(sire_scaffold_order_book,
          path = "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_bionano_sire_corrected.tsv",
          col_names = FALSE)
