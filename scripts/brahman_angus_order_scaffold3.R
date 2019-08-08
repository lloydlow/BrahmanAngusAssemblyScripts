#------------------------------------------------------
# Program name: brahman_angus_order_scaffold3.R
# Objective: read csv of scaffold_order_book for dam and sire
#           and evaluate proportion in unplaced to make
#         unplaced file. This is input to take_unplaced.sh
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(readr)
library(dplyr)

load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/order_scaffolds_to_chr_level/scaff_order_DF.RData")

##### Dam #####
dam_scaffold_order_book <- read_tsv("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book.txt",
                                    col_names = FALSE)
names(dam_scaffold_order_book) <- c("object","object_beg","object_end","part_number","component_type",
                                         "component_id","component_beg","component_end","orientation")

#check expected bases placed in chr
dam_scaffold_order_book_bases <- dam_scaffold_order_book %>% filter(component_type == "A") %>% 
  mutate(bases = as.numeric(component_end) - as.numeric(component_beg) + 1)

sum(dam_scaffold_order_book_bases$bases)/sum(dam_scaff_order_DF$clength)

sum(dam_scaffold_order_book_bases$bases)
#old number: 2533462597 ,new number 2603447850

#does my bases in scaffold tally with expected?
#from expected
f1_dam_salsa_No_Ns_rls <- read_tsv("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/fasta_file_order_checkBases_N/f1_dam_salsa_No_Ns.rls",col_names = FALSE)
names(f1_dam_salsa_No_Ns_rls) <- c("component_id","nothing","gap","length","perc_gap")

f1_dam_salsa_No_Ns_rls <- f1_dam_salsa_No_Ns_rls %>%
  select(component_id,gap,length)

#from order book
dam_scaffold_order_book_bases_by_scaff <- dam_scaffold_order_book_bases %>% group_by(component_id) %>% 
  summarise(expected = sum(bases))

#merging expected with observed
dam_scaff_len_obs_exp_merge <- merge(f1_dam_salsa_No_Ns_rls,dam_scaffold_order_book_bases_by_scaff)

#which scaff is inconsistent?
logic <- dam_scaff_len_obs_exp_merge$length == dam_scaff_len_obs_exp_merge$expected
dam_scaff_len_obs_exp_merge %>% filter(length != expected) %>% mutate(difference = length - expected)

logic <- f1_dam_salsa_No_Ns_rls$component_id %in% dam_scaff_len_obs_exp_merge$component_id
unplaced <- paste0(f1_dam_salsa_No_Ns_rls$component_id[!logic],".fa")
write(unplaced,"unplaced_file")

##### Sire #####
sire_scaffold_order_book <- read_tsv("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_sire.txt",
                                     col_names = FALSE)
names(sire_scaffold_order_book) <- c("object","object_beg","object_end","part_number","component_type",
                                     "component_id","component_beg","component_end","orientation")

#check expected bases placed in chr
sire_scaffold_order_book_bases <- sire_scaffold_order_book %>% filter(component_type == "A") %>% 
  mutate(bases = as.numeric(component_end) - as.numeric(component_beg) + 1)

sum(sire_scaffold_order_book_bases$bases)/sum(sire_scaff_order_DF$clength)

sum(sire_scaffold_order_book_bases$bases)
#old with scaffold_320 problem: 2484314111 new: 2482800784

#does my bases in scaffold tally with expected?
#from expected
f1_sire_salsa_No_Ns_rls <- read_tsv("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/fasta_file_order_checkBases_N/f1_sire_salsa_No_Ns.rls",col_names = FALSE)
names(f1_sire_salsa_No_Ns_rls) <- c("component_id","nothing","gap","length","perc_gap")

f1_sire_salsa_No_Ns_rls <- f1_sire_salsa_No_Ns_rls %>%
  select(component_id,gap,length)

#from order book
sire_scaffold_order_book_bases_by_scaff <- sire_scaffold_order_book_bases %>% group_by(component_id) %>% 
  summarise(expected = sum(bases))

#merging expected with observed
sire_scaff_len_obs_exp_merge <- merge(f1_sire_salsa_No_Ns_rls,sire_scaffold_order_book_bases_by_scaff)

#which scaff is inconsistent?
logic <- sire_scaff_len_obs_exp_merge$length == sire_scaff_len_obs_exp_merge$expected
sire_scaff_len_obs_exp_merge %>% filter(length != expected) %>% mutate(difference = length - expected)

logic <- f1_sire_salsa_No_Ns_rls$component_id %in% sire_scaff_len_obs_exp_merge$component_id
unplaced <- paste0(f1_sire_salsa_No_Ns_rls$component_id[!logic],".fa")
write(unplaced,"unplaced_file_sire")
