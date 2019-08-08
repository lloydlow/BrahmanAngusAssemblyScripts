#------------------------------------------------------
# Program name: brahman_angus_bionano_order_scaffold2.R
# Objective: after making big table now reodering scaffold
#           to best fit HD, rc, align Dom 
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(readr)
library(dplyr)

load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/order_bionano_scaffolds_to_chr_level/scaff_order_DF.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/order_bionano_scaffolds_to_chr_level/f1_dam_bionano_HD_rc_Probes_tab.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/order_bionano_scaffolds_to_chr_level/f1_sire_bionano_HD_rc_Probes_tab.RData")

##### Dam #####
dam_scaff_order_DF_1 <- dam_scaff_order_DF %>% filter(Ref_name == 1) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
dam_scaff_order_DF_1 <- dam_scaff_order_DF %>% filter(Ref_name == 1) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 1) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#write out
dam_scaff_order_DF_1 <- dam_scaff_order_DF_1 %>% select(component_id,clength,orientation)
write_tsv(dam_scaff_order_DF_1,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/dam_scaff_order_DF_1.tsv")

#Super-Scaffold_1563
dam_scaff_order_DF_1_plot1 <- f1_dam_bionano_HDProbes_tab %>% filter(scaffold == "Super-Scaffold_1563")
plot(dam_scaff_order_DF_1_plot1$chr_pos,dam_scaff_order_DF_1_plot1$align_pos)
dam_scaff_order_DF_1_plot1 <- f1_dam_bionano_fasta_rcmap_tab %>% filter(scaffold == "Super-Scaffold_1563")
plot(dam_scaff_order_DF_1_plot1$chr_pos,dam_scaff_order_DF_1_plot1$align_pos)

dam_scaff_order_DF_2 <- dam_scaff_order_DF %>% filter(Ref_name == 2) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
dam_scaff_order_DF_2 <- dam_scaff_order_DF %>% filter(Ref_name == 2) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 2) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
dam_scaff_order_DF_2 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
dam_scaff_order_DF_2 <- dam_scaff_order_DF_2 %>% select(component_id,clength,orientation)
write_tsv(dam_scaff_order_DF_2,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/dam_scaff_order_DF_2.tsv")

dam_scaff_order_DF_3 <- dam_scaff_order_DF %>% filter(Ref_name == 3) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
dam_scaff_order_DF_3 <- dam_scaff_order_DF %>% filter(Ref_name == 3) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 3) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
dam_scaff_order_DF_3 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
dam_scaff_order_DF_3 <- dam_scaff_order_DF_3 %>% select(component_id,clength,orientation)
write_tsv(dam_scaff_order_DF_3,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/dam_scaff_order_DF_3.tsv")

dam_scaff_order_DF_4 <- dam_scaff_order_DF %>% filter(Ref_name == 4) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
dam_scaff_order_DF_4 <- dam_scaff_order_DF %>% filter(Ref_name == 4) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 4) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
dam_scaff_order_DF_4 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
dam_scaff_order_DF_4 <- dam_scaff_order_DF_4 %>% select(component_id,clength,orientation)
write_tsv(dam_scaff_order_DF_4,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/dam_scaff_order_DF_4.tsv")

dam_scaff_order_DF_5 <- dam_scaff_order_DF %>% filter(Ref_name == 5) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
dam_scaff_order_DF_5 <- dam_scaff_order_DF %>% filter(Ref_name == 5) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 5) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
dam_scaff_order_DF_5 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
dam_scaff_order_DF_5 <- dam_scaff_order_DF_5 %>% select(component_id,clength,orientation)
write_tsv(dam_scaff_order_DF_5,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/dam_scaff_order_DF_5.tsv")

dam_scaff_order_DF_6 <- dam_scaff_order_DF %>% filter(Ref_name == 6) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
dam_scaff_order_DF_6 <- dam_scaff_order_DF %>% filter(Ref_name == 6) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 6) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
dam_scaff_order_DF_6 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
dam_scaff_order_DF_6 <- dam_scaff_order_DF_6 %>% select(component_id,clength,orientation)
write_tsv(dam_scaff_order_DF_6,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/dam_scaff_order_DF_6.tsv")

dam_scaff_order_DF_7 <- dam_scaff_order_DF %>% filter(Ref_name == 7) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
dam_scaff_order_DF_7 <- dam_scaff_order_DF %>% filter(Ref_name == 7) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 7) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
dam_scaff_order_DF_7 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
dam_scaff_order_DF_7 <- dam_scaff_order_DF_7 %>% select(component_id,clength,orientation)
write_tsv(dam_scaff_order_DF_7,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/dam_scaff_order_DF_7.tsv")

dam_scaff_order_DF_8 <- dam_scaff_order_DF %>% filter(Ref_name == 8) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
dam_scaff_order_DF_8 <- dam_scaff_order_DF %>% filter(Ref_name == 8) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 8) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
dam_scaff_order_DF_8 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
dam_scaff_order_DF_8 <- dam_scaff_order_DF_8 %>% select(component_id,clength,orientation)
write_tsv(dam_scaff_order_DF_8,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/dam_scaff_order_DF_8.tsv")

dam_scaff_order_DF_9 <- dam_scaff_order_DF %>% filter(Ref_name == 9) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
dam_scaff_order_DF_9 <- dam_scaff_order_DF %>% filter(Ref_name == 9) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 9) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
dam_scaff_order_DF_9 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
dam_scaff_order_DF_9 <- dam_scaff_order_DF_9 %>% select(component_id,clength,orientation)
write_tsv(dam_scaff_order_DF_9,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/dam_scaff_order_DF_9.tsv")

dam_scaff_order_DF_10 <- dam_scaff_order_DF %>% filter(Ref_name == 10) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
dam_scaff_order_DF_10 <- dam_scaff_order_DF %>% filter(Ref_name == 10) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 10) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
dam_scaff_order_DF_10 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
dam_scaff_order_DF_10 <- dam_scaff_order_DF_10 %>% select(component_id,clength,orientation)
write_tsv(dam_scaff_order_DF_10,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/dam_scaff_order_DF_10.tsv")

dam_scaff_order_DF_11 <- dam_scaff_order_DF %>% filter(Ref_name == 11) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
dam_scaff_order_DF_11 <- dam_scaff_order_DF %>% filter(Ref_name == 11) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 11) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
dam_scaff_order_DF_11 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
dam_scaff_order_DF_11 <- dam_scaff_order_DF_11 %>% select(component_id,clength,orientation)
write_tsv(dam_scaff_order_DF_11,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/dam_scaff_order_DF_11.tsv")

dam_scaff_order_DF_12 <- dam_scaff_order_DF %>% filter(Ref_name == 12) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
dam_scaff_order_DF_12 <- dam_scaff_order_DF %>% filter(Ref_name == 12) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 12) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
dam_scaff_order_DF_12 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
dam_scaff_order_DF_12 <- dam_scaff_order_DF_12 %>% select(component_id,clength,orientation)
write_tsv(dam_scaff_order_DF_12,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/dam_scaff_order_DF_12.tsv")

dam_scaff_order_DF_13 <- dam_scaff_order_DF %>% filter(Ref_name == 13) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
dam_scaff_order_DF_13 <- dam_scaff_order_DF %>% filter(Ref_name == 13) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 13) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
dam_scaff_order_DF_13 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
dam_scaff_order_DF_13 <- dam_scaff_order_DF_13 %>% select(component_id,clength,orientation)
write_tsv(dam_scaff_order_DF_13,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/dam_scaff_order_DF_13.tsv")

dam_scaff_order_DF_14 <- dam_scaff_order_DF %>% filter(Ref_name == 14) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
dam_scaff_order_DF_14 <- dam_scaff_order_DF %>% filter(Ref_name == 14) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 14) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
dam_scaff_order_DF_14 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
dam_scaff_order_DF_14 <- dam_scaff_order_DF_14 %>% select(component_id,clength,orientation)
write_tsv(dam_scaff_order_DF_14,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/dam_scaff_order_DF_14.tsv")

dam_scaff_order_DF_15 <- dam_scaff_order_DF %>% filter(Ref_name == 15) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
dam_scaff_order_DF_15 <- dam_scaff_order_DF %>% filter(Ref_name == 15) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 15) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
dam_scaff_order_DF_15 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
dam_scaff_order_DF_15 <- dam_scaff_order_DF_15 %>% select(component_id,clength,orientation)
write_tsv(dam_scaff_order_DF_15,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/dam_scaff_order_DF_15.tsv")

dam_scaff_order_DF_16 <- dam_scaff_order_DF %>% filter(Ref_name == 16) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
dam_scaff_order_DF_16 <- dam_scaff_order_DF %>% filter(Ref_name == 16) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 16) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
dam_scaff_order_DF_16 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
dam_scaff_order_DF_16 <- dam_scaff_order_DF_16 %>% select(component_id,clength,orientation)
write_tsv(dam_scaff_order_DF_16,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/dam_scaff_order_DF_16.tsv")

dam_scaff_order_DF_17 <- dam_scaff_order_DF %>% filter(Ref_name == 17) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
dam_scaff_order_DF_17 <- dam_scaff_order_DF %>% filter(Ref_name == 17) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 17) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
dam_scaff_order_DF_17 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
dam_scaff_order_DF_17 <- dam_scaff_order_DF_17 %>% select(component_id,clength,orientation)
write_tsv(dam_scaff_order_DF_17,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/dam_scaff_order_DF_17.tsv")

dam_scaff_order_DF_18 <- dam_scaff_order_DF %>% filter(Ref_name == 18) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
dam_scaff_order_DF_18 <- dam_scaff_order_DF %>% filter(Ref_name == 18) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 18) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
dam_scaff_order_DF_18 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
dam_scaff_order_DF_18 <- dam_scaff_order_DF_18 %>% select(component_id,clength,orientation)
write_tsv(dam_scaff_order_DF_18,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/dam_scaff_order_DF_18.tsv")

dam_scaff_order_DF_19 <- dam_scaff_order_DF %>% filter(Ref_name == 19) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
dam_scaff_order_DF_19 <- dam_scaff_order_DF %>% filter(Ref_name == 19) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 19) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
dam_scaff_order_DF_19 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
dam_scaff_order_DF_19 <- dam_scaff_order_DF_19 %>% select(component_id,clength,orientation)
write_tsv(dam_scaff_order_DF_19,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/dam_scaff_order_DF_19.tsv")

dam_scaff_order_DF_20 <- dam_scaff_order_DF %>% filter(Ref_name == 20) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
dam_scaff_order_DF_20 <- dam_scaff_order_DF %>% filter(Ref_name == 20) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 20) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
dam_scaff_order_DF_20 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
dam_scaff_order_DF_20 <- dam_scaff_order_DF_20 %>% select(component_id,clength,orientation)
write_tsv(dam_scaff_order_DF_20,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/dam_scaff_order_DF_20.tsv")

dam_scaff_order_DF_21 <- dam_scaff_order_DF %>% filter(Ref_name == 21) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
dam_scaff_order_DF_21 <- dam_scaff_order_DF %>% filter(Ref_name == 21) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 21) %>% filter(rc_n >= 6) %>% arrange(Ref_start)
#special rc probes == 6 bcos of hiC
#check proportion
dam_scaff_order_DF_21 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
dam_scaff_order_DF_21 <- dam_scaff_order_DF_21 %>% select(component_id,clength,orientation)
write_tsv(dam_scaff_order_DF_21,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/dam_scaff_order_DF_21.tsv")

dam_scaff_order_DF_22 <- dam_scaff_order_DF %>% filter(Ref_name == 22) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
dam_scaff_order_DF_22 <- dam_scaff_order_DF %>% filter(Ref_name == 22) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 22) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
dam_scaff_order_DF_22 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
dam_scaff_order_DF_22 <- dam_scaff_order_DF_22 %>% select(component_id,clength,orientation)
write_tsv(dam_scaff_order_DF_22,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/dam_scaff_order_DF_22.tsv")

dam_scaff_order_DF_23 <- dam_scaff_order_DF %>% filter(Ref_name == 23) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
dam_scaff_order_DF_23 <- dam_scaff_order_DF %>% filter(Ref_name == 23) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 23) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
dam_scaff_order_DF_23 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
dam_scaff_order_DF_23 <- dam_scaff_order_DF_23 %>% select(component_id,clength,orientation)
write_tsv(dam_scaff_order_DF_23,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/dam_scaff_order_DF_23.tsv")

dam_scaff_order_DF_24 <- dam_scaff_order_DF %>% filter(Ref_name == 24) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
dam_scaff_order_DF_24 <- dam_scaff_order_DF %>% filter(Ref_name == 24) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 24) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
dam_scaff_order_DF_24 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
dam_scaff_order_DF_24 <- dam_scaff_order_DF_24 %>% select(component_id,clength,orientation)
write_tsv(dam_scaff_order_DF_24,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/dam_scaff_order_DF_24.tsv")

dam_scaff_order_DF_25 <- dam_scaff_order_DF %>% filter(Ref_name == 25) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
dam_scaff_order_DF_25 <- dam_scaff_order_DF %>% filter(Ref_name == 25) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 25) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
dam_scaff_order_DF_25 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
dam_scaff_order_DF_25 <- dam_scaff_order_DF_25 %>% select(component_id,clength,orientation)
write_tsv(dam_scaff_order_DF_25,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/dam_scaff_order_DF_25.tsv")

dam_scaff_order_DF_26 <- dam_scaff_order_DF %>% filter(Ref_name == 26) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
dam_scaff_order_DF_26 <- dam_scaff_order_DF %>% filter(Ref_name == 26) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 26) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
dam_scaff_order_DF_26 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
dam_scaff_order_DF_26 <- dam_scaff_order_DF_26 %>% select(component_id,clength,orientation)
write_tsv(dam_scaff_order_DF_26,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/dam_scaff_order_DF_26.tsv")

dam_scaff_order_DF_27 <- dam_scaff_order_DF %>% filter(Ref_name == 27) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
dam_scaff_order_DF_27 <- dam_scaff_order_DF %>% filter(Ref_name == 27) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 27) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
dam_scaff_order_DF_27 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
dam_scaff_order_DF_27 <- dam_scaff_order_DF_27 %>% select(component_id,clength,orientation)
write_tsv(dam_scaff_order_DF_27,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/dam_scaff_order_DF_27.tsv")

dam_scaff_order_DF_28 <- dam_scaff_order_DF %>% filter(Ref_name == 28) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
dam_scaff_order_DF_28 <- dam_scaff_order_DF %>% filter(Ref_name == 28) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 28) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
dam_scaff_order_DF_28 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
dam_scaff_order_DF_28 <- dam_scaff_order_DF_28 %>% select(component_id,clength,orientation)
write_tsv(dam_scaff_order_DF_28,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/dam_scaff_order_DF_28.tsv")

dam_scaff_order_DF_29 <- dam_scaff_order_DF %>% filter(Ref_name == 29) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
dam_scaff_order_DF_29 <- dam_scaff_order_DF %>% filter(Ref_name == 29) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 29) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
dam_scaff_order_DF_29 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
dam_scaff_order_DF_29 <- dam_scaff_order_DF_29 %>% select(component_id,clength,orientation)
write_tsv(dam_scaff_order_DF_29,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/dam_scaff_order_DF_29.tsv")
