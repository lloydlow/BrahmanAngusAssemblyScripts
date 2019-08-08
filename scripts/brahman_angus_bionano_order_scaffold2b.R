#------------------------------------------------------
# Program name: brahman_angus_bionano_order_scaffold2b.R
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

##### Sire #####
sire_scaff_order_DF_1 <- sire_scaff_order_DF %>% filter(Ref_name == 1) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
sire_scaff_order_DF_1 <- sire_scaff_order_DF %>% filter(Ref_name == 1) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 1) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
sire_scaff_order_DF_1 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
sire_scaff_order_DF_1 <- sire_scaff_order_DF_1 %>% select(component_id,clength,orientation)
write_tsv(sire_scaff_order_DF_1,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/sire_scaff_order_DF_1.tsv")

sire_scaff_order_DF_2 <- sire_scaff_order_DF %>% filter(Ref_name == 2) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
sire_scaff_order_DF_2 <- sire_scaff_order_DF %>% filter(Ref_name == 2) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 2) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
sire_scaff_order_DF_2 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
sire_scaff_order_DF_2 <- sire_scaff_order_DF_2 %>% select(component_id,clength,orientation)
write_tsv(sire_scaff_order_DF_2,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/sire_scaff_order_DF_2.tsv")

sire_scaff_order_DF_3 <- sire_scaff_order_DF %>% filter(Ref_name == 3) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
sire_scaff_order_DF_3 <- sire_scaff_order_DF %>% filter(Ref_name == 3) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 3) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
sire_scaff_order_DF_3 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
sire_scaff_order_DF_3 <- sire_scaff_order_DF_3 %>% select(component_id,clength,orientation)
write_tsv(sire_scaff_order_DF_3,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/sire_scaff_order_DF_3.tsv")

sire_scaff_order_DF_4 <- sire_scaff_order_DF %>% filter(Ref_name == 4) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
sire_scaff_order_DF_4 <- sire_scaff_order_DF %>% filter(Ref_name == 4) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 4) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
sire_scaff_order_DF_4 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
sire_scaff_order_DF_4 <- sire_scaff_order_DF_4 %>% select(component_id,clength,orientation)
write_tsv(sire_scaff_order_DF_4,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/sire_scaff_order_DF_4.tsv")

sire_scaff_order_DF_5 <- sire_scaff_order_DF %>% filter(Ref_name == 5) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
sire_scaff_order_DF_5 <- sire_scaff_order_DF %>% filter(Ref_name == 5) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 5) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
sire_scaff_order_DF_5 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
sire_scaff_order_DF_5 <- sire_scaff_order_DF_5 %>% select(component_id,clength,orientation)
write_tsv(sire_scaff_order_DF_5,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/sire_scaff_order_DF_5.tsv")

sire_scaff_order_DF_6 <- sire_scaff_order_DF %>% filter(Ref_name == 6) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
sire_scaff_order_DF_6 <- sire_scaff_order_DF %>% filter(Ref_name == 6) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 6) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
sire_scaff_order_DF_6 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
sire_scaff_order_DF_6 <- sire_scaff_order_DF_6 %>% select(component_id,clength,orientation)
write_tsv(sire_scaff_order_DF_6,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/sire_scaff_order_DF_6.tsv")

sire_scaff_order_DF_7 <- sire_scaff_order_DF %>% filter(Ref_name == 7) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
sire_scaff_order_DF_7 <- sire_scaff_order_DF %>% filter(Ref_name == 7) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 7) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
sire_scaff_order_DF_7 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
sire_scaff_order_DF_7 <- sire_scaff_order_DF_7 %>% select(component_id,clength,orientation)
write_tsv(sire_scaff_order_DF_7,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/sire_scaff_order_DF_7.tsv")

sire_scaff_order_DF_8 <- sire_scaff_order_DF %>% filter(Ref_name == 8) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
sire_scaff_order_DF_8 <- sire_scaff_order_DF %>% filter(Ref_name == 8) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 8) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
sire_scaff_order_DF_8 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
sire_scaff_order_DF_8 <- sire_scaff_order_DF_8 %>% select(component_id,clength,orientation)
write_tsv(sire_scaff_order_DF_8,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/sire_scaff_order_DF_8.tsv")

sire_scaff_order_DF_9 <- sire_scaff_order_DF %>% filter(Ref_name == 9) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
sire_scaff_order_DF_9 <- sire_scaff_order_DF %>% filter(Ref_name == 9) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 9) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
sire_scaff_order_DF_9 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
sire_scaff_order_DF_9 <- sire_scaff_order_DF_9 %>% select(component_id,clength,orientation)
write_tsv(sire_scaff_order_DF_9,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/sire_scaff_order_DF_9.tsv")

sire_scaff_order_DF_10 <- sire_scaff_order_DF %>% filter(Ref_name == 10) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
sire_scaff_order_DF_10 <- sire_scaff_order_DF %>% filter(Ref_name == 10) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 10) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
sire_scaff_order_DF_10 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
sire_scaff_order_DF_10 <- sire_scaff_order_DF_10 %>% select(component_id,clength,orientation)
write_tsv(sire_scaff_order_DF_10,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/sire_scaff_order_DF_10.tsv")

sire_scaff_order_DF_11 <- sire_scaff_order_DF %>% filter(Ref_name == 11) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
sire_scaff_order_DF_11 <- sire_scaff_order_DF %>% filter(Ref_name == 11) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 11) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
sire_scaff_order_DF_11 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
sire_scaff_order_DF_11 <- sire_scaff_order_DF_11 %>% select(component_id,clength,orientation)
write_tsv(sire_scaff_order_DF_11,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/sire_scaff_order_DF_11.tsv")

sire_scaff_order_DF_12 <- sire_scaff_order_DF %>% filter(Ref_name == 12) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
sire_scaff_order_DF_12 <- sire_scaff_order_DF %>% filter(Ref_name == 12) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 12) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
sire_scaff_order_DF_12 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
sire_scaff_order_DF_12 <- sire_scaff_order_DF_12 %>% select(component_id,clength,orientation)
write_tsv(sire_scaff_order_DF_12,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/sire_scaff_order_DF_12.tsv")

sire_scaff_order_DF_13 <- sire_scaff_order_DF %>% filter(Ref_name == 13) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
sire_scaff_order_DF_13 <- sire_scaff_order_DF %>% filter(Ref_name == 13) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 13) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
sire_scaff_order_DF_13 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
sire_scaff_order_DF_13 <- sire_scaff_order_DF_13 %>% select(component_id,clength,orientation)
write_tsv(sire_scaff_order_DF_13,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/sire_scaff_order_DF_13.tsv")

sire_scaff_order_DF_14 <- sire_scaff_order_DF %>% filter(Ref_name == 14) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
sire_scaff_order_DF_14 <- sire_scaff_order_DF %>% filter(Ref_name == 14) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 14) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
sire_scaff_order_DF_14 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
sire_scaff_order_DF_14 <- sire_scaff_order_DF_14 %>% select(component_id,clength,orientation)
write_tsv(sire_scaff_order_DF_14,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/sire_scaff_order_DF_14.tsv")

sire_scaff_order_DF_15 <- sire_scaff_order_DF %>% filter(Ref_name == 15) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
sire_scaff_order_DF_15 <- sire_scaff_order_DF %>% filter(Ref_name == 15) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 15) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
sire_scaff_order_DF_15 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
sire_scaff_order_DF_15 <- sire_scaff_order_DF_15 %>% select(component_id,clength,orientation)
write_tsv(sire_scaff_order_DF_15,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/sire_scaff_order_DF_15.tsv")

sire_scaff_order_DF_16 <- sire_scaff_order_DF %>% filter(Ref_name == 16) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
sire_scaff_order_DF_16 <- sire_scaff_order_DF %>% filter(Ref_name == 16) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 16) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
sire_scaff_order_DF_16 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
sire_scaff_order_DF_16 <- sire_scaff_order_DF_16 %>% select(component_id,clength,orientation)
write_tsv(sire_scaff_order_DF_16,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/sire_scaff_order_DF_16.tsv")

sire_scaff_order_DF_17 <- sire_scaff_order_DF %>% filter(Ref_name == 17) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
sire_scaff_order_DF_17 <- sire_scaff_order_DF %>% filter(Ref_name == 17) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 17) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
sire_scaff_order_DF_17 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
sire_scaff_order_DF_17 <- sire_scaff_order_DF_17 %>% select(component_id,clength,orientation)
write_tsv(sire_scaff_order_DF_17,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/sire_scaff_order_DF_17.tsv")

sire_scaff_order_DF_18 <- sire_scaff_order_DF %>% filter(Ref_name == 18) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
sire_scaff_order_DF_18 <- sire_scaff_order_DF %>% filter(Ref_name == 18) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 18) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
sire_scaff_order_DF_18 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
sire_scaff_order_DF_18 <- sire_scaff_order_DF_18 %>% select(component_id,clength,orientation)
write_tsv(sire_scaff_order_DF_18,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/sire_scaff_order_DF_18.tsv")

sire_scaff_order_DF_19 <- sire_scaff_order_DF %>% filter(Ref_name == 19) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
sire_scaff_order_DF_19 <- sire_scaff_order_DF %>% filter(Ref_name == 19) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 19) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
sire_scaff_order_DF_19 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
sire_scaff_order_DF_19 <- sire_scaff_order_DF_19 %>% select(component_id,clength,orientation)
write_tsv(sire_scaff_order_DF_19,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/sire_scaff_order_DF_19.tsv")

sire_scaff_order_DF_20 <- sire_scaff_order_DF %>% filter(Ref_name == 20) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
sire_scaff_order_DF_20 <- sire_scaff_order_DF %>% filter(Ref_name == 20) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 20) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
sire_scaff_order_DF_20 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
sire_scaff_order_DF_20 <- sire_scaff_order_DF_20 %>% select(component_id,clength,orientation)
write_tsv(sire_scaff_order_DF_20,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/sire_scaff_order_DF_20.tsv")

sire_scaff_order_DF_21 <- sire_scaff_order_DF %>% filter(Ref_name == 21) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
sire_scaff_order_DF_21 <- sire_scaff_order_DF %>% filter(Ref_name == 21) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 21) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
sire_scaff_order_DF_21 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
sire_scaff_order_DF_21 <- sire_scaff_order_DF_21 %>% select(component_id,clength,orientation)
write_tsv(sire_scaff_order_DF_21,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/sire_scaff_order_DF_21.tsv")

sire_scaff_order_DF_22 <- sire_scaff_order_DF %>% filter(Ref_name == 22) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
sire_scaff_order_DF_22 <- sire_scaff_order_DF %>% filter(Ref_name == 22) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 22) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
sire_scaff_order_DF_22 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
sire_scaff_order_DF_22 <- sire_scaff_order_DF_22 %>% select(component_id,clength,orientation)
write_tsv(sire_scaff_order_DF_22,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/sire_scaff_order_DF_22.tsv")

sire_scaff_order_DF_23 <- sire_scaff_order_DF %>% filter(Ref_name == 23) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
sire_scaff_order_DF_23 <- sire_scaff_order_DF %>% filter(Ref_name == 23) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 23) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
sire_scaff_order_DF_23 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
sire_scaff_order_DF_23 <- sire_scaff_order_DF_23 %>% select(component_id,clength,orientation)
write_tsv(sire_scaff_order_DF_23,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/sire_scaff_order_DF_23.tsv")

sire_scaff_order_DF_24 <- sire_scaff_order_DF %>% filter(Ref_name == 24) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
sire_scaff_order_DF_24 <- sire_scaff_order_DF %>% filter(Ref_name == 24) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 24) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
sire_scaff_order_DF_24 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
sire_scaff_order_DF_24 <- sire_scaff_order_DF_24 %>% select(component_id,clength,orientation)
write_tsv(sire_scaff_order_DF_24,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/sire_scaff_order_DF_24.tsv")

sire_scaff_order_DF_25 <- sire_scaff_order_DF %>% filter(Ref_name == 25) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
sire_scaff_order_DF_25 <- sire_scaff_order_DF %>% filter(Ref_name == 25) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 25) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
sire_scaff_order_DF_25 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
sire_scaff_order_DF_25 <- sire_scaff_order_DF_25 %>% select(component_id,clength,orientation)
write_tsv(sire_scaff_order_DF_25,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/sire_scaff_order_DF_25.tsv")

sire_scaff_order_DF_26 <- sire_scaff_order_DF %>% filter(Ref_name == 26) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
sire_scaff_order_DF_26 <- sire_scaff_order_DF %>% filter(Ref_name == 26) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 26) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
sire_scaff_order_DF_26 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
sire_scaff_order_DF_26 <- sire_scaff_order_DF_26 %>% select(component_id,clength,orientation)
write_tsv(sire_scaff_order_DF_26,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/sire_scaff_order_DF_26.tsv")

sire_scaff_order_DF_27 <- sire_scaff_order_DF %>% filter(Ref_name == 27) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
sire_scaff_order_DF_27 <- sire_scaff_order_DF %>% filter(Ref_name == 27) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 27) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
sire_scaff_order_DF_27 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
sire_scaff_order_DF_27 <- sire_scaff_order_DF_27 %>% select(component_id,clength,orientation)
write_tsv(sire_scaff_order_DF_27,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/sire_scaff_order_DF_27.tsv")

sire_scaff_order_DF_28 <- sire_scaff_order_DF %>% filter(Ref_name == 28) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
sire_scaff_order_DF_28 <- sire_scaff_order_DF %>% filter(Ref_name == 28) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 28) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
sire_scaff_order_DF_28 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
sire_scaff_order_DF_28 <- sire_scaff_order_DF_28 %>% select(component_id,clength,orientation)
write_tsv(sire_scaff_order_DF_28,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/sire_scaff_order_DF_28.tsv")

sire_scaff_order_DF_29 <- sire_scaff_order_DF %>% filter(Ref_name == 29) %>% arrange(desc(scaff_len)) 
#test various conditions to automatically arrange scaffolds
sire_scaff_order_DF_29 <- sire_scaff_order_DF %>% filter(Ref_name == 29) %>% arrange(desc(scaff_len)) %>%
  filter(chr_rc == 29) %>% filter(rc_n >= 10) %>% arrange(Ref_start)

#check proportion
sire_scaff_order_DF_29 %>% group_by(object) %>% summarise(new_proportion = mean(proportion_ref)) %>% 
  summarise(total_proportion = sum(new_proportion))

#write out
sire_scaff_order_DF_29 <- sire_scaff_order_DF_29 %>% select(component_id,clength,orientation)
write_tsv(sire_scaff_order_DF_29,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_book_bionano_archived/sire_scaff_order_DF_29.tsv")



