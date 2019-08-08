#------------------------------------------------------
# Program name: brahman_angus_order_scaffold2.R
# Objective: after making big table now reodering scaffold
#           to best fit HD, rc, align Dom 
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(readr)
library(dplyr)

load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/order_scaffolds_to_chr_level/scaff_order_DF.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/order_scaffolds_to_chr_level/f1_dam_salsa_HD_rc_Probes_tab.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/order_scaffolds_to_chr_level/f1_sire_salsa_HD_rc_Probes_tab.RData")
#load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/brahman_Angus_mashmap/dam_agp_clean_assembly_to_salsa.RData")

##### Dam #####
#write_tsv(dam_scaff_order_DF,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/tig_to_salsa_contig_correspondence/dam_scaff_order_DF.tsv")
#below is a switch to a DF that has canu tigs included as a column
dam_scaff_order_canu_DF <- read_tsv("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/tig_to_salsa_contig_correspondence/dam_scaff_order_canu_DF.txt")
dam_scaff_order_DF <- dam_scaff_order_canu_DF

dam_scaff_order_DF_1 <- dam_scaff_order_DF %>% filter(Ref_name == 1) %>% arrange(desc(scaff_len)) %>% 
  select(18,1:17)
#scaffold_715
dam_scaff_order_DF_1_plot1 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_715")
plot(dam_scaff_order_DF_1_plot1$chr_pos,dam_scaff_order_DF_1_plot1$align_pos)
dam_scaff_order_DF_1_plot1 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_715")
plot(dam_scaff_order_DF_1_plot1$chr_pos,dam_scaff_order_DF_1_plot1$align_pos)

dam_scaff_order_DF_2 <- dam_scaff_order_DF %>% filter(Ref_name == 2) %>% arrange(desc(scaff_len)) %>%
  select(18,1:17)
#scaffold_608 scaffold_306
dam_scaff_order_DF_2_plot1 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_608")
plot(dam_scaff_order_DF_2_plot1$chr_pos,dam_scaff_order_DF_2_plot1$align_pos)
dam_scaff_order_DF_2_plot1 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_608")
plot(dam_scaff_order_DF_2_plot1$chr_pos,dam_scaff_order_DF_2_plot1$align_pos)

dam_scaff_order_DF_2_plot2 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_306")
plot(dam_scaff_order_DF_2_plot2$chr_pos,dam_scaff_order_DF_2_plot2$align_pos)
dam_scaff_order_DF_2_plot2 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_306")
plot(dam_scaff_order_DF_2_plot2$chr_pos,dam_scaff_order_DF_2_plot2$align_pos)
#note: scaffold_306 should be stiched in btween scaffold_608 based on hd and rc map

dam_scaff_order_DF_3 <- dam_scaff_order_DF %>% filter(Ref_name == 3) %>% arrange(desc(scaff_len)) %>%
  select(18,1:17)
#scaffold_1035 scaffold_745
dam_scaff_order_DF_3_plot1 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_1035")
plot(dam_scaff_order_DF_3_plot1$chr_pos,dam_scaff_order_DF_3_plot1$align_pos)
dam_scaff_order_DF_3_plot1 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_1035")
plot(dam_scaff_order_DF_3_plot1$chr_pos,dam_scaff_order_DF_3_plot1$align_pos)

dam_scaff_order_DF_3_plot2 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_745")
plot(dam_scaff_order_DF_3_plot2$chr_pos,dam_scaff_order_DF_3_plot2$align_pos)
dam_scaff_order_DF_3_plot2 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_745")
plot(dam_scaff_order_DF_3_plot2$chr_pos,dam_scaff_order_DF_3_plot2$align_pos)
#note: one inversion in scaffold_1035

dam_scaff_order_DF_4 <- dam_scaff_order_DF %>% filter(Ref_name == 4) %>% arrange(desc(scaff_len)) %>%
  select(18,1:17)
#scaffold_686 scaffold_592
dam_scaff_order_DF_4_plot1 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_686")
plot(dam_scaff_order_DF_4_plot1$chr_pos,dam_scaff_order_DF_4_plot1$align_pos)
dam_scaff_order_DF_4_plot1 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_686")
plot(dam_scaff_order_DF_4_plot1$chr_pos,dam_scaff_order_DF_4_plot1$align_pos)

dam_scaff_order_DF_4_plot2 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_592")
plot(dam_scaff_order_DF_4_plot2$chr_pos,dam_scaff_order_DF_4_plot2$align_pos)
dam_scaff_order_DF_4_plot2 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_592")
plot(dam_scaff_order_DF_4_plot2$chr_pos,dam_scaff_order_DF_4_plot2$align_pos)

dam_scaff_order_DF_5 <- dam_scaff_order_DF %>% filter(Ref_name == 5) %>% arrange(desc(scaff_len)) %>%
  select(18,1:17)
#scaffold_136 scaffold_681
dam_scaff_order_DF_5_plot1 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_136")
plot(dam_scaff_order_DF_5_plot1$chr_pos,dam_scaff_order_DF_5_plot1$align_pos)
dam_scaff_order_DF_5_plot1 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_136")
plot(dam_scaff_order_DF_5_plot1$chr_pos,dam_scaff_order_DF_5_plot1$align_pos)

dam_scaff_order_DF_5_plot2 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_681")
plot(dam_scaff_order_DF_5_plot2$chr_pos,dam_scaff_order_DF_5_plot2$align_pos)
dam_scaff_order_DF_5_plot2 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_681")
plot(dam_scaff_order_DF_5_plot2$chr_pos,dam_scaff_order_DF_5_plot2$align_pos)
#note: 2 rearrangements that could potentially be fixed

dam_scaff_order_DF_6 <- dam_scaff_order_DF %>% filter(Ref_name == 6) %>% arrange(desc(scaff_len)) %>%
  select(18,1:17)
#scaffold_1269 scaffold_37 
dam_scaff_order_DF_6_plot1 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_1269")
plot(dam_scaff_order_DF_6_plot1$chr_pos,dam_scaff_order_DF_6_plot1$align_pos)
dam_scaff_order_DF_6_plot1 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_1269")
plot(dam_scaff_order_DF_6_plot1$chr_pos,dam_scaff_order_DF_6_plot1$align_pos)

dam_scaff_order_DF_6_plot2 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_37")
plot(dam_scaff_order_DF_6_plot2$chr_pos,dam_scaff_order_DF_6_plot2$align_pos)
dam_scaff_order_DF_6_plot2 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_37")
plot(dam_scaff_order_DF_6_plot2$chr_pos,dam_scaff_order_DF_6_plot2$align_pos)

dam_scaff_order_DF_7 <- dam_scaff_order_DF %>% filter(Ref_name == 7) %>% arrange(desc(scaff_len)) %>%
  select(18,1:17)
#scaffold_1185 scaffold_584 
dam_scaff_order_DF_7_plot1 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_1185")
plot(dam_scaff_order_DF_7_plot1$chr_pos,dam_scaff_order_DF_7_plot1$align_pos)
dam_scaff_order_DF_7_plot1 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_1185")
plot(dam_scaff_order_DF_7_plot1$chr_pos,dam_scaff_order_DF_7_plot1$align_pos)

dam_scaff_order_DF_7_plot2 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_584")
plot(dam_scaff_order_DF_7_plot2$chr_pos,dam_scaff_order_DF_7_plot2$align_pos)
dam_scaff_order_DF_7_plot2 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_584")
plot(dam_scaff_order_DF_7_plot2$chr_pos,dam_scaff_order_DF_7_plot2$align_pos)

dam_scaff_order_DF_8 <- dam_scaff_order_DF %>% filter(Ref_name == 8) %>% arrange(desc(scaff_len)) %>%
  select(18,1:17)
#scaffold_204 scaffold_1164
dam_scaff_order_DF_8_plot1 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_204")
plot(dam_scaff_order_DF_8_plot1$chr_pos,dam_scaff_order_DF_8_plot1$align_pos)
dam_scaff_order_DF_8_plot1 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_204")
plot(dam_scaff_order_DF_8_plot1$chr_pos,dam_scaff_order_DF_8_plot1$align_pos)

dam_scaff_order_DF_8_plot2 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_1164")
plot(dam_scaff_order_DF_8_plot2$chr_pos,dam_scaff_order_DF_8_plot2$align_pos)
dam_scaff_order_DF_8_plot2 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_1164")
plot(dam_scaff_order_DF_8_plot2$chr_pos,dam_scaff_order_DF_8_plot2$align_pos)

dam_scaff_order_DF_9 <- dam_scaff_order_DF %>% filter(Ref_name == 9) %>% arrange(desc(scaff_len)) %>%
  select(18,1:17)
#scaffold_69 scaffold_36
dam_scaff_order_DF_9_plot1 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_69")
plot(dam_scaff_order_DF_9_plot1$chr_pos,dam_scaff_order_DF_9_plot1$align_pos)
dam_scaff_order_DF_9_plot1 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_69")
plot(dam_scaff_order_DF_9_plot1$chr_pos,dam_scaff_order_DF_9_plot1$align_pos)

dam_scaff_order_DF_9_plot2 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_36")
plot(dam_scaff_order_DF_9_plot2$chr_pos,dam_scaff_order_DF_9_plot2$align_pos)
dam_scaff_order_DF_9_plot2 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_36")
plot(dam_scaff_order_DF_9_plot2$chr_pos,dam_scaff_order_DF_9_plot2$align_pos)

dam_scaff_order_DF_10 <- dam_scaff_order_DF %>% filter(Ref_name == 10) %>% arrange(desc(scaff_len)) %>%
  select(18,1:17)
#scaffold_65
dam_scaff_order_DF_10_plot1 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_65")
plot(dam_scaff_order_DF_10_plot1$chr_pos,dam_scaff_order_DF_10_plot1$align_pos)
dam_scaff_order_DF_10_plot1 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_65")
plot(dam_scaff_order_DF_10_plot1$chr_pos,dam_scaff_order_DF_10_plot1$align_pos)

dam_scaff_order_DF_11 <- dam_scaff_order_DF %>% filter(Ref_name == 11) %>% arrange(desc(scaff_len)) %>%
  select(18,1:17)
#scaffold_1263 scaffold_344
dam_scaff_order_DF_11_plot1 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_1263")
plot(dam_scaff_order_DF_11_plot1$chr_pos,dam_scaff_order_DF_11_plot1$align_pos)
dam_scaff_order_DF_11_plot1 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_1263")
plot(dam_scaff_order_DF_11_plot1$chr_pos,dam_scaff_order_DF_11_plot1$align_pos)

dam_scaff_order_DF_11_plot2 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_344")
plot(dam_scaff_order_DF_11_plot2$chr_pos,dam_scaff_order_DF_11_plot2$align_pos)
dam_scaff_order_DF_11_plot2 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_344")
plot(dam_scaff_order_DF_11_plot2$chr_pos,dam_scaff_order_DF_11_plot2$align_pos)

dam_scaff_order_DF_12 <- dam_scaff_order_DF %>% filter(Ref_name == 12) %>% arrange(desc(scaff_len)) %>%
  select(18,1:17)
#scaffold_671 scaffold_186
dam_scaff_order_DF_12_plot1 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_671")
plot(dam_scaff_order_DF_12_plot1$chr_pos,dam_scaff_order_DF_12_plot1$align_pos)
dam_scaff_order_DF_12_plot1 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_671")
plot(dam_scaff_order_DF_12_plot1$chr_pos,dam_scaff_order_DF_12_plot1$align_pos)

dam_scaff_order_DF_12_plot2 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_186")
plot(dam_scaff_order_DF_12_plot2$chr_pos,dam_scaff_order_DF_12_plot2$align_pos)
dam_scaff_order_DF_12_plot2 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_186")
plot(dam_scaff_order_DF_12_plot2$chr_pos,dam_scaff_order_DF_12_plot2$align_pos)

dam_scaff_order_DF_13 <- dam_scaff_order_DF %>% filter(Ref_name == 13) %>% arrange(desc(scaff_len)) %>%
  select(18,1:17)
#scaffold_700 scaffold_781 scaffold_539 partial
dam_scaff_order_DF_13_plot1 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_700")
plot(dam_scaff_order_DF_13_plot1$chr_pos,dam_scaff_order_DF_13_plot1$align_pos)
dam_scaff_order_DF_13_plot1 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_700")
plot(dam_scaff_order_DF_13_plot1$chr_pos,dam_scaff_order_DF_13_plot1$align_pos)

dam_scaff_order_DF_13_plot2 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_781")
plot(dam_scaff_order_DF_13_plot2$chr_pos,dam_scaff_order_DF_13_plot2$align_pos)
dam_scaff_order_DF_13_plot2 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_781")
plot(dam_scaff_order_DF_13_plot2$chr_pos,dam_scaff_order_DF_13_plot2$align_pos)

dam_scaff_order_DF_14 <- dam_scaff_order_DF %>% filter(Ref_name == 14) %>% arrange(desc(scaff_len)) %>%
  select(18,1:17)
#scaffold_400
dam_scaff_order_DF_14_plot1 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_400")
plot(dam_scaff_order_DF_14_plot1$chr_pos,dam_scaff_order_DF_14_plot1$align_pos)
dam_scaff_order_DF_14_plot1 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_400")
plot(dam_scaff_order_DF_14_plot1$chr_pos,dam_scaff_order_DF_14_plot1$align_pos)

dam_scaff_order_DF_15 <- dam_scaff_order_DF %>% filter(Ref_name == 15) %>% arrange(desc(scaff_len)) %>%
  select(18,1:17)
#scaffold_539 scaffold_901 scaffold_67
dam_scaff_order_DF_15_plot1 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_539")
plot(dam_scaff_order_DF_15_plot1$chr_pos,dam_scaff_order_DF_15_plot1$align_pos)
dam_scaff_order_DF_15_plot1 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_539")
plot(dam_scaff_order_DF_15_plot1$chr_pos,dam_scaff_order_DF_15_plot1$align_pos)

dam_scaff_order_DF_15_plot2 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_901")
plot(dam_scaff_order_DF_15_plot2$chr_pos,dam_scaff_order_DF_15_plot2$align_pos)
dam_scaff_order_DF_15_plot2 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_901")
plot(dam_scaff_order_DF_15_plot2$chr_pos,dam_scaff_order_DF_15_plot2$align_pos)

dam_scaff_order_DF_15_plot3 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_67")
plot(dam_scaff_order_DF_15_plot3$chr_pos,dam_scaff_order_DF_15_plot3$align_pos)
dam_scaff_order_DF_15_plot3 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_67")
plot(dam_scaff_order_DF_15_plot3$chr_pos,dam_scaff_order_DF_15_plot3$align_pos)

dam_scaff_order_DF_16 <- dam_scaff_order_DF %>% filter(Ref_name == 16) %>% arrange(desc(scaff_len))  %>%
  select(18,1:17)
#scaffold_1052 scaffold_1311 scaffold_727
dam_scaff_order_DF_16_plot1 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_1052")
plot(dam_scaff_order_DF_16_plot1$chr_pos,dam_scaff_order_DF_16_plot1$align_pos)
dam_scaff_order_DF_16_plot1 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_1052")
plot(dam_scaff_order_DF_16_plot1$chr_pos,dam_scaff_order_DF_16_plot1$align_pos)

dam_scaff_order_DF_16_plot2 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_1311")
plot(dam_scaff_order_DF_16_plot2$chr_pos,dam_scaff_order_DF_16_plot2$align_pos)
dam_scaff_order_DF_16_plot2 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_1311")
plot(dam_scaff_order_DF_16_plot2$chr_pos,dam_scaff_order_DF_16_plot2$align_pos)

dam_scaff_order_DF_16_plot3 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_727")
plot(dam_scaff_order_DF_16_plot3$chr_pos,dam_scaff_order_DF_16_plot3$align_pos)
dam_scaff_order_DF_16_plot3 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_727")
plot(dam_scaff_order_DF_16_plot3$chr_pos,dam_scaff_order_DF_16_plot3$align_pos)

dam_scaff_order_DF_17 <- dam_scaff_order_DF %>% filter(Ref_name == 17) %>% arrange(desc(scaff_len))  %>%
  select(18,1:17)
#scaffold_969
dam_scaff_order_DF_17_plot1 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_969")
plot(dam_scaff_order_DF_17_plot1$chr_pos,dam_scaff_order_DF_17_plot1$align_pos)
dam_scaff_order_DF_17_plot1 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_969")
plot(dam_scaff_order_DF_17_plot1$chr_pos,dam_scaff_order_DF_17_plot1$align_pos)

dam_scaff_order_DF_18 <- dam_scaff_order_DF %>% filter(Ref_name == 18) %>% arrange(desc(scaff_len)) %>%
  select(18,1:17)
#scaffold_516 scaffold_313 scaffold_120 scaffold_91
dam_scaff_order_DF_18_plot1 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_516")
plot(dam_scaff_order_DF_18_plot1$chr_pos,dam_scaff_order_DF_18_plot1$align_pos)
dam_scaff_order_DF_18_plot1 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_516")
plot(dam_scaff_order_DF_18_plot1$chr_pos,dam_scaff_order_DF_18_plot1$align_pos)

dam_scaff_order_DF_18_plot2 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_313")
plot(dam_scaff_order_DF_18_plot2$chr_pos,dam_scaff_order_DF_18_plot2$align_pos)
dam_scaff_order_DF_18_plot2 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_313")
plot(dam_scaff_order_DF_18_plot2$chr_pos,dam_scaff_order_DF_18_plot2$align_pos)

dam_scaff_order_DF_18_plot3 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_120")
plot(dam_scaff_order_DF_18_plot3$chr_pos,dam_scaff_order_DF_18_plot3$align_pos)
dam_scaff_order_DF_18_plot3 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_120")
plot(dam_scaff_order_DF_18_plot3$chr_pos,dam_scaff_order_DF_18_plot3$align_pos)

dam_scaff_order_DF_18_plot4 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_91")
plot(dam_scaff_order_DF_18_plot4$chr_pos,dam_scaff_order_DF_18_plot4$align_pos)
dam_scaff_order_DF_18_plot4 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_91")
plot(dam_scaff_order_DF_18_plot4$chr_pos,dam_scaff_order_DF_18_plot4$align_pos)
#note: complicated rearrangements here but I wont break contigs within scaffold in this case

dam_scaff_order_DF_19 <- dam_scaff_order_DF %>% filter(Ref_name == 19) %>% arrange(desc(scaff_len)) %>%
  select(18,1:17)
#scaffold_623
dam_scaff_order_DF_19_plot1 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_623")
plot(dam_scaff_order_DF_19_plot1$chr_pos,dam_scaff_order_DF_19_plot1$align_pos)
dam_scaff_order_DF_19_plot1 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_623")
plot(dam_scaff_order_DF_19_plot1$chr_pos,dam_scaff_order_DF_19_plot1$align_pos)

dam_scaff_order_DF_20 <- dam_scaff_order_DF %>% filter(Ref_name == 20) %>% arrange(desc(scaff_len)) %>%
  select(18,1:17)
#scaffold_919
dam_scaff_order_DF_20_plot1 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_919")
plot(dam_scaff_order_DF_20_plot1$chr_pos,dam_scaff_order_DF_20_plot1$align_pos)
dam_scaff_order_DF_20_plot1 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_919")
plot(dam_scaff_order_DF_20_plot1$chr_pos,dam_scaff_order_DF_20_plot1$align_pos)

dam_scaff_order_DF_21 <- dam_scaff_order_DF %>% filter(Ref_name == 21) %>% arrange(desc(scaff_len)) %>%
  select(18,1:17)
#scaffold_5 scaffold_985
dam_scaff_order_DF_21_plot1 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_5")
plot(dam_scaff_order_DF_21_plot1$chr_pos,dam_scaff_order_DF_21_plot1$align_pos)
dam_scaff_order_DF_21_plot1 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_5")
plot(dam_scaff_order_DF_21_plot1$chr_pos,dam_scaff_order_DF_21_plot1$align_pos)

dam_scaff_order_DF_21_plot2 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_985")
plot(dam_scaff_order_DF_21_plot2$chr_pos,dam_scaff_order_DF_21_plot2$align_pos)
dam_scaff_order_DF_21_plot2 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_985")
plot(dam_scaff_order_DF_21_plot2$chr_pos,dam_scaff_order_DF_21_plot2$align_pos)

dam_scaff_order_DF_22 <- dam_scaff_order_DF %>% filter(Ref_name == 22) %>% arrange(desc(scaff_len)) %>%
  select(18,1:17)
#scaffold_442
dam_scaff_order_DF_22_plot1 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_442")
plot(dam_scaff_order_DF_22_plot1$chr_pos,dam_scaff_order_DF_22_plot1$align_pos)
dam_scaff_order_DF_22_plot1 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_442")
plot(dam_scaff_order_DF_22_plot1$chr_pos,dam_scaff_order_DF_22_plot1$align_pos)

dam_scaff_order_DF_23 <- dam_scaff_order_DF %>% filter(Ref_name == 23) %>% arrange(desc(scaff_len)) %>%
  select(18,1:17)
#scaffold_921 scaffold_513
dam_scaff_order_DF_23_plot1 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_921")
plot(dam_scaff_order_DF_23_plot1$chr_pos,dam_scaff_order_DF_23_plot1$align_pos)
dam_scaff_order_DF_23_plot1 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_921")
plot(dam_scaff_order_DF_23_plot1$chr_pos,dam_scaff_order_DF_23_plot1$align_pos)

dam_scaff_order_DF_23_plot2 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_513")
plot(dam_scaff_order_DF_23_plot2$chr_pos,dam_scaff_order_DF_23_plot2$align_pos)
dam_scaff_order_DF_23_plot2 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_513")
plot(dam_scaff_order_DF_23_plot2$chr_pos,dam_scaff_order_DF_23_plot2$align_pos)

dam_scaff_order_DF_24 <- dam_scaff_order_DF %>% filter(Ref_name == 24) %>% arrange(desc(scaff_len)) %>%
  select(18,1:17)
#scaffold_815
dam_scaff_order_DF_24_plot1 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_815")
plot(dam_scaff_order_DF_24_plot1$chr_pos,dam_scaff_order_DF_24_plot1$align_pos)
dam_scaff_order_DF_24_plot1 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_815")
plot(dam_scaff_order_DF_24_plot1$chr_pos,dam_scaff_order_DF_24_plot1$align_pos)

dam_scaff_order_DF_25 <- dam_scaff_order_DF %>% filter(Ref_name == 25) %>% arrange(desc(scaff_len))  %>%
  select(18,1:17)
#scaffold_280
dam_scaff_order_DF_25_plot1 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_280")
plot(dam_scaff_order_DF_25_plot1$chr_pos,dam_scaff_order_DF_25_plot1$align_pos)
dam_scaff_order_DF_25_plot1 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_280")
plot(dam_scaff_order_DF_25_plot1$chr_pos,dam_scaff_order_DF_25_plot1$align_pos)

dam_scaff_order_DF_26 <- dam_scaff_order_DF %>% filter(Ref_name == 26) %>% arrange(desc(scaff_len)) %>%
  select(18,1:17)
#scaffold_701
dam_scaff_order_DF_26_plot1 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_701")
plot(dam_scaff_order_DF_26_plot1$chr_pos,dam_scaff_order_DF_26_plot1$align_pos)
dam_scaff_order_DF_26_plot1 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_701")
plot(dam_scaff_order_DF_26_plot1$chr_pos,dam_scaff_order_DF_26_plot1$align_pos)

dam_scaff_order_DF_27 <- dam_scaff_order_DF %>% filter(Ref_name == 27) %>% arrange(desc(scaff_len)) %>%
  select(18,1:17)
#scaffold_1283
dam_scaff_order_DF_27_plot1 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_1283")
plot(dam_scaff_order_DF_27_plot1$chr_pos,dam_scaff_order_DF_27_plot1$align_pos)
dam_scaff_order_DF_27_plot1 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_1283")
plot(dam_scaff_order_DF_27_plot1$chr_pos,dam_scaff_order_DF_27_plot1$align_pos)

dam_scaff_order_DF_28 <- dam_scaff_order_DF %>% filter(Ref_name == 28) %>% arrange(desc(scaff_len)) %>%
  select(18,1:17)
#scaffold_474
dam_scaff_order_DF_28_plot1 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_474")
plot(dam_scaff_order_DF_28_plot1$chr_pos,dam_scaff_order_DF_28_plot1$align_pos)
dam_scaff_order_DF_28_plot1 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_474")
plot(dam_scaff_order_DF_28_plot1$chr_pos,dam_scaff_order_DF_28_plot1$align_pos)

dam_scaff_order_DF_29 <- dam_scaff_order_DF %>% filter(Ref_name == 29) %>% arrange(desc(scaff_len)) %>%
  select(18,1:17)
#scaffold_343
dam_scaff_order_DF_29_plot1 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_343")
plot(dam_scaff_order_DF_29_plot1$chr_pos,dam_scaff_order_DF_29_plot1$align_pos)
dam_scaff_order_DF_29_plot1 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_343")
plot(dam_scaff_order_DF_29_plot1$chr_pos,dam_scaff_order_DF_29_plot1$align_pos)

dam_scaff_order_DF_29_plot2 <- f1_dam_salsa_HDProbes_tab %>% filter(scaffold == "")
plot(dam_scaff_order_DF_29_plot2$chr_pos,dam_scaff_order_DF_29_plot2$align_pos)
dam_scaff_order_DF_29_plot2 <- f1_dam_salsa_fasta_rcmap_tab %>% filter(scaffold == "")
plot(dam_scaff_order_DF_29_plot2$chr_pos,dam_scaff_order_DF_29_plot2$align_pos)
