#------------------------------------------------------
# Program name: brahman_angus_order_scaffold2b.R
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

##### Sire #####
#write_tsv(sire_scaff_order_DF,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/tig_to_salsa_contig_correspondence/sire_scaff_order_DF.tsv")
#below is a switch to a DF that has canu tigs included as a column
sire_scaff_order_canu_DF <- read_tsv("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/tig_to_salsa_contig_correspondence/sire_scaff_order_canu_DF.txt")
sire_scaff_order_DF <- sire_scaff_order_canu_DF

sire_scaff_order_DF_1 <- sire_scaff_order_DF %>% filter(Ref_name == 1) %>% arrange(desc(scaff_len)) %>% 
  select(18,1:17)
#scaffold_96 scaffold_1146 scaffold_1151
sire_scaff_order_DF_1_plot1 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_96")
plot(sire_scaff_order_DF_1_plot1$align_pos,sire_scaff_order_DF_1_plot1$chr_pos)
sire_scaff_order_DF_1_plot1 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_96")
plot(sire_scaff_order_DF_1_plot1$align_pos,sire_scaff_order_DF_1_plot1$chr_pos)
abline(v=c(6427245,25882560,98806698,140985929))

sire_scaff_order_DF_1_plot2 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_1146")
plot(sire_scaff_order_DF_1_plot2$align_pos,sire_scaff_order_DF_1_plot2$chr_pos)
sire_scaff_order_DF_1_plot2 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_1146")
plot(sire_scaff_order_DF_1_plot2$align_pos,sire_scaff_order_DF_1_plot2$chr_pos)

sire_scaff_order_DF_1_plot3 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_1151")
plot(sire_scaff_order_DF_1_plot3$align_pos,sire_scaff_order_DF_1_plot3$chr_pos)
sire_scaff_order_DF_1_plot3 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_1151")
plot(sire_scaff_order_DF_1_plot3$align_pos,sire_scaff_order_DF_1_plot3$chr_pos)

sire_scaff_order_DF_2 <- sire_scaff_order_DF %>% filter(Ref_name == 2) %>% arrange(desc(scaff_len)) %>% 
  select(18,1:17)
#scaffold_520
sire_scaff_order_DF_2_plot1 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_520")
plot(sire_scaff_order_DF_2_plot1$align_pos,sire_scaff_order_DF_2_plot1$chr_pos)
sire_scaff_order_DF_2_plot1 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_520")
plot(sire_scaff_order_DF_2_plot1$align_pos,sire_scaff_order_DF_2_plot1$chr_pos)
abline(v=c(12790660,48282210,134183855,219048332,253622839))

sire_scaff_order_DF_3 <- sire_scaff_order_DF %>% filter(Ref_name == 3) %>% arrange(desc(scaff_len)) %>% 
  select(18,1:17)
#scaffold_824 scaffold_475 scaffold_794
sire_scaff_order_DF_3_plot1 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_824")
plot(sire_scaff_order_DF_3_plot1$align_pos,sire_scaff_order_DF_3_plot1$chr_pos)
sire_scaff_order_DF_3_plot1 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_824")
plot(sire_scaff_order_DF_3_plot1$align_pos,sire_scaff_order_DF_3_plot1$chr_pos)

sire_scaff_order_DF_3_plot2 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_475")
plot(sire_scaff_order_DF_3_plot2$align_pos,sire_scaff_order_DF_3_plot2$chr_pos)
sire_scaff_order_DF_3_plot2 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_475")
plot(sire_scaff_order_DF_3_plot2$align_pos,sire_scaff_order_DF_3_plot2$chr_pos)

sire_scaff_order_DF_3_plot3 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_794")
plot(sire_scaff_order_DF_3_plot3$align_pos,sire_scaff_order_DF_3_plot3$chr_pos)
sire_scaff_order_DF_3_plot3 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_794")
plot(sire_scaff_order_DF_3_plot3$align_pos,sire_scaff_order_DF_3_plot3$chr_pos)

sire_scaff_order_DF_4 <- sire_scaff_order_DF %>% filter(Ref_name == 4) %>% arrange(desc(scaff_len)) %>% 
  select(18,1:17)
#scaffold_1470 scaffold_174
sire_scaff_order_DF_4_plot1 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_1470")
plot(sire_scaff_order_DF_4_plot1$align_pos,sire_scaff_order_DF_4_plot1$chr_pos)
sire_scaff_order_DF_4_plot1 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_1470")
plot(sire_scaff_order_DF_4_plot1$align_pos,sire_scaff_order_DF_4_plot1$chr_pos)

sire_scaff_order_DF_4_plot2 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_174")
plot(sire_scaff_order_DF_4_plot2$align_pos,sire_scaff_order_DF_4_plot2$chr_pos)
sire_scaff_order_DF_4_plot2 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_174")
plot(sire_scaff_order_DF_4_plot2$align_pos,sire_scaff_order_DF_4_plot2$chr_pos)

sire_scaff_order_DF_5 <- sire_scaff_order_DF %>% filter(Ref_name == 5) %>% arrange(desc(scaff_len)) %>% 
  select(18,1:17)
#
sire_scaff_order_DF_5_plot1 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "")
plot(sire_scaff_order_DF_5_plot1$align_pos,sire_scaff_order_DF_5_plot1$chr_pos)
sire_scaff_order_DF_5_plot1 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "")
plot(sire_scaff_order_DF_5_plot1$align_pos,sire_scaff_order_DF_5_plot1$chr_pos)
#note: chr5 contigs got scaffolded together in scaffold_520

sire_scaff_order_DF_6 <- sire_scaff_order_DF %>% filter(Ref_name == 6) %>% arrange(desc(scaff_len)) %>% 
  select(18,1:17)
#scaffold_1426 scaffold_1141 scaffold_350 scaffold_1165
sire_scaff_order_DF_6_plot1 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_1426")
plot(sire_scaff_order_DF_6_plot1$align_pos,sire_scaff_order_DF_6_plot1$chr_pos)
sire_scaff_order_DF_6_plot1 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_1426")
plot(sire_scaff_order_DF_6_plot1$align_pos,sire_scaff_order_DF_6_plot1$chr_pos)

sire_scaff_order_DF_6_plot2 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_1141")
plot(sire_scaff_order_DF_6_plot2$align_pos,sire_scaff_order_DF_6_plot2$chr_pos)
sire_scaff_order_DF_6_plot2 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_1141")
plot(sire_scaff_order_DF_6_plot2$align_pos,sire_scaff_order_DF_6_plot2$chr_pos)

sire_scaff_order_DF_6_plot3 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_350")
plot(sire_scaff_order_DF_6_plot3$align_pos,sire_scaff_order_DF_6_plot3$chr_pos)
sire_scaff_order_DF_6_plot3 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_350")
plot(sire_scaff_order_DF_6_plot3$align_pos,sire_scaff_order_DF_6_plot3$chr_pos)

sire_scaff_order_DF_6_plot4 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_1165")
plot(sire_scaff_order_DF_6_plot4$align_pos,sire_scaff_order_DF_6_plot4$chr_pos)
sire_scaff_order_DF_6_plot4 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_1165")
plot(sire_scaff_order_DF_6_plot4$align_pos,sire_scaff_order_DF_6_plot4$chr_pos)

sire_scaff_order_DF_7 <- sire_scaff_order_DF %>% filter(Ref_name == 7) %>% arrange(desc(scaff_len)) %>% 
  select(18,1:17)
#scaffold_906
sire_scaff_order_DF_7_plot1 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_906")
plot(sire_scaff_order_DF_7_plot1$align_pos,sire_scaff_order_DF_7_plot1$chr_pos)
sire_scaff_order_DF_7_plot1 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_906")
plot(sire_scaff_order_DF_7_plot1$align_pos,sire_scaff_order_DF_7_plot1$chr_pos)
abline(v=c(9605195,12412605,14386151,21946165,67018291,80238494,91446336,97940944,108621554))

sire_scaff_order_DF_8 <- sire_scaff_order_DF %>% filter(Ref_name == 8) %>% arrange(desc(scaff_len)) %>% 
  select(18,1:17)
#scaffold_1314 scaffold_783 scaffold_1039 partial
sire_scaff_order_DF_8_plot1 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_1314")
plot(sire_scaff_order_DF_8_plot1$align_pos,sire_scaff_order_DF_8_plot1$chr_pos)
sire_scaff_order_DF_8_plot1 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_1314")
plot(sire_scaff_order_DF_8_plot1$align_pos,sire_scaff_order_DF_8_plot1$chr_pos)

sire_scaff_order_DF_8_plot2 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_783")
plot(sire_scaff_order_DF_8_plot2$align_pos,sire_scaff_order_DF_8_plot2$chr_pos)
sire_scaff_order_DF_8_plot2 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_783")
plot(sire_scaff_order_DF_8_plot2$align_pos,sire_scaff_order_DF_8_plot2$chr_pos)
#note: scaffold_1039 mapped to chr 12 but has part of contig_81 that belongs to chr8

sire_scaff_order_DF_9 <- sire_scaff_order_DF %>% filter(Ref_name == 9) %>% arrange(desc(scaff_len)) %>% 
  select(18,1:17)
#scaffold_1032
sire_scaff_order_DF_9_plot1 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_1032")
plot(sire_scaff_order_DF_9_plot1$align_pos,sire_scaff_order_DF_9_plot1$chr_pos)
sire_scaff_order_DF_9_plot1 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_1032")
plot(sire_scaff_order_DF_9_plot1$align_pos,sire_scaff_order_DF_9_plot1$chr_pos)
abline(v=c(5213419,7691168,9956846,10301955,13194950,39498261,88562625,121061233,142725547))

sire_scaff_order_DF_10 <- sire_scaff_order_DF %>% filter(Ref_name == 10) %>% arrange(desc(scaff_len)) %>% 
  select(18,1:17)
#scaffold_1045 scaffold_1202 scaffold_740
sire_scaff_order_DF_10_plot1 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_1045")
plot(sire_scaff_order_DF_10_plot1$align_pos,sire_scaff_order_DF_10_plot1$chr_pos)
sire_scaff_order_DF_10_plot1 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_1045")
plot(sire_scaff_order_DF_10_plot1$align_pos,sire_scaff_order_DF_10_plot1$chr_pos)

sire_scaff_order_DF_10_plot2 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_1202")
plot(sire_scaff_order_DF_10_plot2$align_pos,sire_scaff_order_DF_10_plot2$chr_pos)
sire_scaff_order_DF_10_plot2 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_1202")
plot(sire_scaff_order_DF_10_plot2$align_pos,sire_scaff_order_DF_10_plot2$chr_pos)

sire_scaff_order_DF_10_plot3 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_740")
plot(sire_scaff_order_DF_10_plot3$align_pos,sire_scaff_order_DF_10_plot3$chr_pos)
sire_scaff_order_DF_10_plot3 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_740")
plot(sire_scaff_order_DF_10_plot3$align_pos,sire_scaff_order_DF_10_plot3$chr_pos)

sire_scaff_order_DF_11 <- sire_scaff_order_DF %>% filter(Ref_name == 11) %>% arrange(desc(scaff_len)) %>% 
  select(18,1:17)
#scaffold_930
sire_scaff_order_DF_11_plot1 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_930")
plot(sire_scaff_order_DF_11_plot1$align_pos,sire_scaff_order_DF_11_plot1$chr_pos)
sire_scaff_order_DF_11_plot1 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_930")
plot(sire_scaff_order_DF_11_plot1$align_pos,sire_scaff_order_DF_11_plot1$chr_pos)

sire_scaff_order_DF_12 <- sire_scaff_order_DF %>% filter(Ref_name == 12) %>% arrange(desc(scaff_len)) %>% 
  select(18,1:17)
#scaffold_1039 
sire_scaff_order_DF_12_plot1 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "")
plot(sire_scaff_order_DF_12_plot1$align_pos,sire_scaff_order_DF_12_plot1$chr_pos)
sire_scaff_order_DF_12_plot1 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "")
plot(sire_scaff_order_DF_12_plot1$align_pos,sire_scaff_order_DF_12_plot1$chr_pos)

sire_scaff_order_DF_12_plot2 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "")
plot(sire_scaff_order_DF_12_plot2$align_pos,sire_scaff_order_DF_12_plot2$chr_pos)
sire_scaff_order_DF_12_plot2 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "")
plot(sire_scaff_order_DF_12_plot2$align_pos,sire_scaff_order_DF_12_plot2$chr_pos)

sire_scaff_order_DF_12_plot3 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "")
plot(sire_scaff_order_DF_12_plot3$align_pos,sire_scaff_order_DF_12_plot3$chr_pos)
sire_scaff_order_DF_12_plot3 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "")
plot(sire_scaff_order_DF_12_plot3$align_pos,sire_scaff_order_DF_12_plot3$chr_pos)

sire_scaff_order_DF_12_plot4 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "")
plot(sire_scaff_order_DF_12_plot4$align_pos,sire_scaff_order_DF_12_plot4$chr_pos)
sire_scaff_order_DF_12_plot4 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "")
plot(sire_scaff_order_DF_12_plot4$align_pos,sire_scaff_order_DF_12_plot4$chr_pos)

sire_scaff_order_DF_13 <- sire_scaff_order_DF %>% filter(Ref_name == 13) %>% arrange(desc(scaff_len)) %>% 
  select(18,1:17)
#scaffold_920
sire_scaff_order_DF_13_plot1 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_920")
plot(sire_scaff_order_DF_13_plot1$align_pos,sire_scaff_order_DF_13_plot1$chr_pos)
sire_scaff_order_DF_13_plot1 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_920")
plot(sire_scaff_order_DF_13_plot1$align_pos,sire_scaff_order_DF_13_plot1$chr_pos)

sire_scaff_order_DF_14 <- sire_scaff_order_DF %>% filter(Ref_name == 14) %>% arrange(desc(scaff_len)) %>% 
  select(18,1:17)
#scaffold_1267
sire_scaff_order_DF_14_plot1 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_1267")
plot(sire_scaff_order_DF_14_plot1$align_pos,sire_scaff_order_DF_14_plot1$chr_pos)
sire_scaff_order_DF_14_plot1 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_1267")
plot(sire_scaff_order_DF_14_plot1$align_pos,sire_scaff_order_DF_14_plot1$chr_pos)

sire_scaff_order_DF_15 <- sire_scaff_order_DF %>% filter(Ref_name == 15) %>% arrange(desc(scaff_len)) %>% 
  select(18,1:17)
#scaffold_1277 scaffold_664
sire_scaff_order_DF_15_plot1 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_1277")
plot(sire_scaff_order_DF_15_plot1$align_pos,sire_scaff_order_DF_15_plot1$chr_pos)
sire_scaff_order_DF_15_plot1 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_1277")
plot(sire_scaff_order_DF_15_plot1$align_pos,sire_scaff_order_DF_15_plot1$chr_pos)

sire_scaff_order_DF_15_plot2 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_664")
plot(sire_scaff_order_DF_15_plot2$align_pos,sire_scaff_order_DF_15_plot2$chr_pos)
sire_scaff_order_DF_15_plot2 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_664")
plot(sire_scaff_order_DF_15_plot2$align_pos,sire_scaff_order_DF_15_plot2$chr_pos)

sire_scaff_order_DF_16 <- sire_scaff_order_DF %>% filter(Ref_name == 16) %>% arrange(desc(scaff_len)) %>% 
  select(18,1:17)
#scaffold_451 scaffold_694
sire_scaff_order_DF_16_plot1 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_451")
plot(sire_scaff_order_DF_16_plot1$align_pos,sire_scaff_order_DF_16_plot1$chr_pos)
sire_scaff_order_DF_16_plot1 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_451")
plot(sire_scaff_order_DF_16_plot1$align_pos,sire_scaff_order_DF_16_plot1$chr_pos)

sire_scaff_order_DF_16_plot2 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_694")
plot(sire_scaff_order_DF_16_plot2$align_pos,sire_scaff_order_DF_16_plot2$chr_pos)
sire_scaff_order_DF_16_plot2 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_694")
plot(sire_scaff_order_DF_16_plot2$align_pos,sire_scaff_order_DF_16_plot2$chr_pos)

sire_scaff_order_DF_17 <- sire_scaff_order_DF %>% filter(Ref_name == 17) %>% arrange(desc(scaff_len)) %>% 
  select(18,1:17)
#scaffold_1161 scaffold_320
sire_scaff_order_DF_17_plot1 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_1161")
plot(sire_scaff_order_DF_17_plot1$align_pos,sire_scaff_order_DF_17_plot1$chr_pos)
sire_scaff_order_DF_17_plot1 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_1161")
plot(sire_scaff_order_DF_17_plot1$align_pos,sire_scaff_order_DF_17_plot1$chr_pos)

sire_scaff_order_DF_17_plot2 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_320")
plot(sire_scaff_order_DF_17_plot2$align_pos,sire_scaff_order_DF_17_plot2$chr_pos)
sire_scaff_order_DF_17_plot2 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_320")
plot(sire_scaff_order_DF_17_plot2$align_pos,sire_scaff_order_DF_17_plot2$chr_pos)

sire_scaff_order_DF_18 <- sire_scaff_order_DF %>% filter(Ref_name == 18) %>% arrange(desc(scaff_len)) %>% 
  select(18,1:17)
#scaffold_1215 scaffold_1383 scaffold_1164
sire_scaff_order_DF_18_plot1 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_1215")
plot(sire_scaff_order_DF_18_plot1$align_pos,sire_scaff_order_DF_18_plot1$chr_pos)
sire_scaff_order_DF_18_plot1 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_1215")
plot(sire_scaff_order_DF_18_plot1$align_pos,sire_scaff_order_DF_18_plot1$chr_pos)

sire_scaff_order_DF_18_plot2 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_1383")
plot(sire_scaff_order_DF_18_plot2$align_pos,sire_scaff_order_DF_18_plot2$chr_pos)
sire_scaff_order_DF_18_plot2 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_1383")
plot(sire_scaff_order_DF_18_plot2$align_pos,sire_scaff_order_DF_18_plot2$chr_pos)

sire_scaff_order_DF_18_plot3 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_1164")
plot(sire_scaff_order_DF_18_plot3$align_pos,sire_scaff_order_DF_18_plot3$chr_pos)
sire_scaff_order_DF_18_plot3 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_1164")
plot(sire_scaff_order_DF_18_plot3$align_pos,sire_scaff_order_DF_18_plot3$chr_pos)

sire_scaff_order_DF_19 <- sire_scaff_order_DF %>% filter(Ref_name == 19) %>% arrange(desc(scaff_len)) %>% 
  select(18,1:17)
#scaffold_1089
sire_scaff_order_DF_19_plot1 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_1089")
plot(sire_scaff_order_DF_19_plot1$align_pos,sire_scaff_order_DF_19_plot1$chr_pos)
sire_scaff_order_DF_19_plot1 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_1089")
plot(sire_scaff_order_DF_19_plot1$align_pos,sire_scaff_order_DF_19_plot1$chr_pos)

sire_scaff_order_DF_20 <- sire_scaff_order_DF %>% filter(Ref_name == 20) %>% arrange(desc(scaff_len)) %>% 
  select(18,1:17)
#scaffold_68
sire_scaff_order_DF_20_plot1 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_68")
plot(sire_scaff_order_DF_20_plot1$align_pos,sire_scaff_order_DF_20_plot1$chr_pos)
sire_scaff_order_DF_20_plot1 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_68")
plot(sire_scaff_order_DF_20_plot1$align_pos,sire_scaff_order_DF_20_plot1$chr_pos)
#note: potentially can join scaffold_856 in between contigs of scaffold_68

sire_scaff_order_DF_21 <- sire_scaff_order_DF %>% filter(Ref_name == 21) %>% arrange(desc(scaff_len)) %>% 
  select(18,1:17)
#scaffold_644
sire_scaff_order_DF_21_plot1 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_644")
plot(sire_scaff_order_DF_21_plot1$align_pos,sire_scaff_order_DF_21_plot1$chr_pos)
sire_scaff_order_DF_21_plot1 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_644")
plot(sire_scaff_order_DF_21_plot1$align_pos,sire_scaff_order_DF_21_plot1$chr_pos)

sire_scaff_order_DF_22 <- sire_scaff_order_DF %>% filter(Ref_name == 22) %>% arrange(desc(scaff_len)) %>% 
  select(18,1:17)
#scaffold_667 scaffold_1285
sire_scaff_order_DF_22_plot1 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_667")
plot(sire_scaff_order_DF_22_plot1$align_pos,sire_scaff_order_DF_22_plot1$chr_pos)
sire_scaff_order_DF_22_plot1 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_667")
plot(sire_scaff_order_DF_22_plot1$align_pos,sire_scaff_order_DF_22_plot1$chr_pos)

sire_scaff_order_DF_22_plot2 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_1285")
plot(sire_scaff_order_DF_22_plot2$align_pos,sire_scaff_order_DF_22_plot2$chr_pos)
sire_scaff_order_DF_22_plot2 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_1285")
plot(sire_scaff_order_DF_22_plot2$align_pos,sire_scaff_order_DF_22_plot2$chr_pos)

sire_scaff_order_DF_23 <- sire_scaff_order_DF %>% filter(Ref_name == 23) %>% arrange(desc(scaff_len)) %>% 
  select(18,1:17)
#scaffold_381 scaffold_20 scaffold_68 partial
sire_scaff_order_DF_23_plot1 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_381")
plot(sire_scaff_order_DF_23_plot1$align_pos,sire_scaff_order_DF_23_plot1$chr_pos)
sire_scaff_order_DF_23_plot1 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_381")
plot(sire_scaff_order_DF_23_plot1$align_pos,sire_scaff_order_DF_23_plot1$chr_pos)

sire_scaff_order_DF_23_plot2 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_20")
plot(sire_scaff_order_DF_23_plot2$align_pos,sire_scaff_order_DF_23_plot2$chr_pos)
sire_scaff_order_DF_23_plot2 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_20")
plot(sire_scaff_order_DF_23_plot2$align_pos,sire_scaff_order_DF_23_plot2$chr_pos)

sire_scaff_order_DF_24 <- sire_scaff_order_DF %>% filter(Ref_name == 24) %>% arrange(desc(scaff_len)) %>% 
  select(18,1:17)
#scaffold_1397
sire_scaff_order_DF_24_plot1 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_1397")
plot(sire_scaff_order_DF_24_plot1$align_pos,sire_scaff_order_DF_24_plot1$chr_pos)
sire_scaff_order_DF_24_plot1 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_1397")
plot(sire_scaff_order_DF_24_plot1$align_pos,sire_scaff_order_DF_24_plot1$chr_pos)

sire_scaff_order_DF_25 <- sire_scaff_order_DF %>% filter(Ref_name == 25) %>% arrange(desc(scaff_len)) %>% 
  select(18,1:17)
#scaffold_476
sire_scaff_order_DF_25_plot1 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_476")
plot(sire_scaff_order_DF_25_plot1$align_pos,sire_scaff_order_DF_25_plot1$chr_pos)
sire_scaff_order_DF_25_plot1 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_476")
plot(sire_scaff_order_DF_25_plot1$align_pos,sire_scaff_order_DF_25_plot1$chr_pos)

sire_scaff_order_DF_26 <- sire_scaff_order_DF %>% filter(Ref_name == 26) %>% arrange(desc(scaff_len)) %>% 
  select(18,1:17)
#scaffold_1091
sire_scaff_order_DF_26_plot1 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_1091")
plot(sire_scaff_order_DF_26_plot1$align_pos,sire_scaff_order_DF_26_plot1$chr_pos)
sire_scaff_order_DF_26_plot1 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_1091")
plot(sire_scaff_order_DF_26_plot1$align_pos,sire_scaff_order_DF_26_plot1$chr_pos)

sire_scaff_order_DF_27 <- sire_scaff_order_DF %>% filter(Ref_name == 27) %>% arrange(desc(scaff_len)) %>% 
  select(18,1:17)
#scaffold_1128 
sire_scaff_order_DF_27_plot1 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_1128")
plot(sire_scaff_order_DF_27_plot1$align_pos,sire_scaff_order_DF_27_plot1$chr_pos)
sire_scaff_order_DF_27_plot1 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_1128")
plot(sire_scaff_order_DF_27_plot1$align_pos,sire_scaff_order_DF_27_plot1$chr_pos)

sire_scaff_order_DF_28 <- sire_scaff_order_DF %>% filter(Ref_name == 28) %>% arrange(desc(scaff_len)) %>% 
  select(18,1:17)
#scaffold_410 scaffold_249 scaffold_1032 partial
sire_scaff_order_DF_28_plot1 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_410")
plot(sire_scaff_order_DF_28_plot1$align_pos,sire_scaff_order_DF_28_plot1$chr_pos)
sire_scaff_order_DF_28_plot1 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_410")
plot(sire_scaff_order_DF_28_plot1$align_pos,sire_scaff_order_DF_28_plot1$chr_pos)

sire_scaff_order_DF_28_plot2 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_249")
plot(sire_scaff_order_DF_28_plot2$align_pos,sire_scaff_order_DF_28_plot2$chr_pos)
sire_scaff_order_DF_28_plot2 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_249")
plot(sire_scaff_order_DF_28_plot2$align_pos,sire_scaff_order_DF_28_plot2$chr_pos)

sire_scaff_order_DF_29 <- sire_scaff_order_DF %>% filter(Ref_name == 29) %>% arrange(desc(scaff_len)) %>% 
  select(18,1:17)
#scaffold_190 scaffold_790
sire_scaff_order_DF_29_plot1 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_190")
plot(sire_scaff_order_DF_29_plot1$align_pos,sire_scaff_order_DF_29_plot1$chr_pos)
sire_scaff_order_DF_29_plot1 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_190")
plot(sire_scaff_order_DF_29_plot1$align_pos,sire_scaff_order_DF_29_plot1$chr_pos)

sire_scaff_order_DF_29_plot2 <- f1_sire_salsa_HDProbes_tab %>% filter(scaffold == "scaffold_790")
plot(sire_scaff_order_DF_29_plot2$align_pos,sire_scaff_order_DF_29_plot2$chr_pos)
sire_scaff_order_DF_29_plot2 <- f1_sire_salsa_fasta_rcmap_tab %>% filter(scaffold == "scaffold_790")
plot(sire_scaff_order_DF_29_plot2$align_pos,sire_scaff_order_DF_29_plot2$chr_pos)

