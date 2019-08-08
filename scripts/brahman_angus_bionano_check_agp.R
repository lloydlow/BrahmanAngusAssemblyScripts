#------------------------------------------------------
# Program name: brahman_angus_bionano_check_agp.R
# Objective: check agp file from the start of 
#           scaffolding of the various versions bionano
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(readr)
library(dplyr)

# reading brahman-selected bionano scaffolds agp
dir1 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/agp/"
path1 <- paste0(dir1,"EXP_REFINEFINAL1_bppAdjust_cmap_bostaurus_brahma_fasta_NGScontigs_HYBRID_SCAFFOLD.agp")

dam_scaffolds_FINAL_agp <- read_tsv(path1,col_names = FALSE, comment = "#")


dam_scaffolds_FINAL_agp_modi <- dam_scaffolds_FINAL_agp %>% select(X1:X9)
names(dam_scaffolds_FINAL_agp_modi) <- c("object","object_beg","object_end","part_number","component_type",
                                         "component_id","component_beg","component_end","orientation")

#no gaps in df
dam_scaffolds_FINAL_agp_modi <- dam_scaffolds_FINAL_agp_modi %>% filter(component_type == "W")

#changing to numeric for input start and stop
dam_scaffolds_FINAL_agp_modi$component_end <- as.numeric(dam_scaffolds_FINAL_agp_modi$component_end)
dam_scaffolds_FINAL_agp_modi$component_beg <- as.numeric(dam_scaffolds_FINAL_agp_modi$component_beg)

#make length variables calc for final scaffold and begin contig input
dam_scaffolds_FINAL_agp_modi$slength <- dam_scaffolds_FINAL_agp_modi$object_end - 
  dam_scaffolds_FINAL_agp_modi$object_beg +1

dam_scaffolds_FINAL_agp_modi$clength <- dam_scaffolds_FINAL_agp_modi$component_end - 
  dam_scaffolds_FINAL_agp_modi$component_beg +1

#slength matches clength row by row #total length in dam 2678769087

#try to summarise count by object and sum of length
dam_scaffolds_FINAL_agp_modi_1st_summ <- dam_scaffolds_FINAL_agp_modi %>% 
  select(object,slength) %>% group_by(object) %>%
  summarise(count = n(), total_length = sum(slength)) %>% 
  arrange(desc(total_length))

dam_scaffolds_FINAL_agp_modi_1st_summ_multiContig <- dam_scaffolds_FINAL_agp_modi_1st_summ %>%
  filter(count != 1)

#68% bases in 84 scaffolds that was put together with more than 1 contig
sum(dam_scaffolds_FINAL_agp_modi_1st_summ_multiContig$total_length)/
  sum(dam_scaffolds_FINAL_agp_modi_1st_summ$total_length)

summary(dam_scaffolds_FINAL_agp_modi_1st_summ$total_length)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#60    26638    38053  1979874    53571 97845705 

#scaffold is NOT uniq
test_uniq_scaff <- unique(dam_scaffolds_FINAL_agp_modi$object)
length(test_uniq_scaff)

scaff_dup <- dam_scaffolds_FINAL_agp_modi$object[duplicated(dam_scaffolds_FINAL_agp_modi$object)]
scaff_dup_uniq <- unique(scaff_dup)

logic_dup <- dam_scaffolds_FINAL_agp_modi$object %in% scaff_dup_uniq
dam_scaffolds_FINAL_agp_modi_multiContig <- dam_scaffolds_FINAL_agp_modi[logic_dup,]
dam_scaffolds_FINAL_agp_modi_multiContig <- dam_scaffolds_FINAL_agp_modi_multiContig %>%
  arrange(object,object_beg)

#save(dam_scaffolds_FINAL_agp_modi,dam_scaffolds_FINAL_agp_modi_1st_summ,
#     dam_scaffolds_FINAL_agp_modi_1st_summ_multiContig,dam_scaffolds_FINAL_agp_modi_multiContig,
#     file = "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/brahman_Angus_mashmap/dam_agp_clean_assembly_to_salsa.RData")

#this is after 20180920 
dam_contig_scaff <- dam_scaffolds_FINAL_agp_modi %>% select(component_id,object,clength,object_beg,object_end,orientation)
save(dam_contig_scaff,
     file = "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/order_bionano_scaffolds_to_chr_level/dam_contig_scaff.RData")

############sire####################################################################
# reading angus-selected bionano scaffolds agp
dir2 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/agp/"
path2 <- paste0(dir2,"EXP_REFINEFINAL1_bppAdjust_cmap_bostaurus_angus_fasta_NGScontigs_HYBRID_SCAFFOLD.agp")

sire_scaffolds_FINAL_agp <- read_tsv(path2,col_names = FALSE, comment = "#")

sire_scaffolds_FINAL_agp_modi <- sire_scaffolds_FINAL_agp %>% select(X1:X9)
names(sire_scaffolds_FINAL_agp_modi) <- c("object","object_beg","object_end","part_number","component_type",
                                          "component_id","component_beg","component_end","orientation")

#no gaps in df
sire_scaffolds_FINAL_agp_modi <- sire_scaffolds_FINAL_agp_modi %>% filter(component_type == "W")

#changing to numeric for input start and stop
sire_scaffolds_FINAL_agp_modi$component_end <- as.numeric(sire_scaffolds_FINAL_agp_modi$component_end)
sire_scaffolds_FINAL_agp_modi$component_beg <- as.numeric(sire_scaffolds_FINAL_agp_modi$component_beg)

#make length variables calc for final scaffold and begin contig input
sire_scaffolds_FINAL_agp_modi$slength <- sire_scaffolds_FINAL_agp_modi$object_end - 
  sire_scaffolds_FINAL_agp_modi$object_beg +1

sire_scaffolds_FINAL_agp_modi$clength <- sire_scaffolds_FINAL_agp_modi$component_end - 
  sire_scaffolds_FINAL_agp_modi$component_beg +1

#slength matches clength row by row #total length in sire

#try to summarise count by object and sum of length
sire_scaffolds_FINAL_agp_modi_1st_summ <- sire_scaffolds_FINAL_agp_modi %>% 
  select(object,slength) %>% group_by(object) %>%
  summarise(count = n(), total_length = sum(slength)) %>% 
  arrange(desc(total_length))

sire_scaffolds_FINAL_agp_modi_1st_summ_multiContig <- sire_scaffolds_FINAL_agp_modi_1st_summ %>%
  filter(count != 1)

#97% bases in 82 scaffolds that was put together with more than 1 contig
sum(sire_scaffolds_FINAL_agp_modi_1st_summ_multiContig$total_length)/
  sum(sire_scaffolds_FINAL_agp_modi_1st_summ$total_length)

summary(sire_scaffolds_FINAL_agp_modi_1st_summ_multiContig$total_length)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#95828    273347   3074763  16159584  24755263 106258578 

#scaffold is NOT uniq
test_uniq_scaff <- unique(sire_scaffolds_FINAL_agp_modi$object)
length(test_uniq_scaff)

scaff_dup <- sire_scaffolds_FINAL_agp_modi$object[duplicated(sire_scaffolds_FINAL_agp_modi$object)]
scaff_dup_uniq <- unique(scaff_dup)

logic_dup <- sire_scaffolds_FINAL_agp_modi$object %in% scaff_dup_uniq
sire_scaffolds_FINAL_agp_modi_multiContig <- sire_scaffolds_FINAL_agp_modi[logic_dup,]
sire_scaffolds_FINAL_agp_modi_multiContig <- sire_scaffolds_FINAL_agp_modi_multiContig %>%
  arrange(object,object_beg)

# save(sire_scaffolds_FINAL_agp_modi,sire_scaffolds_FINAL_agp_modi_1st_summ,
#      sire_scaffolds_FINAL_agp_modi_1st_summ_multiContig,sire_scaffolds_FINAL_agp_modi_multiContig,
#      file = "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/brahman_Angus_mashmap/sire_agp_clean_assembly_to_salsa.RData")

#this is after 20180905 
sire_contig_scaff <- sire_scaffolds_FINAL_agp_modi %>% select(component_id,object,clength,object_beg,object_end,orientation)
save(sire_contig_scaff,
     file = "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/order_bionano_scaffolds_to_chr_level/sire_contig_scaff.RData")
