#------------------------------------------------------
# Program name: brahman_angus_mashmap.R
# Objective: check mashmap of cleaned assembly contigs 
#           mapping to latest UCDv25
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(readr)
library(dplyr)

load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/brahman_Angus_mashmap/dam_agp_clean_assembly_to_salsa.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/brahman_Angus_mashmap/sire_agp_clean_assembly_to_salsa.RData")

# reading dam_assembly_cleaned_mash_UCDv25.mashmap
dam_mashmap <- read_delim("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/brahman_Angus_mashmap/script/dam_assembly_cleaned_mash_UCDv25.mashmap",
                          " ", col_names = FALSE)
names(dam_mashmap) <- c("Query_name","Query_length","Query_start","Query_end","Orientation",
                        "Ref_name","Ref_length","Ref_start","Ref_end","Percentid")

dam_mashmap_modi <- dam_mashmap %>% arrange(Ref_name,Ref_start) %>% mutate(alignedL = Ref_end - Ref_start +1) %>%
  mutate(alignProportion = alignedL/Query_length)

length(unique(dam_mashmap_modi$Query_name))

dam_mashmap_modi_chr1_to_X <- dam_mashmap_modi %>% filter(Ref_name == "1" | Ref_name == "2" | Ref_name == "3" | 
                                                            Ref_name == "4" | Ref_name == "5" | Ref_name == "6" | 
                                                            Ref_name == "7" | Ref_name == "8" | Ref_name == "9" |
                                                            Ref_name == "10" | Ref_name == "11" | Ref_name == "12" |
                                                            Ref_name == "13" | Ref_name == "14" | Ref_name == "15" |
                                                            Ref_name == "16" | Ref_name == "17" | Ref_name == "18" |
                                                            Ref_name == "19" | Ref_name == "20" | Ref_name == "21" |
                                                            Ref_name == "22" | Ref_name == "23" | Ref_name == "24" |
                                                            Ref_name == "25" | Ref_name == "26" | Ref_name == "27" |
                                                            Ref_name == "28" | Ref_name == "29" | Ref_name == "X")

#dam_mashmap_modi_length_notmatched <- dam_mashmap_modi[!dam_mashmap_modi$alignedL == dam_mashmap_modi$Query_length,]

object <- c()
for (i in 1:nrow(dam_mashmap_modi_chr1_to_X)){
 logic_selc <- dam_scaffolds_FINAL_agp_modi$component_id %in% dam_mashmap_modi_chr1_to_X$Query_name[i]
 object <- c(object,dam_scaffolds_FINAL_agp_modi$object[logic_selc])
}
dam_mashmap_modi_chr1_to_X$object <- object

#group scaffolds mapped by mashmap
dam_scaffold_count_mapped <- dam_mashmap_modi_chr1_to_X %>% group_by(Ref_name,object) %>% summarise(count = n())

dam_scaffold_count_afterSalsa <- dam_scaffolds_FINAL_agp_modi %>% group_by(object) %>% summarise(count = n())

#loop thro mapped count and check tally with straight from salsa
count_Salsa <- c()
for (i in 1:nrow(dam_scaffold_count_mapped)){
  logic_selc2 <- dam_scaffold_count_afterSalsa$object %in% dam_scaffold_count_mapped$object[i]
  count_Salsa <- c(count_Salsa,dam_scaffold_count_afterSalsa$count[logic_selc2])
}
dam_scaffold_count_mapped$count_Salsa <- count_Salsa

#scaffolded but not mapped to the same place
dam_scaffold_count_mapped %>% filter(count != count_Salsa)

#looking for totally missing dam contigs and their lengths
logic_selc4 <- dam_scaffolds_FINAL_agp_modi$component_id %in% dam_mashmap_modi$Query_name
dam_scaffolds_missing_contigs <- dam_scaffolds_FINAL_agp_modi[!logic_selc4,]
dam_scaffolds_missing_contigs <- dam_scaffolds_missing_contigs %>% arrange(slength)

dam_scaffolds_missing_contigs %>% select(slength,object) %>% group_by(object) %>% 
  summarise(count = n(),total_length = sum(slength)) %>% arrange(desc(total_length))

sum(dam_scaffolds_missing_contigs$slength)/1e6
#0.491877 ~ 0.5Mb

#############sire######################################################
# reading sire_assembly_cleaned_mash_UCDv25.mashmap
sire_mashmap <- read_delim("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/brahman_Angus_mashmap/script/sire_assembly_cleaned_mash_UCDv25.mashmap",
                          " ", col_names = FALSE)
names(sire_mashmap) <- c("Query_name","Query_length","Query_start","Query_end","Orientation",
                        "Ref_name","Ref_length","Ref_start","Ref_end","Percentid")

sire_mashmap_modi <- sire_mashmap %>% arrange(Ref_name,Ref_start) %>% mutate(alignedL = Ref_end - Ref_start +1) %>%
  mutate(alignProportion = alignedL/Query_length)

length(unique(sire_mashmap_modi$Query_name))

sire_mashmap_modi_chr1_to_X <- sire_mashmap_modi %>% filter(Ref_name == "1" | Ref_name == "2" | Ref_name == "3" | 
                                                            Ref_name == "4" | Ref_name == "5" | Ref_name == "6" | 
                                                            Ref_name == "7" | Ref_name == "8" | Ref_name == "9" |
                                                            Ref_name == "10" | Ref_name == "11" | Ref_name == "12" |
                                                            Ref_name == "13" | Ref_name == "14" | Ref_name == "15" |
                                                            Ref_name == "16" | Ref_name == "17" | Ref_name == "18" |
                                                            Ref_name == "19" | Ref_name == "20" | Ref_name == "21" |
                                                            Ref_name == "22" | Ref_name == "23" | Ref_name == "24" |
                                                            Ref_name == "25" | Ref_name == "26" | Ref_name == "27" |
                                                            Ref_name == "28" | Ref_name == "29" | Ref_name == "X")

#sire_mashmap_modi_length_notmatched <- sire_mashmap_modi[!sire_mashmap_modi$alignedL == sire_mashmap_modi$Query_length,]

object <- c()
for (i in 1:nrow(sire_mashmap_modi_chr1_to_X)){
  logic_selc <- sire_scaffolds_FINAL_agp_modi$component_id %in% sire_mashmap_modi_chr1_to_X$Query_name[i]
  object <- c(object,sire_scaffolds_FINAL_agp_modi$object[logic_selc])
}
sire_mashmap_modi_chr1_to_X$object <- object

#group scaffolds mapped by mashmap
sire_scaffold_count_mapped <- sire_mashmap_modi_chr1_to_X %>% group_by(Ref_name,object) %>% summarise(count = n())

sire_scaffold_count_afterSalsa <- sire_scaffolds_FINAL_agp_modi %>% group_by(object) %>% summarise(count = n())

#loop thro mapped count and check tally with straight from salsa
count_Salsa <- c()
for (i in 1:nrow(sire_scaffold_count_mapped)){
  logic_selc2 <- sire_scaffold_count_afterSalsa$object %in% sire_scaffold_count_mapped$object[i]
  count_Salsa <- c(count_Salsa,sire_scaffold_count_afterSalsa$count[logic_selc2])
}
sire_scaffold_count_mapped$count_Salsa <- count_Salsa

#scaffolded but not mapped to the same place
sire_scaffold_count_mapped %>% filter(count != count_Salsa)

#looking for totally missing sire contigs and their lengths
logic_selc3 <- sire_scaffolds_FINAL_agp_modi$component_id %in% sire_mashmap_modi$Query_name
sire_scaffolds_missing_contigs <- sire_scaffolds_FINAL_agp_modi[!logic_selc3,]
sire_scaffolds_missing_contigs <- sire_scaffolds_missing_contigs %>% arrange(slength)

sire_scaffolds_missing_contigs %>% select(slength,object) %>% group_by(object) %>% 
  summarise(count = n(),total_length = sum(slength)) %>% arrange(desc(total_length))

sum(sire_scaffolds_missing_contigs$slength)/1e6
#3.7 Mb

#any missing contig in sire that maps to human chr X or Y
sire_scaffolds_missing_contigs$component_id %in% c("contig_445","contig_1544","contig_752","contig_751",
                                      "contig_1060","contig_514","contig_1076","contig_1075",
                                      "contig_331","contig_1080","contig_656","contig_808",
                                      "contig_1216","contig_1358","contig_1056","contig_194")

#any missing contig in sire that maps to pig chr X or Y
sire_scaffolds_missing_contigs$component_id %in% c("contig_1163","contig_373","contig_981",
                                                   "contig_977","contig_1743","contig_751",
                                                   "contig_1060","contig_1076","contig_331",
                                                   "contig_1080","contig_1248","contig_1095",
                                                   "contig_110","contig_808","contig_802",
                                                   "contig_874","contig_1216","contig_1358",
                                                   "contig_1591","contig_194")

#any missing contig in sire that maps to CM001061


#for cynthia
save(dam_mashmap_modi_chr1_to_X,sire_mashmap_modi_chr1_to_X,dam_mashmap_modi,sire_mashmap_modi,
     file="/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/brahman_Angus_mashmap/contig_scaff_UCDv25.RData")
