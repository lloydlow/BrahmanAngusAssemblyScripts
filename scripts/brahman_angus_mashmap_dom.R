#------------------------------------------------------
# Program name: brahman_angus_mashmap_dom.R
# Objective: checking discrepancy in contig mapped to
#           Dominette starting with contig_232 brahman
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(readr)
library(dplyr)

# reading dam_contig_232_vs_cattleUCD.mashmap
dam_contig_232_vs_cattleUCD_mashmap <- read_delim("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/mashmap_vs_dominette/dam_contig_232_vs_cattleUCD.mashmap",
                          " ", col_names = FALSE)
names(dam_contig_232_vs_cattleUCD_mashmap) <- c("Query_name","Query_length","Query_start","Query_end","Orientation",
                        "Ref_name","Ref_length","Ref_start","Ref_end","Percentid")

dam_contig_232_vs_cattleUCD_mashmap_modi <- dam_contig_232_vs_cattleUCD_mashmap %>% arrange(Ref_name,Ref_start) %>% 
  mutate(alignedL = Ref_end - Ref_start +1) %>%
  mutate(alignProportion = alignedL/Query_length) %>% arrange(Query_start) %>% filter(Percentid > 98)

#length(unique(dam_mashmap_modi$Query_name))

# reading sire_contig_81_vs_cattleUCD.mashmap
sire_contig_81_vs_cattleUCD_mashmap <- read_delim("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/mashmap_vs_dominette/sire_contig_81_vs_cattleUCD.mashmap",
                                                  " ", col_names = FALSE)
names(sire_contig_81_vs_cattleUCD_mashmap) <- c("Query_name","Query_length","Query_start","Query_end","Orientation",
                                                "Ref_name","Ref_length","Ref_start","Ref_end","Percentid")

sire_contig_81_vs_cattleUCD_mashmap_modi <- sire_contig_81_vs_cattleUCD_mashmap %>% arrange(Ref_name,Ref_start) %>% 
  mutate(alignedL = Ref_end - Ref_start +1) %>%
  mutate(alignProportion = alignedL/Query_length) %>% arrange(Query_start) %>% filter(Percentid > 98)

# reading sire_contig_802_vs_cattleUCD.mashmap
sire_contig_802_vs_cattleUCD_mashmap <- read_delim("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/mashmap_vs_dominette/sire_contig_802_vs_cattleUCD.mashmap",
                                                  " ", col_names = FALSE)
names(sire_contig_802_vs_cattleUCD_mashmap) <- c("Query_name","Query_length","Query_start","Query_end","Orientation",
                                                "Ref_name","Ref_length","Ref_start","Ref_end","Percentid")

sire_contig_802_vs_cattleUCD_mashmap_modi <- sire_contig_802_vs_cattleUCD_mashmap %>% arrange(Ref_name,Ref_start) %>% 
  mutate(alignedL = Ref_end - Ref_start +1) %>%
  mutate(alignProportion = alignedL/Query_length) %>% arrange(Query_start) %>% filter(Percentid > 98)

# reading sire_contig_1722_vs_cattleUCD.mashmap
sire_contig_1722_vs_cattleUCD_mashmap <- read_delim("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/mashmap_vs_dominette/sire_contig_1722_vs_cattleUCD.mashmap",
                                                   " ", col_names = FALSE)
names(sire_contig_1722_vs_cattleUCD_mashmap) <- c("Query_name","Query_length","Query_start","Query_end","Orientation",
                                                 "Ref_name","Ref_length","Ref_start","Ref_end","Percentid")

sire_contig_1722_vs_cattleUCD_mashmap_modi <- sire_contig_1722_vs_cattleUCD_mashmap %>% arrange(Ref_name,Ref_start) %>% 
  mutate(alignedL = Ref_end - Ref_start +1) %>%
  mutate(alignProportion = alignedL/Query_length) %>% arrange(Query_start) %>% filter(Percentid > 98)
