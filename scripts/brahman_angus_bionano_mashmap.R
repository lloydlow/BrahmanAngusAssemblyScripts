#------------------------------------------------------
# Program name: brahman_angus_bionano_mashmap.R
# Objective: check mashmap of CANU haplo resolved contigs 
#           to SALSA broken assembly.cleaned.fasta
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(readr)
library(dplyr)

# reading bostaurus_angus_vs_sire_cleaned_assembly.mashmap
bostaurus_angus_vs_sire_cleaned_assembly <- 
  read_delim("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/mashmap_bionano_stuff/bostaurus_angus_vs_sire_cleaned_assembly.mashmap",
             " ", col_names = FALSE)
names(bostaurus_angus_vs_sire_cleaned_assembly) <- c("Query_name","Query_length","Query_start","Query_end","Orientation",
                        "Ref_name","Ref_length","Ref_start","Ref_end","Percentid")

bostaurus_angus_vs_sire_cleaned_assembly_modi <- bostaurus_angus_vs_sire_cleaned_assembly %>% 
  arrange(Ref_name,Ref_start) %>% mutate(alignedL = Ref_end - Ref_start +1) %>%
  mutate(alignProportion = alignedL/Query_length)

#is *tigs uniq? Ans is No.
length(unique(bostaurus_angus_vs_sire_cleaned_assembly_modi$Query_name))
#1739
length(bostaurus_angus_vs_sire_cleaned_assembly_modi$Query_name)
#1827

#is contig uniq? Ans is No.
length(unique(bostaurus_angus_vs_sire_cleaned_assembly_modi$Ref_name))
#1603
length(bostaurus_angus_vs_sire_cleaned_assembly_modi$Ref_name)
#1827

#how many the contig lined up with tigs at 0.99 proportion align length
bostaurus_angus_vs_sire_cleaned_assembly_modi_correspond1 <- bostaurus_angus_vs_sire_cleaned_assembly_modi %>%
  filter(alignProportion > 0.98)

#1484 *tigs are mapped to contigs here
bostaurus_angus_vs_sire_cleaned_assembly_modi_correspond1_group <- bostaurus_angus_vs_sire_cleaned_assembly_modi_correspond1 %>% 
  group_by(Query_name) %>% summarise(count = n())

#analyse the remaining *tigs mapping
bostaurus_angus_vs_sire_cleaned_assembly_modi_correspond2 <- bostaurus_angus_vs_sire_cleaned_assembly_modi %>%
  filter(alignProportion <= 0.98)

##### sire #####
#this isnt done yet because the previous mashmap file didnt finish running
#and the purpose of this script is correspondence of *tigs to salsa contigs
#which is already solved by looking at sequence lengths
