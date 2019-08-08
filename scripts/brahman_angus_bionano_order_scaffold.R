#------------------------------------------------------
# Program name: brahman_angus_bionano_order_scaffold.R
# Objective: make a big table to see where scaffold belong
#           and order them accordingly
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(readr)
library(dplyr)

#table to build
#contig scaffold contig_len ungap_scaff_len chr_HD HD_n chr_rc rc_n chr_Dom Dom_align_len Dom_len

is.integer0 <- function(x){
  is.integer(x) && length(x) == 0L
}

# START scaffolds alignment to Dominette #
##### Dam #####
dam_mashmap <- read_delim("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/brahman_Angus_mashmap/script/bostaurus_brahma_bionano_NCBI_full_mash_UCDv25.mashmap",
                          " ", col_names = FALSE)
names(dam_mashmap) <- c("Query_name","Query_length","Query_start","Query_end","Orientation",
                        "Ref_name","Ref_length","Ref_start","Ref_end","Percentid")

dam_mashmap <- dam_mashmap %>% mutate(query_align = Query_end - Query_start, proportion_ref = query_align/Ref_length) %>%
  select(Query_name,Ref_name,Ref_length,Ref_start,Ref_end,query_align,proportion_ref)

##### Sire #####
sire_mashmap <- read_delim("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/brahman_Angus_mashmap/script/bostaurus_angus_bionano_NCBI_full_mash_UCDv25.mashmap",
                           " ", col_names = FALSE)
names(sire_mashmap) <- c("Query_name","Query_length","Query_start","Query_end","Orientation",
                         "Ref_name","Ref_length","Ref_start","Ref_end","Percentid")

sire_mashmap <- sire_mashmap %>% mutate(query_align = Query_end - Query_start, proportion_ref = query_align/Ref_length) %>%
  select(Query_name,Ref_name,Ref_length,Ref_start,Ref_end,query_align, proportion_ref)
# END scaffolds alignment to Dominette #

# START contig mapping to scaffold #
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/order_bionano_scaffolds_to_chr_level/dam_contig_scaff.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/order_bionano_scaffolds_to_chr_level/sire_contig_scaff.RData")
# END contig mapping to scaffold #

# START HD and rc probes #
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/order_bionano_scaffolds_to_chr_level/f1_dam_bionano_fasta_rcmap_tab_by_chrScaff_mostprobes.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/order_bionano_scaffolds_to_chr_level/f1_sire_bionano_fasta_rcmap_tab_by_chrScaff_mostprobes.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/order_bionano_scaffolds_to_chr_level/f1_dam_bionano_HD_rc_Probes_tab.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/order_bionano_scaffolds_to_chr_level/f1_sire_bionano_HD_rc_Probes_tab.RData")

# END HD and rc probes #

#to get ungapped scaffold length
dam_scaffold_length <- dam_contig_scaff %>% group_by(object) %>% summarise(scaff_len = sum(clength))
sire_scaffold_length <- sire_contig_scaff %>% group_by(object) %>% summarise(scaff_len = sum(clength))

##### Dam #####
#new var to cbind with *_contig_scaff
scaff_len <- c() #scaffold_length var
chr_HD <- c() #probes HD 
HD_n <- c()
chr_rc <- c() #probes rc
rc_n <- c()
Ref_name <- c() #mashmap vars
Ref_length <- c()
Ref_start <- c()
Ref_end <- c()
query_align <- c()
proportion_ref <- c()

for (i in 1:nrow(dam_contig_scaff)){
  #scaffold_length var
  logic <- dam_scaffold_length$object %in% dam_contig_scaff$object[i]
  holder <- dam_scaffold_length$scaff_len[logic]
  if (!is.integer0(holder)) {scaff_len <- c(scaff_len,holder)} else {scaff_len <- c(scaff_len,NA)}
  
  #probes HD 
  logic <- f1_dam_bionano_HDProbes_tab_by_chrScaff_mostprobes$scaffold %in% dam_contig_scaff$object[i]
  holder <- f1_dam_bionano_HDProbes_tab_by_chrScaff_mostprobes$chromosome[logic]
  if (!is.integer0(holder)) {chr_HD <- c(chr_HD,holder)} else {chr_HD <- c(chr_HD,NA)}
  holder <- f1_dam_bionano_HDProbes_tab_by_chrScaff_mostprobes$n[logic]
  if (!is.integer0(holder)) {HD_n <- c(HD_n,holder)} else {HD_n <- c(HD_n,NA)}
  
  #probes rc
  logic <- f1_dam_bionano_fasta_rcmap_tab_by_chrScaff_mostprobes$scaffold %in% dam_contig_scaff$object[i]
  holder <- f1_dam_bionano_fasta_rcmap_tab_by_chrScaff_mostprobes$chromosome[logic]
  if (!is.integer0(holder)) {chr_rc <- c(chr_rc,holder)} else {chr_rc <- c(chr_rc,NA)}
  holder <- f1_dam_bionano_fasta_rcmap_tab_by_chrScaff_mostprobes$n[logic]
  if (!is.integer0(holder)) {rc_n <- c(rc_n,holder)} else {rc_n <- c(rc_n,NA)}
  
  #mashmap var
  logic <- dam_mashmap$Query_name %in% dam_contig_scaff$object[i]
  holder <- dam_mashmap$Ref_name[logic]
  if (!identical(holder, character(0))) {Ref_name <- c(Ref_name,holder)} else {Ref_name <- c(Ref_name,NA)}
  holder <- dam_mashmap$Ref_length[logic]
  if (!is.integer0(holder)) {Ref_length <- c(Ref_length,holder)} else {Ref_length <- c(Ref_length,NA)}
  holder <- dam_mashmap$Ref_start[logic]
  if (!is.integer0(holder)) {Ref_start <- c(Ref_start,holder)} else {Ref_start <- c(Ref_start,NA)}
  holder <- dam_mashmap$Ref_end[logic]
  if (!is.integer0(holder)) {Ref_end <- c(Ref_end,holder)} else {Ref_end <- c(Ref_end,NA)}
  holder <- dam_mashmap$query_align[logic]
  if (!is.integer0(holder)) {query_align <- c(query_align,holder)} else {query_align <- c(query_align,NA)}
  holder <- dam_mashmap$proportion_ref[logic]
  if (!identical(holder, numeric(0))) {proportion_ref <- c(proportion_ref,holder)} else {proportion_ref <- c(proportion_ref,NA)}
}

#cbind all vars to make big table
dam_scaff_order_DF <- dam_contig_scaff %>% mutate(scaff_len = scaff_len,chr_HD = chr_HD, HD_n = HD_n,chr_rc = chr_rc,
                                                  rc_n = rc_n, Ref_name = Ref_name, Ref_length = Ref_length, 
                                                  Ref_start = Ref_start, Ref_end = Ref_end, query_align = query_align,
                                                  proportion_ref = proportion_ref)

##### Sire #####
#new var to cbind with *_contig_scaff
scaff_len <- c() #scaffold_length var
chr_HD <- c() #probes HD 
HD_n <- c()
chr_rc <- c() #probes rc
rc_n <- c()
Ref_name <- c() #mashmap vars
Ref_length <- c()
Ref_start <- c()
Ref_end <- c()
query_align <- c()
proportion_ref <- c()

for (i in 1:nrow(sire_contig_scaff)){
  #scaffold_length var
  logic <- sire_scaffold_length$object %in% sire_contig_scaff$object[i]
  holder <- sire_scaffold_length$scaff_len[logic]
  if (!is.integer0(holder)) {scaff_len <- c(scaff_len,holder)} else {scaff_len <- c(scaff_len,NA)}
  
  #probes HD 
  logic <- f1_sire_bionano_HDProbes_tab_by_chrScaff_mostprobes$scaffold %in% sire_contig_scaff$object[i]
  holder <- f1_sire_bionano_HDProbes_tab_by_chrScaff_mostprobes$chromosome[logic]
  if (!is.integer0(holder)) {chr_HD <- c(chr_HD,holder)} else {chr_HD <- c(chr_HD,NA)}
  holder <- f1_sire_bionano_HDProbes_tab_by_chrScaff_mostprobes$n[logic]
  if (!is.integer0(holder)) {HD_n <- c(HD_n,holder)} else {HD_n <- c(HD_n,NA)}
  
  #probes rc
  logic <- f1_sire_bionano_fasta_rcmap_tab_by_chrScaff_mostprobes$scaffold %in% sire_contig_scaff$object[i]
  holder <- f1_sire_bionano_fasta_rcmap_tab_by_chrScaff_mostprobes$chromosome[logic]
  if (!is.integer0(holder)) {chr_rc <- c(chr_rc,holder)} else {chr_rc <- c(chr_rc,NA)}
  holder <- f1_sire_bionano_fasta_rcmap_tab_by_chrScaff_mostprobes$n[logic]
  if (!is.integer0(holder)) {rc_n <- c(rc_n,holder)} else {rc_n <- c(rc_n,NA)}
  
  #mashmap var
  logic <- sire_mashmap$Query_name %in% sire_contig_scaff$object[i]
  holder <- sire_mashmap$Ref_name[logic]
  if (!identical(holder, character(0))) {Ref_name <- c(Ref_name,holder)} else {Ref_name <- c(Ref_name,NA)}
  holder <- sire_mashmap$Ref_length[logic]
  if (!is.integer0(holder)) {Ref_length <- c(Ref_length,holder)} else {Ref_length <- c(Ref_length,NA)}
  holder <- sire_mashmap$Ref_start[logic]
  if (!is.integer0(holder)) {Ref_start <- c(Ref_start,holder)} else {Ref_start <- c(Ref_start,NA)}
  holder <- sire_mashmap$Ref_end[logic]
  if (!is.integer0(holder)) {Ref_end <- c(Ref_end,holder)} else {Ref_end <- c(Ref_end,NA)}
  holder <- sire_mashmap$query_align[logic]
  if (!is.integer0(holder)) {query_align <- c(query_align,holder)} else {query_align <- c(query_align,NA)}
  holder <- sire_mashmap$proportion_ref[logic]
  if (!identical(holder, numeric(0))) {proportion_ref <- c(proportion_ref,holder)} else {proportion_ref <- c(proportion_ref,NA)}
}

#cbind all vars to make big table
sire_scaff_order_DF <- sire_contig_scaff %>% mutate(scaff_len = scaff_len,chr_HD = chr_HD, HD_n = HD_n,chr_rc = chr_rc,
                                                    rc_n = rc_n, Ref_name = Ref_name, Ref_length = Ref_length, 
                                                    Ref_start = Ref_start, Ref_end = Ref_end, query_align = query_align,
                                                    proportion_ref = proportion_ref)

save(dam_scaff_order_DF,sire_scaff_order_DF,
     file = "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/order_bionano_scaffolds_to_chr_level/scaff_order_DF.RData")
