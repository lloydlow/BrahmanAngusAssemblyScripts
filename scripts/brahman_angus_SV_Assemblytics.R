#------------------------------------------------------
# Program name: brahman_angus_SV_Assemblytics.R
# Objective: analyse the categories of SV, (i) diff in Angus
#           (ii) diff in Brahman, (iii) diff in both Angus Brahman
#           (iv) same in Brahman Angus, Hereford is the different one
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(readr)
library(dplyr)
library(ggplot2)
library(VennDiagram)
library(tidyr)
library(easyGgplot2)
library(stringr)

#angus
# path to assemblytics folder
dir1 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/genomeAlignmentToCattle/assemblytics_results/ARSUCD_VS_ANGUS.Assemblytics_results/"

# reading ARSUCD_VS_ANGUS.Assemblytics_structural_variants.bed
# this is ungapped ARSUCD as ref, Angus ungapped ref as query
path1 <- paste0(dir1,"ARSUCD_VS_ANGUS.Assemblytics_structural_variants.bed")

ARSUCD_VS_ANGUS <- read_tsv(path1)

ARSUCD_VS_ANGUS$uniqueName <- paste(ARSUCD_VS_ANGUS$reference,ARSUCD_VS_ANGUS$ref_start,ARSUCD_VS_ANGUS$ref_stop,ARSUCD_VS_ANGUS$size,
                                    ARSUCD_VS_ANGUS$strand,ARSUCD_VS_ANGUS$type,sep = "_")

table(ARSUCD_VS_ANGUS$method,ARSUCD_VS_ANGUS$type)

# Deletion Insertion Repeat_contraction Repeat_expansion Tandem_contraction Tandem_expansion
# between_alignments     3078      3401                769              627               1283             1381
# within_alignment        838       888                  0                0                  0                0

# determine duplicates i.e. exactly same REF coordinate called for same SV
n_occur_ARSUCD_VS_ANGUS <- data.frame(table(ARSUCD_VS_ANGUS$uniqueName))
n_occur_ARSUCD_VS_ANGUS[n_occur_ARSUCD_VS_ANGUS$Freq > 1,]

ARSUCD_VS_ANGUS[ARSUCD_VS_ANGUS$uniqueName %in% n_occur_ARSUCD_VS_ANGUS$Var1[n_occur_ARSUCD_VS_ANGUS$Freq > 1],]

# reference   ref_start ref_stop ID       size strand type    ref_gap_size query_gap_size query_coordinates    method  uniqueName      
# <chr>           <dbl>    <dbl> <chr>   <dbl> <chr>  <chr>          <dbl>          <dbl> <chr>                <chr>   <chr>           
# 1 12:1-21495…     54606    54744 Assemb…   138 +      Deleti…          138              0 12:1-4790906:22858-… betwee… 12:1-21495996_5…
# 2 12:1-21495…     54606    54744 Assemb…   138 +      Deleti…          138              0 12:1-4790906:31695-… betwee… 12:1-21495996_5…
# 3 17:7115260…    492848   493225 Assemb…    75 +      Tandem…         -377           -452 17:25313034-7106450… betwee… 17:71152603-731…
# 4 17:7115260…    492848   493225 Assemb…    75 +      Tandem…         -377           -452 17:71065006-7257845… betwee… 17:71152603-731…
# 5 18:1-57154…  49256213 49256514 Assemb…   267 +      Deleti…          301             34 18:14821707-3485755… betwee… 18:1-57154076_4…
# 6 18:1-57154…  49256213 49256514 Assemb…   267 +      Deleti…          301             34 18:14821707-3485755… betwee… 18:1-57154076_4…

dupid_ARSUCD_VS_ANGUS <- ARSUCD_VS_ANGUS[ARSUCD_VS_ANGUS$uniqueName %in% n_occur_ARSUCD_VS_ANGUS$Var1[n_occur_ARSUCD_VS_ANGUS$Freq > 1],]$ID

dupid_ARSUCD_VS_ANGUS <- dupid_ARSUCD_VS_ANGUS[seq(1,length(dupid_ARSUCD_VS_ANGUS),by=2)]

#remove duplicate hits to same REF
ARSUCD_VS_ANGUS <- ARSUCD_VS_ANGUS[!ARSUCD_VS_ANGUS$ID %in% dupid_ARSUCD_VS_ANGUS,]

#split ref and query back to chr:start:stop for easy fasta extraction later
ARSUCD_VS_ANGUS <- ARSUCD_VS_ANGUS %>% 
  separate(reference, c("chr","start","end"),extra = "drop",remove = FALSE,convert = TRUE) %>%
  separate(query_coordinates, c("query_chr","query_start","query_end","cut_start","cut_end","query_strand"),
           extra = "drop",remove = FALSE,convert = TRUE)

ARSUCD_VS_ANGUS$query_strand <- "+"

#giving ref proper coordinates (ARSUCD1.2)
ARSUCD_VS_ANGUS$start <- ARSUCD_VS_ANGUS$start + ARSUCD_VS_ANGUS$ref_start - 1
ARSUCD_VS_ANGUS$end <- ARSUCD_VS_ANGUS$ref_stop - ARSUCD_VS_ANGUS$ref_start + ARSUCD_VS_ANGUS$start

#giving query proper coordinates
ARSUCD_VS_ANGUS$query_start <- ARSUCD_VS_ANGUS$query_start + ARSUCD_VS_ANGUS$cut_start - 1
ARSUCD_VS_ANGUS$query_end <- ARSUCD_VS_ANGUS$cut_end - ARSUCD_VS_ANGUS$cut_start + ARSUCD_VS_ANGUS$query_start

#for paper, total size of Angus affected SV
sum(ARSUCD_VS_ANGUS$size)

#chr X/Y filter out from both query and ref
ARSUCD_VS_ANGUS <- ARSUCD_VS_ANGUS %>% filter(chr != "X") %>% filter(query_chr != "Y")

#filter table into 6 categories
#Deletion
ARSUCD_VS_ANGUS_Deletion <- ARSUCD_VS_ANGUS %>% filter(type == "Deletion")

#Insertion
ARSUCD_VS_ANGUS_Insertion <- ARSUCD_VS_ANGUS %>% filter(type == "Insertion")
  
#Repeat_contraction 
ARSUCD_VS_ANGUS_Repeat_contraction <- ARSUCD_VS_ANGUS %>% filter(type == "Repeat_contraction")
  
#Repeat_expansion
ARSUCD_VS_ANGUS_Repeat_expansion <- ARSUCD_VS_ANGUS %>% filter(type == "Repeat_expansion")
  
#Tandem_contraction
ARSUCD_VS_ANGUS_Tandem_contraction <- ARSUCD_VS_ANGUS %>% filter(type == "Tandem_contraction")

#Tandem_expansion
ARSUCD_VS_ANGUS_Tandem_expansion <- ARSUCD_VS_ANGUS %>% filter(type == "Tandem_expansion")

#brahman
# path to assemblytics folder
dir1 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/genomeAlignmentToCattle/assemblytics_results/ARSUCD_VS_BRAHMAN.Assemblytics_results/"

# reading ARSUCD_VS_BRAHMAN.Assemblytics_structural_variants.bed
# this is ungapped ARSUCD as ref, Brahman ungapped ref as query
path1 <- paste0(dir1,"ARSUCD_VS_BRAHMAN.Assemblytics_structural_variants.bed")

ARSUCD_VS_BRAHMAN <- read_tsv(path1)

ARSUCD_VS_BRAHMAN$uniqueName <- paste(ARSUCD_VS_BRAHMAN$reference,ARSUCD_VS_BRAHMAN$ref_start,ARSUCD_VS_BRAHMAN$ref_stop,ARSUCD_VS_BRAHMAN$size,
                                      ARSUCD_VS_BRAHMAN$strand,ARSUCD_VS_BRAHMAN$type,sep = "_")

table(ARSUCD_VS_BRAHMAN$method,ARSUCD_VS_BRAHMAN$type)

# Deletion Insertion Repeat_contraction Repeat_expansion Tandem_contraction Tandem_expansion
# between_alignments     6861      7545               2069             1682               2008             2291
# within_alignment       2156      2302                  0                0                  0                0

# determine duplicates i.e. exactly same REF coordinate called for same SV
n_occur_ARSUCD_VS_BRAHMAN <- data.frame(table(ARSUCD_VS_BRAHMAN$uniqueName))
n_occur_ARSUCD_VS_BRAHMAN[n_occur_ARSUCD_VS_BRAHMAN$Freq > 1,]

ARSUCD_VS_BRAHMAN[ARSUCD_VS_BRAHMAN$uniqueName %in% n_occur_ARSUCD_VS_BRAHMAN$Var1[n_occur_ARSUCD_VS_BRAHMAN$Freq > 1],]

# reference  ref_start ref_stop ID       size strand type   ref_gap_size query_gap_size query_coordinates  method uniqueName  
# <chr>          <dbl>    <dbl> <chr>   <dbl> <chr>  <chr>         <dbl>          <dbl> <chr>              <chr>  <chr>       
# 1 X:9022757…     21147    21201 Assemb…    54 +      Delet…           54              0 X:53349294-535251… withi… X:90227579-…
# 2 X:9022757…     21147    21201 Assemb…    54 +      Delet…           54              0 X:53525658-536626… withi… X:90227579-…
# 3 12:1-2149…     84453    90037 Assemb…  5900 +      Tande…        -5584            316 12:1-7194098:3106… betwe… 12:1-214959…
# 4 12:1-2149…     85826    85940 Assemb…   114 +      Delet…          114              0 12:1-7194098:3277… betwe… 12:1-214959…
# 5 12:1-2149…     84453    90037 Assemb…  5900 +      Tande…        -5584            316 12:1-7194098:3687… betwe… 12:1-214959…
# 6 12:1-2149…     85826    85940 Assemb…   114 +      Delet…          114              0 12:1-7194098:3858… betwe… 12:1-214959…
# 7 24:1-6216…  61908433 61908732 Assemb…   299 +      Delet…          299              0 24:448299-6211282… betwe… 24:1-621618…
# 8 24:1-6216…  61908433 61908732 Assemb…   299 +      Delet…          299              0 24:62113326-62414… betwe… 24:1-621618…
# 9 X:3224075…    272065   272201 Assemb…   136 +      Delet…          136              0 X:113589379-11364… betwe… X:32240750-…
# 10 X:3224075…    272065   272201 Assemb…   136 +      Delet…          136              0 X:113643604-11729… betwe… X:32240750-…

dupid_ARSUCD_VS_BRAHMAN <- ARSUCD_VS_BRAHMAN[ARSUCD_VS_BRAHMAN$uniqueName %in% n_occur_ARSUCD_VS_BRAHMAN$Var1[n_occur_ARSUCD_VS_BRAHMAN$Freq > 1],]$ID

#special treament due to same SV called different types, can be solved by arrange ref, then by type, the take and skip by 2
dupid_ARSUCD_VS_BRAHMAN <- dupid_ARSUCD_VS_BRAHMAN[c(1,3,4,7,9)]

#remove duplicate hits to same REF
ARSUCD_VS_BRAHMAN <- ARSUCD_VS_BRAHMAN[!ARSUCD_VS_BRAHMAN$ID %in% dupid_ARSUCD_VS_BRAHMAN,]

#split ref and query back to chr:start:stop for easy fasta extraction later
ARSUCD_VS_BRAHMAN <- ARSUCD_VS_BRAHMAN %>% 
  separate(reference, c("chr","start","end"),extra = "drop",remove = FALSE,convert = TRUE) %>%
  separate(query_coordinates, c("query_chr","query_start","query_end","cut_start","cut_end","query_strand"),
           extra = "drop",remove = FALSE,convert = TRUE)

ARSUCD_VS_BRAHMAN$query_strand <- "+"

#giving ref proper coordinates (ARSUCD1.2)
ARSUCD_VS_BRAHMAN$start <- ARSUCD_VS_BRAHMAN$start + ARSUCD_VS_BRAHMAN$ref_start - 1
ARSUCD_VS_BRAHMAN$end <- ARSUCD_VS_BRAHMAN$ref_stop - ARSUCD_VS_BRAHMAN$ref_start + ARSUCD_VS_BRAHMAN$start

#giving query proper coordinates
ARSUCD_VS_BRAHMAN$query_start <- ARSUCD_VS_BRAHMAN$query_start + ARSUCD_VS_BRAHMAN$cut_start - 1
ARSUCD_VS_BRAHMAN$query_end <- ARSUCD_VS_BRAHMAN$cut_end - ARSUCD_VS_BRAHMAN$cut_start + ARSUCD_VS_BRAHMAN$query_start

#for paper, size Brahman affected SV
sum(ARSUCD_VS_BRAHMAN$size)

#chr X/Y filter out from both query and ref
ARSUCD_VS_BRAHMAN <- ARSUCD_VS_BRAHMAN %>% filter(chr != "X") %>% filter(query_chr != "X")

#filter table into 6 categories
#Deletion
ARSUCD_VS_BRAHMAN_Deletion <- ARSUCD_VS_BRAHMAN %>% filter(type == "Deletion")

#Insertion
ARSUCD_VS_BRAHMAN_Insertion <- ARSUCD_VS_BRAHMAN %>% filter(type == "Insertion")

#Repeat_contraction 
ARSUCD_VS_BRAHMAN_Repeat_contraction <- ARSUCD_VS_BRAHMAN %>% filter(type == "Repeat_contraction")

#Repeat_expansion
ARSUCD_VS_BRAHMAN_Repeat_expansion <- ARSUCD_VS_BRAHMAN %>% filter(type == "Repeat_expansion")

#Tandem_contraction
ARSUCD_VS_BRAHMAN_Tandem_contraction <- ARSUCD_VS_BRAHMAN %>% filter(type == "Tandem_contraction")

#Tandem_expansion
ARSUCD_VS_BRAHMAN_Tandem_expansion <- ARSUCD_VS_BRAHMAN %>% filter(type == "Tandem_expansion")

#####################################################################################################
# Crossing both Brahman and Angus tables
# Deletion
length(intersect(ARSUCD_VS_ANGUS_Deletion$uniqueName,ARSUCD_VS_BRAHMAN_Deletion$uniqueName))

venn.diagram(
  x = list(ARSUCD_VS_ANGUS_Deletion$uniqueName, ARSUCD_VS_BRAHMAN_Deletion$uniqueName),cat.pos = 0,
  category.names = c("Angus" , "Brahman"),fill = c('turquoise', 'red'),cex=2,cat.cex=1.5,
  filename = 'AngusVsBrahmanSVDeletion_venn.png')

# extract just the Brahman SV with no intersection to Angus
ARSUCD_VS_ANGUS_Deletion_BrahmanOnly <- setdiff(ARSUCD_VS_BRAHMAN_Deletion$uniqueName,ARSUCD_VS_ANGUS_Deletion$uniqueName)

ARSUCD_VS_ANGUS_Deletion_BrahmanOnly_DF <- ARSUCD_VS_BRAHMAN_Deletion[ARSUCD_VS_BRAHMAN_Deletion$uniqueName %in% ARSUCD_VS_ANGUS_Deletion_BrahmanOnly,]

# extract just the Angus SV with no intersection to Brahman
ARSUCD_VS_ANGUS_Deletion_AngusOnly <- setdiff(ARSUCD_VS_ANGUS_Deletion$uniqueName,ARSUCD_VS_BRAHMAN_Deletion$uniqueName)

ARSUCD_VS_ANGUS_Deletion_AngusOnly_DF <- ARSUCD_VS_ANGUS_Deletion[ARSUCD_VS_ANGUS_Deletion$uniqueName %in% ARSUCD_VS_ANGUS_Deletion_AngusOnly,]

# histograms to compare Angus- and Brahman-specific Deletion
# ggplot2.histogram(data=ARSUCD_VS_ANGUS_Deletion_BrahmanOnly_DF, xName= 'start', xtitle="Position",
#                   groupName='chr', legendPosition="right",
#                   faceting=TRUE, facetingVarNames="chr",
#                   binwidth = 0.1e6,yShowTitle=FALSE,yShowTickLabel=FALSE,
#                   hideAxisTicks=TRUE) 
# 
# ggplot2.histogram(data=ARSUCD_VS_ANGUS_Deletion_AngusOnly_DF, xName= 'start', xtitle="Position",
#                   groupName='chr', legendPosition="right",
#                   faceting=TRUE, facetingVarNames="chr",
#                   binwidth = 0.1e6,yShowTitle=FALSE,yShowTickLabel=FALSE,
#                   hideAxisTicks=TRUE) 

# Insertion
length(intersect(ARSUCD_VS_ANGUS_Insertion$uniqueName,ARSUCD_VS_BRAHMAN_Insertion$uniqueName))

venn.diagram(
  x = list(ARSUCD_VS_ANGUS_Insertion$uniqueName, ARSUCD_VS_BRAHMAN_Insertion$uniqueName),cat.pos = 0,
  category.names = c("Angus" , "Brahman"),fill = c('turquoise', 'red'),cex=2,cat.cex=1.5,
  filename = 'AngusVsBrahmanSVInsertion_venn.png')

# extract just the Brahman SV with no intersection to Angus
ARSUCD_VS_ANGUS_Insertion_BrahmanOnly <- setdiff(ARSUCD_VS_BRAHMAN_Insertion$uniqueName,ARSUCD_VS_ANGUS_Insertion$uniqueName)

ARSUCD_VS_ANGUS_Insertion_BrahmanOnly_DF <- ARSUCD_VS_BRAHMAN_Insertion[ARSUCD_VS_BRAHMAN_Insertion$uniqueName %in% ARSUCD_VS_ANGUS_Insertion_BrahmanOnly,]

# extract just the Angus SV with no intersection to Brahman
ARSUCD_VS_ANGUS_Insertion_AngusOnly <- setdiff(ARSUCD_VS_ANGUS_Insertion$uniqueName,ARSUCD_VS_BRAHMAN_Insertion$uniqueName)

ARSUCD_VS_ANGUS_Insertion_AngusOnly_DF <- ARSUCD_VS_ANGUS_Insertion[ARSUCD_VS_ANGUS_Insertion$uniqueName %in% ARSUCD_VS_ANGUS_Insertion_AngusOnly,]

# Repeat_contraction
length(intersect(ARSUCD_VS_ANGUS_Repeat_contraction$uniqueName,ARSUCD_VS_BRAHMAN_Repeat_contraction$uniqueName))

venn.diagram(
  x = list(ARSUCD_VS_ANGUS_Repeat_contraction$uniqueName, ARSUCD_VS_BRAHMAN_Repeat_contraction$uniqueName),cat.pos = 0,
  category.names = c("Angus" , "Brahman"),fill = c('turquoise', 'red'),cex=2,cat.cex=1.5,
  filename = 'AngusVsBrahmanSVRepeat_contraction_venn.png')

# extract just the Brahman SV with no intersection to Angus
ARSUCD_VS_ANGUS_Repeat_contraction_BrahmanOnly <- setdiff(ARSUCD_VS_BRAHMAN_Repeat_contraction$uniqueName,ARSUCD_VS_ANGUS_Repeat_contraction$uniqueName)

ARSUCD_VS_ANGUS_Repeat_contraction_BrahmanOnly_DF <- ARSUCD_VS_BRAHMAN_Repeat_contraction[ARSUCD_VS_BRAHMAN_Repeat_contraction$uniqueName %in% ARSUCD_VS_ANGUS_Repeat_contraction_BrahmanOnly,]

# extract just the Angus SV with no intersection to Brahman
ARSUCD_VS_ANGUS_Repeat_contraction_AngusOnly <- setdiff(ARSUCD_VS_ANGUS_Repeat_contraction$uniqueName,ARSUCD_VS_BRAHMAN_Repeat_contraction$uniqueName)

ARSUCD_VS_ANGUS_Repeat_contraction_AngusOnly_DF <- ARSUCD_VS_ANGUS_Repeat_contraction[ARSUCD_VS_ANGUS_Repeat_contraction$uniqueName %in% ARSUCD_VS_ANGUS_Repeat_contraction_AngusOnly,]

# Repeat_expansion
length(intersect(ARSUCD_VS_ANGUS_Repeat_expansion$uniqueName,ARSUCD_VS_BRAHMAN_Repeat_expansion$uniqueName))

venn.diagram(
  x = list(ARSUCD_VS_ANGUS_Repeat_expansion$uniqueName, ARSUCD_VS_BRAHMAN_Repeat_expansion$uniqueName),cat.pos = 0,
  category.names = c("Angus" , "Brahman"),fill = c('turquoise', 'red'),cex=2,cat.cex=1.5,
  filename = 'AngusVsBrahmanSVRepeat_expansion_venn.png')

# extract just the Brahman SV with no intersection to Angus
ARSUCD_VS_ANGUS_Repeat_expansion_BrahmanOnly <- setdiff(ARSUCD_VS_BRAHMAN_Repeat_expansion$uniqueName,ARSUCD_VS_ANGUS_Repeat_expansion$uniqueName)

ARSUCD_VS_ANGUS_Repeat_expansion_BrahmanOnly_DF <- ARSUCD_VS_BRAHMAN_Repeat_expansion[ARSUCD_VS_BRAHMAN_Repeat_expansion$uniqueName %in% ARSUCD_VS_ANGUS_Repeat_expansion_BrahmanOnly,]

# extract just the Angus SV with no intersection to Brahman
ARSUCD_VS_ANGUS_Repeat_expansion_AngusOnly <- setdiff(ARSUCD_VS_ANGUS_Repeat_expansion$uniqueName,ARSUCD_VS_BRAHMAN_Repeat_expansion$uniqueName)

ARSUCD_VS_ANGUS_Repeat_expansion_AngusOnly_DF <- ARSUCD_VS_ANGUS_Repeat_expansion[ARSUCD_VS_ANGUS_Repeat_expansion$uniqueName %in% ARSUCD_VS_ANGUS_Repeat_expansion_AngusOnly,]

# Tandem_contraction
length(intersect(ARSUCD_VS_ANGUS_Tandem_contraction$uniqueName,ARSUCD_VS_BRAHMAN_Tandem_contraction$uniqueName))

venn.diagram(
  x = list(ARSUCD_VS_ANGUS_Tandem_contraction$uniqueName, ARSUCD_VS_BRAHMAN_Tandem_contraction$uniqueName),cat.pos = 0,
  category.names = c("Angus" , "Brahman"),fill = c('turquoise', 'red'),cex=2,cat.cex=1.5,
  filename = 'AngusVsBrahmanSVTandem_contraction_venn.png')

# extract just the Brahman SV with no intersection to Angus
ARSUCD_VS_ANGUS_Tandem_contraction_BrahmanOnly <- setdiff(ARSUCD_VS_BRAHMAN_Tandem_contraction$uniqueName,ARSUCD_VS_ANGUS_Tandem_contraction$uniqueName)

ARSUCD_VS_ANGUS_Tandem_contraction_BrahmanOnly_DF <- ARSUCD_VS_BRAHMAN_Tandem_contraction[ARSUCD_VS_BRAHMAN_Tandem_contraction$uniqueName %in% ARSUCD_VS_ANGUS_Tandem_contraction_BrahmanOnly,]

# extract just the Angus SV with no intersection to Brahman
ARSUCD_VS_ANGUS_Tandem_contraction_AngusOnly <- setdiff(ARSUCD_VS_ANGUS_Tandem_contraction$uniqueName,ARSUCD_VS_BRAHMAN_Tandem_contraction$uniqueName)

ARSUCD_VS_ANGUS_Tandem_contraction_AngusOnly_DF <- ARSUCD_VS_ANGUS_Tandem_contraction[ARSUCD_VS_ANGUS_Tandem_contraction$uniqueName %in% ARSUCD_VS_ANGUS_Tandem_contraction_AngusOnly,]

# Tandem_expansion
length(intersect(ARSUCD_VS_ANGUS_Tandem_expansion$uniqueName,ARSUCD_VS_BRAHMAN_Tandem_expansion$uniqueName))

venn.diagram(
  x = list(ARSUCD_VS_ANGUS_Tandem_expansion$uniqueName, ARSUCD_VS_BRAHMAN_Tandem_expansion$uniqueName),cat.pos = 0,
  category.names = c("Angus" , "Brahman"),fill = c('turquoise', 'red'),cex=2,cat.cex=1.5,
  filename = 'AngusVsBrahmanSVTandem_expansion_venn.png')

# extract just the Brahman SV with no intersection to Angus
ARSUCD_VS_ANGUS_Tandem_expansion_BrahmanOnly <- setdiff(ARSUCD_VS_BRAHMAN_Tandem_expansion$uniqueName,ARSUCD_VS_ANGUS_Tandem_expansion$uniqueName)

ARSUCD_VS_ANGUS_Tandem_expansion_BrahmanOnly_DF <- ARSUCD_VS_BRAHMAN_Tandem_expansion[ARSUCD_VS_BRAHMAN_Tandem_expansion$uniqueName %in% ARSUCD_VS_ANGUS_Tandem_expansion_BrahmanOnly,]

# extract just the Angus SV with no intersection to Brahman
ARSUCD_VS_ANGUS_Tandem_expansion_AngusOnly <- setdiff(ARSUCD_VS_ANGUS_Tandem_expansion$uniqueName,ARSUCD_VS_BRAHMAN_Tandem_expansion$uniqueName)

ARSUCD_VS_ANGUS_Tandem_expansion_AngusOnly_DF <- ARSUCD_VS_ANGUS_Tandem_expansion[ARSUCD_VS_ANGUS_Tandem_expansion$uniqueName %in% ARSUCD_VS_ANGUS_Tandem_expansion_AngusOnly,]

#Save the 6 SV types as Rdata
save(ARSUCD_VS_ANGUS_Deletion_BrahmanOnly_DF,ARSUCD_VS_ANGUS_Deletion_AngusOnly_DF,
     ARSUCD_VS_ANGUS_Insertion_BrahmanOnly_DF,ARSUCD_VS_ANGUS_Insertion_AngusOnly_DF,
     ARSUCD_VS_ANGUS_Repeat_contraction_BrahmanOnly_DF,ARSUCD_VS_ANGUS_Repeat_contraction_AngusOnly_DF,
     ARSUCD_VS_ANGUS_Repeat_expansion_BrahmanOnly_DF,ARSUCD_VS_ANGUS_Repeat_expansion_AngusOnly_DF,
     ARSUCD_VS_ANGUS_Tandem_contraction_BrahmanOnly_DF,ARSUCD_VS_ANGUS_Tandem_contraction_AngusOnly_DF,
     ARSUCD_VS_ANGUS_Tandem_expansion_BrahmanOnly_DF,ARSUCD_VS_ANGUS_Tandem_expansion_AngusOnly_DF, 
     file = "SV6typesBrahmanAngus.RData")

