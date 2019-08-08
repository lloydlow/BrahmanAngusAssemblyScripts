#------------------------------------------------------
# Program name: brahman_angus_CNVR_liftover_v1.R
# Objective: analyse Derek CNVR liftover that will give 
#     common arsucd coor     
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(UpSetR)
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(GenomicFeatures)
library(ComplexHeatmap)
library(circlize)

#reproducing Derek's upset results
dir1 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/CopyNumberVariation/20190530/cnvr_data/"

upset_angus <- scan(paste0(dir1,"angus.upset.list"),sep="\n")
upset_brahman <- scan(paste0(dir1,"brahman.upset.list"),sep="\n")
upset_arsucd <- scan(paste0(dir1,"arsucd.upset.list"),sep="\n")

listInput <- list(angus=upset_angus,brahman=upset_brahman,hereford=upset_arsucd)

upset(fromList(listInput), order.by = "freq")

############################# Start from lifover modified results, without sex chr and unplaced ###########################
rm(list=ls())

# path to CNV results
dir2 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/CopyNumberVariation/20190530/cnvr_data/results_liftover/"

# angus
path1 <- paste0(dir2,"angus_arsucd_modi_coor")

angus.cnvrs.mapped <- read_tsv(path1,col_names = FALSE, col_types = "cddcdd")
colnames(angus.cnvrs.mapped) <- c("chr","start","end","species","unsure1","unsure2")

#filter for size less than 1 mil
angus.cnvrs.mapped_modi <- angus.cnvrs.mapped %>% mutate(size=end - start) %>% filter(size < 1e6)
angus.cnvrs.mapped_modi$unique_name <- paste(angus.cnvrs.mapped_modi$chr,angus.cnvrs.mapped_modi$start,angus.cnvrs.mapped_modi$end,sep = "_")

angus.cnvrs.mapped_modi$strand <- rep("*",nrow(angus.cnvrs.mapped_modi))

angus.cnvrs.mapped_modi <- angus.cnvrs.mapped_modi %>% filter(chr != "tig00020276_arrow_arrow_40739300_45952718")

cnv_interval_angus <- makeGRangesFromDataFrame(angus.cnvrs.mapped_modi, keep.extra.columns = TRUE, 
                                               seqnames.field="chr", start.field="start", 
                                               end.field="end", strand.field="strand")

#brahman
path2 <- paste0(dir2,"brahman_arsucd_modi_coor")

brahman.cnvrs.mapped <- read_tsv(path2,col_names = FALSE, col_types = "cddcdd")
colnames(brahman.cnvrs.mapped) <- c("chr","start","end","species","unsure1","unsure2")

#filter for size less than 1 mil
brahman.cnvrs.mapped_modi <- brahman.cnvrs.mapped %>% mutate(size=end - start) %>% filter(size < 1e6)
brahman.cnvrs.mapped_modi$unique_name <- paste(brahman.cnvrs.mapped_modi$chr,brahman.cnvrs.mapped_modi$start,brahman.cnvrs.mapped_modi$end,sep = "_")

brahman.cnvrs.mapped_modi$strand <- rep("*",nrow(brahman.cnvrs.mapped_modi))

brahman.cnvrs.mapped_modi <- brahman.cnvrs.mapped_modi %>% filter(chr != "tig00000831_arrow_arrow_obj") %>% 
  filter(chr != "tig00001951_arrow_arrow_obj") %>% filter(chr != "tig00002091_arrow_arrow_obj") %>%
  filter(chr != "X")

cnv_interval_brahman <- makeGRangesFromDataFrame(brahman.cnvrs.mapped_modi, keep.extra.columns = TRUE, 
                                                 seqnames.field="chr", start.field="start", 
                                                 end.field="end", strand.field="strand")

#arsucd
path3 <- paste0(dir2,"arsucd.cnvrs_regions.bed")

arsucd.cnvrs.mapped <- read_tsv(path3,col_names = FALSE, col_types = "cddcdd")
colnames(arsucd.cnvrs.mapped) <- c("chr","start","end","species","unsure1","unsure2")

#filter for size less than 1 mil
arsucd.cnvrs.mapped_modi <- arsucd.cnvrs.mapped %>% mutate(size=end - start) %>% filter(size < 1e6)
arsucd.cnvrs.mapped_modi$unique_name <- paste(arsucd.cnvrs.mapped_modi$chr,arsucd.cnvrs.mapped_modi$start,arsucd.cnvrs.mapped_modi$end,sep = "_")

arsucd.cnvrs.mapped_modi$strand <- rep("*",nrow(arsucd.cnvrs.mapped_modi))

arsucd.cnvrs.mapped_modi <- arsucd.cnvrs.mapped_modi %>% mutate(selc = as.numeric(chr)) %>%
  filter(!is.na(selc)) %>% dplyr::select(chr:strand)

cnv_interval_arsucd <- makeGRangesFromDataFrame(arsucd.cnvrs.mapped_modi, keep.extra.columns = TRUE, 
                                                seqnames.field="chr", start.field="start", 
                                                end.field="end", strand.field="strand")

#combine cnv granges as list
listInput_complxhtmap <- list(angus=cnv_interval_angus,brahman=cnv_interval_brahman,hereford=cnv_interval_arsucd)

#make UpSet plot using the function from ComplexHeatmap
m = make_comb_mat(listInput_complxhtmap)

set_size(m)

comb_size(m)

UpSet(m)

tiff(filename = "FigFinal_Upset_basepair_resolution.tiff",width = 500,height = 300)
UpSet(m, pt_size = unit(5, "mm"), lwd = 3, comb_col = c("red", "blue", "black")[comb_degree(m)])
dev.off()

#On average, 0.5% of each cattle genome was covered by CNV regions (CNVRs) 
#mean((set_size(m)/2.7e9)*100)
#[1] 0.51061

#The majority of CNVRs (at least 76% from each assembly) were found to be unique to one assembly
# (comb_size(m)[1:3]/set_size(m))*100
# angus brahman hereford
# 100      010      001 
# 76.88974 82.24537 87.84588

# Angus vs Brahman
# comb_size(m)[4]
# 110 
# 1345463 
# Angus vs Hereford
# comb_size(m)[5]
# 101 
# 988764 

#region of Angus intersect with Brahman only
Angus_vs_Brahman_intersect_gr <- GenomicRanges::intersect(cnv_interval_angus,cnv_interval_brahman, ignore.strand = TRUE)

Angus_vs_Brahman_intersect_only_gr <- GenomicRanges::setdiff(Angus_vs_Brahman_intersect_gr,cnv_interval_arsucd, ignore.strand = TRUE)

Angus_vs_Brahman_intersect_only_df <- as(Angus_vs_Brahman_intersect_only_gr, "data.frame")

write_csv(Angus_vs_Brahman_intersect_only_df,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/CopyNumberVariation/20190530/suspicious_overlap_btwn_brahman_angus_only/Angus_vs_Brahman_intersect_only_df.csv")

# > sessionInfo()
# R version 3.5.3 (2019-03-11)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS Mojave 10.14.5
# 
# Matrix products: default
# BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_AU.UTF-8/en_AU.UTF-8/en_AU.UTF-8/C/en_AU.UTF-8/en_AU.UTF-8
# 
# attached base packages:
#   [1] grid      stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] circlize_0.4.6         ComplexHeatmap_2.1.0   GenomicFeatures_1.34.8 AnnotationDbi_1.44.0   Biobase_2.42.0        
# [6] GenomicRanges_1.34.0   GenomeInfoDb_1.18.2    IRanges_2.16.0         S4Vectors_0.20.1       BiocGenerics_0.28.0   
# [11] tidyr_0.8.3            readr_1.3.1            dplyr_0.8.1            ggplot2_3.1.1          UpSetR_1.4.0          
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.1                  lattice_0.20-38             prettyunits_1.0.2           png_0.1-7                  
# [5] Rsamtools_1.34.1            Biostrings_2.50.2           assertthat_0.2.1            digest_0.6.19              
# [9] R6_2.4.0                    plyr_1.8.4                  RSQLite_2.1.1               httr_1.4.0                 
# [13] pillar_1.4.1                GlobalOptions_0.1.0         zlibbioc_1.28.0             rlang_0.3.4                
# [17] progress_1.2.2              lazyeval_0.2.2              rstudioapi_0.10             blob_1.1.1                 
# [21] GetoptLong_0.1.7            Matrix_1.2-17               BiocParallel_1.16.6         stringr_1.4.0              
# [25] RCurl_1.95-4.12             bit_1.1-14                  biomaRt_2.38.0              munsell_0.5.0              
# [29] DelayedArray_0.8.0          compiler_3.5.3              rtracklayer_1.42.2          pkgconfig_2.0.2            
# [33] shape_1.4.4                 tidyselect_0.2.5            SummarizedExperiment_1.12.0 tibble_2.1.3               
# [37] gridExtra_2.3               GenomeInfoDbData_1.2.0      matrixStats_0.54.0          XML_3.98-1.20              
# [41] crayon_1.3.4                withr_2.1.2                 GenomicAlignments_1.18.1    bitops_1.0-6               
# [45] gtable_0.3.0                DBI_1.0.0                   magrittr_1.5                scales_1.0.0               
# [49] stringi_1.4.3               XVector_0.22.0              RColorBrewer_1.1-2          rjson_0.2.20               
# [53] tools_3.5.3                 bit64_0.9-7                 glue_1.3.1                  purrr_0.3.2                
# [57] hms_0.4.2                   yaml_2.2.0                  clue_0.3-57                 colorspace_1.4-1           
# [61] cluster_2.0.9               memoise_1.1.0         

#####extra
############################# Start from lifover results ###########################
# # path to CNV results
# dir2 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/CopyNumberVariation/20190530/cnvr_data/results_liftover/"
# 
# # angus
# path1 <- paste0(dir2,"angus_arsucd_coor")
# 
# angus.cnvrs.mapped <- read_tsv(path1,col_names = FALSE, col_types = "cddcdd")
# colnames(angus.cnvrs.mapped) <- c("chr","start","end","species","unsure1","unsure2")
# 
# #filter for size less than 1 mil
# angus.cnvrs.mapped_modi <- angus.cnvrs.mapped %>% mutate(size=end - start) %>% filter(size < 1e6)
# angus.cnvrs.mapped_modi$unique_name <- paste(angus.cnvrs.mapped_modi$chr,angus.cnvrs.mapped_modi$start,angus.cnvrs.mapped_modi$end,sep = "_")
# 
# angus.cnvrs.mapped_modi$strand <- rep("*",nrow(angus.cnvrs.mapped_modi))
# 
# cnv_interval_angus <- makeGRangesFromDataFrame(angus.cnvrs.mapped_modi, keep.extra.columns = TRUE, 
#                                                 seqnames.field="chr", start.field="start", 
#                                                 end.field="end", strand.field="strand")
# 
# #brahman
# path2 <- paste0(dir2,"brahman_arsucd_coor")
# 
# brahman.cnvrs.mapped <- read_tsv(path2,col_names = FALSE, col_types = "cddcdd")
# colnames(brahman.cnvrs.mapped) <- c("chr","start","end","species","unsure1","unsure2")
# 
# #filter for size less than 1 mil
# brahman.cnvrs.mapped_modi <- brahman.cnvrs.mapped %>% mutate(size=end - start) %>% filter(size < 1e6)
# brahman.cnvrs.mapped_modi$unique_name <- paste(brahman.cnvrs.mapped_modi$chr,brahman.cnvrs.mapped_modi$start,brahman.cnvrs.mapped_modi$end,sep = "_")
# 
# brahman.cnvrs.mapped_modi$strand <- rep("*",nrow(brahman.cnvrs.mapped_modi))
# 
# cnv_interval_brahman <- makeGRangesFromDataFrame(brahman.cnvrs.mapped_modi, keep.extra.columns = TRUE, 
#                                                  seqnames.field="chr", start.field="start", 
#                                                  end.field="end", strand.field="strand")
# 
# #arsucd
# path3 <- paste0(dir2,"arsucd.cnvrs_regions.bed")
# 
# arsucd.cnvrs.mapped <- read_tsv(path3,col_names = FALSE, col_types = "cddcdd")
# colnames(arsucd.cnvrs.mapped) <- c("chr","start","end","species","unsure1","unsure2")
# 
# #filter for size less than 1 mil
# arsucd.cnvrs.mapped_modi <- arsucd.cnvrs.mapped %>% mutate(size=end - start) %>% filter(size < 1e6)
# arsucd.cnvrs.mapped_modi$unique_name <- paste(arsucd.cnvrs.mapped_modi$chr,arsucd.cnvrs.mapped_modi$start,arsucd.cnvrs.mapped_modi$end,sep = "_")
# 
# arsucd.cnvrs.mapped_modi$strand <- rep("*",nrow(arsucd.cnvrs.mapped_modi))
# 
# cnv_interval_arsucd <- makeGRangesFromDataFrame(arsucd.cnvrs.mapped_modi, keep.extra.columns = TRUE, 
#                                                 seqnames.field="chr", start.field="start", 
#                                                 end.field="end", strand.field="strand")
# 
# #combine cnv granges as list
# listInput_complxhtmap <- list(angus=cnv_interval_angus,brahman=cnv_interval_brahman,hereford=cnv_interval_arsucd)
# 
# #make UpSet plot using the function from ComplexHeatmap
# # install.packages("remotes")
# # remotes::install_github("jokergoo/ComplexHeatmap")
# 
# m = make_comb_mat(listInput_complxhtmap)
# 
# set_size(m)
# 
# comb_size(m)
# 
# UpSet(m)
# 
# tiff(filename = "FigFinal_Upset_basepair_resolution.tiff",width = 500,height = 300)
# UpSet(m, pt_size = unit(5, "mm"), lwd = 3, comb_col = c("red", "blue", "black")[comb_degree(m)])
# dev.off()
