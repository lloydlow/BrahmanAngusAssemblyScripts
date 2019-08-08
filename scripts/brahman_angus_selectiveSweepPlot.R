#------------------------------------------------------
# Program name: brahman_angus_selectiveSweepPlot.R
# Objective: zooming in details on candidate sweep
#           
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
require(scales)

#chr 7
# load("Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_7.RData")
load("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Brahman_SevenBreeds/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_7.RData")

# for_plot_chr7 <- combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_7 %>% filter(POS >= 66100000) %>% 
#   filter(POS < 66200000) %>% select(POS,an_ALT_proportion,br_ALT_proportion)
for_plot_chr7 <- combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_2_7 %>% filter(POS >= 66100000) %>% 
  filter(POS < 66200000) %>% dplyr::select(POS,an_ALT_proportion,br_ALT_proportion)

#rename var
names(for_plot_chr7) <- c("Position","Angus","Brahman")

g <- ggplot(for_plot_chr7, aes(x=Position, y = value, color = breed))
g <- g + geom_point(aes(y = Angus, col = "Angus")) + 
  geom_point(aes(y = Brahman, col = "Brahman")) 
g <- g + ylab("Proportion of alternate allele") + xlab("Position on Brahman chr 7")
g <- g + scale_x_continuous(labels = comma)
g

#chr 19
load("Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_19.RData")

for_plot_chr19 <- combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_19 %>% filter(POS >= 27700000) %>% 
  filter(POS < 27800000) %>% select(POS,an_ALT_proportion,br_ALT_proportion)

#rename var
names(for_plot_chr19) <- c("Position","Angus","Brahman")

g <- ggplot(for_plot_chr19, aes(x=Position, y = value, color = breed))
g <- g + geom_point(aes(y = Angus, col = "Angus")) + 
  geom_point(aes(y = Brahman, col = "Brahman")) 
g <- g + ylab("Proportion of alternate allele") + xlab("Position on Brahman chr 19")
g <- g + scale_x_continuous(labels = comma)
g

#chr 24
load("Assembly_version/final_to_correct_20180905/Annotation/EBI/combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_24.RData")

for_plot_chr24 <- combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_24 %>% filter(POS >= 400000) %>% 
  filter(POS < 800000) %>% select(POS,an_ALT_proportion,br_ALT_proportion)

#rename var
names(for_plot_chr24) <- c("Position","Angus","Brahman")

g <- ggplot(for_plot_chr24, aes(x=Position, y = value, color = breed))
g <- g + geom_point(aes(y = Angus, col = "Angus")) + 
  geom_point(aes(y = Brahman, col = "Brahman")) 
g <- g + ylab("Proportion of alternate allele") + xlab("Position on Brahman chr 24")
g <- g + scale_x_continuous(labels = comma)
g



