#------------------------------------------------------
# Program name: brahman_angus_gapCount_ungappedContigLength_Figure.R
# Objective: All gap count and ungapped contig length comparisons
#           for Brahman and Angus paper
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(readr)
library(dplyr)
library(ggplot2)
library(easyGgplot2)

#calculate N50 function #revised with arranging the df
N50 <- function(assembly_df){
  assembly_df <- assembly_df %>% arrange(desc(length))
  vect <- c()
  for (i in 1:nrow(assembly_df)){
    N50length <- 0.5*sum(as.numeric(assembly_df$length))
    vect <- c(vect,assembly_df$length[i])
    cumsumlength <- sum(as.numeric(vect))
    if (cumsumlength >= N50length){
      print(assembly_df$length[i])
      break
    } 
  }
}

########## checking gaps by chromosome by species #############
# path to all coor results
dir1 <- "/Users/lloyd/Documents/lloyd_2017/Research/RiverBuffalo/buffalo_NextGenAssembly/buffalo_NextGenAssembly/gap_genome_analysis/species/gap_coor/"

# reading bostaurus_angus_bionano_NCBI_full_corrected_gapfill_arrow_fil_chronly.coor
path1 <- paste0(dir1,"bostaurus_angus_bionano_NCBI_full_corrected_gapfill_arrow_fil_chronly.coor")

cattle_angus_coor <- read_tsv(path1,col_names = FALSE)
names(cattle_angus_coor) <- c("chromosome","start","stop")

cattle_angus_coor$species <- 'UOA_Angus_1'

cattle_angus_coor$chromosome[cattle_angus_coor$chromosome == "Y"] <- "X/Y"

# reading bostaurus_brahma_bionano_NCBI_full_corrected_gapfill_arrow_fil_withM_chronly.coor
path2 <- paste0(dir1,"bostaurus_brahma_bionano_NCBI_full_corrected_gapfill_arrow_fil_withM_chronly.coor")

cattle_brahman_coor <- read_tsv(path2,col_names = FALSE)
names(cattle_brahman_coor) <- c("chromosome","start","stop")

cattle_brahman_coor$species <- 'UOA_Brahman_1'

cattle_brahman_coor$chromosome[cattle_brahman_coor$chromosome == "X"] <- "X/Y"

# reading cattle_arsucd_chr_only.coor
path3 <- paste0(dir1,"cattle_arsucd_chr_only.coor")

cattle_arsucd12_coor <- read_tsv(path3,col_names = FALSE)
names(cattle_arsucd12_coor) <- c("chromosome","start","stop")

cattle_arsucd12_coor$species <- 'ARS-UCD1.2'

cattle_arsucd12_coor$chromosome[cattle_arsucd12_coor$chromosome == "X"] <- "X/Y"

# reading water_buffalo_20180219_gapf_noMito_arrowRename4_pilon_chr_only.coor
path4 <- paste0(dir1,"water_buffalo_20180219_gapf_noMito_arrowRename4_pilon_chr_only.coor")

buffalo_coor <- read_tsv(path4,col_names = FALSE,col_types = list(col_character(), col_double(),col_double()))
names(buffalo_coor) <- c("chromosome","start","stop")

buffalo_coor$species <- 'UOA_WB_1'

buffalo_coor$chromosome[buffalo_coor$chromosome == "X"] <- "X/Y"

# reading human_chr_only.coor
path5 <- paste0(dir1,"human_chr_only.coor")

human_coor <- read_tsv(path5,col_names = FALSE,col_types = list(col_character(), col_double(),col_double()))
names(human_coor) <- c("chromosome","start","stop")

human_coor$species <- 'GRCh38'

human_coor$chromosome[human_coor$chromosome == "X"] <- "X/Y"

all_spp_gap <- rbind(cattle_angus_coor,cattle_brahman_coor,cattle_arsucd12_coor,buffalo_coor,human_coor)

all_spp_gap_group <- all_spp_gap %>% dplyr::group_by(species,chromosome) %>% dplyr::summarise(n =n())

#order chr
ordered <- c("1","2","3","4","5","6","7","8","9","10","11",
             "12","13","14","15","16","17","18","19","20","21",
             "22","23","24","25","26","27","28","29","X/Y")
all_spp_gap_group$chromosome <- factor(all_spp_gap_group$chromosome, 
                                       levels = ordered)

#order spp
order_spp <- c("UOA_Angus_1","UOA_Brahman_1","ARS-UCD1.2","UOA_WB_1","GRCh38")
all_spp_gap_group$species <- factor(all_spp_gap_group$species, 
                                    levels = order_spp)

tiff(filename = "FigFinal_NumberOfGapsBySpeciesBrahmanAngus.tiff",width = 800, height = 600)
g <- ggplot(data = all_spp_gap_group, aes(x = species, y = n, fill = chromosome)) + 
  geom_bar(stat = "identity") + ylab("Number of gaps") + xlab("Genome")
g <- g + ggtitle("Number of gaps per chromosome by genome assembly")
g <- g + theme(plot.title = element_text(hjust = 0.5),axis.text=element_text(size=13),
               axis.title=element_text(size=17))
g
dev.off()

# cattle_arsucd12_coor$chromosome <- factor(cattle_arsucd12_coor$chromosome,
#                                           levels = ordered)
#
# ggplot2.histogram(data=cattle_arsucd12_coor, xName= 'start', xtitle="Position",
#                   groupName='chromosome', legendPosition="right",
#                   faceting=TRUE, facetingVarNames="chromosome",
#                   binwidth = 0.1e6,yShowTitle=FALSE,yShowTickLabel=FALSE,
#                   hideAxisTicks=TRUE)
#
# cattle_brahman_coor$chromosome <- factor(cattle_brahman_coor$chromosome,levels = ordered)
# 
# ggplot2.histogram(data=cattle_brahman_coor, xName= 'start', xtitle="Position",
#                   groupName='chromosome', legendPosition="right",
#                   faceting=TRUE, facetingVarNames="chromosome",
#                   binwidth = 0.1e6,yShowTitle=FALSE,yShowTickLabel=FALSE,
#                   hideAxisTicks=TRUE)
# 
# cattle_angus_coor$chromosome <- factor(cattle_angus_coor$chromosome,levels = ordered)
# 
# ggplot2.histogram(data=cattle_angus_coor, xName= 'start', xtitle="Position",
#                   groupName='chromosome', legendPosition="right",
#                   faceting=TRUE, facetingVarNames="chromosome",
#                   binwidth = 0.1e6,yShowTitle=FALSE,yShowTickLabel=FALSE,
#                   hideAxisTicks=TRUE)

########## checking ungapped contig length in various assemblies #############

# path to all gapLength results
dir2 <- "/Users/lloyd/Documents/lloyd_2017/Research/RiverBuffalo/buffalo_NextGenAssembly/buffalo_NextGenAssembly/reordering_HiRise_result/fasta_file_order_checkBases_N/"

# reading bostaurus_angus_bionano_NCBI_full_corrected_gapfill_arrow_fil_chronly_ungapped_No_Ns.rls
path1 <- paste0(dir2,"bostaurus_angus_bionano_NCBI_full_corrected_gapfill_arrow_fil_chronly_ungapped_No_Ns.rls")

angus_chr_only_ungapped_No_Ns <- read_tsv(path1,col_names = FALSE)
names(angus_chr_only_ungapped_No_Ns) <- c("scaffold","nothing","gap","length","perc_gap")

angus_chr_only_ungapped_No_Ns <- angus_chr_only_ungapped_No_Ns %>%
  dplyr::select(scaffold,gap,length) %>% arrange(desc(length))

# reading bostaurus_brahma_bionano_NCBI_full_corrected_gapfill_arrow_fil_withM_chronly_ungapped_No_Ns.rls
path2 <- paste0(dir2,"bostaurus_brahma_bionano_NCBI_full_corrected_gapfill_arrow_fil_withM_chronly_ungapped_No_Ns.rls")

brahman_chr_only_ungapped_No_Ns <- read_tsv(path2,col_names = FALSE)
names(brahman_chr_only_ungapped_No_Ns) <- c("scaffold","nothing","gap","length","perc_gap")

brahman_chr_only_ungapped_No_Ns <- brahman_chr_only_ungapped_No_Ns %>%
  dplyr::select(scaffold,gap,length) %>% arrange(desc(length))

# reading cattle_arsucd_chr_only_ungapped_No_Ns.rls
path3 <- paste0(dir2,"cattle_arsucd_chr_only_ungapped_No_Ns.rls")

cattle_arsucd_chr_only_ungapped_No_Ns <- read_tsv(path3,col_names = FALSE)
names(cattle_arsucd_chr_only_ungapped_No_Ns) <- c("scaffold","nothing","gap","length","perc_gap")

cattle_arsucd_chr_only_ungapped_No_Ns <- cattle_arsucd_chr_only_ungapped_No_Ns %>%
  dplyr::select(scaffold,gap,length) %>% arrange(desc(length))

# reading water_buffalo_20180219_gapf_noMito_arrowRename4_pilon_chr_only_ungapped_No_Ns.rls
path4 <- paste0(dir2,"water_buffalo_20180219_gapf_noMito_arrowRename4_pilon_chr_only_ungapped_No_Ns.rls")

water_buffalo_20180219_gapf_noMito_arrowRename4_pilon_chr_only_ungapped_No_Ns <- read_tsv(path4,col_names = FALSE)
names(water_buffalo_20180219_gapf_noMito_arrowRename4_pilon_chr_only_ungapped_No_Ns) <- c("scaffold","nothing","gap","length","perc_gap")

water_buffalo_20180219_gapf_noMito_arrowRename4_pilon_chr_only_ungapped_No_Ns <- water_buffalo_20180219_gapf_noMito_arrowRename4_pilon_chr_only_ungapped_No_Ns %>%
  dplyr::select(scaffold,gap,length) %>% arrange(desc(length))

#give each DF a species variable
water_buffalo_20180219_gapf_noMito_arrowRename4_pilon_chr_only_ungapped_No_Ns$species <- 'UOA_WB_1'
cattle_arsucd_chr_only_ungapped_No_Ns$species <- 'ARS-UCD1.2'
angus_chr_only_ungapped_No_Ns$species <- 'UOA_Angus_1'
brahman_chr_only_ungapped_No_Ns$species <- 'UOA_Brahman_1'

spp_df <- rbind(angus_chr_only_ungapped_No_Ns,
                brahman_chr_only_ungapped_No_Ns,
                water_buffalo_20180219_gapf_noMito_arrowRename4_pilon_chr_only_ungapped_No_Ns,
                cattle_arsucd_chr_only_ungapped_No_Ns)

#order spp
order_spp <- c("UOA_Angus_1","UOA_Brahman_1","ARS-UCD1.2","UOA_WB_1")
spp_df$species <- factor(spp_df$species, levels = order_spp)

#density plot
ggplot(spp_df, aes(log10(length), fill = species)) + geom_density(alpha = 0.2) +
  xlab(expression(log[10]*" ungapped contig length")) + ggtitle("Density plots of ungapped contig length distribution") + 
  theme(plot.title = element_text(hjust = 0.5))

tiff(filename = "UngappedContigLengthBrahmanAngus.tiff",width = 600, height = 600, compression = 'none')
g <- ggplot(spp_df, aes(x=as.factor(species),y=(length/1e6)))
g <- g + geom_boxplot(fill="slateblue", alpha=0.5)
g <- g + xlab("Genome") + ylab("Un-gapped contig length (Mbp)")
g <- g + ggtitle("Un-gapped contig length distribution by genome assembly") 
g <-  g + theme(plot.title = element_text(hjust = 0.5),axis.text=element_text(size=15),axis.title=element_text(size=17))
g
dev.off()
