library(tidyverse)
library(rGMAP)
## to install rGMAP
# library(devtools)
# install_github("ningbioinfostruggling/rGMAP")

#Angus

## first need to load in the HiC matrix index from the HiC-Pro pipeline.
## to get the bin num of you desire chromosome
## can be found in phoenix:"/data/biohub/20180511_Ning_HiCow/Angus_splited_out/hic_results/matrix/Angus/raw/40000/"

dir1 <-"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/TADs/input_files_plot_TAD/angus/"
path1 <- paste0(dir1,"Angus_40000_abs.bed")

#abin <- './splited_analysis/Angus_40000_abs.bed'
abin <- path1
abinindex <- read_delim(path1, delim = '\t', col_names = F)

#names abinindex chr start end bin_number

chr <- "chr1"
startbin <- as.integer(as.vector(filter(abinindex, X1 == chr)[1,4]))
endbin <- as.integer(as.vector(filter(abinindex, X1 == chr) %>% .[nrow(.),4]))

## then get the HiC matrix with specific chr
## file can be found in phoenix:"/data/biohub/20180511_Ning_HiCow/Angus_splited_out/hic_results/matrix/Angus/iced/40000/"
path2 <- paste0(dir1,"Angus_40000.matrix")

amat <- read_delim(path2, delim = '\t', col_names = F) %>%
  filter(X1 >= startbin & X1 <= endbin & X2 >= startbin & X2 <= endbin)

#names amat bin1 bin2 normalised_count

## call TADs
resl <- 40000 ## resolution
atad <- rGMAP(amat, index_file = abin, resl = resl)

startcor <- 152000000 ## the start coordinate to plot
endcor <- 157000000 ## the end to plot

## get the bin number to plot from the coordinates
sb <- startcor/resl - startbin + 1
eb <- endcor/resl - startbin + 1

## get the plot
ap <- plotdom(amat, atad$hierTads, start_bin = sb, end_bin = eb, cthr = 20, resl = resl)

## plot

## CNVs <- x/10^6  here can be the CNV start and end
## CNVe <- y/10^6

tiff("Angus_polled_LL_plot_wLine.tiff")
g <- ap$p2  +  ggtitle(paste("TADs and sub-TADs of Angus at 40 kb resolution"))
g <- g + theme_bw() + xlab("Chromosome 1 coordinate (Mb)") + ylab("Count")
g <- g + geom_vline(xintercept = 154851456/1e6, linetype = "dotted", color = "blue", lwd = 1.2)
g <- g + theme(plot.title = element_text(hjust = 0.5))
g
dev.off()

# ap$p2  +  ggtitle(paste("TADs and sub-TADs of Angus at 40 kb resolution")) +  
#   theme_bw() + xlab("Chromosome 1 coordinate (Mb)") + geom_vline(xintercept = 154852238/1e6, linetype = "dotted", color = "blue", lwd = 1.2)
## +  geom_vline(xintercept = c(CNVs, CNVe), linetype = "dotted" )

#Brahman

## first need to load in the HiC matrix index from the HiC-Pro pipeline.
## to get the bin num of you desire chromosome
## can be found in phoenix:"/data/biohub/20180511_Ning_HiCow/Angus_splited_out/hic_results/matrix/Angus/raw/40000/"

dir1 <-"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/TADs/input_files_plot_TAD/brahman/"
path1 <- paste0(dir1,"Brahman_40000_abs.bed")

#abin <- './splited_analysis/Angus_40000_abs.bed'
abin <- path1
abinindex <- read_delim(path1, delim = '\t', col_names = F)

#names abinindex chr start end bin_number

chr <- "chr1"
startbin <- as.integer(as.vector(filter(abinindex, X1 == chr)[1,4]))
endbin <- as.integer(as.vector(filter(abinindex, X1 == chr) %>% .[nrow(.),4]))

## then get the HiC matrix with specific chr
## file can be found in phoenix:"/data/biohub/20180511_Ning_HiCow/Angus_splited_out/hic_results/matrix/Angus/iced/40000/"
path2 <- paste0(dir1,"Brahman_40000_iced.matrix")

amat <- read_delim(path2, delim = '\t', col_names = F) %>%
  filter(X1 >= startbin & X1 <= endbin & X2 >= startbin & X2 <= endbin)

#names amat bin1 bin2 normalised_count

## call TADs
resl <- 40000 ## resolution
atad <- rGMAP(amat, index_file = abin, resl = resl)

startcor <- 152000000 ## the start coordinate to plot
endcor <- 157000000 ## the end to plot

## get the bin number to plot from the coordinates
sb <- startcor/resl - startbin + 1
eb <- endcor/resl - startbin + 1

## get the plot
ap <- plotdom(amat, atad$hierTads, start_bin = sb, end_bin = eb, cthr = 20, resl = resl)

## plot

## CNVs <- x/10^6  here can be the CNV start and end
## CNVe <- y/10^6

tiff("Brahman_polled_LL_plot_wLine.tiff")
g <- ap$p2 + ggtitle(paste("TADs and sub-TADs of Brahman at 40k resolution"))
g <- g + theme_bw() + xlab("Chromosome 1 coordinate (Mb)") + ylab("Count")
g <- g + geom_vline(xintercept = 154748363/1e6, linetype = "dotted", color = "blue", lwd = 1.2)
g <- g + theme(plot.title = element_text(hjust = 0.5))
g
dev.off()

# ap$p2 + 
#   ggtitle(paste("TADs and sub-TADs of Brahman at 40k resolution")) +
#   theme_bw() +
#   xlab("Coordinate (Mb)") + geom_vline(xintercept = 154749145/1e6, linetype = "dotted", color = "blue", lwd = 1.2)
## +  geom_vline(xintercept = c(CNVs, CNVe), linetype = "dotted" )
