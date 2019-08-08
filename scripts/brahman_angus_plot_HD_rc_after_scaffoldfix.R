#------------------------------------------------------
# Program name: brahman_angus_plot_HD_rc_after_scaffoldfix.R
# Objective: after reodering scaffolds now check how they best fit 
#          each chr HD and rc
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(readr)
library(dplyr)
library(ggplot2)

##### Dam #####
# this is after first round fixing; reading HD probes
dir1 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/alignAndOrderSnp_result/"
#path1 <- paste0(dir1,"dam_hd_rc_dom_v2.HDProbes.tab")
#this is after bionano cut fix
path1 <- paste0(dir1,"bostaurus_brahma_bionano_NCBI_full_corrected.HDProbes.tab")

dam_hd_rc_dom_v2_HDProbes <- read_tsv(path1,col_names = FALSE)
names(dam_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

#test
#dam_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 1) %>% filter(scaffold == "1")

#create a vec to store chr1 to chr29
chr_vec <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20",
             "21","22","23","24","25","26","27","28","29")

#loop thro each chr and ggplot it
for (i in 1:length(chr_vec)){
  DF <- dam_hd_rc_dom_v2_HDProbes %>% filter(chromosome == as.numeric(chr_vec[i])) %>% filter(scaffold == chr_vec[i])
  #plot(DF$align_pos,DF$chr_pos)
  dom_chr <- paste0("HD chr",chr_vec[i])
  brahman_chr <- paste0("Brahman chr",chr_vec[i])
  brahman_chr_pdf <- paste0("Brahman chr",chr_vec[i],"_HD_.pdf")
  g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab(brahman_chr) + ylab(dom_chr)
  ggsave(brahman_chr_pdf,device = "pdf")
}

# this is after first round fixing; reading HD probes
# path2 <- paste0(dir1,"dam_hd_rc_dom_v2.rcmap.tab")
#this is after bionano cut fix
path2 <- paste0(dir1,"bostaurus_brahma_bionano_NCBI_full_corrected.rcmap.tab")

dam_hd_rc_dom_v2_rcmap <- read_tsv(path2,col_names = FALSE)
names(dam_hd_rc_dom_v2_rcmap) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

#loop thro each chr and ggplot it
for (i in 1:length(chr_vec)){
  DF <- dam_hd_rc_dom_v2_rcmap %>% filter(chromosome == as.numeric(chr_vec[i])) %>% filter(scaffold == chr_vec[i])
  dom_chr <- paste0("rc chr",chr_vec[i])
  brahman_chr <- paste0("Brahman chr",chr_vec[i])
  brahman_chr_pdf <- paste0("Brahman chr",chr_vec[i],"_rcmap_.pdf")
  g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab(brahman_chr) + ylab(dom_chr)
  ggsave(brahman_chr_pdf,device = "pdf")
}

##### Sire #####
# this is after first round fixing; reading HD probes
dir1 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/alignAndOrderSnp_result/"
# path1 <- paste0(dir1,"sire_hd_rc_dom_v2.HDProbes.tab")
#this is after bionano cut fix
path1 <- paste0(dir1,"bostaurus_angus_bionano_NCBI_full_corrected.HDProbes.tab")

sire_hd_rc_dom_v2_HDProbes <- read_tsv(path1,col_names = FALSE)
names(sire_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

#test
#sire_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 1) %>% filter(scaffold == "1")

#create a vec to store chr1 to chr29
chr_vec <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20",
             "21","22","23","24","25","26","27","28","29")

#loop thro each chr and ggplot it
for (i in 1:length(chr_vec)){
  DF <- sire_hd_rc_dom_v2_HDProbes %>% filter(chromosome == as.numeric(chr_vec[i])) %>% filter(scaffold == chr_vec[i])
  #plot(DF$align_pos,DF$chr_pos)
  dom_chr <- paste0("HD chr",chr_vec[i])
  brahman_chr <- paste0("Angus chr",chr_vec[i])
  brahman_chr_pdf <- paste0("Angus chr",chr_vec[i],"_HD_.pdf")
  g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab(brahman_chr) + ylab(dom_chr)
  ggsave(brahman_chr_pdf,device = "pdf")
}

# this is after first round fixing; reading HD probes
# path2 <- paste0(dir1,"sire_hd_rc_dom_v2.rcmap.tab")
#this is after bionano cut fix
path2 <- paste0(dir1,"bostaurus_angus_bionano_NCBI_full_corrected.rcmap.tab")

sire_hd_rc_dom_v2_rcmap <- read_tsv(path2,col_names = FALSE)
names(sire_hd_rc_dom_v2_rcmap) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

#test reading agp to draw lines where joins occur in the scaffolds
scaffold_order_book <- read_tsv("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/agp/scaffold_order_book_sire.agp.tsv",
                                     col_names = FALSE)
names(scaffold_order_book) <- c("object","object_beg","object_end","part_number","component_type","component_id","component_beg","component_end","orientation")

#loop thro each chr and ggplot it
for (i in 1:length(chr_vec)){
  DF <- sire_hd_rc_dom_v2_rcmap %>% filter(chromosome == as.numeric(chr_vec[i])) %>% filter(scaffold == chr_vec[i])
  
  #vline
  # vline <- scaffold_order_book %>% filter(object == as.numeric(chr_vec[i])) %>% select(object_beg)
  # vline_vec <- vline$object_beg
  
  dom_chr <- paste0("rc chr",chr_vec[i])
  brahman_chr <- paste0("Angus chr",chr_vec[i])
  brahman_chr_pdf <- paste0("Angus chr",chr_vec[i],"_rcmap_.pdf")
  # g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab(brahman_chr) + ylab(dom_chr) + 
  #   geom_vline(xintercept = vline_vec, linetype="dotted")
  g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab(brahman_chr) + ylab(dom_chr) 
  ggsave(brahman_chr_pdf,device = "pdf")
}

##### Dominette #####
# this is after first round fixing; reading HD probes
dir1 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/alignAndOrderSnp_result/"
path1 <- paste0(dir1,"ARS-UCD1.0.25.HDProbes.tab")

dom_hd_rc_dom_v2_HDProbes <- read_tsv(path1,col_names = FALSE)
names(dom_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

#test
#dom_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 1) %>% filter(scaffold == "1")

#create a vec to store chr1 to chr29
chr_vec <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20",
             "21","22","23","24","25","26","27","28","29")

#loop thro each chr and ggplot it
for (i in 1:length(chr_vec)){
  DF <- dom_hd_rc_dom_v2_HDProbes %>% filter(chromosome == as.numeric(chr_vec[i])) %>% filter(scaffold == chr_vec[i])
  #plot(DF$align_pos,DF$chr_pos)
  dom_chr <- paste0("HD chr",chr_vec[i])
  brahman_chr <- paste0("Dominette chr",chr_vec[i])
  brahman_chr_pdf <- paste0("Dominette chr",chr_vec[i],"_HD_.pdf")
  g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab(brahman_chr) + ylab(dom_chr)
  ggsave(brahman_chr_pdf,device = "pdf")
}

# this is after first round fixing; reading HD probes
path2 <- paste0(dir1,"ARS-UCD1.0.25.rcmap.tab")

dom_hd_rc_dom_v2_rcmap <- read_tsv(path2,col_names = FALSE)
names(dom_hd_rc_dom_v2_rcmap) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

#loop thro each chr and ggplot it
for (i in 1:length(chr_vec)){
  DF <- dom_hd_rc_dom_v2_rcmap %>% filter(chromosome == as.numeric(chr_vec[i])) %>% filter(scaffold == chr_vec[i])
  dom_chr <- paste0("rc chr",chr_vec[i])
  brahman_chr <- paste0("Dominette chr",chr_vec[i])
  brahman_chr_pdf <- paste0("Dominette chr",chr_vec[i],"_rcmap_.pdf")
  g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab(brahman_chr) + ylab(dom_chr)
  ggsave(brahman_chr_pdf,device = "pdf")
}

