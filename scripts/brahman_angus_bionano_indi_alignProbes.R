#------------------------------------------------------
# Program name: brahman_angus_bionano_indi_alignProbes.R
# Objective: this is to test arrangement of chr by probes
#           one by one
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------


library(readr)
library(dplyr)
library(ggplot2)

##### Dam #####
# this is after first round fixing; reading rc probes
dir1 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/alignAndOrderSnp_result/dam_chr/"
path1 <- paste0(dir1,"1.rcmap.tab")

dam_hd_rc_dom_v2_HDProbes <- read_tsv(path1,col_names = FALSE)
names(dam_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- dam_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 1) %>% filter(scaffold == 1)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("brahman_chr") + ylab("HDprobe")
g

path2 <- paste0(dir1,"2.rcmap.tab")

dam_hd_rc_dom_v2_HDProbes <- read_tsv(path2,col_names = FALSE)
names(dam_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- dam_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 2) %>% filter(scaffold == 2)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("brahman_chr") + ylab("HDprobe")
g

path3 <- paste0(dir1,"3.rcmap.tab")

dam_hd_rc_dom_v2_HDProbes <- read_tsv(path3,col_names = FALSE)
names(dam_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- dam_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 3) %>% filter(scaffold == 3)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("brahman_chr") + ylab("HDprobe") +
  geom_vline(xintercept=c(3483180,3627691,12963347,21015764,21015764+1917829,21015764+1917829+6791776,21015764+86060879,21015764+1917829+97845705), color = "blue")
g

path4 <- paste0(dir1,"4.rcmap.tab")

dam_hd_rc_dom_v2_HDProbes <- read_tsv(path4,col_names = FALSE)
names(dam_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- dam_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 4) %>% filter(scaffold == 4)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("brahman_chr") + ylab("HDprobe") +
  geom_vline(xintercept=c(6929177,6929177+572272,6929177+572272+4120589,6929177+572272+17595171), color = "blue")
g

path5 <- paste0(dir1,"5.rcmap.tab")

dam_hd_rc_dom_v2_HDProbes <- read_tsv(path5,col_names = FALSE)
names(dam_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- dam_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 5) %>% filter(scaffold == 5)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("brahman_chr") + ylab("HDprobe") +
  geom_vline(xintercept=c(24091847,24091847+1175520,24091847+11636309,24091847+11636309+859718,
                          24091847+11636309+859718+53711729,24091847+11636309+859718+53711729+13159056,
                          24091847+11636309+859718+53711729+13159056+15125639), color = "blue")
g

path6 <- paste0(dir1,"6.rcmap.tab")

dam_hd_rc_dom_v2_HDProbes <- read_tsv(path6,col_names = FALSE)
names(dam_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- dam_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 6) %>% filter(scaffold == 6)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("brahman_chr") + ylab("HDprobe") +
  geom_vline(xintercept=c(1161469+4936806+29117862+49688823, 1161469+4936806+29117862+49688823+31945074), color = "blue")
g

path7 <- paste0(dir1,"7.rcmap.tab")

dam_hd_rc_dom_v2_HDProbes <- read_tsv(path7,col_names = FALSE)
names(dam_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- dam_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 7) %>% filter(scaffold == 7)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("brahman_chr") + ylab("HDprobe") +
  geom_vline(xintercept=c(31743416), color = "blue")
g

path8 <- paste0(dir1,"8.rcmap.tab")

dam_hd_rc_dom_v2_HDProbes <- read_tsv(path8,col_names = FALSE)
names(dam_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- dam_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 8) %>% filter(scaffold == 8)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("brahman_chr") + ylab("HDprobe") +
  geom_vline(xintercept=c(19338071+39527988+11229811,19338071+39527988+11229811+16545126), color = "blue")
g

path9 <- paste0(dir1,"9.rcmap.tab")

dam_hd_rc_dom_v2_HDProbes <- read_tsv(path9,col_names = FALSE)
names(dam_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- dam_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 9) %>% filter(scaffold == 9)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("brahman_chr") + ylab("HDprobe") +
  geom_vline(xintercept=c(53370553), color = "blue")
g

path10 <- paste0(dir1,"10.rcmap.tab")

dam_hd_rc_dom_v2_HDProbes <- read_tsv(path10,col_names = FALSE)
names(dam_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- dam_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 10) %>% filter(scaffold == 10)
#93210416
g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("brahman_chr10") + ylab("rcprobe") +
  geom_vline(xintercept=cumsum(c(16358552,73796,507344,630400,193765,476635,1843297,4890602,61294,29439823,22902169,15655560,6176548,3937383)), color = "blue")
g

path11 <- paste0(dir1,"11.rcmap.tab")

dam_hd_rc_dom_v2_HDProbes <- read_tsv(path11,col_names = FALSE)
names(dam_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- dam_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 11) %>% filter(scaffold == 11)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("brahman_chr") + ylab("HDprobe") +
  geom_vline(xintercept=c(17563123), color = "blue")
g

path12 <- paste0(dir1,"12.rcmap.tab")

dam_hd_rc_dom_v2_HDProbes <- read_tsv(path12,col_names = FALSE)
names(dam_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- dam_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 12) %>% filter(scaffold == 12)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("brahman_chr") + ylab("HDprobe") +
  geom_vline(xintercept=c(7193955,7193955+43039549), color = "blue")
g

path13 <- paste0(dir1,"13.rcmap.tab")

dam_hd_rc_dom_v2_HDProbes <- read_tsv(path13,col_names = FALSE)
names(dam_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- dam_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 13) %>% filter(scaffold == 13)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("brahman_chr") + ylab("HDprobe") +
  geom_vline(xintercept=c(1005886), color = "blue")
g

path14 <- paste0(dir1,"14.rcmap.tab")

dam_hd_rc_dom_v2_HDProbes <- read_tsv(path14,col_names = FALSE)
names(dam_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- dam_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 14) %>% filter(scaffold == 14)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("brahman_chr") + ylab("HDprobe") +
  geom_vline(xintercept=c(14557452+23445275+21581435+21532498), color = "blue")
g

path15 <- paste0(dir1,"15.rcmap.tab")

dam_hd_rc_dom_v2_HDProbes <- read_tsv(path15,col_names = FALSE)
names(dam_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- dam_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 15) %>% filter(scaffold == 15)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("brahman_chr") + ylab("HDprobe") +
  geom_vline(xintercept=c(19575447,19575447+33167117,19575447+33167117+31513481), color = "blue")
g

path16 <- paste0(dir1,"16.rcmap.tab")

dam_hd_rc_dom_v2_HDProbes <- read_tsv(path16,col_names = FALSE)
names(dam_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- dam_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 16) %>% filter(scaffold == 16)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("brahman_chr") + ylab("HDprobe") +
  geom_vline(xintercept=c(6314326), color = "blue")
g

path17 <- paste0(dir1,"17.rcmap.tab")

dam_hd_rc_dom_v2_HDProbes <- read_tsv(path17,col_names = FALSE)
names(dam_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- dam_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 17) %>% filter(scaffold == 17)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("brahman_chr") + ylab("HDprobe") +
  geom_vline(xintercept=c(13950234), color = "blue")
g

path18 <- paste0(dir1,"18.rcmap.tab")

dam_hd_rc_dom_v2_HDProbes <- read_tsv(path18,col_names = FALSE)
names(dam_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- dam_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 18) %>% filter(scaffold == 18)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("brahman_chr") + ylab("HDprobe") +
  geom_vline(xintercept=c(46554213,46554213+770152), color = "blue")
g

path19 <- paste0(dir1,"19.rcmap.tab")

dam_hd_rc_dom_v2_HDProbes <- read_tsv(path19,col_names = FALSE)
names(dam_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- dam_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 19) %>% filter(scaffold == 19)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("brahman_chr") + ylab("HDprobe") +
  geom_vline(xintercept=c(7221563), color = "blue")
g

path20 <- paste0(dir1,"20.rcmap.tab")

dam_hd_rc_dom_v2_HDProbes <- read_tsv(path20,col_names = FALSE)
names(dam_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- dam_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 20) %>% filter(scaffold == 20)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("brahman_chr") + ylab("HDprobe") +
  geom_vline(xintercept=c(44969044), color = "blue")
g

path21 <- paste0(dir1,"21.rcmap.tab")

dam_hd_rc_dom_v2_HDProbes <- read_tsv(path21,col_names = FALSE)
names(dam_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- dam_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 21) %>% filter(scaffold == 21)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("brahman_chr") + ylab("HDprobe") +
  geom_vline(xintercept=c(2555379), color = "blue")
g

path22 <- paste0(dir1,"22.rcmap.tab")

dam_hd_rc_dom_v2_HDProbes <- read_tsv(path22,col_names = FALSE)
names(dam_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- dam_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 22) %>% filter(scaffold == 22)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("brahman_chr") + ylab("HDprobe") +
  geom_vline(xintercept=c(6879691), color = "blue")
g

path23 <- paste0(dir1,"23.rcmap.tab")

dam_hd_rc_dom_v2_HDProbes <- read_tsv(path23,col_names = FALSE)
names(dam_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- dam_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 23) %>% filter(scaffold == 23)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("brahman_chr") + ylab("HDprobe") +
  geom_vline(xintercept=c(35185741), color = "blue")
g

path24 <- paste0(dir1,"24.rcmap.tab")

dam_hd_rc_dom_v2_HDProbes <- read_tsv(path24,col_names = FALSE)
names(dam_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- dam_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 24) %>% filter(scaffold == 24)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("brahman_chr") + ylab("HDprobe") +
  geom_vline(xintercept=c(301221), color = "blue")
g

path25 <- paste0(dir1,"25.rcmap.tab")

dam_hd_rc_dom_v2_HDProbes <- read_tsv(path25,col_names = FALSE)
names(dam_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- dam_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 25) %>% filter(scaffold == 25)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("brahman_chr") + ylab("HDprobe") +
  geom_vline(xintercept=c(183242), color = "blue")
g

path26 <- paste0(dir1,"26.rcmap.tab")

dam_hd_rc_dom_v2_HDProbes <- read_tsv(path26,col_names = FALSE)
names(dam_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- dam_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 26) %>% filter(scaffold == 26)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("brahman_chr") + ylab("HDprobe") +
  geom_vline(xintercept=c(3659221), color = "blue")
g

path27 <- paste0(dir1,"27.rcmap.tab")

dam_hd_rc_dom_v2_HDProbes <- read_tsv(path27,col_names = FALSE)
names(dam_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- dam_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 27) %>% filter(scaffold == 27)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("brahman_chr") + ylab("HDprobe") +
  geom_vline(xintercept=c(6550254), color = "blue")
g

path28 <- paste0(dir1,"28.rcmap.tab")

dam_hd_rc_dom_v2_HDProbes <- read_tsv(path28,col_names = FALSE)
names(dam_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- dam_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 28) %>% filter(scaffold == 28)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("brahman_chr") + ylab("HDprobe") +
  geom_vline(xintercept=c(24278565), color = "blue")
g

path29 <- paste0(dir1,"29.rcmap.tab")

dam_hd_rc_dom_v2_HDProbes <- read_tsv(path29,col_names = FALSE)
names(dam_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- dam_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 29) %>% filter(scaffold == 29)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("brahman_chr") + ylab("HDprobe") +
  geom_vline(xintercept=c(5999048), color = "blue")
g

##### Sire #####
# this is after first round fixing; reading rc probes
dir1 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/alignAndOrderSnp_result/sire_chr/"
path1 <- paste0(dir1,"1.rcmap.tab")

sire_hd_rc_dom_v2_HDProbes <- read_tsv(path1,col_names = FALSE)
names(sire_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- sire_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 1) %>% filter(scaffold == 1)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("sire_chr") + ylab("rcprobe") +
  geom_vline(xintercept=c(48793504,48793504+84326795,48793504+84326795+818583,
                          48793504+84326795+818583+19454815,
                          48793504+84326795+818583+19454815+3130680,
                          48793504+84326795+818583+19454815+3130680+3737369,
                          48793504+84326795+818583+19454815+3130680+3737369+1726305), color = "blue")
g

path2 <- paste0(dir1,"2.rcmap.tab")

sire_hd_rc_dom_v2_HDProbes <- read_tsv(path2,col_names = FALSE)
names(sire_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- sire_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 2) %>% filter(scaffold == 2)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("sire_chr") + ylab("rcprobe") +
  geom_vline(xintercept=c(1), color = "blue")
g

path3 <- paste0(dir1,"3.rcmap.tab")

sire_hd_rc_dom_v2_HDProbes <- read_tsv(path3,col_names = FALSE)
names(sire_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- sire_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 3) %>% filter(scaffold == 3)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("sire_chr") + ylab("rcprobe") +
  geom_vline(xintercept=c(1), color = "blue")
g

path4 <- paste0(dir1,"4.rcmap.tab")

sire_hd_rc_dom_v2_HDProbes <- read_tsv(path4,col_names = FALSE)
names(sire_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- sire_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 4) %>% filter(scaffold == 4)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("sire_chr") + ylab("rcprobe") +
  geom_vline(xintercept=c(1), color = "blue")
g

path5 <- paste0(dir1,"5.rcmap.tab")

sire_hd_rc_dom_v2_HDProbes <- read_tsv(path5,col_names = FALSE)
names(sire_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- sire_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 5) %>% filter(scaffold == 5)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("sire_chr") + ylab("rcprobe") +
  geom_vline(xintercept=c(1), color = "blue")
g

path6 <- paste0(dir1,"6.rcmap.tab")

sire_hd_rc_dom_v2_HDProbes <- read_tsv(path6,col_names = FALSE)
names(sire_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- sire_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 6) %>% filter(scaffold == 6)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("sire_chr") + ylab("rcprobe") +
  geom_vline(xintercept=c(1), color = "blue")
g

path7 <- paste0(dir1,"7.rcmap.tab")

sire_hd_rc_dom_v2_HDProbes <- read_tsv(path7,col_names = FALSE)
names(sire_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- sire_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 7) %>% filter(scaffold == 7)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("sire_chr") + ylab("rcprobe") +
  geom_vline(xintercept=c(18816456,18816456+11207342,18816456+11207342+13219703,
                          18816456+11207342+13219703+9605195,
                          18816456+11207342+13219703+9605195+2806910,
                          18816456+11207342+13219703+9605195+2806910+45071626,
                          18816456+11207342+13219703+9605195+2806910+45071626+7559514), color = "blue")
g

path8 <- paste0(dir1,"8.rcmap.tab")

sire_hd_rc_dom_v2_HDProbes <- read_tsv(path8,col_names = FALSE)
names(sire_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- sire_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 8) %>% filter(scaffold == 8)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("sire_chr") + ylab("rcprobe") +
  geom_vline(xintercept=c(19795968,19795968+2703322,19795968+2703322+8713867,
                          19795968+2703322+8713867+13660270,
                          19795968+2703322+8713867+13660270+15505350,
                          19795968+2703322+8713867+13660270+15505350+599653,
                          19795968+2703322+8713867+13660270+15505350+599653+5742948,
                          19795968+2703322+8713867+13660270+15505350+599653+5742948+5508537,
                          19795968+2703322+8713867+13660270+15505350+599653+5742948+5508537+34962082), color = "blue")
g

path9 <- paste0(dir1,"9.rcmap.tab")

sire_hd_rc_dom_v2_HDProbes <- read_tsv(path9,col_names = FALSE)
names(sire_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- sire_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 9) %>% filter(scaffold == 9)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("sire_chr") + ylab("rcprobe") +
  geom_vline(xintercept=c(54161935,54161935+48959662), color = "blue")
g

path10 <- paste0(dir1,"10.rcmap.tab")

sire_hd_rc_dom_v2_HDProbes <- read_tsv(path10,col_names = FALSE)
names(sire_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- sire_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 10) %>% filter(scaffold == 10)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("sire_chr") + ylab("rcprobe") +
  geom_vline(xintercept=c(4928945,4928945+95520976), color = "blue")
g

path11 <- paste0(dir1,"11.rcmap.tab")

sire_hd_rc_dom_v2_HDProbes <- read_tsv(path11,col_names = FALSE)
names(sire_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- sire_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 11) %>% filter(scaffold == 11)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("sire_chr") + ylab("rcprobe") +
  geom_vline(xintercept=c(12203522), color = "blue")
g

path12 <- paste0(dir1,"12.rcmap.tab")

sire_hd_rc_dom_v2_HDProbes <- read_tsv(path12,col_names = FALSE)
names(sire_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- sire_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 12) %>% filter(scaffold == 12)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("sire_chr") + ylab("rcprobe") +
  geom_vline(xintercept=c(59797708,59797708+6246420), color = "blue")
g

path13 <- paste0(dir1,"13.rcmap.tab")

sire_hd_rc_dom_v2_HDProbes <- read_tsv(path13,col_names = FALSE)
names(sire_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- sire_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 13) %>% filter(scaffold == 13)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("sire_chr") + ylab("rcprobe") +
  geom_vline(xintercept=c(1), color = "blue")
g

path14 <- paste0(dir1,"14.rcmap.tab")

sire_hd_rc_dom_v2_HDProbes <- read_tsv(path14,col_names = FALSE)
names(sire_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- sire_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 14) %>% filter(scaffold == 14)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("sire_chr") + ylab("rcprobe") +
  geom_vline(xintercept=c(1), color = "blue")
g

path15 <- paste0(dir1,"15.rcmap.tab")

sire_hd_rc_dom_v2_HDProbes <- read_tsv(path15,col_names = FALSE)
names(sire_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- sire_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 15) %>% filter(scaffold == 15)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("sire_chr") + ylab("rcprobe") +
  geom_vline(xintercept=c(1), color = "blue")
g


path16 <- paste0(dir1,"16.rcmap.tab")

sire_hd_rc_dom_v2_HDProbes <- read_tsv(path16,col_names = FALSE)
names(sire_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- sire_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 16) %>% filter(scaffold == 16)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("sire_chr") + ylab("rcprobe") +
  geom_vline(xintercept=c(6279414), color = "blue")
g

path17 <- paste0(dir1,"17.rcmap.tab")

sire_hd_rc_dom_v2_HDProbes <- read_tsv(path17,col_names = FALSE)
names(sire_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- sire_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 17) %>% filter(scaffold == 17)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("sire_chr") + ylab("rcprobe") +
  geom_vline(xintercept=c(1), color = "blue")
g

path18 <- paste0(dir1,"18.rcmap.tab")

sire_hd_rc_dom_v2_HDProbes <- read_tsv(path18,col_names = FALSE)
names(sire_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- sire_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 18) %>% filter(scaffold == 18)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("sire_chr") + ylab("rcprobe") +
  geom_vline(xintercept=c(29439230,29439230+8276155,29439230+8276155+6296727,
                          29439230+8276155+6296727+6307394,
                          29439230+8276155+6296727+6307394+8057580,
                          29439230+8276155+6296727+6307394+8057580+898024), color = "blue")
g

path19 <- paste0(dir1,"19.rcmap.tab")

sire_hd_rc_dom_v2_HDProbes <- read_tsv(path19,col_names = FALSE)
names(sire_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- sire_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 19) %>% filter(scaffold == 19)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("sire_chr") + ylab("rcprobe") +
  geom_vline(xintercept=c(1), color = "blue")
g

path20 <- paste0(dir1,"20.rcmap.tab")

sire_hd_rc_dom_v2_HDProbes <- read_tsv(path20,col_names = FALSE)
names(sire_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- sire_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 20) %>% filter(scaffold == 20)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("sire_chr") + ylab("rcprobe") +
  geom_vline(xintercept=c(12805794+19889672+4900778+2125641+18390672+2343122), color = "blue")
g

path21 <- paste0(dir1,"21.rcmap.tab")

sire_hd_rc_dom_v2_HDProbes <- read_tsv(path21,col_names = FALSE)
names(sire_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- sire_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 21) %>% filter(scaffold == 21)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("sire_chr") + ylab("rcprobe") +
  geom_vline(xintercept=c(1), color = "blue")
g

path22 <- paste0(dir1,"22.rcmap.tab")

sire_hd_rc_dom_v2_HDProbes <- read_tsv(path22,col_names = FALSE)
names(sire_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- sire_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 22) %>% filter(scaffold == 22)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("sire_chr") + ylab("rcprobe") +
  geom_vline(xintercept=c(1), color = "blue")
g

path23 <- paste0(dir1,"23.rcmap.tab")

sire_hd_rc_dom_v2_HDProbes <- read_tsv(path23,col_names = FALSE)
names(sire_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- sire_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 23) %>% filter(scaffold == 23)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("sire_chr") + ylab("rcprobe") +
  geom_vline(xintercept=c(8084886,8084886+34250757,8084886+34250757+8574807), color = "blue")
g

path24 <- paste0(dir1,"24.rcmap.tab")

sire_hd_rc_dom_v2_HDProbes <- read_tsv(path24,col_names = FALSE)
names(sire_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- sire_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 24) %>% filter(scaffold == 24)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("sire_chr") + ylab("rcprobe") +
  geom_vline(xintercept=c(1), color = "blue")
g

path25 <- paste0(dir1,"25.rcmap.tab")

sire_hd_rc_dom_v2_HDProbes <- read_tsv(path25,col_names = FALSE)
names(sire_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- sire_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 25) %>% filter(scaffold == 25)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("sire_chr") + ylab("rcprobe") +
  geom_vline(xintercept=c(1), color = "blue")
g

path26 <- paste0(dir1,"26.rcmap.tab")

sire_hd_rc_dom_v2_HDProbes <- read_tsv(path26,col_names = FALSE)
names(sire_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- sire_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 26) %>% filter(scaffold == 26)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("sire_chr") + ylab("rcprobe") +
  geom_vline(xintercept=c(1), color = "blue")
g

path27 <- paste0(dir1,"27.rcmap.tab")

sire_hd_rc_dom_v2_HDProbes <- read_tsv(path27,col_names = FALSE)
names(sire_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- sire_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 27) %>% filter(scaffold == 27)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("sire_chr") + ylab("rcprobe") +
  geom_vline(xintercept=c(1), color = "blue")
g

path28 <- paste0(dir1,"28.rcmap.tab")

sire_hd_rc_dom_v2_HDProbes <- read_tsv(path28,col_names = FALSE)
names(sire_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- sire_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 28) %>% filter(scaffold == 28)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("sire_chr") + ylab("rcprobe") +
  geom_vline(xintercept=c(2889604,2889604+4191532,2889604+4191532+9509999), color = "blue")
g

path29 <- paste0(dir1,"29.rcmap.tab")

sire_hd_rc_dom_v2_HDProbes <- read_tsv(path29,col_names = FALSE)
names(sire_hd_rc_dom_v2_HDProbes) <- c("marker","scaffold","align_pos","orientation","chromosome","chr_pos")

DF <- sire_hd_rc_dom_v2_HDProbes %>% filter(chromosome == 29) %>% filter(scaffold == 29)

g <- ggplot(DF, aes(x=align_pos,y=chr_pos)) + geom_point() + xlab("sire_chr") + ylab("rcprobe") +
  geom_vline(xintercept=c(5544007), color = "blue")
g
