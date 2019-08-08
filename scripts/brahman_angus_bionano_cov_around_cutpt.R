#------------------------------------------------------
# Program name: brahman_angus_bionano_cov_around_cutpt.R
# Objective: check cov Illumina to determine cutpt
#           
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(readr)
library(dplyr)
library(ggplot2)

#Dam
dir1 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/conflict_resolution/cutpt_for_looping/coverage/dam/"

#create a list to read all *cov file
file_list <- list.files(path = dir1, pattern = "*.cov")
breakpt <- read_tsv("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/conflict_resolution/cutpt_for_looping/cutpt_brahma.tsv",
                    col_names = FALSE)
colnames(breakpt) <- c("input","startseq","breakpt","start","end")

afterbreakptcorrect <- read_tsv("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/conflict_resolution/cutpt_for_looping/cutpt_brahma_breakpt_modi.tsv",
                                col_names = FALSE)
colnames(afterbreakptcorrect) <- c("input","startseq","breakpt","start","end","change")

for (i in 1:length(file_list)){
  file <- assign(file_list[i], read_tsv(paste0(dir1,file_list[i]),col_names = FALSE))
  colnames(file) <- c("input","pos","cov")
  brahman_breakpt_pdf <- paste0("Brahman_",file_list[i],".png")
  intercept <- breakpt$breakpt[as.integer(gsub("_.*","",file_list[i]))]
  newintercept <- afterbreakptcorrect$breakpt[as.integer(gsub("_.*","",file_list[i]))]
  
  title <- gsub("^[0-9]+_","",file_list[i])
  title <- gsub("_.*","",title)
  g <- ggplot(file, aes(x=pos,y=cov)) + geom_line() + xlab("position") + ylab("coverage") +
    geom_vline(xintercept=intercept, color = "red") + 
    geom_vline(xintercept=newintercept, color = "blue") +
    geom_hline(yintercept=c(32.5,65), color = "blue") + 
    ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
  g
  ggsave(brahman_breakpt_pdf,device = "png")
}

#Individual breakpt check
#this is for default graphics plot
for (i in 1:length(file_list)){
  file <- assign(file_list[i], read_tsv(paste0(dir1,file_list[i]),col_names = FALSE))
  colnames(file) <- c("input","pos","cov")
  plot(file$pos,file$cov,type = "l",xlab = "Position",
       ylab = "Coverage",ylim = c(0,range(ill_cov$cov)[2]))
  abline(v=breakpt$breakpt[as.integer(gsub("_.*","",file_list[i]))],col = "red")
}

#1_tig00000424_arrow_arrow_42251924_42451924.cov
#2_tig00000424_arrow_arrow_42259539_42459539.cov
file <- read_tsv(paste0(dir1,"1_tig00000424_arrow_arrow_42251924_42451924.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(42351924,42359539), col="red")
abline(v=42346923,col = "blue")
abline(v=42349538,col = "green")
#leave as is

#3_tig00000512_arrow_arrow_19308477_19508477.cov
#4_tig00000512_arrow_arrow_19308537_19508537.cov
file <- read_tsv(paste0(dir1,"3_tig00000512_arrow_arrow_19308477_19508477.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(19408477,19408537), col="red")
abline(v=19408477, col="blue")
#leave as is

#5_tig00000573_arrow_arrow_8661778_8861778.cov
#6_tig00000573_arrow_arrow_8675991_8875991.cov
file <- read_tsv(paste0(dir1,"5_tig00000573_arrow_arrow_8661778_8861778.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(8761778,8775991), col="red")
abline(v=c(8751778,8780991), col="blue")
#changed both

#7_tig00000925_arrow_arrow_16866328_17066328.cov
file <- read_tsv(paste0(dir1,"7_tig00000925_arrow_arrow_16866328_17066328.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(16966328), col="red")
abline(v=c(17006327), col="blue")

#8_tig00002026_arrow_arrow_94521_294521.cov
#9_tig00002026_arrow_arrow_125156_325156.cov
file <- read_tsv(paste0(dir1,"8_tig00002026_arrow_arrow_94521_294521.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(194521,225156), col="red")
abline(v=c(189520), col="blue")

#10_tig00003119_arrow_arrow_248010_448010.cov
#11_tig00003119_arrow_arrow_262143_462143.cov
#12_tig00003119_arrow_arrow_274201_474201.cov
file <- read_tsv(paste0(dir1,"10_tig00003119_arrow_arrow_248010_448010.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(348010,362143,374201), col="red")
abline(v=c(357143), col="blue")
#changed 362143 to 357143

#13_tig00003246_arrow_arrow_9862_209862.cov
file <- read_tsv(paste0(dir1,"13_tig00003246_arrow_arrow_9862_209862.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(109862), col="red")
abline(v=c(159862), col="blue")
#changed 109862 to 159862

#14_tig00003407_arrow_arrow_1_159250.cov
file <- read_tsv(paste0(dir1,"14_tig00003407_arrow_arrow_1_159250.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(59250), col="red")
abline(v=c(114999), col="blue")
#changed 59250 to 114999

#15_tig00011620_arrow_arrow_29339823_29539823.cov
file <- read_tsv(paste0(dir1,"15_tig00011620_arrow_arrow_29339823_29539823.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(29439823), col="red")
abline(v=c(29439823), col="blue")
#leave as is

#16_tig00011626_arrow_arrow_19304922_19504922.cov
file <- read_tsv(paste0(dir1,"16_tig00011626_arrow_arrow_19304922_19504922.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(19404922), col="red")
abline(v=c(19400000), col="blue")
#changed from 19404922 to 19400000

#17_tig00011633_arrow_arrow_1896832_2096832.cov
file <- read_tsv(paste0(dir1,"17_tig00011633_arrow_arrow_1896832_2096832.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(1996832), col="red")
abline(v=c(1996832), col="blue")
#leave as is

#18_tig00011639_arrow_arrow_1_191084.cov
#19_tig00011639_arrow_arrow_19146_219146.cov
file <- read_tsv(paste0(dir1,"18_tig00011639_arrow_arrow_1_191084.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(91084,119146), col="red")
abline(v=c(144999), col="blue")
#leave as is

#20_tig00011653_arrow_arrow_3531345_3731345.cov
file <- read_tsv(paste0(dir1,"20_tig00011653_arrow_arrow_3531345_3731345.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(3631345), col="red")
abline(v=c(3626345), col="blue")
#change from 3631345 to 3626345

#21_tig00011723_arrow_arrow_43816_243816.cov
file <- read_tsv(paste0(dir1,"21_tig00011723_arrow_arrow_43816_243816.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(143816), col="red")
abline(v=c(128816), col="blue")
#leave as is

#22_tig00011745_arrow_arrow_530400_730400.cov
#23_tig00011745_arrow_arrow_541949_741949.cov
file <- read_tsv(paste0(dir1,"22_tig00011745_arrow_arrow_530400_730400.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(630400,641949), col="red")
abline(v=c(656949), col="blue")
#leave as is

#24_tig00011795_arrow_arrow_73016_273016.cov
file <- read_tsv(paste0(dir1,"24_tig00011795_arrow_arrow_73016_273016.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(173016), col="red")
abline(v=c(113016), col="blue")
#leave as is

#25_tig00011797_arrow_arrow_42496_242496.cov
#26_tig00011797_arrow_arrow_49185_249185.cov
file <- read_tsv(paste0(dir1,"25_tig00011797_arrow_arrow_42496_242496.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(142496,149185), col="red")
abline(v=c(165000), col="blue")
#leave as is

#27_tig00011826_arrow_arrow_1_140847.cov
file <- read_tsv(paste0(dir1,"27_tig00011826_arrow_arrow_1_140847.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(40847), col="red")
abline(v=c(14999), col="blue")
#leave as is

#28_tig00011833_arrow_arrow_5756939_5956939.cov
file <- read_tsv(paste0(dir1,"28_tig00011833_arrow_arrow_5756939_5956939.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(5856939), col="red")
abline(v=c(5826939,5861938), col="blue")
#leave as is

#29_tig00011838_arrow_arrow_1_153571.cov
file <- read_tsv(paste0(dir1,"29_tig00011838_arrow_arrow_1_153571.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(53571), col="red")
abline(v=c(69999), col="blue")
#changed from 53571 to 69999

#30_tig00011886_arrow_arrow_1_144014.cov
file <- read_tsv(paste0(dir1,"30_tig00011886_arrow_arrow_1_144014.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(44014), col="red")
abline(v=c(44014), col="blue")
#leave as is

#31_tig00011901_arrow_arrow_1_165007.cov
file <- read_tsv(paste0(dir1,"31_tig00011901_arrow_arrow_1_165007.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(65007), col="red")
abline(v=c(44014), col="blue")
#leave as is

#32_tig00011907_arrow_arrow_1_162789.cov
file <- read_tsv(paste0(dir1,"32_tig00011907_arrow_arrow_1_162789.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(62789), col="red")
abline(v=c(44014), col="blue")
#leave as is

#33_tig00478914_arrow_arrow_61564303_61764303.cov
file <- read_tsv(paste0(dir1,"33_tig00478914_arrow_arrow_61564303_61764303.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(61664303,61702992), col="red")
abline(v=c(44014), col="blue")
#leave as is

#35_tig00479033_arrow_arrow_4803034_5003034.cov
#36_tig00479033_arrow_arrow_4841920_5041920.cov
file <- read_tsv(paste0(dir1,"35_tig00479033_arrow_arrow_4803034_5003034.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(4903034,4941920), col="red")
abline(v=c(44014), col="blue")
#leave as is


#Sire
dir2 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/conflict_resolution/cutpt_for_looping/coverage/sire/"

#create a list to read all *cov file
file_list <- list.files(path = dir2, pattern = "*.cov")
breakpt <- read_tsv("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/conflict_resolution/cutpt_for_looping/cutpt_angus.tsv",
                    col_names = FALSE)
colnames(breakpt) <- c("input","startseq","breakpt","start","end")

afterbreakptcorrect <- read_tsv("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/conflict_resolution/cutpt_for_looping/cutpt_angus_breakpt_modi.tsv",
                                col_names = FALSE)
colnames(afterbreakptcorrect) <- c("input","startseq","breakpt","start","end","change")

for (i in 1:length(file_list)){
  file <- assign(file_list[i], read_tsv(paste0(dir2,file_list[i]),col_names = FALSE))
  colnames(file) <- c("input","pos","cov")
  angus_breakpt_pdf <- paste0("Angus_",file_list[i],".png")
  intercept <- breakpt$breakpt[as.integer(gsub("_.*","",file_list[i]))]
  newintercept <- afterbreakptcorrect$breakpt[as.integer(gsub("_.*","",file_list[i]))]
  
  title <- gsub("^[0-9]+_","",file_list[i])
  title <- gsub("_.*","",title)
  g <- ggplot(file, aes(x=pos,y=cov)) + geom_line() + xlab("position") + ylab("coverage") +
    geom_vline(xintercept=intercept, color = "red") + 
    geom_vline(xintercept=newintercept, color = "blue") +
    geom_hline(yintercept=c(32.5,65), color = "blue") + 
    ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
  g
  ggsave(angus_breakpt_pdf,device = "png")
}

#1_tig00000427_arrow_arrow_37646845_37846845.cov
file <- read_tsv(paste0(dir2,"1_tig00000427_arrow_arrow_37646845_37846845.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(37746845,37792551,37822567), col="red")
abline(v=c(37751845,37796845,37816845), col="blue")
#changed as above

#4_tig00000717_arrow_arrow_10804473_11004473.cov
file <- read_tsv(paste0(dir2,"4_tig00000717_arrow_arrow_10804473_11004473.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(10904473), col="red")
abline(v=c(10919472), col="blue")
#changed from 10904473 to 10919472

#5_tig00000717_arrow_arrow_26348453_26548453.cov
file <- read_tsv(paste0(dir2,"5_tig00000717_arrow_arrow_26348453_26548453.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(26448453), col="red")
abline(v=c(10919472), col="blue")
#leave as is

#6_tig00001028_arrow_arrow_6278041_6478041.cov
file <- read_tsv(paste0(dir2,"6_tig00001028_arrow_arrow_6278041_6478041.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(6378041), col="red")
abline(v=c(6368041), col="blue")
#leave as is

#7_tig00001643_arrow_arrow_219172_419172.cov
#8_tig00001643_arrow_arrow_225245_425245.cov
file <- read_tsv(paste0(dir2,"7_tig00001643_arrow_arrow_219172_419172.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(319172,325245), col="red")
abline(v=c(319172,325245), col="blue")
#leave as is

#9_tig00002180_arrow_arrow_1_172692.cov
#10_tig00002180_arrow_arrow_72693_195963.cov
file <- read_tsv(paste0(dir2,"9_tig00002180_arrow_arrow_1_172692.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(72692,95963), col="red")
abline(v=c(64999), col="blue")
#changed 72692 to 64999

#11_tig00002359_arrow_arrow_1_195167.cov
#15_tig00002359_arrow_arrow_11214_211214.cov
file <- read_tsv(paste0(dir2,"11_tig00002359_arrow_arrow_1_195167.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(95167,99179,103192,107203,111214), col="red")
abline(v=c(64999), col="blue")
#leave as is

#16_tig00020269_arrow_arrow_85582436_85782436.cov
file <- read_tsv(paste0(dir2,"16_tig00020269_arrow_arrow_85582436_85782436.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(85682436,85762408), col="red")
abline(v=c(64999), col="blue")
#leave as is

#18_tig00020275_arrow_arrow_1_197865.cov
file <- read_tsv(paste0(dir2,"18_tig00020275_arrow_arrow_1_197865.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(97865), col="red")
abline(v=c(64999), col="blue")
#leave as is

#19_tig00020280_arrow_arrow_1_156605.cov
file <- read_tsv(paste0(dir2,"19_tig00020280_arrow_arrow_1_156605.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(56605), col="red")
abline(v=c(55000), col="blue")
#changed from 56605 to 55000

#20_tig00020313_arrow_arrow_15079554_15279554.cov
file <- read_tsv(paste0(dir2,"20_tig00020313_arrow_arrow_15079554_15279554.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(15179554), col="red")
abline(v=c(55000), col="blue")
#leave as is

#21_tig00020332_arrow_arrow_4690891_4890891.cov
file <- read_tsv(paste0(dir2,"21_tig00020332_arrow_arrow_4690891_4890891.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(4790891), col="red")
abline(v=c(4790890), col="blue")
#leave as is

#22_tig00020357_arrow_arrow_449466_649466.cov
file <- read_tsv(paste0(dir2,"22_tig00020357_arrow_arrow_449466_649466.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(549466), col="red")
abline(v=c(4790890), col="blue")
#leave as is

#23_tig00020407_arrow_arrow_42106_242106.cov
#24_tig00020407_arrow_arrow_114451_314451.cov
file <- read_tsv(paste0(dir2,"23_tig00020407_arrow_arrow_42106_242106.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(142106,214451), col="red")
abline(v=c(4790890), col="blue")
#leave as is

#25_tig00020408_arrow_arrow_640955_840955.cov
file <- read_tsv(paste0(dir2,"25_tig00020408_arrow_arrow_640955_840955.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(740955), col="red")
abline(v=c(750954), col="blue")
#leave as is

#26_tig00020445_arrow_arrow_16170_216170.cov
file <- read_tsv(paste0(dir2,"26_tig00020445_arrow_arrow_16170_216170.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(116170,173534,203344), col="red")
abline(v=c(750954), col="blue")
#leave as is

#29_tig00516280_arrow_arrow_4202_204202.cov
file <- read_tsv(paste0(dir2,"29_tig00516280_arrow_arrow_4202_204202.cov"),col_names = FALSE)
colnames(file) <- c("input","pos","cov")
plot(file$pos,file$cov,type = "l",xlab = "Position",
     ylab = "Coverage",xlim = c())
abline(v=c(104202), col="red")
abline(v=c(79202,149201), col="blue")
#leave as is
#This is the one that didn't match my chimeric contig break based on dominette by a large dist


