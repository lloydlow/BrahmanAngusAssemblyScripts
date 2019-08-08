#------------------------------------------------------
# Program name: brahman_angus_IsoPhase_allelic_imbalance.R
# Objective: analyse any allelic imbalance in general, not 
#         specific to any gene
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(epade)
library(tidyr)
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(ggbiplot)
library(preprocessCore)

# path to isophase results
dir1 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/isoseq/IsoPhase/April2019_IsoPhase_with_Lloyd/"

#Isoform not separated #evaled_isophase.demux_hap_count.txt
path1 <- paste0(dir1,"evaled_isophase.demux_hap_count.txt")

evaled_isophase.demux_hap_count <- read_tsv(path1)

#proportion of brahman allele modi
evaled_isophase.demux_hap_count_modi <- evaled_isophase.demux_hap_count %>% 
  mutate(heart_prop = heart_p0/(heart_p0+heart_p1)) %>%
  mutate(liver_prop = liver_p0/(liver_p0+liver_p1)) %>%
  mutate(kidney_prop = kidney_p0/(kidney_p0+kidney_p1)) %>%
  mutate(brain_prop = brain_p0/(brain_p0+brain_p1)) %>%
  mutate(lung_prop = lung_p0/(lung_p0+lung_p1)) %>%
  mutate(muscle_prop = muscle_p0/(muscle_p0+muscle_p1)) %>%
  mutate(placenta_prop = placenta_p0/(placenta_p0+placenta_p1))

#total transcript count modi
evaled_isophase.demux_hap_count_modi <- evaled_isophase.demux_hap_count_modi %>% 
  mutate(heart_total_transcript = heart_p0+heart_p1) %>%
  mutate(liver_total_transcript = liver_p0+liver_p1) %>%
  mutate(kidney_total_transcript = kidney_p0+kidney_p1) %>%
  mutate(brain_total_transcript = brain_p0+brain_p1) %>%
  mutate(lung_total_transcript = lung_p0+lung_p1) %>%
  mutate(muscle_total_transcript = muscle_p0+muscle_p1) %>%
  mutate(placenta_total_transcript = placenta_p0+placenta_p1)

#log fold change
# evaled_isophase.demux_hap_count_modi <- evaled_isophase.demux_hap_count_modi %>% 
#   mutate(heart_log_fold_change = log2(heart_p0/heart_p1)) %>%
#   mutate(liver_log_fold_change = log2(liver_p0/liver_p1)) %>%
#   mutate(kidney_log_fold_change = log2(kidney_p0/kidney_p1)) %>%
#   mutate(brain_log_fold_change = log2(brain_p0/brain_p1)) %>%
#   mutate(lung_log_fold_change = log2(lung_p0/lung_p1)) %>%
#   mutate(muscle_log_fold_change = log2(muscle_p0/muscle_p1)) %>%
#   mutate(placenta_log_fold_change = log2(placenta_p0/placenta_p1))

################################################################################################
#####Normalization
#by lib size, TPM
evaled_isophase.demux_hap_count_TPM <- evaled_isophase.demux_hap_count %>%
  mutate(heart_p0_new = (heart_p0/(sum(heart_p0)+sum(heart_p1)))*1e6) %>%
  mutate(heart_p1_new = (heart_p1/(sum(heart_p0)+sum(heart_p1)))*1e6) %>%
  mutate(liver_p0_new = (liver_p0/(sum(liver_p0)+sum(liver_p1)))*1e6) %>%
  mutate(liver_p1_new = (liver_p1/(sum(liver_p0)+sum(liver_p1)))*1e6) %>%
  mutate(kidney_p0_new = (kidney_p0/(sum(kidney_p0)+sum(kidney_p1)))*1e6) %>%
  mutate(kidney_p1_new = (kidney_p1/(sum(kidney_p0)+sum(kidney_p1)))*1e6) %>%
  mutate(brain_p0_new = (brain_p0/(sum(brain_p0)+sum(brain_p1)))*1e6) %>%
  mutate(brain_p1_new = (brain_p1/(sum(brain_p0)+sum(brain_p1)))*1e6) %>%
  mutate(lung_p0_new = (lung_p0/(sum(lung_p0)+sum(lung_p1)))*1e6) %>%
  mutate(lung_p1_new = (lung_p1/(sum(lung_p0)+sum(lung_p1)))*1e6) %>%
  mutate(muscle_p0_new = (muscle_p0/(sum(muscle_p0)+sum(muscle_p1)))*1e6) %>%
  mutate(muscle_p1_new = (muscle_p1/(sum(muscle_p0)+sum(muscle_p1)))*1e6) %>%
  mutate(placenta_p0_new = (placenta_p0/(sum(placenta_p0)+sum(placenta_p1)))*1e6) %>%
  mutate(placenta_p1_new = (placenta_p1/(sum(placenta_p0)+sum(placenta_p1)))*1e6)

#remove original count columns, replace with and rename the lib size TPM normalized columns
ori_col_names <- colnames(evaled_isophase.demux_hap_count_TPM)[4:17]
evaled_isophase.demux_hap_count_TPM <- evaled_isophase.demux_hap_count_TPM %>% select(locus:phase1,heart_p0_new:placenta_p1_new)
colnames(evaled_isophase.demux_hap_count_TPM)[4:17] <- ori_col_names

#total count per tissue
total_count_all_tissues <- evaled_isophase.demux_hap_count_modi %>%
  summarise(total_count_heart = sum(heart_total_transcript),total_count_liver = sum(liver_total_transcript),
            total_count_kidney = sum(kidney_total_transcript),total_count_brain = sum(brain_total_transcript),
            total_count_lung = sum(lung_total_transcript),total_count_muscle = sum(muscle_total_transcript),total_count_placenta = sum(placenta_total_transcript))

#quantile normalization
mat <- as.matrix(evaled_isophase.demux_hap_count_modi[,4:17])
evaled_isophase.demux_hap_count_modi_Qnorm <- normalize.quantiles(mat)

colnames(evaled_isophase.demux_hap_count_modi_Qnorm) <- c("heart_p0", "heart_p1","liver_p0", "liver_p1","kidney_p0", "kidney_p1",
  "brain_p0","brain_p1","lung_p0","lung_p1","muscle_p0","muscle_p1",
  "placenta_p0","placenta_p1")

evaled_isophase.demux_hap_count_Qnorm <- cbind(evaled_isophase.demux_hap_count_modi[,1:3],evaled_isophase.demux_hap_count_modi_Qnorm)
evaled_isophase.demux_hap_count_Qnorm <- as.tbl(evaled_isophase.demux_hap_count_Qnorm)

################################################################################################
#####boxplots
#boxplot to see that the distribution of counts vary across samples per allele
boxplot(evaled_isophase.demux_hap_count_modi$heart_p0,
        evaled_isophase.demux_hap_count_modi$heart_p1,
        evaled_isophase.demux_hap_count_modi$liver_p0,
        evaled_isophase.demux_hap_count_modi$liver_p1,
        evaled_isophase.demux_hap_count_modi$kidney_p0,
        evaled_isophase.demux_hap_count_modi$kidney_p1,
        evaled_isophase.demux_hap_count_modi$brain_p0,
        evaled_isophase.demux_hap_count_modi$brain_p1,
        evaled_isophase.demux_hap_count_modi$lung_p0,
        evaled_isophase.demux_hap_count_modi$lung_p1,
        evaled_isophase.demux_hap_count_modi$muscle_p0,
        evaled_isophase.demux_hap_count_modi$muscle_p1,
        evaled_isophase.demux_hap_count_modi$placenta_p0,
        evaled_isophase.demux_hap_count_modi$placenta_p1,
        outline=FALSE,names = c("heart_B", "heart_A","liver_B", "liver_A","kidney_B", "kidney_A",
                                "brain_B","brain_A","lung_B","lung_A","muscle_B","muscle_A",
                                "placenta_B","placenta_A"))

#boxplot to see that the distribution of counts vary across samples
# boxplot(evaled_isophase.demux_hap_count_modi$heart_total_transcript,
#         evaled_isophase.demux_hap_count_modi$liver_total_transcript,
#         evaled_isophase.demux_hap_count_modi$kidney_total_transcript,
#         evaled_isophase.demux_hap_count_modi$brain_total_transcript,
#         evaled_isophase.demux_hap_count_modi$lung_total_transcript,
#         evaled_isophase.demux_hap_count_modi$muscle_total_transcript,
#         evaled_isophase.demux_hap_count_modi$placenta_total_transcript,
#         outline=FALSE,names = c("heart", "liver", "kidney", "brain","lung","muscle","placenta"))

#boxplot to see that the distribution of counts vary across samples per allele after TPM
boxplot(evaled_isophase.demux_hap_count_TPM$heart_p0,
        evaled_isophase.demux_hap_count_TPM$heart_p1,
        evaled_isophase.demux_hap_count_TPM$liver_p0,
        evaled_isophase.demux_hap_count_TPM$liver_p1,
        evaled_isophase.demux_hap_count_TPM$kidney_p0,
        evaled_isophase.demux_hap_count_TPM$kidney_p1,
        evaled_isophase.demux_hap_count_TPM$brain_p0,
        evaled_isophase.demux_hap_count_TPM$brain_p1,
        evaled_isophase.demux_hap_count_TPM$lung_p0,
        evaled_isophase.demux_hap_count_TPM$lung_p1,
        evaled_isophase.demux_hap_count_TPM$muscle_p0,
        evaled_isophase.demux_hap_count_TPM$muscle_p1,
        evaled_isophase.demux_hap_count_TPM$placenta_p0,
        evaled_isophase.demux_hap_count_TPM$placenta_p1,
        outline=FALSE, names = c("heart_B", "heart_A","liver_B", "liver_A","kidney_B", "kidney_A",
                                 "brain_B","brain_A","lung_B","lung_A","muscle_B","muscle_A",
                                 "placenta_B","placenta_A"))

#boxplot to see that the distribution of counts vary across samples per allele after Qnorm
boxplot(evaled_isophase.demux_hap_count_Qnorm$heart_p0,
        evaled_isophase.demux_hap_count_Qnorm$heart_p1,
        evaled_isophase.demux_hap_count_Qnorm$liver_p0,
        evaled_isophase.demux_hap_count_Qnorm$liver_p1,
        evaled_isophase.demux_hap_count_Qnorm$kidney_p0,
        evaled_isophase.demux_hap_count_Qnorm$kidney_p1,
        evaled_isophase.demux_hap_count_Qnorm$brain_p0,
        evaled_isophase.demux_hap_count_Qnorm$brain_p1,
        evaled_isophase.demux_hap_count_Qnorm$lung_p0,
        evaled_isophase.demux_hap_count_Qnorm$lung_p1,
        evaled_isophase.demux_hap_count_Qnorm$muscle_p0,
        evaled_isophase.demux_hap_count_Qnorm$muscle_p1,
        evaled_isophase.demux_hap_count_Qnorm$placenta_p0,
        evaled_isophase.demux_hap_count_Qnorm$placenta_p1,
        outline=FALSE, names = c("heart_B", "heart_A","liver_B", "liver_A","kidney_B", "kidney_A",
                                 "brain_B","brain_A","lung_B","lung_A","muscle_B","muscle_A",
                                 "placenta_B","placenta_A"))

################################################################################################
#####Proportion and total transcript count
#proportion of brahman allele TPM
evaled_isophase.demux_hap_count_TPM <- evaled_isophase.demux_hap_count_TPM %>% 
  mutate(heart_prop = heart_p0/(heart_p0+heart_p1)) %>%
  mutate(liver_prop = liver_p0/(liver_p0+liver_p1)) %>%
  mutate(kidney_prop = kidney_p0/(kidney_p0+kidney_p1)) %>%
  mutate(brain_prop = brain_p0/(brain_p0+brain_p1)) %>%
  mutate(lung_prop = lung_p0/(lung_p0+lung_p1)) %>%
  mutate(muscle_prop = muscle_p0/(muscle_p0+muscle_p1)) %>%
  mutate(placenta_prop = placenta_p0/(placenta_p0+placenta_p1))

#total transcript count TPM
evaled_isophase.demux_hap_count_TPM <- evaled_isophase.demux_hap_count_TPM %>% 
  mutate(heart_total_transcript = heart_p0+heart_p1) %>%
  mutate(liver_total_transcript = liver_p0+liver_p1) %>%
  mutate(kidney_total_transcript = kidney_p0+kidney_p1) %>%
  mutate(brain_total_transcript = brain_p0+brain_p1) %>%
  mutate(lung_total_transcript = lung_p0+lung_p1) %>%
  mutate(muscle_total_transcript = muscle_p0+muscle_p1) %>%
  mutate(placenta_total_transcript = placenta_p0+placenta_p1)

save(evaled_isophase.demux_hap_count_TPM,file = "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/isoseq/IsoPhase/April2019_IsoPhase_with_Lloyd/evaled_isophase.demux_hap_count_TPM.RData")

#log fold change TPM
# evaled_isophase.demux_hap_count_TPM <- evaled_isophase.demux_hap_count_TPM %>% 
#   mutate(heart_log_fold_change = log2(heart_p0/heart_p1)) %>%
#   mutate(liver_log_fold_change = log2(liver_p0/liver_p1)) %>%
#   mutate(kidney_log_fold_change = log2(kidney_p0/kidney_p1)) %>%
#   mutate(brain_log_fold_change = log2(brain_p0/brain_p1)) %>%
#   mutate(lung_log_fold_change = log2(lung_p0/lung_p1)) %>%
#   mutate(muscle_log_fold_change = log2(muscle_p0/muscle_p1)) %>%
#   mutate(placenta_log_fold_change = log2(placenta_p0/placenta_p1))

#proportion of brahman allele Qnorm
evaled_isophase.demux_hap_count_Qnorm <- evaled_isophase.demux_hap_count_Qnorm %>% 
  mutate(heart_prop = heart_p0/(heart_p0+heart_p1)) %>%
  mutate(liver_prop = liver_p0/(liver_p0+liver_p1)) %>%
  mutate(kidney_prop = kidney_p0/(kidney_p0+kidney_p1)) %>%
  mutate(brain_prop = brain_p0/(brain_p0+brain_p1)) %>%
  mutate(lung_prop = lung_p0/(lung_p0+lung_p1)) %>%
  mutate(muscle_prop = muscle_p0/(muscle_p0+muscle_p1)) %>%
  mutate(placenta_prop = placenta_p0/(placenta_p0+placenta_p1))

#total transcript count Qnorm
evaled_isophase.demux_hap_count_Qnorm <- evaled_isophase.demux_hap_count_Qnorm %>% 
  mutate(heart_total_transcript = heart_p0+heart_p1) %>%
  mutate(liver_total_transcript = liver_p0+liver_p1) %>%
  mutate(kidney_total_transcript = kidney_p0+kidney_p1) %>%
  mutate(brain_total_transcript = brain_p0+brain_p1) %>%
  mutate(lung_total_transcript = lung_p0+lung_p1) %>%
  mutate(muscle_total_transcript = muscle_p0+muscle_p1) %>%
  mutate(placenta_total_transcript = placenta_p0+placenta_p1)

#log fold change Qnorm
# evaled_isophase.demux_hap_count_Qnorm <- evaled_isophase.demux_hap_count_Qnorm %>% 
#   mutate(heart_log_fold_change = log2(heart_p0/heart_p1)) %>%
#   mutate(liver_log_fold_change = log2(liver_p0/liver_p1)) %>%
#   mutate(kidney_log_fold_change = log2(kidney_p0/kidney_p1)) %>%
#   mutate(brain_log_fold_change = log2(brain_p0/brain_p1)) %>%
#   mutate(lung_log_fold_change = log2(lung_p0/lung_p1)) %>%
#   mutate(muscle_log_fold_change = log2(muscle_p0/muscle_p1)) %>%
#   mutate(placenta_log_fold_change = log2(placenta_p0/placenta_p1))

################################################################################################
#####modi
#heart
plot(density(evaled_isophase.demux_hap_count_modi$heart_prop,na.rm = TRUE))

shapiro.test(sample(evaled_isophase.demux_hap_count_modi$heart_prop,5000))

qqnorm(evaled_isophase.demux_hap_count_modi$heart_prop);qqline(evaled_isophase.demux_hap_count_modi$heart_prop, col = 2)

#liver
plot(density(evaled_isophase.demux_hap_count_modi$liver_prop,na.rm = TRUE))

shapiro.test(sample(evaled_isophase.demux_hap_count_modi$liver_prop,5000))

qqnorm(evaled_isophase.demux_hap_count_modi$liver_prop);qqline(evaled_isophase.demux_hap_count_modi$liver_prop, col = 2)

#kidney
plot(density(evaled_isophase.demux_hap_count_modi$kidney_prop,na.rm = TRUE))

shapiro.test(sample(evaled_isophase.demux_hap_count_modi$kidney_prop,5000))

qqnorm(evaled_isophase.demux_hap_count_modi$kidney_prop);qqline(evaled_isophase.demux_hap_count_modi$kidney_prop, col = 2)

#brain
plot(density(evaled_isophase.demux_hap_count_modi$brain_prop,na.rm = TRUE))

shapiro.test(sample(evaled_isophase.demux_hap_count_modi$brain_prop,5000))

qqnorm(evaled_isophase.demux_hap_count_modi$brain_prop);qqline(evaled_isophase.demux_hap_count_modi$brain_prop, col = 2)

#lung
plot(density(evaled_isophase.demux_hap_count_modi$lung_prop,na.rm = TRUE))

shapiro.test(sample(evaled_isophase.demux_hap_count_modi$lung_prop,5000))

qqnorm(evaled_isophase.demux_hap_count_modi$lung_prop);qqline(evaled_isophase.demux_hap_count_modi$lung_prop, col = 2)

#muscle
plot(density(evaled_isophase.demux_hap_count_modi$muscle_prop,na.rm = TRUE))

shapiro.test(sample(evaled_isophase.demux_hap_count_modi$muscle_prop,5000))

qqnorm(evaled_isophase.demux_hap_count_modi$muscle_prop);qqline(evaled_isophase.demux_hap_count_modi$muscle_prop, col = 2)

#placenta
plot(density(evaled_isophase.demux_hap_count_modi$placenta_prop,na.rm = TRUE))

shapiro.test(sample(evaled_isophase.demux_hap_count_modi$placenta_prop,5000))

qqnorm(evaled_isophase.demux_hap_count_modi$placenta_prop);qqline(evaled_isophase.demux_hap_count_modi$placenta_prop, col = 2)

#multiplot of proportion density
par(mfrow=c(3,3))
plot(density(evaled_isophase.demux_hap_count_modi$heart_prop,na.rm = TRUE),xlab="",main="heart")
plot(density(evaled_isophase.demux_hap_count_modi$liver_prop,na.rm = TRUE),xlab="",main="liver")
plot(density(evaled_isophase.demux_hap_count_modi$kidney_prop,na.rm = TRUE),xlab="",main="kidney")
plot(density(evaled_isophase.demux_hap_count_modi$brain_prop,na.rm = TRUE),xlab="",main="brain")
plot(density(evaled_isophase.demux_hap_count_modi$lung_prop,na.rm = TRUE),xlab="",main="lung")
plot(density(evaled_isophase.demux_hap_count_modi$muscle_prop,na.rm = TRUE),xlab="",main="muscle")
plot(density(evaled_isophase.demux_hap_count_modi$placenta_prop,na.rm = TRUE),xlab="",main="placenta")
par(mfrow=c(1,1))

par(mfrow=c(3,3))
hist(evaled_isophase.demux_hap_count_modi$heart_prop,xlab="",main="heart")
hist(evaled_isophase.demux_hap_count_modi$liver_prop,xlab="",main="liver")
hist(evaled_isophase.demux_hap_count_modi$kidney_prop,xlab="",main="kidney")
hist(evaled_isophase.demux_hap_count_modi$brain_prop,xlab="",main="brain")
hist(evaled_isophase.demux_hap_count_modi$lung_prop,xlab="",main="lung")
hist(evaled_isophase.demux_hap_count_modi$muscle_prop,xlab="",main="muscle")
hist(evaled_isophase.demux_hap_count_modi$placenta_prop,xlab="",main="placenta")
par(mfrow=c(1,1))

##### summary of proportion values distribution #####
summary(evaled_isophase.demux_hap_count_modi$heart_prop)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.0000  0.4286  0.5000  0.5083  0.5964  1.0000     540 
summary(evaled_isophase.demux_hap_count_modi$liver_prop)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.0000  0.2000  0.5000  0.5102  1.0000  1.0000    2861 
summary(evaled_isophase.demux_hap_count_modi$kidney_prop)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.0000  0.4286  0.5000  0.5045  0.5833  1.0000     385 
summary(evaled_isophase.demux_hap_count_modi$brain_prop)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.0000  0.4364  0.5000  0.5087  0.5771  1.0000     110 
summary(evaled_isophase.demux_hap_count_modi$lung_prop)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.0000  0.3333  0.5000  0.5066  0.6667  1.0000    1167 
summary(evaled_isophase.demux_hap_count_modi$muscle_prop)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.0000  0.3333  0.5000  0.5187  0.7500  1.0000    1522 
summary(evaled_isophase.demux_hap_count_modi$placenta_prop)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.0000  0.2000  0.5000  0.5106  0.8571  1.0000    2506 


################################################################################################
#####TPM
#heart
plot(density(evaled_isophase.demux_hap_count_TPM$heart_prop,na.rm = TRUE))

shapiro.test(sample(evaled_isophase.demux_hap_count_TPM$heart_prop,5000))

qqnorm(evaled_isophase.demux_hap_count_TPM$heart_prop);qqline(evaled_isophase.demux_hap_count_TPM$heart_prop, col = 2)

#liver
plot(density(evaled_isophase.demux_hap_count_TPM$liver_prop,na.rm = TRUE))

shapiro.test(sample(evaled_isophase.demux_hap_count_TPM$liver_prop,5000))

qqnorm(evaled_isophase.demux_hap_count_TPM$liver_prop);qqline(evaled_isophase.demux_hap_count_TPM$liver_prop, col = 2)

#kidney
plot(density(evaled_isophase.demux_hap_count_TPM$kidney_prop,na.rm = TRUE))

shapiro.test(sample(evaled_isophase.demux_hap_count_TPM$kidney_prop,5000))

qqnorm(evaled_isophase.demux_hap_count_TPM$kidney_prop);qqline(evaled_isophase.demux_hap_count_TPM$kidney_prop, col = 2)

#brain
plot(density(evaled_isophase.demux_hap_count_TPM$brain_prop,na.rm = TRUE))

shapiro.test(sample(evaled_isophase.demux_hap_count_TPM$brain_prop,5000))

qqnorm(evaled_isophase.demux_hap_count_TPM$brain_prop);qqline(evaled_isophase.demux_hap_count_TPM$brain_prop, col = 2)

#lung
plot(density(evaled_isophase.demux_hap_count_TPM$lung_prop,na.rm = TRUE))

shapiro.test(sample(evaled_isophase.demux_hap_count_TPM$lung_prop,5000))

qqnorm(evaled_isophase.demux_hap_count_TPM$lung_prop);qqline(evaled_isophase.demux_hap_count_TPM$lung_prop, col = 2)

#muscle
plot(density(evaled_isophase.demux_hap_count_TPM$muscle_prop,na.rm = TRUE))

shapiro.test(sample(evaled_isophase.demux_hap_count_TPM$muscle_prop,5000))

qqnorm(evaled_isophase.demux_hap_count_TPM$muscle_prop);qqline(evaled_isophase.demux_hap_count_TPM$muscle_prop, col = 2)

#placenta
plot(density(evaled_isophase.demux_hap_count_TPM$placenta_prop,na.rm = TRUE))

shapiro.test(sample(evaled_isophase.demux_hap_count_TPM$placenta_prop,5000))

qqnorm(evaled_isophase.demux_hap_count_TPM$placenta_prop);qqline(evaled_isophase.demux_hap_count_TPM$placenta_prop, col = 2)

#multiplot of proportion density
par(mfrow=c(3,3))
plot(density(evaled_isophase.demux_hap_count_TPM$heart_prop,na.rm = TRUE),xlab="",main="heart")
plot(density(evaled_isophase.demux_hap_count_TPM$liver_prop,na.rm = TRUE),xlab="",main="liver")
plot(density(evaled_isophase.demux_hap_count_TPM$kidney_prop,na.rm = TRUE),xlab="",main="kidney")
plot(density(evaled_isophase.demux_hap_count_TPM$brain_prop,na.rm = TRUE),xlab="",main="brain")
plot(density(evaled_isophase.demux_hap_count_TPM$lung_prop,na.rm = TRUE),xlab="",main="lung")
plot(density(evaled_isophase.demux_hap_count_TPM$muscle_prop,na.rm = TRUE),xlab="",main="muscle")
plot(density(evaled_isophase.demux_hap_count_TPM$placenta_prop,na.rm = TRUE),xlab="",main="placenta")
par(mfrow=c(1,1))

par(mfrow=c(3,3))
hist(evaled_isophase.demux_hap_count_TPM$heart_prop,xlab="",main="heart")
hist(evaled_isophase.demux_hap_count_TPM$liver_prop,xlab="",main="liver")
hist(evaled_isophase.demux_hap_count_TPM$kidney_prop,xlab="",main="kidney")
hist(evaled_isophase.demux_hap_count_TPM$brain_prop,xlab="",main="brain")
hist(evaled_isophase.demux_hap_count_TPM$lung_prop,xlab="",main="lung")
hist(evaled_isophase.demux_hap_count_TPM$muscle_prop,xlab="",main="muscle")
hist(evaled_isophase.demux_hap_count_TPM$placenta_prop,xlab="",main="placenta")
par(mfrow=c(1,1))

#for publication
tiff(filename = "FigFinal_hist_isoseq_prop.tiff",width = 800,height = 800)
par(mfrow=c(3,3))
hist(evaled_isophase.demux_hap_count_TPM$heart_prop,xlab="",main="heart",breaks = seq(0,1,0.05),ylim = c(0,1200))
hist(evaled_isophase.demux_hap_count_TPM$liver_prop,xlab="",main="liver",breaks = seq(0,1,0.05),ylim = c(0,1200))
hist(evaled_isophase.demux_hap_count_TPM$kidney_prop,xlab="",main="kidney",breaks = seq(0,1,0.05),ylim = c(0,1200))
hist(evaled_isophase.demux_hap_count_TPM$brain_prop,xlab="",main="brain",breaks = seq(0,1,0.05),ylim = c(0,1200))
hist(evaled_isophase.demux_hap_count_TPM$lung_prop,xlab="",main="lung",breaks = seq(0,1,0.05),ylim = c(0,1200))
hist(evaled_isophase.demux_hap_count_TPM$muscle_prop,xlab="",main="muscle",breaks = seq(0,1,0.05),ylim = c(0,1200))
hist(evaled_isophase.demux_hap_count_TPM$placenta_prop,xlab="",main="placenta",breaks = seq(0,1,0.05),ylim = c(0,1200))
par(mfrow=c(1,1))
dev.off()

#ggplot violin for colours
heart_cat <- rep("heart",length(evaled_isophase.demux_hap_count_TPM$heart_prop))
liver_cat <- rep("liver",length(evaled_isophase.demux_hap_count_TPM$liver_prop))
kidney_cat <- rep("kidney",length(evaled_isophase.demux_hap_count_TPM$kidney_prop))
brain_cat <- rep("brain",length(evaled_isophase.demux_hap_count_TPM$brain_prop))
lung_cat <- rep("lung",length(evaled_isophase.demux_hap_count_TPM$lung_prop))
muscle_cat <- rep("muscle",length(evaled_isophase.demux_hap_count_TPM$muscle_prop))
placenta_cat <- rep("placenta",length(evaled_isophase.demux_hap_count_TPM$placenta_prop))

all_tissue_cat <- c(heart_cat,liver_cat,kidney_cat,brain_cat,lung_cat,muscle_cat,placenta_cat)

evaled_isophase.demux_hap_count_TPM_violinplot <- 
  c(evaled_isophase.demux_hap_count_TPM$heart_prop,evaled_isophase.demux_hap_count_TPM$liver_prop,
      evaled_isophase.demux_hap_count_TPM$kidney_prop,evaled_isophase.demux_hap_count_TPM$brain_prop,
      evaled_isophase.demux_hap_count_TPM$lung_prop,evaled_isophase.demux_hap_count_TPM$muscle_prop,
      evaled_isophase.demux_hap_count_TPM$placenta_prop)

evaled_isophase.demux_hap_count_TPM_violinplot_df <- 
  data.frame("Proportion" = evaled_isophase.demux_hap_count_TPM_violinplot, "Tissue" = all_tissue_cat)

tiff(filename = "FigFinal_hist_isoseq_prop_violin.tiff",width = 400,height = 250)
bp <- ggplot(evaled_isophase.demux_hap_count_TPM_violinplot_df, aes(x=Tissue, y=Proportion, group=Tissue))
bp <- bp + geom_violin(aes(fill=Tissue)) + guides(fill=FALSE)
bp <- bp + coord_flip()
bp <- bp + ylab("Proportion of Brahman allele") + theme_bw()
bp <- bp + theme(plot.title = element_text(hjust = 0.5),
                 axis.text.x = element_text(angle = 0, hjust = 0.5,color = "black"),
                 axis.text.y = element_text(color = "black")) 
# bp <- bp + scale_y_log10(labels = scales::comma,breaks=c(0,5000,10000,15000,20000,25000,30000,40000,50000,60000,80000,100000,120000,150000,250000,500000))
# bp <- bp + facet_grid(. ~ family) + theme_bw(base_size = 13.5)
# bp <- bp + theme(strip.text.x = element_text(size = 20),element_line(colour = "black"))
# bp <- bp + ylab("Length of reference matched to repeat (bp)") + xlab("Assembly")
bp
dev.off()

##### summary of proportion values distribution #####
summary(evaled_isophase.demux_hap_count_TPM$heart_prop)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.0000  0.4286  0.5000  0.5083  0.5964  1.0000     540 
summary(evaled_isophase.demux_hap_count_TPM$liver_prop)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.0000  0.2000  0.5000  0.5102  1.0000  1.0000    2861 
summary(evaled_isophase.demux_hap_count_TPM$kidney_prop)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.0000  0.4286  0.5000  0.5045  0.5833  1.0000     385 
summary(evaled_isophase.demux_hap_count_TPM$brain_prop)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.0000  0.4364  0.5000  0.5087  0.5771  1.0000     110 
summary(evaled_isophase.demux_hap_count_TPM$lung_prop)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.0000  0.3333  0.5000  0.5066  0.6667  1.0000    1167 
summary(evaled_isophase.demux_hap_count_TPM$muscle_prop)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.0000  0.3333  0.5000  0.5187  0.7500  1.0000    1522 
summary(evaled_isophase.demux_hap_count_TPM$placenta_prop)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.0000  0.2000  0.5000  0.5106  0.8571  1.0000    2506

################################################################################################
#####Qnorm
#heart
plot(density(evaled_isophase.demux_hap_count_Qnorm$heart_prop,na.rm = TRUE))

shapiro.test(sample(evaled_isophase.demux_hap_count_Qnorm$heart_prop,5000))

qqnorm(evaled_isophase.demux_hap_count_Qnorm$heart_prop);qqline(evaled_isophase.demux_hap_count_Qnorm$heart_prop, col = 2)

#liver
plot(density(evaled_isophase.demux_hap_count_Qnorm$liver_prop,na.rm = TRUE))

shapiro.test(sample(evaled_isophase.demux_hap_count_Qnorm$liver_prop,5000))

qqnorm(evaled_isophase.demux_hap_count_Qnorm$liver_prop);qqline(evaled_isophase.demux_hap_count_Qnorm$liver_prop, col = 2)

#kidney
plot(density(evaled_isophase.demux_hap_count_Qnorm$kidney_prop,na.rm = TRUE))

shapiro.test(sample(evaled_isophase.demux_hap_count_Qnorm$kidney_prop,5000))

qqnorm(evaled_isophase.demux_hap_count_Qnorm$kidney_prop);qqline(evaled_isophase.demux_hap_count_Qnorm$kidney_prop, col = 2)

#brain
plot(density(evaled_isophase.demux_hap_count_Qnorm$brain_prop,na.rm = TRUE))

shapiro.test(sample(evaled_isophase.demux_hap_count_Qnorm$brain_prop,5000))

qqnorm(evaled_isophase.demux_hap_count_Qnorm$brain_prop);qqline(evaled_isophase.demux_hap_count_Qnorm$brain_prop, col = 2)

#lung
plot(density(evaled_isophase.demux_hap_count_Qnorm$lung_prop,na.rm = TRUE))

shapiro.test(sample(evaled_isophase.demux_hap_count_Qnorm$lung_prop,5000))

qqnorm(evaled_isophase.demux_hap_count_Qnorm$lung_prop);qqline(evaled_isophase.demux_hap_count_Qnorm$lung_prop, col = 2)

#muscle
plot(density(evaled_isophase.demux_hap_count_Qnorm$muscle_prop,na.rm = TRUE))

shapiro.test(sample(evaled_isophase.demux_hap_count_Qnorm$muscle_prop,5000))

qqnorm(evaled_isophase.demux_hap_count_Qnorm$muscle_prop);qqline(evaled_isophase.demux_hap_count_Qnorm$muscle_prop, col = 2)

#placenta
plot(density(evaled_isophase.demux_hap_count_Qnorm$placenta_prop,na.rm = TRUE))

shapiro.test(sample(evaled_isophase.demux_hap_count_Qnorm$placenta_prop,5000))

qqnorm(evaled_isophase.demux_hap_count_Qnorm$placenta_prop);qqline(evaled_isophase.demux_hap_count_Qnorm$placenta_prop, col = 2)

#multiplot of proportion density
par(mfrow=c(3,3))
plot(density(evaled_isophase.demux_hap_count_Qnorm$heart_prop,na.rm = TRUE),xlab="",main="heart")
plot(density(evaled_isophase.demux_hap_count_Qnorm$liver_prop,na.rm = TRUE),xlab="",main="liver")
plot(density(evaled_isophase.demux_hap_count_Qnorm$kidney_prop,na.rm = TRUE),xlab="",main="kidney")
plot(density(evaled_isophase.demux_hap_count_Qnorm$brain_prop,na.rm = TRUE),xlab="",main="brain")
plot(density(evaled_isophase.demux_hap_count_Qnorm$lung_prop,na.rm = TRUE),xlab="",main="lung")
plot(density(evaled_isophase.demux_hap_count_Qnorm$muscle_prop,na.rm = TRUE),xlab="",main="muscle")
plot(density(evaled_isophase.demux_hap_count_Qnorm$placenta_prop,na.rm = TRUE),xlab="",main="placenta")
par(mfrow=c(1,1))

par(mfrow=c(3,3))
hist(evaled_isophase.demux_hap_count_Qnorm$heart_prop,xlab="",main="heart")
hist(evaled_isophase.demux_hap_count_Qnorm$liver_prop,xlab="",main="liver")
hist(evaled_isophase.demux_hap_count_Qnorm$kidney_prop,xlab="",main="kidney")
hist(evaled_isophase.demux_hap_count_Qnorm$brain_prop,xlab="",main="brain")
hist(evaled_isophase.demux_hap_count_Qnorm$lung_prop,xlab="",main="lung")
hist(evaled_isophase.demux_hap_count_Qnorm$muscle_prop,xlab="",main="muscle")
hist(evaled_isophase.demux_hap_count_Qnorm$placenta_prop,xlab="",main="placenta")
par(mfrow=c(1,1))

##### summary of proportion values distribution #####
summary(evaled_isophase.demux_hap_count_Qnorm$heart_prop)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.01342 0.43318 0.50000 0.50150 0.57399 0.98246 
summary(evaled_isophase.demux_hap_count_Qnorm$liver_prop)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.04741 0.48120 0.48120 0.49095 0.49660 0.97642 
summary(evaled_isophase.demux_hap_count_Qnorm$kidney_prop)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.02273 0.42190 0.50000 0.50015 0.57965 0.99649 
summary(evaled_isophase.demux_hap_count_Qnorm$brain_prop)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.0000  0.4248  0.5000  0.5026  0.5775  1.0000     110 
summary(evaled_isophase.demux_hap_count_Qnorm$lung_prop)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0243  0.3996  0.4815  0.4968  0.6065  0.9836 
summary(evaled_isophase.demux_hap_count_Qnorm$muscle_prop)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.01946 0.40317 0.47368 0.49473 0.62103 0.96805 
summary(evaled_isophase.demux_hap_count_Qnorm$placenta_prop)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.01002 0.49145 0.50000 0.50118 0.50554 0.96931 

################################################################################################
#####allelic imbalance? high in Brahman, TPM
#heart
evaled_isophase.demux_hap_count_TPM_heart <- evaled_isophase.demux_hap_count_TPM %>% filter(heart_prop >= 0.85) %>%
  select(locus:phase1,heart_p0,heart_p1,heart_prop) %>% arrange(desc(heart_p0))

#liver
evaled_isophase.demux_hap_count_TPM_liver <- evaled_isophase.demux_hap_count_TPM %>% filter(liver_prop >= 0.85) %>%
  select(locus:phase1,liver_p0,liver_p1,liver_prop) %>% arrange(desc(liver_p0))

#kidney
evaled_isophase.demux_hap_count_TPM_kidney <- evaled_isophase.demux_hap_count_TPM %>% filter(kidney_prop >= 0.85) %>%
  select(locus:phase1,kidney_p0,kidney_p1,kidney_prop) %>% arrange(desc(kidney_p0))

#brain
evaled_isophase.demux_hap_count_TPM_brain <- evaled_isophase.demux_hap_count_TPM %>% filter(brain_prop >= 0.85) %>%
  select(locus:phase1,brain_p0,brain_p1,brain_prop) %>% arrange(desc(brain_p0))

#lung
evaled_isophase.demux_hap_count_TPM_lung <- evaled_isophase.demux_hap_count_TPM %>% filter(lung_prop >= 0.85) %>%
  select(locus:phase1,lung_p0,lung_p1,lung_prop) %>% arrange(desc(lung_p0))

#muscle
evaled_isophase.demux_hap_count_TPM_muscle <- evaled_isophase.demux_hap_count_TPM %>% filter(muscle_prop >= 0.85) %>%
  select(locus:phase1,muscle_p0,muscle_p1,muscle_prop) %>% arrange(desc(muscle_p0))

#placenta
evaled_isophase.demux_hap_count_TPM_placenta <- evaled_isophase.demux_hap_count_TPM %>% filter(placenta_prop >= 0.85) %>%
  select(locus:phase1,placenta_p0,placenta_p1,placenta_prop) %>% arrange(desc(placenta_p0))

#####Top 10 high in Brahman in 7 tissues
#heart
evaled_isophase.demux_hap_count_TPM_heart_top10 <- head(evaled_isophase.demux_hap_count_TPM_heart$locus,10)
tissue_evaled_isophase.demux_hap_count_TPM_heart_top10 <- rep("heart",10)

evaled_isophase.demux_hap_count_TPM_heart_top10_df <- 
  cbind(evaled_isophase.demux_hap_count_TPM_heart_top10,tissue_evaled_isophase.demux_hap_count_TPM_heart_top10)

#liver
evaled_isophase.demux_hap_count_TPM_liver_top10 <- head(evaled_isophase.demux_hap_count_TPM_liver$locus,10)
tissue_evaled_isophase.demux_hap_count_TPM_liver_top10 <- rep("liver",10)

evaled_isophase.demux_hap_count_TPM_liver_top10_df <- 
  cbind(evaled_isophase.demux_hap_count_TPM_liver_top10,tissue_evaled_isophase.demux_hap_count_TPM_liver_top10)

#kidney
evaled_isophase.demux_hap_count_TPM_kidney_top10 <- head(evaled_isophase.demux_hap_count_TPM_kidney$locus,10)
tissue_evaled_isophase.demux_hap_count_TPM_kidney_top10 <- rep("kidney",10)

evaled_isophase.demux_hap_count_TPM_kidney_top10_df <- 
  cbind(evaled_isophase.demux_hap_count_TPM_kidney_top10,tissue_evaled_isophase.demux_hap_count_TPM_kidney_top10)

#brain
evaled_isophase.demux_hap_count_TPM_brain_top10 <- head(evaled_isophase.demux_hap_count_TPM_brain$locus,10)
tissue_evaled_isophase.demux_hap_count_TPM_brain_top10 <- rep("brain",10)

evaled_isophase.demux_hap_count_TPM_brain_top10_df <- 
  cbind(evaled_isophase.demux_hap_count_TPM_brain_top10,tissue_evaled_isophase.demux_hap_count_TPM_brain_top10)

#lung
evaled_isophase.demux_hap_count_TPM_lung_top10 <- head(evaled_isophase.demux_hap_count_TPM_lung$locus,10)
tissue_evaled_isophase.demux_hap_count_TPM_lung_top10 <- rep("lung",10)

evaled_isophase.demux_hap_count_TPM_lung_top10_df <- 
  cbind(evaled_isophase.demux_hap_count_TPM_lung_top10,tissue_evaled_isophase.demux_hap_count_TPM_lung_top10)

#muscle
evaled_isophase.demux_hap_count_TPM_muscle_top10 <- head(evaled_isophase.demux_hap_count_TPM_muscle$locus,10)
tissue_evaled_isophase.demux_hap_count_TPM_muscle_top10 <- rep("muscle",10)

evaled_isophase.demux_hap_count_TPM_muscle_top10_df <- 
  cbind(evaled_isophase.demux_hap_count_TPM_muscle_top10,tissue_evaled_isophase.demux_hap_count_TPM_muscle_top10)

#placenta
evaled_isophase.demux_hap_count_TPM_placenta_top10 <- head(evaled_isophase.demux_hap_count_TPM_placenta$locus,10)
tissue_evaled_isophase.demux_hap_count_TPM_placenta_top10 <- rep("placenta",10)

evaled_isophase.demux_hap_count_TPM_placenta_top10_df <- 
  cbind(evaled_isophase.demux_hap_count_TPM_placenta_top10,tissue_evaled_isophase.demux_hap_count_TPM_placenta_top10)

#combine all top 10 brahman high
top10_Brahman_high <- 
  rbind(evaled_isophase.demux_hap_count_TPM_heart_top10_df,evaled_isophase.demux_hap_count_TPM_liver_top10_df,
      evaled_isophase.demux_hap_count_TPM_kidney_top10_df,evaled_isophase.demux_hap_count_TPM_brain_top10_df,
      evaled_isophase.demux_hap_count_TPM_lung_top10_df,evaled_isophase.demux_hap_count_TPM_muscle_top10_df,
      evaled_isophase.demux_hap_count_TPM_placenta_top10_df)

colnames(top10_Brahman_high) <- c("locus","tissue")

top10_Brahman_high_groupby_locus <- as.tbl(as.data.frame(top10_Brahman_high))

top10_Brahman_high_groupby_locus_df <- as.data.frame.matrix(table(top10_Brahman_high_groupby_locus$locus,
                                                                     top10_Brahman_high_groupby_locus$tissue))

top10_Brahman_high_groupby_locus_df$sum_all_tissue <- rowSums(top10_Brahman_high_groupby_locus_df)

################################################################################################
################################################################################################
#####allelic imbalance? high in Angus, TPM
#heart
evaled_isophase.demux_hap_count_TPM_low_heart <- evaled_isophase.demux_hap_count_TPM %>% filter(heart_prop <= 0.15) %>%
  select(locus:phase1,heart_p0,heart_p1,heart_prop) %>% arrange(desc(heart_p1))

#liver
evaled_isophase.demux_hap_count_TPM_low_liver <- evaled_isophase.demux_hap_count_TPM %>% filter(liver_prop <= 0.15) %>%
  select(locus:phase1,liver_p0,liver_p1,liver_prop) %>% arrange(desc(liver_p1))

#kidney
evaled_isophase.demux_hap_count_TPM_low_kidney <- evaled_isophase.demux_hap_count_TPM %>% filter(kidney_prop <= 0.15) %>%
  select(locus:phase1,kidney_p0,kidney_p1,kidney_prop) %>% arrange(desc(kidney_p1))

#brain
evaled_isophase.demux_hap_count_TPM_low_brain <- evaled_isophase.demux_hap_count_TPM %>% filter(brain_prop <= 0.15) %>%
  select(locus:phase1,brain_p0,brain_p1,brain_prop) %>% arrange(desc(brain_p1))

#lung
evaled_isophase.demux_hap_count_TPM_low_lung <- evaled_isophase.demux_hap_count_TPM %>% filter(lung_prop <= 0.15) %>%
  select(locus:phase1,lung_p0,lung_p1,lung_prop) %>% arrange(desc(lung_p1))

#muscle
evaled_isophase.demux_hap_count_TPM_low_muscle <- evaled_isophase.demux_hap_count_TPM %>% filter(muscle_prop <= 0.15) %>%
  select(locus:phase1,muscle_p0,muscle_p1,muscle_prop) %>% arrange(desc(muscle_p1))

#placenta
evaled_isophase.demux_hap_count_TPM_low_placenta <- evaled_isophase.demux_hap_count_TPM %>% filter(placenta_prop <= 0.15) %>%
  select(locus:phase1,placenta_p0,placenta_p1,placenta_prop) %>% arrange(desc(placenta_p1))

#####Top 10 high in Angus in 7 tissues
#heart
evaled_isophase.demux_hap_count_TPM_low_heart_top10 <- head(evaled_isophase.demux_hap_count_TPM_low_heart$locus,10)
tissue_evaled_isophase.demux_hap_count_TPM_low_heart_top10 <- rep("heart",10)

evaled_isophase.demux_hap_count_TPM_low_heart_top10_df <- 
  cbind(evaled_isophase.demux_hap_count_TPM_low_heart_top10,tissue_evaled_isophase.demux_hap_count_TPM_low_heart_top10)

#liver
evaled_isophase.demux_hap_count_TPM_low_liver_top10 <- head(evaled_isophase.demux_hap_count_TPM_low_liver$locus,10)
tissue_evaled_isophase.demux_hap_count_TPM_low_liver_top10 <- rep("liver",10)

evaled_isophase.demux_hap_count_TPM_low_liver_top10_df <- 
  cbind(evaled_isophase.demux_hap_count_TPM_low_liver_top10,tissue_evaled_isophase.demux_hap_count_TPM_low_liver_top10)

#kidney
evaled_isophase.demux_hap_count_TPM_low_kidney_top10 <- head(evaled_isophase.demux_hap_count_TPM_low_kidney$locus,10)
tissue_evaled_isophase.demux_hap_count_TPM_low_kidney_top10 <- rep("kidney",10)

evaled_isophase.demux_hap_count_TPM_low_kidney_top10_df <- 
  cbind(evaled_isophase.demux_hap_count_TPM_low_kidney_top10,tissue_evaled_isophase.demux_hap_count_TPM_low_kidney_top10)

#brain
evaled_isophase.demux_hap_count_TPM_low_brain_top10 <- head(evaled_isophase.demux_hap_count_TPM_low_brain$locus,10)
tissue_evaled_isophase.demux_hap_count_TPM_low_brain_top10 <- rep("brain",10)

evaled_isophase.demux_hap_count_TPM_low_brain_top10_df <- 
  cbind(evaled_isophase.demux_hap_count_TPM_low_brain_top10,tissue_evaled_isophase.demux_hap_count_TPM_low_brain_top10)

#lung
evaled_isophase.demux_hap_count_TPM_low_lung_top10 <- head(evaled_isophase.demux_hap_count_TPM_low_lung$locus,10)
tissue_evaled_isophase.demux_hap_count_TPM_low_lung_top10 <- rep("lung",10)

evaled_isophase.demux_hap_count_TPM_low_lung_top10_df <- 
  cbind(evaled_isophase.demux_hap_count_TPM_low_lung_top10,tissue_evaled_isophase.demux_hap_count_TPM_low_lung_top10)

#muscle
evaled_isophase.demux_hap_count_TPM_low_muscle_top10 <- head(evaled_isophase.demux_hap_count_TPM_low_muscle$locus,10)
tissue_evaled_isophase.demux_hap_count_TPM_low_muscle_top10 <- rep("muscle",10)

evaled_isophase.demux_hap_count_TPM_low_muscle_top10_df <- 
  cbind(evaled_isophase.demux_hap_count_TPM_low_muscle_top10,tissue_evaled_isophase.demux_hap_count_TPM_low_muscle_top10)

#placenta
evaled_isophase.demux_hap_count_TPM_low_placenta_top10 <- head(evaled_isophase.demux_hap_count_TPM_low_placenta$locus,10)
tissue_evaled_isophase.demux_hap_count_TPM_low_placenta_top10 <- rep("placenta",10)

evaled_isophase.demux_hap_count_TPM_low_placenta_top10_df <- 
  cbind(evaled_isophase.demux_hap_count_TPM_low_placenta_top10,tissue_evaled_isophase.demux_hap_count_TPM_low_placenta_top10)

#combine all top 10 Angus high
top10_Angus_high <- 
  rbind(evaled_isophase.demux_hap_count_TPM_low_heart_top10_df,evaled_isophase.demux_hap_count_TPM_low_liver_top10_df,
        evaled_isophase.demux_hap_count_TPM_low_kidney_top10_df,evaled_isophase.demux_hap_count_TPM_low_brain_top10_df,
        evaled_isophase.demux_hap_count_TPM_low_lung_top10_df,evaled_isophase.demux_hap_count_TPM_low_muscle_top10_df,
        evaled_isophase.demux_hap_count_TPM_low_placenta_top10_df)

colnames(top10_Angus_high) <- c("locus","tissue")

top10_Angus_high_groupby_locus <- as.tbl(as.data.frame(top10_Angus_high))

top10_Angus_high_groupby_locus_df <- as.data.frame.matrix(table(top10_Angus_high_groupby_locus$locus,
                                                                  top10_Angus_high_groupby_locus$tissue))

top10_Angus_high_groupby_locus_df$sum_all_tissue <- rowSums(top10_Angus_high_groupby_locus_df)

################################################################################################

