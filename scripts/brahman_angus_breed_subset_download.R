#------------------------------------------------------
# Program name: brahman_angus_breed_subset_download.R
# Objective: subset the runInfo table on breed for the
#           purpose of download
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(readr)
library(dplyr)

breed <- read_tsv("runInfo/MBCDPv29_10xWGS_PublicLinksByBreed26RunInfo.txt")

# Getting each of the 7 breeds table out
# Brahman
breed_brahman <- breed %>% filter(breed == "Brahman") %>% arrange(isolate) %>% 
  select(Run:tissue)

write.table(breed_brahman,file = "runInfo/breed_brahman.tsv",sep="\t", row.names = FALSE,
            col.names = FALSE)

# Angus
breed_angus <- breed %>% filter(breed == "Angus") %>% arrange(isolate) %>% 
  select(Run:tissue)

write.table(breed_angus,file = "runInfo/breed_angus.tsv",sep="\t", row.names = FALSE,
            col.names = FALSE)

# RedAngus
breed_redangus <- breed %>% filter(breed == "Red Angus") %>% arrange(isolate) %>% 
  select(Run:tissue)

write.table(breed_redangus,file = "runInfo/breed_redangus.tsv",sep="\t", row.names = FALSE,
            col.names = FALSE)

# Gelbvieh
breed_gelbvieh <- breed %>% filter(breed == "Gelbvieh") %>% arrange(isolate) %>% 
  select(Run:tissue)

write.table(breed_gelbvieh,file = "runInfo/breed_gelbvieh.tsv",sep="\t", row.names = FALSE,
            col.names = FALSE)

# Hereford
breed_hereford <- breed %>% filter(breed == "Hereford") %>% arrange(isolate) %>% 
  select(Run:tissue)

write.table(breed_hereford,file = "runInfo/breed_hereford.tsv",sep="\t", row.names = FALSE,
            col.names = FALSE)

# Shorthorn
breed_shorthorn <- breed %>% filter(breed == "Shorthorn") %>% arrange(isolate) %>% 
  select(Run:tissue)

write.table(breed_shorthorn,file = "runInfo/breed_shorthorn.tsv",sep="\t", row.names = FALSE,
            col.names = FALSE)

# Simmental
breed_simmental <- breed %>% filter(breed == "Simmental") %>% arrange(isolate) %>% 
  select(Run:tissue)

write.table(breed_simmental,file = "runInfo/breed_simmental.tsv",sep="\t", row.names = FALSE,
            col.names = FALSE)

#### Writing a function to generate all breeds tsv ####
uniq <- unique(breed$breed)

tsv_generator <- function(breed_obj){
  for (i in 1:length(breed_obj)){
    breed_DF <- breed %>% filter(breed == breed_obj[i]) %>% arrange(isolate) %>% 
      select(Run:tissue)
    breed_DF$breed <- gsub("![:alnum:]","",breed_DF$breed)
    breed_DF$isolate <- gsub("![:alnum:]","",breed_DF$isolate)
    breed_name <- gsub("![:alnum:]","",breed_obj[i])
    filename <- paste0("runInfo/breed_",breed_name,".tsv")
    write.table(breed_DF,file = filename,sep="\t", row.names = FALSE,col.names = FALSE)
  }
}

#After learning the new 49 Brahman seq available, I decided to download them as a backup
#for selective sweep analysis
MooreBrahman <- read_tsv("steveMooreBrahman_runinfo/PRJNA432125_runinfo.txt")

#filter for those with cov more than 10 based on genome size of 2.8 Gb
#even though brahman is probably closer to 2.7 Gb bcos I allow for collapsed repeats
MooreBrahman_10cov <- MooreBrahman %>% mutate(cov = MBases/2.8e3) %>% 
  filter(cov >= 10) %>% select(Run:isolate)

write.table(MooreBrahman_10cov,file = "steveMooreBrahman_runinfo/MooreBrahman.tsv",sep="\t", row.names = FALSE,
            col.names = FALSE)
