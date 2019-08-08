#------------------------------------------------------
# Program name: brahman_angus_convert_selectionRegion_to_genotype.R
# Objective: instead of looking at 0/1 etc, change it to
#           genotype across chr position
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(readr)
library(dplyr)
library(tidyr)

#right now, use in conjuction with brahman_angus_vcf_score_alt_allele_BrahmanRef.R to get the per chr
#combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split

#assume df in this format 
#CHROM POS REF ALT <an1_GT> <an2_GT> ... all_geno_var
#show with chr 7
df1 <- combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split %>% filter(POS >= 66100000) %>% filter(POS < 66200000) %>% 
  dplyr::select(CHROM,POS,REF,ALT,an1_GT,an2_GT,an3_GT,an4_GT,an5_GT,an6_GT,br1_GT,br2_GT,br3_GT,br4_GT,br5_GT)

#another with chr 19
df1 <- combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split %>% filter(POS >= 27700000) %>% filter(POS < 27800000) %>% 
  select(CHROM,POS,REF,ALT,an1_GT,an2_GT,an3_GT,an4_GT,an5_GT,an6_GT,br1_GT,br2_GT,br3_GT,br4_GT,br5_GT)

#another with chr 24
df1 <- combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split %>% filter(POS >= 400000) %>% filter(POS < 800000) %>% 
  select(CHROM,POS,REF,ALT,an1_GT,an2_GT,an3_GT,an4_GT,an5_GT,an6_GT,br1_GT,br2_GT,br3_GT,br4_GT,br5_GT)

#with all SevenBreeds data
#show with chr 7
df1 <- combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split %>% filter(POS >= 66100000) %>% filter(POS < 66200000) %>% 
  dplyr::select(CHROM,POS,REF,ALT,an1_GT,an2_GT,an3_GT,an4_GT,an5_GT,an6_GT,br1_GT,br2_GT,br3_GT,br4_GT,br5_GT,
         ge1_GT,ge2_GT,ge3_GT,ge4_GT,ge5_GT,ge6_GT,he1_GT,he2_GT,he3_GT,he4_GT,he5_GT,he6_GT,
         re2_GT,re3_GT,re4_GT,re5_GT,re6_GT,sh1_GT,sh2_GT,sh3_GT,sh4_GT,sh5_GT,
         si1_GT,si2_GT,si3_GT,si4_GT,si5_GT)

#show PLAG1 on chr 14 
df1 <- combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split %>% filter(POS >= 23633375-1e5) %>% filter(POS < 23685278+1e5) %>% 
  select(CHROM,POS,REF,ALT,an1_GT,an2_GT,an3_GT,an4_GT,an5_GT,an6_GT,br1_GT,br2_GT,br3_GT,br4_GT,br5_GT,
         ge1_GT,ge2_GT,ge3_GT,ge4_GT,ge5_GT,ge6_GT,he1_GT,he2_GT,he3_GT,he4_GT,he5_GT,he6_GT,
         re2_GT,re3_GT,re4_GT,re5_GT,re6_GT,sh1_GT,sh2_GT,sh3_GT,sh4_GT,sh5_GT,
         si1_GT,si2_GT,si3_GT,si4_GT,si5_GT)

#show IGF1R with chr 21
df1 <- combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split %>% filter(POS >= 7800000) %>% filter(POS < 7900000) %>% 
  dplyr::select(CHROM,POS,REF,ALT,an1_GT,an2_GT,an3_GT,an4_GT,an5_GT,an6_GT,br1_GT,br2_GT,br3_GT,br4_GT,br5_GT,
                ge1_GT,ge2_GT,ge3_GT,ge4_GT,ge5_GT,ge6_GT,he1_GT,he2_GT,he3_GT,he4_GT,he5_GT,he6_GT,
                re2_GT,re3_GT,re4_GT,re5_GT,re6_GT,sh1_GT,sh2_GT,sh3_GT,sh4_GT,sh5_GT,
                si1_GT,si2_GT,si3_GT,si4_GT,si5_GT)

#the three 200-kb selection interval
#chr8
df1 <- combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split %>% filter(POS >= 57500000) %>% filter(POS < 57700000) %>% 
  dplyr::select(CHROM,POS,REF,ALT,an1_GT,an2_GT,an3_GT,an4_GT,an5_GT,an6_GT,br1_GT,br2_GT,br3_GT,br4_GT,br5_GT,
                ge1_GT,ge2_GT,ge3_GT,ge4_GT,ge5_GT,ge6_GT,he1_GT,he2_GT,he3_GT,he4_GT,he5_GT,he6_GT,
                re2_GT,re3_GT,re4_GT,re5_GT,re6_GT,sh1_GT,sh2_GT,sh3_GT,sh4_GT,sh5_GT,
                si1_GT,si2_GT,si3_GT,si4_GT,si5_GT)

#chr11
df1 <- combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split %>% filter(POS >= 45700000) %>% filter(POS < 45900000) %>% 
  dplyr::select(CHROM,POS,REF,ALT,an1_GT,an2_GT,an3_GT,an4_GT,an5_GT,an6_GT,br1_GT,br2_GT,br3_GT,br4_GT,br5_GT,
                ge1_GT,ge2_GT,ge3_GT,ge4_GT,ge5_GT,ge6_GT,he1_GT,he2_GT,he3_GT,he4_GT,he5_GT,he6_GT,
                re2_GT,re3_GT,re4_GT,re5_GT,re6_GT,sh1_GT,sh2_GT,sh3_GT,sh4_GT,sh5_GT,
                si1_GT,si2_GT,si3_GT,si4_GT,si5_GT)

#loop thro each row
for (i in 1:nrow(df1)){
  # geno0 <- df1$REF[i]
  geno0 <- "."
  #separate(an1, c("an1_GT","an1_AD","an1_DP","an1_GQ","an1_PL"),sep = ":",extra = "drop", fill = "right")
  alt <- unlist(strsplit(df1$ALT[i], ","))
  geno1 <- alt[1]
  geno2 <- alt[2]
  geno3 <- alt[3]
  
  #check genotype status in var such as an1_GT and replace with the right geno<no>
  #"0/0","0/1","0/2","0/3","1/1","1/2","1/3","2/2","2/3","3/3"
  for (j in 5:ncol(df1)){
    #check if genotype column for that row in particular cell, starting with row 1, col 5
    #match any of the 10 genotypes. It must match at least one.
    if (df1[i,j] == "0/0"){
      df1[i,j] <- paste0(geno0,geno0)
    } else if (df1[i,j] == "0/1") {
      df1[i,j] <- paste0(geno0,geno1)
    } else if (df1[i,j] == "0/2") {
      df1[i,j] <- paste0(geno0,geno2)
    } else if (df1[i,j] == "0/3") {
      df1[i,j] <- paste0(geno0,geno3)
    } else if (df1[i,j] == "1/1") {
      df1[i,j] <- paste0(geno1,geno1)
    } else if (df1[i,j] == "1/2") {
      df1[i,j] <- paste0(geno1,geno2)
    } else if (df1[i,j] == "1/3") {
      df1[i,j] <- paste0(geno1,geno3)
    } else if (df1[i,j] == "2/2") {
      df1[i,j] <- paste0(geno2,geno2)
    } else if (df1[i,j] == "2/3") {
      df1[i,j] <- paste0(geno2,geno3)
    } else if (df1[i,j] == "3/3") {
      df1[i,j] <- paste0(geno3,geno3)
    } else {
      df1[i,j] <- "NA"
    }
  }
}

#test whether the code works for chr7 hspa4 region
#t <- df1 %>% filter(POS >= 66100000) %>% filter(POS < 66200000)
t.select <- df1 %>% dplyr::select(POS,REF,an1_GT,an2_GT,an3_GT,an4_GT,an5_GT,an6_GT,br1_GT,br2_GT,br3_GT,br4_GT,br5_GT)
t(t.select)

#with the SevenBreeds data
#t <- df1 %>% filter(POS >= 66100000) %>% filter(POS < 66200000)
t.select <- df1 %>% dplyr::select(POS,REF,an1_GT:si5_GT)
df_transpose <- as.data.frame(t(t.select))

write_tsv(df_transpose,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Brahman_SevenBreeds/chr7_66100000_66200000_BrahmanRef_HSPA4.tsv")

