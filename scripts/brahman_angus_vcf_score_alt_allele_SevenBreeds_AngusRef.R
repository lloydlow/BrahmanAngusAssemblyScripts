#------------------------------------------------------
# Program name: brahman_angus_vcf_score_alt_allele_SevenBreeds_AngusRef.R
# Objective: tested code to detect extended haplotype homozygosity
#           and apply to after annovar results now. This one for Angus ref.
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(readr)
library(dplyr)
library(tidyr)

#need to run subset_brahman_vcf_by_chr.sh to get subset by chr first
#folder: /Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus_SevenBreeds

#get chr no
chr <- gsub("SevenBreedsRefAngus_filteredSNP.Angus_multianno_","",vcfFile)
chr <- gsub(".vcf","",chr)

# path to annotated vcf
dir1 <- 
  "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Angus_SevenBreeds/"

# reading SevenBreedsRefAngus_filteredSNP.Angus_multianno_19.vcf
#path1 <- paste0(dir1,"SevenBreedsRefAngus_filteredSNP.Angus_multianno_19.vcf")

path1 <- paste0(dir1,"SevenBreedsRefAngus_filteredSNP.Angus_multianno_19.vcf")

combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1 <- 
  read_tsv(path1,col_names = FALSE)
names(combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1) <- 
  c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","an1","an2","an3","an4","an5",
    "an6","br1","br2","br3","br4","br5","ge1","ge2","ge3","ge4","ge5","ge6","he1","he2","he3","he4",
    "he5","he6","re2","re3","re4","re5","re6","sh1","sh2","sh3","sh4","sh5","si1","si2","si3","si4",
    "si5")

# combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1 <- combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1 %>%
#   filter(FILTER == "PASS")

#split based on an1-an5, br1-br6 etc
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split <- combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1 %>% 
  separate(an1, c("an1_GT","an1_AD","an1_DP","an1_GQ","an1_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(an2, c("an2_GT","an2_AD","an2_DP","an2_GQ","an2_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(an3, c("an3_GT","an3_AD","an3_DP","an3_GQ","an3_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(an4, c("an4_GT","an4_AD","an4_DP","an4_GQ","an4_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(an5, c("an5_GT","an5_AD","an5_DP","an5_GQ","an5_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(an6, c("an6_GT","an6_AD","an6_DP","an6_GQ","an6_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(br1, c("br1_GT","br1_AD","br1_DP","br1_GQ","br1_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(br2, c("br2_GT","br2_AD","br2_DP","br2_GQ","br2_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(br3, c("br3_GT","br3_AD","br3_DP","br3_GQ","br3_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(br4, c("br4_GT","br4_AD","br4_DP","br4_GQ","br4_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(br5, c("br5_GT","br5_AD","br5_DP","br5_GQ","br5_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(ge1, c("ge1_GT","ge1_AD","ge1_DP","ge1_GQ","ge1_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(ge2, c("ge2_GT","ge2_AD","ge2_DP","ge2_GQ","ge2_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(ge3, c("ge3_GT","ge3_AD","ge3_DP","ge3_GQ","ge3_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(ge4, c("ge4_GT","ge4_AD","ge4_DP","ge4_GQ","ge4_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(ge5, c("ge5_GT","ge5_AD","ge5_DP","ge5_GQ","ge5_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(ge6, c("ge6_GT","ge6_AD","ge6_DP","ge6_GQ","ge6_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(he1, c("he1_GT","he1_AD","he1_DP","he1_GQ","he1_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(he2, c("he2_GT","he2_AD","he2_DP","he2_GQ","he2_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(he3, c("he3_GT","he3_AD","he3_DP","he3_GQ","he3_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(he4, c("he4_GT","he4_AD","he4_DP","he4_GQ","he4_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(he5, c("he5_GT","he5_AD","he5_DP","he5_GQ","he5_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(he6, c("he6_GT","he6_AD","he6_DP","he6_GQ","he6_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(re2, c("re2_GT","re2_AD","re2_DP","re2_GQ","re2_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(re3, c("re3_GT","re3_AD","re3_DP","re3_GQ","re3_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(re4, c("re4_GT","re4_AD","re4_DP","re4_GQ","re4_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(re5, c("re5_GT","re5_AD","re5_DP","re5_GQ","re5_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(re6, c("re6_GT","re6_AD","re6_DP","re6_GQ","re6_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(sh1, c("sh1_GT","sh1_AD","sh1_DP","sh1_GQ","sh1_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(sh2, c("sh2_GT","sh2_AD","sh2_DP","sh2_GQ","sh2_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(sh3, c("sh3_GT","sh3_AD","sh3_DP","sh3_GQ","sh3_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(sh4, c("sh4_GT","sh4_AD","sh4_DP","sh4_GQ","sh4_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(sh5, c("sh5_GT","sh5_AD","sh5_DP","sh5_GQ","sh5_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(si1, c("si1_GT","si1_AD","si1_DP","si1_GQ","si1_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(si2, c("si2_GT","si2_AD","si2_DP","si2_GQ","si2_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(si3, c("si3_GT","si3_AD","si3_DP","si3_GQ","si3_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(si4, c("si4_GT","si4_AD","si4_DP","si4_GQ","si4_PL"),sep = ":",extra = "drop", fill = "right") %>%
  separate(si5, c("si5_GT","si5_AD","si5_DP","si5_GQ","si5_PL"),sep = ":",extra = "drop", fill = "right")

#complete calls for all genotypes
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split <- combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split %>% 
  filter(an1_GT != "./.") %>%
  filter(an2_GT != "./.") %>%
  filter(an3_GT != "./.") %>%
  filter(an4_GT != "./.") %>%
  filter(an5_GT != "./.") %>%
  filter(an6_GT != "./.") %>%
  filter(br1_GT != "./.") %>%
  filter(br2_GT != "./.") %>%
  filter(br3_GT != "./.") %>%
  filter(br4_GT != "./.") %>%
  filter(br5_GT != "./.") %>%
  filter(ge1_GT != "./.") %>%
  filter(ge2_GT != "./.") %>%
  filter(ge3_GT != "./.") %>%
  filter(ge4_GT != "./.") %>%
  filter(ge5_GT != "./.") %>%
  filter(ge6_GT != "./.") %>%
  filter(he1_GT != "./.") %>%
  filter(he2_GT != "./.") %>%
  filter(he3_GT != "./.") %>%
  filter(he4_GT != "./.") %>%
  filter(he5_GT != "./.") %>%
  filter(he6_GT != "./.") %>%
  filter(re2_GT != "./.") %>%
  filter(re3_GT != "./.") %>%
  filter(re4_GT != "./.") %>%
  filter(re5_GT != "./.") %>%
  filter(re6_GT != "./.") %>%
  filter(sh1_GT != "./.") %>%
  filter(sh2_GT != "./.") %>%
  filter(sh3_GT != "./.") %>%
  filter(sh4_GT != "./.") %>%
  filter(sh5_GT != "./.") %>%
  filter(si1_GT != "./.") %>%
  filter(si2_GT != "./.") %>%
  filter(si3_GT != "./.") %>%
  filter(si4_GT != "./.") %>%
  filter(si5_GT != "./.")

#filter for at least 5 reads mapped
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split <- combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split %>% 
  filter(an1_DP > 4) %>%
  filter(an2_DP > 4) %>%
  filter(an3_DP > 4) %>%
  filter(an4_DP > 4) %>%
  filter(an5_DP > 4) %>%
  filter(an6_DP > 4) 
  # filter(br1_DP > 4) %>%
  # filter(br2_DP > 4) %>%
  # filter(br3_DP > 4) %>%
  # filter(br4_DP > 4) %>%
  # filter(br5_DP > 4)

#split GT column into 2 (e.g. "an1_GT_A","an1_GT_B")
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT <- combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split %>%
  separate(an1_GT, c("an1_GT_A","an1_GT_B"),sep = "/",convert = TRUE) %>%
  separate(an2_GT, c("an2_GT_A","an2_GT_B"),sep = "/",convert = TRUE) %>%
  separate(an3_GT, c("an3_GT_A","an3_GT_B"),sep = "/",convert = TRUE) %>%
  separate(an4_GT, c("an4_GT_A","an4_GT_B"),sep = "/",convert = TRUE) %>%
  separate(an5_GT, c("an5_GT_A","an5_GT_B"),sep = "/",convert = TRUE) %>%
  separate(an6_GT, c("an6_GT_A","an6_GT_B"),sep = "/",convert = TRUE) %>%
  separate(br1_GT, c("br1_GT_A","br1_GT_B"),sep = "/",convert = TRUE) %>%
  separate(br2_GT, c("br2_GT_A","br2_GT_B"),sep = "/",convert = TRUE) %>%
  separate(br3_GT, c("br3_GT_A","br3_GT_B"),sep = "/",convert = TRUE) %>%
  separate(br4_GT, c("br4_GT_A","br4_GT_B"),sep = "/",convert = TRUE) %>%
  separate(br5_GT, c("br5_GT_A","br5_GT_B"),sep = "/",convert = TRUE) %>%
  separate(ge1_GT, c("ge1_GT_A","ge1_GT_B"),sep = "/",convert = TRUE) %>%
  separate(ge2_GT, c("ge2_GT_A","ge2_GT_B"),sep = "/",convert = TRUE) %>%
  separate(ge3_GT, c("ge3_GT_A","ge3_GT_B"),sep = "/",convert = TRUE) %>%
  separate(ge4_GT, c("ge4_GT_A","ge4_GT_B"),sep = "/",convert = TRUE) %>%
  separate(ge5_GT, c("ge5_GT_A","ge5_GT_B"),sep = "/",convert = TRUE) %>%
  separate(ge6_GT, c("ge6_GT_A","ge6_GT_B"),sep = "/",convert = TRUE) %>%
  separate(he1_GT, c("he1_GT_A","he1_GT_B"),sep = "/",convert = TRUE) %>%
  separate(he2_GT, c("he2_GT_A","he2_GT_B"),sep = "/",convert = TRUE) %>%
  separate(he3_GT, c("he3_GT_A","he3_GT_B"),sep = "/",convert = TRUE) %>%
  separate(he4_GT, c("he4_GT_A","he4_GT_B"),sep = "/",convert = TRUE) %>%
  separate(he5_GT, c("he5_GT_A","he5_GT_B"),sep = "/",convert = TRUE) %>%
  separate(he6_GT, c("he6_GT_A","he6_GT_B"),sep = "/",convert = TRUE) %>%
  separate(re2_GT, c("re2_GT_A","re2_GT_B"),sep = "/",convert = TRUE) %>%
  separate(re3_GT, c("re3_GT_A","re3_GT_B"),sep = "/",convert = TRUE) %>%
  separate(re4_GT, c("re4_GT_A","re4_GT_B"),sep = "/",convert = TRUE) %>%
  separate(re5_GT, c("re5_GT_A","re5_GT_B"),sep = "/",convert = TRUE) %>%
  separate(re6_GT, c("re6_GT_A","re6_GT_B"),sep = "/",convert = TRUE) %>%
  separate(sh1_GT, c("sh1_GT_A","sh1_GT_B"),sep = "/",convert = TRUE) %>%
  separate(sh2_GT, c("sh2_GT_A","sh2_GT_B"),sep = "/",convert = TRUE) %>%
  separate(sh3_GT, c("sh3_GT_A","sh3_GT_B"),sep = "/",convert = TRUE) %>%
  separate(sh4_GT, c("sh4_GT_A","sh4_GT_B"),sep = "/",convert = TRUE) %>%
  separate(sh5_GT, c("sh5_GT_A","sh5_GT_B"),sep = "/",convert = TRUE) %>%
  separate(si1_GT, c("si1_GT_A","si1_GT_B"),sep = "/",convert = TRUE) %>%
  separate(si2_GT, c("si2_GT_A","si2_GT_B"),sep = "/",convert = TRUE) %>%
  separate(si3_GT, c("si3_GT_A","si3_GT_B"),sep = "/",convert = TRUE) %>%
  separate(si4_GT, c("si4_GT_A","si4_GT_B"),sep = "/",convert = TRUE) %>%
  separate(si5_GT, c("si5_GT_A","si5_GT_B"),sep = "/",convert = TRUE)

#convert non-zero to 1 for allele column
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$an1_GT_A[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$an1_GT_A != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$an1_GT_B[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$an1_GT_B != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$an2_GT_A[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$an2_GT_A != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$an2_GT_B[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$an2_GT_B != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$an3_GT_A[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$an3_GT_A != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$an3_GT_B[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$an3_GT_B != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$an4_GT_A[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$an4_GT_A != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$an4_GT_B[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$an4_GT_B != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$an5_GT_A[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$an5_GT_A != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$an5_GT_B[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$an5_GT_B != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$an6_GT_A[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$an6_GT_A != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$an6_GT_B[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$an6_GT_B != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$br1_GT_A[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$br1_GT_A != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$br1_GT_B[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$br1_GT_B != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$br2_GT_A[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$br2_GT_A != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$br2_GT_B[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$br2_GT_B != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$br3_GT_A[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$br3_GT_A != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$br3_GT_B[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$br3_GT_B != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$br4_GT_A[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$br4_GT_A != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$br4_GT_B[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$br4_GT_B != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$br5_GT_A[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$br5_GT_A != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$br5_GT_B[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$br5_GT_B != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$ge1_GT_A[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$ge1_GT_A != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$ge1_GT_B[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$ge1_GT_B != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$ge2_GT_A[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$ge2_GT_A != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$ge2_GT_B[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$ge2_GT_B != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$ge3_GT_A[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$ge3_GT_A != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$ge3_GT_B[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$ge3_GT_B != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$ge4_GT_A[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$ge4_GT_A != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$ge4_GT_B[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$ge4_GT_B != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$ge5_GT_A[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$ge5_GT_A != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$ge5_GT_B[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$ge5_GT_B != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$ge6_GT_A[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$ge6_GT_A != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$ge6_GT_B[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$ge6_GT_B != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$he1_GT_A[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$he1_GT_A != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$he1_GT_B[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$he1_GT_B != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$he2_GT_A[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$he2_GT_A != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$he2_GT_B[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$he2_GT_B != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$he3_GT_A[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$he3_GT_A != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$he3_GT_B[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$he3_GT_B != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$he4_GT_A[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$he4_GT_A != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$he4_GT_B[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$he4_GT_B != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$he5_GT_A[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$he5_GT_A != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$he5_GT_B[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$he5_GT_B != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$he6_GT_A[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$he6_GT_A != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$he6_GT_B[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$he6_GT_B != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$re2_GT_A[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$re2_GT_A != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$re2_GT_B[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$re2_GT_B != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$re3_GT_A[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$re3_GT_A != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$re3_GT_B[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$re3_GT_B != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$re4_GT_A[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$re4_GT_A != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$re4_GT_B[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$re4_GT_B != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$re5_GT_A[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$re5_GT_A != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$re5_GT_B[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$re5_GT_B != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$re6_GT_A[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$re6_GT_A != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$re6_GT_B[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$re6_GT_B != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$sh1_GT_A[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$sh1_GT_A != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$sh1_GT_B[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$sh1_GT_B != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$sh2_GT_A[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$sh2_GT_A != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$sh2_GT_B[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$sh2_GT_B != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$sh3_GT_A[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$sh3_GT_A != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$sh3_GT_B[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$sh3_GT_B != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$sh4_GT_A[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$sh4_GT_A != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$sh4_GT_B[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$sh4_GT_B != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$sh5_GT_A[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$sh5_GT_A != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$sh5_GT_B[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$sh5_GT_B != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$si1_GT_A[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$si1_GT_A != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$si1_GT_B[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$si1_GT_B != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$si2_GT_A[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$si2_GT_A != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$si2_GT_B[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$si2_GT_B != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$si3_GT_A[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$si3_GT_A != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$si3_GT_B[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$si3_GT_B != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$si4_GT_A[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$si4_GT_A != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$si4_GT_B[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$si4_GT_B != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$si5_GT_A[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$si5_GT_A != 0] <- 1
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$si5_GT_B[combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT$si5_GT_B != 0] <- 1

#create new var for scoring number of ALT allele
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT <- combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT %>%
  mutate(an1_GT_ALT_count = an1_GT_A + an1_GT_B) %>%
  mutate(an2_GT_ALT_count = an2_GT_A + an2_GT_B) %>%
  mutate(an3_GT_ALT_count = an3_GT_A + an3_GT_B) %>%
  mutate(an4_GT_ALT_count = an4_GT_A + an4_GT_B) %>%
  mutate(an5_GT_ALT_count = an5_GT_A + an5_GT_B) %>%
  mutate(an6_GT_ALT_count = an6_GT_A + an6_GT_B) %>%
  mutate(br1_GT_ALT_count = br1_GT_A + br1_GT_B) %>%
  mutate(br2_GT_ALT_count = br2_GT_A + br2_GT_B) %>%
  mutate(br3_GT_ALT_count = br3_GT_A + br3_GT_B) %>%
  mutate(br4_GT_ALT_count = br4_GT_A + br4_GT_B) %>%
  mutate(br5_GT_ALT_count = br5_GT_A + br5_GT_B) %>%
  mutate(ge1_GT_ALT_count = ge1_GT_A + ge1_GT_B) %>%
  mutate(ge2_GT_ALT_count = ge2_GT_A + ge2_GT_B) %>%
  mutate(ge3_GT_ALT_count = ge3_GT_A + ge3_GT_B) %>%
  mutate(ge4_GT_ALT_count = ge4_GT_A + ge4_GT_B) %>%
  mutate(ge5_GT_ALT_count = ge5_GT_A + ge5_GT_B) %>%
  mutate(ge6_GT_ALT_count = ge6_GT_A + ge6_GT_B) %>%
  mutate(he1_GT_ALT_count = he1_GT_A + he1_GT_B) %>%
  mutate(he2_GT_ALT_count = he2_GT_A + he2_GT_B) %>%
  mutate(he3_GT_ALT_count = he3_GT_A + he3_GT_B) %>%
  mutate(he4_GT_ALT_count = he4_GT_A + he4_GT_B) %>%
  mutate(he5_GT_ALT_count = he5_GT_A + he5_GT_B) %>%
  mutate(he6_GT_ALT_count = he6_GT_A + he6_GT_B) %>%
  mutate(re2_GT_ALT_count = re2_GT_A + re2_GT_B) %>%
  mutate(re3_GT_ALT_count = re3_GT_A + re3_GT_B) %>%
  mutate(re4_GT_ALT_count = re4_GT_A + re4_GT_B) %>%
  mutate(re5_GT_ALT_count = re5_GT_A + re5_GT_B) %>%
  mutate(re6_GT_ALT_count = re6_GT_A + re6_GT_B) %>%
  mutate(sh1_GT_ALT_count = sh1_GT_A + sh1_GT_B) %>%
  mutate(sh2_GT_ALT_count = sh2_GT_A + sh2_GT_B) %>%
  mutate(sh3_GT_ALT_count = sh3_GT_A + sh3_GT_B) %>%
  mutate(sh4_GT_ALT_count = sh4_GT_A + sh4_GT_B) %>%
  mutate(sh5_GT_ALT_count = sh5_GT_A + sh5_GT_B) %>%
  mutate(si1_GT_ALT_count = si1_GT_A + si1_GT_B) %>%
  mutate(si2_GT_ALT_count = si2_GT_A + si2_GT_B) %>%
  mutate(si3_GT_ALT_count = si3_GT_A + si3_GT_B) %>%
  mutate(si4_GT_ALT_count = si4_GT_A + si4_GT_B) %>%
  mutate(si5_GT_ALT_count = si5_GT_A + si5_GT_B)

#create breed proportion variable
combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro <- combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT %>%
  rowwise() %>%
  mutate(an_ALT_proportion = sum(an1_GT_ALT_count,an2_GT_ALT_count,an3_GT_ALT_count,an4_GT_ALT_count,an5_GT_ALT_count,an6_GT_ALT_count)/12,
         br_ALT_proportion = sum(br1_GT_ALT_count,br2_GT_ALT_count,br3_GT_ALT_count,br4_GT_ALT_count,br5_GT_ALT_count)/10,
         ge_ALT_proportion = sum(ge1_GT_ALT_count,ge2_GT_ALT_count,ge3_GT_ALT_count,ge4_GT_ALT_count,ge5_GT_ALT_count,ge6_GT_ALT_count)/12,
         he_ALT_proportion = sum(he1_GT_ALT_count,he2_GT_ALT_count,he3_GT_ALT_count,he4_GT_ALT_count,he5_GT_ALT_count,he6_GT_ALT_count)/12,
         re_ALT_proportion = sum(re2_GT_ALT_count,re3_GT_ALT_count,re4_GT_ALT_count,re5_GT_ALT_count,re6_GT_ALT_count)/10,
         sh_ALT_proportion = sum(sh1_GT_ALT_count,sh2_GT_ALT_count,sh3_GT_ALT_count,sh4_GT_ALT_count,sh5_GT_ALT_count)/10,
         si_ALT_proportion = sum(si1_GT_ALT_count,si2_GT_ALT_count,si3_GT_ALT_count,si4_GT_ALT_count,si5_GT_ALT_count)/10)

#create proportion difference var, angus minus brahman
# combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro <- combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro %>%
#   mutate(proportion_diff = an_ALT_proportion - br_ALT_proportion)

#for converting genotypes in selection interval across the genome
assign(paste('combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split', chr, sep=''), 
       combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split)

path2 <- paste0(dir1,'combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split',chr,".RData")

save(list = paste0('combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split', chr),file = 
       path2)

#test different window sizes
window_size <- seq(50e3,150e3,by = 50e3)
for (i in 1:length(window_size)){
  
  #create interval to summarise proportion
  combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro$interval <- 
    with(combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro, cut(POS, seq(1, max(POS), by = window_size[i])))
  
  #use proportion per breed instead of proportion diff
  proportion_diff_df <- combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro %>% 
    group_by(interval) %>% summarise(an_proportion_mean = mean(an_ALT_proportion),
                                     br_proportion_mean= mean(br_ALT_proportion),
                                     ge_proportion_mean= mean(ge_ALT_proportion),
                                     he_proportion_mean= mean(he_ALT_proportion),
                                     re_proportion_mean= mean(re_ALT_proportion),
                                     sh_proportion_mean= mean(sh_ALT_proportion),
                                     si_proportion_mean= mean(si_ALT_proportion),
                                     snp_counted = n())
  
  proportion_diff_df$interval_no <- 1:nrow(proportion_diff_df)
  
  #any consecutive low Brahman ALT allele ref?
  # proportion_diff_df_summ <- proportion_diff_df %>% filter(br_proportion_mean <= 0.1) %>%
  #   filter(an_proportion_mean > 0.4) %>% filter(snp_counted >= 10)
  # 
  # proportion_diff_df_summ$chr <- chr
  
  #write out tables per chr
  #final filtered vcf file with new columns from my split
  assign(paste('combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_', i,"_",chr, sep=''), 
         combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro)
  
  path2 <- paste0(dir1,'combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_',i,"_",chr,".RData")
  
  save(list = paste0('combinedBreedRefBrahman_filteredSNP.Brahman_multianno_subset_1_split_GT_pro_', i,"_",chr),file = 
         path2)
}

