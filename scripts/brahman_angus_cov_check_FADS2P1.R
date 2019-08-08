#------------------------------------------------------
# Program name: brahman_angus_cov_check_FADS2P1.R
# Objective: this script is to loop thro cov file from
#   bam file around the FADS2P1 locus, Brahman ref
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

# path to cov file
dir1 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/FADS2P1_other_indi_copy_number/cov/"

#FADS2P1 indicus specific
# x.lower.limit <- 4743181
# x.upper.limit <- 4796329

# x.lower.limit <- 5095415
# x.upper.limit <- 5140465

x.lower.limit <- 3748952
x.upper.limit <- 3789148

all_cov_filenames <- c("Brahman1.merged.addrg.cov",
                       "Brahman2.merged.addrg.cov",
                       "Brahman3.merged.addrg.cov",
                       "Brahman4.merged.addrg.cov",
                       "Brahman5.merged.addrg.cov",
                       "Angus1.merged.addrg.cov",
                       "Angus2.merged.addrg.cov",
                       "Angus3.merged.addrg.cov",
                       "Angus4.merged.addrg.cov",
                       "Angus5.merged.addrg.cov",
                       "Angus6.merged.addrg.cov",
                       "Gelbvieh1.merged.addrg.cov",
                       "Gelbvieh2.merged.addrg.cov",
                       "Gelbvieh3.merged.addrg.cov",
                       "Gelbvieh4.merged.addrg.cov",
                       "Gelbvieh5.merged.addrg.cov",
                       "Gelbvieh6.merged.addrg.cov",
                       "Hereford1.merged.addrg.cov",
                       "Hereford2.merged.addrg.cov",
                       "Hereford3.merged.addrg.cov",
                       "Hereford4.merged.addrg.cov",
                       "Hereford5.merged.addrg.cov",
                       "Hereford6.merged.addrg.cov",
                       "RedAngus2.merged.addrg.cov",
                       "RedAngus3.merged.addrg.cov",
                       "RedAngus4.merged.addrg.cov",
                       "RedAngus5.merged.addrg.cov",
                       "RedAngus6.merged.addrg.cov",
                       "Shorthorn1.merged.addrg.cov",
                       "Shorthorn2.merged.addrg.cov",
                       "Shorthorn3.merged.addrg.cov",
                       "Shorthorn4.merged.addrg.cov",
                       "Shorthorn5.merged.addrg.cov",
                       "Simmental1.merged.addrg.cov",
                       "Simmental2.merged.addrg.cov",
                       "Simmental3.merged.addrg.cov",
                       "Simmental4.merged.addrg.cov",
                       "Simmental5.merged.addrg.cov")

tiff(filename = "FADS2P1_cov_3748952_3789148.tiff", height = 1000, width = 1000)
par(mfrow=c(7,6))
for (i in 1:length(all_cov_filenames)){
  #print(all_cov_filenames[i])

filename <- all_cov_filenames[i]
path1 <- paste0(dir1,filename)

title <- gsub(".merged.addrg.cov","",filename)

#read in cov file
Ill_cov <- read.delim(path1,header = FALSE, stringsAsFactors = FALSE)

colnames(Ill_cov) <- c("Input","POS","cov")

boolean1 <- which(Ill_cov$POS >= x.lower.limit & Ill_cov$POS <= x.upper.limit)

Ill_cov_region <- Ill_cov[boolean1,]

plot(Ill_cov_region$POS,Ill_cov_region$cov,type = "h",xlab = "Position on chr 15",
     ylab = "Coverage",ylim = c(0,35),
     main = title)
#ylim = c(0,range(Ill_cov_region$cov)[2])
abline(h=5,col="blue",lwd=1.5,lty="dashed")

}
dev.off()
par(mfrow=c(1,1))

#ori 38 animals
# all_cov_filenames <- c("Brahman1.merged.addrg.cov",
#                        "Brahman2.merged.addrg.cov",
#                        "Brahman3.merged.addrg.cov",
#                        "Brahman4.merged.addrg.cov",
#                        "Brahman5.merged.addrg.cov",
#                        "Angus1.merged.addrg.cov",
#                        "Angus2.merged.addrg.cov",
#                        "Angus3.merged.addrg.cov",
#                        "Angus4.merged.addrg.cov",
#                        "Angus5.merged.addrg.cov",
#                        "Angus6.merged.addrg.cov",
#                        "Gelbvieh1.merged.addrg.cov",
#                        "Gelbvieh2.merged.addrg.cov",
#                        "Gelbvieh3.merged.addrg.cov",
#                        "Gelbvieh4.merged.addrg.cov",
#                        "Gelbvieh5.merged.addrg.cov",
#                        "Gelbvieh6.merged.addrg.cov",
#                        "Hereford1.merged.addrg.cov",
#                        "Hereford2.merged.addrg.cov",
#                        "Hereford3.merged.addrg.cov",
#                        "Hereford4.merged.addrg.cov",
#                        "Hereford5.merged.addrg.cov",
#                        "Hereford6.merged.addrg.cov",
#                        "RedAngus2.merged.addrg.cov",
#                        "RedAngus3.merged.addrg.cov",
#                        "RedAngus4.merged.addrg.cov",
#                        "RedAngus5.merged.addrg.cov",
#                        "RedAngus6.merged.addrg.cov",
#                        "Shorthorn1.merged.addrg.cov",
#                        "Shorthorn2.merged.addrg.cov",
#                        "Shorthorn3.merged.addrg.cov",
#                        "Shorthorn4.merged.addrg.cov",
#                        "Shorthorn5.merged.addrg.cov",
#                        "Simmental1.merged.addrg.cov",
#                        "Simmental2.merged.addrg.cov",
#                        "Simmental3.merged.addrg.cov",
#                        "Simmental4.merged.addrg.cov",
#                        "Simmental5.merged.addrg.cov")
