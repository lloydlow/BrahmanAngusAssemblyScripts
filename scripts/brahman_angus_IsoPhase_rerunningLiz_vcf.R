#------------------------------------------------------
# Program name: brahman_angus_IsoPhase_rerunningLiz_vcf.R
# Objective: reanalyse Liz Venn diagram and her RNA editing  
#         call
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(ggplot2)
library(dplyr)
library(VennDiagram)

# path to isophase results
dir1 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/isoseq/IsoPhase/April2019_IsoPhase_with_Lloyd/"

# reading data files
path1 <- paste0(dir1,"evaled.isophase.combined_97_1modi.txt")
path2 <- paste0(dir1,"evaled.isophase.F1RefBrahman_filteredSNP.txt")

# “ref”: Brahman nucleotide at the position
# “alt_Short”:  the alternative base called by the short read (RNA-seq or genomic) data
# “alt_PB”:   this is always the alternative base called by IsoPhase (Iso-Seq, PacBio). If PacBio did not call a SNP here it is “NA”
# “in_Short”:  a “Y” if alt_Short!=’NA’
# “in_PB”:   a “Y” if alt_PB!=’NA’
# “cov_Short”:  the number of short read (RNA-seq or genomic) coverage at this position, based on the VCF file
# “cov_PB”:  the number of FLNC reads covering this position
# “genomic_HP”:  this field is currently un-used. I used it for the maize genome where I marked locations of potential homopolymers, which is very common in maize.

x1 <- read.table(path1,sep='\t',header=T)
x2 <- read.table(path2,sep='\t',header=T)

x1[1:5,]

x1$coord <- paste(x1$chrom, x1$pos, sep=':')
x2$coord <- paste(x2$chrom, x2$pos, sep=':')

x3 <- merge(x1, x2, by="coord", all=TRUE)
x3 <- x3[,c("coord", "dir.x", "ref.x", "alt_Short.x", "alt_Short.y", "alt_PB.x", "cov_PB.x")]
x3 <- distinct(x3, coord, .keep_all=T)

# let's draw a venn diagram
# area1:PB, area2:RNA-seq, area3:genomic
area1 <- sum(!is.na(x3$alt_PB.x))
area2 <- sum(!is.na(x3$alt_Short.x)) # RNA-seq
area3 <- sum(!is.na(x3$alt_Short.y)) # genomic

n12 <- sum(!is.na(x3$alt_PB.x)&!is.na(x3$alt_Short.x))
n23 <- sum(!is.na(x3$alt_Short.x)&!is.na(x3$alt_Short.y))
n13 <- sum(!is.na(x3$alt_PB.x)&!is.na(x3$alt_Short.y))
n123 <- sum(!is.na(x3$alt_PB.x)&!is.na(x3$alt_Short.x)&!is.na(x3$alt_Short.y))

draw.triple.venn(area1, area2, area3, n12, n23, n13, n123, category=c("Iso-Seq", "RNA-seq", "Genomic"), fill=c("lightpink", 'lightblue', 'lightgreen'))

tiff(filename = "FigFinal_IsoPhase_SNP_validation.tiff",width = 300,height = 250)
draw.triple.venn(area1, area2, area3, n12, n23, n13, n123, category=c("Iso-Seq", "RNA-seq", "Genomic"), fill=c("lightpink", 'lightblue', 'lightgreen'))
dev.off()
