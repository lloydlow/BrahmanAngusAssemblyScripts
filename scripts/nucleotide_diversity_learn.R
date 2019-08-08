#Learning nucleotide diversity calc, use in conjunction with nucleotide_diversity_step_by_step.pdf
# install.packages("pegas")
library(pegas)

data(woodmouse)
nuc.div(woodmouse)

## a small extract from data(woddmouse) in sequential format:
cat("5 500",
    paste0("No304     TGTCT",paste(rep("A",495),collapse = ""),collapse = ""),
    paste0("No305     TGTCT",paste(rep("A",495),collapse = ""),collapse = ""),
    paste0("No306     TATTA",paste(rep("A",495),collapse = ""),collapse = ""),
    paste0("No307     CGTCT",paste(rep("A",495),collapse = ""),collapse = ""),
    paste0("No308     CGGCT",paste(rep("A",495),collapse = ""),collapse = ""),
    file = "exdna.txt", sep = "\n")

ex.dna <- read.dna("exdna.txt", format = "sequential")
str(ex.dna)
ex.dna

nuc.div(ex.dna)
