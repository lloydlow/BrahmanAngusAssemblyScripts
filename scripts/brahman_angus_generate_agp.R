#------------------------------------------------------
# Program name: brahman_angus_generate_agp.R
# Objective: when given "object","component_id",
#           "component_beg","component_end","orientation"
#           the script makes agp file
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(readr)
library(dplyr)

#first step is to gather 5 columns in Excel or vim text file and save it as tsv file
#"object","component_id","component_beg","component_end","orientation"

#then read in scaffold order book in tsv
# scaffold_order_book <- read_tsv("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_indi/dam/scaffold_order_book_Aonly.tsv",
#                                      col_names = FALSE)
# names(scaffold_order_book) <- c("object","object_beg","object_end","part_number","component_type","component_id","component_beg","component_end","orientation")

# scaffold_order_book <- read_tsv("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_template.tsv",
#                                 col_names = FALSE)
#dam
# scaffold_order_book <- read_tsv("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_bionano_brahma.txt",
#                                 col_names = FALSE)
# names(scaffold_order_book) <- c("object","component_id","component_beg","component_end","orientation")
#for final agp after bionano correction
# scaffold_order_book <- read_tsv("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_bionano_brahma_corrected.tsv",
#                                 col_names = FALSE)
# names(scaffold_order_book) <- c("object","component_id","component_beg","component_end","orientation")
#sire
# scaffold_order_book <- read_tsv("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_bionano_angus.txt",
#                                  col_names = FALSE)
# names(scaffold_order_book) <- c("object","component_id","component_beg","component_end","orientation")
#for final agp after bionano correction
scaffold_order_book <- read_tsv("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/scaffold_order_bionano_sire_corrected.tsv",
                                col_names = FALSE)
names(scaffold_order_book) <- c("object","component_id","component_beg","component_end","orientation")


object <- c()
object_beg <- c()
object_end <- c()
part_number <- c()
component_type <- c()
component_id <- c()
component_beg <- c()
component_end <- c()
orientation <- c()

for (i in 1:nrow(scaffold_order_book)){
  if (i == 1 || scaffold_order_book$object[i-1] != scaffold_order_book$object[i]){
    #col 1 agp
    holder <- scaffold_order_book$object[i]
    object <- c(object,holder)
    
    #col 2 agp
    holder <- 1
    object_beg <- c(object_beg,holder)
    
    #col 3 agp
    #holder <- scaffold_order_book$object_beg[i] + scaffold_order_book$component_end[i] - scaffold_order_book$component_beg[i]
    holder <- 1 + scaffold_order_book$component_end[i] - scaffold_order_book$component_beg[i]
    object_end <- c(object_end,holder)
    
    #col 4 agp
    holder <- 1
    part_number <- c(part_number,holder)
    
    #col 5 agp
    holder <- "A"
    component_type <- c(component_type,holder)
    
    #col 6 agp
    holder <- scaffold_order_book$component_id[i]
    component_id <- c(component_id,holder)
    
    #col 7 agp
    holder <- scaffold_order_book$component_beg[i]
    component_beg <- c(component_beg,holder)
    
    #col 8 agp
    holder <- scaffold_order_book$component_end[i]
    component_end <- c(component_end,holder)
    
    #col 9 agp
    holder <- scaffold_order_book$orientation[i]
    orientation <- c(orientation,holder)
    
  } else if (scaffold_order_book$object[i-1] == scaffold_order_book$object[i]){
    #col 1 agp
    holder <- scaffold_order_book$object[i]
    object <- c(object,holder)
    holder <- scaffold_order_book$object[i]
    object <- c(object,holder)
    
    #col 2 & 3 agp
    holder <- object_end[length(object_end)] + 1
    object_beg <- c(object_beg,holder)
    
    holder <- object_end[length(object_end)] + 500
    object_end <- c(object_end,holder)
    
    holder <- object_end[length(object_end)] + 1
    object_beg <- c(object_beg,holder)
    
    holder <- object_end[length(object_end)] + scaffold_order_book$component_end[i] - scaffold_order_book$component_beg[i] + 1
    object_end <- c(object_end,holder)
    
    #col 4 agp
    holder <- part_number[length(part_number)] + 1
    part_number <- c(part_number,holder)
    holder <- part_number[length(part_number)] + 1
    part_number <- c(part_number,holder)
    
    #col 5 agp
    holder <- "U"
    component_type <- c(component_type,holder)
    holder <- "A"
    component_type <- c(component_type,holder)
    
    #col 6 agp
    holder <- "500"
    component_id <- c(component_id,holder)
    holder <- scaffold_order_book$component_id[i]
    component_id <- c(component_id,holder)
    
    #col 7 agp
    holder <- "scaffold"
    component_beg <- c(component_beg,holder)
    holder <- scaffold_order_book$component_beg[i]
    component_beg <- c(component_beg,holder)

    #col 8 agp
    holder <- "yes"
    component_end <- c(component_end,holder)
    holder <- scaffold_order_book$component_end[i]
    component_end <- c(component_end,holder)
    
    #col 9 agp
    holder <- "map"
    orientation <- c(orientation,holder)
    holder <- scaffold_order_book$orientation[i]
    orientation <- c(orientation,holder)
    
  }
}

#to join them all
DF <- data.frame(object=object,object_beg=object_beg,object_end=object_end, part_number=part_number,
                 component_type=component_type,component_id=component_id,component_beg=component_beg,
                 component_end=component_end,orientation=orientation,stringsAsFactors = FALSE)

write_tsv(DF,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/file.agp",col_names = FALSE)
