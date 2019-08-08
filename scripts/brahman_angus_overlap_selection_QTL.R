#------------------------------------------------------
# Program name: brahman_angus_overlap_selection_QTL.R
# Objective: analyse the overlap of selective sweep interval 
#   with QTL
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(ggplot2)
library(readr)

#learning pie chart: https://www.displayr.com/how-to-make-a-pie-chart-in-r/ 
dir1 <- "/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/Annotation/EBI_Brahman_SevenBreeds/"

filename <- "overlap_selection_QTL_as_df.csv"

overlap_selection_QTL_as_df <- read_csv(paste0(dir1,filename))

#unique QTL ID only because selection interval may hit same qtl multiple times
overlap_selection_QTL_as_df <- overlap_selection_QTL_as_df %>% distinct(QTL_ID, .keep_all = TRUE)

overlap_selection_QTL_as_df_QTL_type <- overlap_selection_QTL_as_df %>% group_by(QTL_GR.type) %>% summarise(count = n()) %>% 
  mutate(proportion = count/sum(count))

overlap_selection_QTL_as_df_QTL_type$QTL_GR.type <- gsub("_"," ",overlap_selection_QTL_as_df_QTL_type$QTL_GR.type)

names(overlap_selection_QTL_as_df_QTL_type) <- c("QTL_type", "count", "proportion")

overlap_selection_QTL_as_df_QTL_type

#pie chart by QTL type
# Create a basic bar
pie <- ggplot(overlap_selection_QTL_as_df_QTL_type, aes(x="", y=proportion, fill=QTL_type)) + geom_bar(stat="identity", width=1)

# Convert to pie (polar coordinates) and add labels
pie <- pie + coord_polar("y", start=0) + geom_text(aes(label = paste0(round(proportion*100,1), "%")), 
                                                   position = position_stack(vjust = 0.5), size =3.5)

# Add color scale (hex colors)
#pie = pie + scale_fill_manual(values=c("#55DDE0", "#33658A", "#2F4858", "#F6AE2D", "#F26419", "#999999")) 

# Remove labels and add title
pie <- pie + labs(x = NULL, y = NULL, fill = "QTL type", title = "Overlap of Selective Sweep Intervals with QTL")

# Tidy up the theme
pie <- pie + theme_classic() + theme(axis.line = element_blank(),
                                    axis.text = element_blank(),
                                    axis.ticks = element_blank(),
                                    plot.title = element_text(hjust = 0.5))
pie

#merging QTL and association
overlap_selection_QTL_as_df$QTL_GR.type <- gsub("_"," ",overlap_selection_QTL_as_df$QTL_GR.type)
overlap_selection_QTL_as_df$QTL_GR.type <- gsub(" Association","",overlap_selection_QTL_as_df$QTL_GR.type)
overlap_selection_QTL_as_df$QTL_GR.type <- gsub(" QTL","",overlap_selection_QTL_as_df$QTL_GR.type)

overlap_selection_QTL_as_df_QTL_type <- overlap_selection_QTL_as_df %>% group_by(QTL_GR.type) %>% summarise(count = n()) %>% 
  mutate(proportion = count/sum(count))

names(overlap_selection_QTL_as_df_QTL_type) <- c("QTL_type", "count", "proportion")

overlap_selection_QTL_as_df_QTL_type

tiff(filename = "FinalFig_Overlap_Selective_Sweep_QTL.tiff",height = 400, width = 400)
#pie chart by QTL type
# Create a basic bar
pie <- ggplot(overlap_selection_QTL_as_df_QTL_type, aes(x="", y=proportion, fill=QTL_type)) + geom_bar(stat="identity", width=1)

# Convert to pie (polar coordinates) and add labels
pie <- pie + coord_polar("y", start=0) + geom_text(aes(label = paste0(round(proportion*100,1), "%")), 
                                                   position = position_stack(vjust = 0.5), size =3.5)

# Add color scale (hex colors)
#pie = pie + scale_fill_manual(values=c("#55DDE0", "#33658A", "#2F4858", "#F6AE2D", "#F26419", "#999999")) 

# Remove labels and add title
pie <- pie + labs(x = NULL, y = NULL, fill = "QTL type", title = "Overlap of Selective Sweep Intervals with QTL")

# Tidy up the theme
pie <- pie + theme_classic() + theme(axis.line = element_blank(),
                                     axis.text = element_blank(),
                                     axis.ticks = element_blank(),
                                     plot.title = element_text(hjust = 0.5))
pie
dev.off()
