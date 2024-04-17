library(ggplot2)
library(reshape2)
library(vegan)
library(grid)
library(ggh4x)

#Generating barplot of all OTUs from all samples

#Read in data
data <- read.csv(file = "barplot_03_input.txt", header=T, sep="\t")

#convert data frame from a "wide" format to a "long" format
melted_data = melt(data, id = c("Sample","Source"))

#Keep samples ordered the way they're ordered in the data
melted_data$Sample <- factor(melted_data$Sample,levels=unique(melted_data$Sample))
melted_data$Source <- factor(melted_data$Source,levels=unique(melted_data$Source))

#Generating plot
ggplot() +
  geom_bar(data = melted_data,aes(x = Sample, y = value, fill = forcats::fct_rev(variable)), stat = "identity", width = 1) +
  theme_bw() +
  facet_nested(. ~ Source, scales = "free", space = "free") +
  theme(strip.text.x = element_text(angle = 90,size=5),
        panel.spacing=unit(0,"lines"),
        panel.background = element_blank(),
        axis.text.x=element_text(angle = 90, hjust = 0.95, vjust = 0.2),
        axis.ticks.x=element_blank())
