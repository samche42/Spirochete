library(ggplot2)
library(reshape2)

data <- read.csv("chemotype_boxplot.txt", header=TRUE, sep="\t", row.names=1)

plot <- ggplot(data, aes(x=Chemotype, y=OTU6, group=Chemotype)) + geom_boxplot(aes(fill=Chemotype))+ facet_grid(. ~ Source)
