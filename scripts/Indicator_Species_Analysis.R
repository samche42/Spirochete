#Taken from https://jkzorz.github.io/2019/07/02/Indicator-species-analysis.html

library(indicspecies)
data = read.csv("KEGG_data.txt", header= TRUE, sep= "\t")

abund = data[,4:ncol(data)] #Collect only numerical data, starts in 4th column
host = data$Broad_host

inv = multipatt(abund, host, func = "r.g", control = how(nperm=1000),duleg = TRUE) #Run analysis. The duleg parameter prevents group site combinations being assessed

summary(inv)
