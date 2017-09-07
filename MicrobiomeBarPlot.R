library(ggplot2)
library(reshape2)

setwd("~/Desktop/OysterMicrobiome")

## Read in OTU taxonomy data.  This matches OTUs down to a genus level label.
OTULabels <-read.table("OysterMicrobiomeTaxonomy.csv", header=TRUE, sep=",",as.is=T)

## Convert genus level labels to a single row vector
Genus <- t(OTULabels$Genus)

## Read in data representing relative abundance values for each OTU by oyster sample.
## The absolute sequence counts by OTU are normalized by total sequences per sample.
## The data and design factors here are for the oyster microbiome samples.
#oyster <- read.table("oyster", header = TRUE, sep="\t", as.is=T)
oysterAbundances <- read.table("oysterAbundances", header = TRUE, sep="\t", as.is=T)

## Read in data representing relative abundance values for each OTU by oyster sample.
## The absolute sequence counts by OTU are normalized by total sequences per sample.
## The data and design factors here are for the feedstock microbiome samples.
feed <- read.table("feedstock", header = TRUE, sep="\t", as.is=T)
feedAbundances <- read.table("feedAbundances", header = TRUE, sep="\t", as.is=T)

## Merge the experimental design factors and genu labels with the relative OTU abundance data.
## Choose this merge command set for the feedstock data.
RabundData <- merge(feed,feedAbundances,all.Group=TRUE)
colnames(RabundData)[6:4014] <- Genus
RabundData <- RabundData[,c(1,3,6:7)]

## Merge the experimental design factors and genu labels with the relative OTU abundance data.
## Choose this merge command set for the oyster data.
RabundData <- merge(oyster,oysterAbundances,all.Group=TRUE)
colnames(RabundData)[9:4017] <- Genus
RabundData <- RabundData[,c(1,9:19)]

## Use this command to store community data as a plot object.
RabundPlot <- melt(RabundData, id.vars = "WeekFromStart", variable.name = "Genus")
RabundPlot <- RabundPlot[-c(1:17),]

## Plot the relative abundance of each taxa for a given plot object.
p <- ggplot(RabundPlot, aes(x = WeekFromStart, y = value, fill = Genus)) 
p + geom_bar(stat = "identity") + labs(x="Week from start",y="Relative abundance", title="Relative abundance of taxa by sample")
