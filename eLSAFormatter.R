##############
## Read in data
## There is a separate feedstock and oyster file.
## The first columns contain important experimental information.
setwd("~/Desktop/OysterGut/OysterMicrobiome")

## Read in data representing relative abundance values for each OTU by oyster sample.
## The absolute sequence counts by OTU are normalized by total sequences per sample.
oyster <- read.table("oyster", header = TRUE, sep="\t", as.is=T)
oysterAbundances <- read.table("oysterAbundances", header = TRUE, sep="\t", as.is=T)

## Merge the experimental design factors with the relative OTU abundance data.
oysters <- merge(oyster,oysterAbundances,all.Group=TRUE)
oysters <- oysters[,c(1,4,9:48)]

## Read in data representing relative abundance values for each OTU by feedstock sample.
## The absolute sequence counts by OTU are normalized by total sequences per sample.
feedstock <- read.table("feedstock", header = TRUE, sep="\t", as.is=T)
feedAbundances <- read.table("feedAbundances", header = TRUE, sep="\t", as.is=T)

## Merge the experimental design factors with the relative OTU abundance data.
feed <- merge(feedstock,feedAbundances,all.Group=TRUE)

microbiome <- rbind(feedAbundances,oysterAbundances,by="Group")
