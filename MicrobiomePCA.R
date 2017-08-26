library(ggplot2)
library(ggfortify)

setwd("~/Desktop/OysterGut/Analysis")

## Read in relative abundance per OTUS data for oyster microbiome study samples
RabundData <- read.table("oystersRabund", header = TRUE, sep="\t", as.is=T)
## Ignore metadata columns.
## For the fecal pellet and stomach samples this is columns 1 through 8.
## For the feedstock samples this is columns 1 through 5.
RabundOTUs <- RabundData[ -c(1:8)]
## Run PCA on relative abundance data.
RabundOTUs.pca <- prcomp(RabundOTUs)
## Calculate the percent variance attributed to each principal component.
RabundOTUsSTDEV <- RabundOTUs.pca$sdev
var <- RabundOTUsSTDEV^2
var.percent <- var/sum(var) * 100
dev.new()
## Generate histogram of the percent variance attributed to each principal component.
barplot(var.percent, xlab="Principal Component", ylab="Percent Variance", names.arg=1:length(var.percent), las=1, ylim=c(0,max(var.percent)), col="gray")
## Plot the two dominant principal components from the data.
autoplot(prcomp(RabundOTUs), data=RabundData, 
         xlab=paste('PCA1, percent variation: ',round(var.percent[1], digits=3)), 
         ylab=paste('PCA2, percent variation: ',round(var.percent[2], digits=3)),
         main='PCA2 vs. PCA1 Oyster stomach samples and fecal pellet communities',
         shape=TRUE, label=FALSE, loadings.label.size = 3, colour='FeedType')
