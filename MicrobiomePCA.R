library(ggplot2)
library(ggfortify)
library(rgl)

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
         main='PCA2 vs. PCA1 Oyster stomach and fecal pellet communities',
         loadings=FALSE, loadings.label=FALSE, 
         shape=TRUE, label=FALSE, loadings.label.size = 3, colour='FeedType')
plot3d(pc$scores[,1:3], col=iris$Species)
summary(prcomp(RabundOTUs))

head(unclass(RabundOTUs.pca$rotation)[,1:4])
## Take the top four principal components from your data.
RabundOTUsPC <- princomp(RabundOTUs[,1:4])
## Convert feed type from a name to a numerical scale for later plotting.
FeedNum <- as.numeric(as.factor(RabundData$FeedType))
## Plot the three largest principal components.
plot3d(RabundOTUsPC$scores[,1:3], col=FeedNum)
## Include a legend for the scale for feed type.
legend3d("topleft", col=1:3, legend = c("CHAE","ISO","TET"), pch = 20, bty='n', cex=.75)
summary(RabundOTUsPC)
