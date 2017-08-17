library(glm2)
library(lme4)
library(plyr)
library(car)
library(bestglm)
library(codyn)
library(corrgram)
library(igraph)
library(corrplot)
##############
## Read in data
## There is a separate feedstock and oyster file.
## The first columns contain important experimental information.
setwd("~/Desktop/OysterGut/Analysis")

## Read in an OTU abundance by sample data .shared file.
## The first three columns are removed as they contain labels unnecessary for later analysis.
Top30 <- read.table("stability.feed.CHAE.label.Top30.shared.tsv", header = TRUE, sep="\t", as.is=T)
Top30 <- Top30[ -c(1:3)]
M <- cor(Top30)

# Tally the total OTU count across all samples being studied.
# This is to resize the nodes in a network graph according to their OTU abundance.
TotOTUs <- colSums(Top30)

# Generate correlation matrix from OTU abundance data.
# Nodes are sized by the OTU abundance across samples.
# Links are sized by the strength of the Pearson correlation coefficient across samples.
mat=cor((Top30))
mat[abs(mat)<0.9]=0
network=graph_from_adjacency_matrix(mat, weighted=T, mode="undirected", diag=F)
fine = 50 # this will adjust the resolving power.
palette = colorRampPalette(c('red','blue'))
#this gives you the colors you want for every point
graphCol = palette(fine)[as.numeric(cut(mat,breaks = fine))]

# Plot network map in a circle with the largest OTUs starting at a 3 o'clock position
# spiraling counter-clockwise in decreasing total OTU abundance.
# Label nodes with OTU names.
# Show links weighted by correlation coefficients between OTU abundances across samples.
# Color links blue for positive correlations between OTU abundances across samples, and red for negative correlations.
plot(network, 
     main="Chaetoceros feedstock communities weeks 9-12",
     sub="abs(r)>0.9 p<0.05 \n nodes sized by total OTU abundance (97% similarity cutoff) \n link weight sized by correlation(blue = positive, red = negative)",
     vertex.label.cex = 0.7, 
     layout=layout_in_circle(network), 
     vertex.size=10^-4*(TotOTUs), 
     edge.width=abs(10*E(network)$weight), 
     edge.color=ifelse(E(network)$weight > 0, "blue","red"), 
     edge.label.cex = 0.7)

## This generates a correlation matrix of the OTU abundance data across
## a designated range of samples.  Correlations are plotted using a piechart
## representation and are only displayed if the correlation's
## significance is below the 5% threshold.
cor.mtest <- function(mat, conf.level = 0.95){
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      tmp <- cor.test(mat[,i], mat[,j], conf.level = conf.level)
      p.mat[i,j] <- p.mat[j,i] <- tmp$p.value
      lowCI.mat[i,j] <- lowCI.mat[j,i] <- tmp$conf.int[1]
      uppCI.mat[i,j] <- uppCI.mat[j,i] <- tmp$conf.int[2]
    }
  }
  return(list(p.mat, lowCI.mat, uppCI.mat))
}
res1 <- cor.mtest(Top30,0.95)
res2 <- cor.mtest(Top30,0.99)
## specialized the insignificant value according to the significant level
## plot matrix
corrplot(M, p.mat = res1[[1]], sig.level=0.05, insig = "blank", tl.cex=0.6, method="pie")
