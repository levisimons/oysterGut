## This script uses MANOVA to check for significant experimental factors
## in shaping the microbial community composition found in the feedstock,
## stomach, and fecal pellet communities of the oyster microbiome project.
library(glm2)
library(lme4)
library(dplyr)
library(car)
library(bestglm)
library(codyn)
library(corrgram)
library(dummies)
library(heatmap3)
##############
## Read in data
## There is a separate feedstock and oyster file.
## The first columns contain important experimental information.
setwd("~/Desktop/OysterGut/Analysis")

## Read in data representing relative abundance values for each OTU by oyster sample.
## The absolute sequence counts by OTU are normalized by total sequences per sample.
oysters <- read.table("oystersRabund", header = TRUE, sep="\t", as.is=T)

## Ignore the oyster gut community sample rows.
oysters <- oysters[-c(6,22,31,48,64,74,86,93),]

## Treat the relative abundances of the selected OTUs as a variable.
TopOysterOTUs <- as.matrix(oysters[,c(9:18)])
Marimononas <- as.matrix(oysters[,c(72,3267,4023)])

## Test the signifance of experimental variables on the selected OTU data.
OysterManova <- manova(TopOysterOTUs ~ OysterID*FeedType, data=oysters)
summary(OysterManova)

## Read in data representing relative abundance values for each OTU by feedstock sample.
## The absolute sequence counts by OTU are normalized by total sequences per sample.
feed <- read.table("feedRabund", header = TRUE, sep="\t", as.is=T)

## Treat the relative abundances of the selected OTUs as a variable.
TopFeedOTUs <- as.matrix(feed[,c(6:15)])

## Test the signifance of experimental variables on the selected OTU data.
FeedManova <- manova(TopFeedOTUs ~ FeedType, data=feed)
summary(FeedManova)
