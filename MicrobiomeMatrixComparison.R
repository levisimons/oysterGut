library(glm2)
library(lme4)
library(plyr)
library(car)
library(bestglm)
library(codyn)
library(corrgram)
library(igraph)
library(corrplot)
library(ade4)
library(covmat)
##############
## Read in data
## There is a separate feedstock and oyster file.
## The first columns contain important experimental information.
setwd("~/Desktop/OysterGut/Analysis")

## Read in the first .corr file.
## This file contains four columns.
## Columns 1 and 2 are the OTU labels.
## Column 3 is the correlation between the abundances of those two OTUS are observed across a set of samples.
## Column 4 is the p-value of the correlation.
## If the p-value of a particular correlation is above 0.05 then it is disregarded.
M1 <- read.table("stability.ISOFeed.0.03.Top30.pearson.otu.corr", header = TRUE, sep="\t", as.is=T)
M1 <- ifelse(M1$Significance>0.05,NA,M1$pearsonCoef)

## Read in the second .corr file.
## This file contains four columns.
## Columns 1 and 2 are the OTU labels.
## Column 3 is the correlation between the abundances of those two OTUS are observed across a set of samples.
## Column 4 is the p-value of the correlation.
## If the p-value of a particular correlation is above 0.05 then it is disregarded.
M2 <- read.table("stability.TETFecalControl.0.03.Top30.pearson.otu.corr", header = TRUE, sep="\t", as.is=T)
M2 <- ifelse(M2$Significance>0.05,NA,M2$pearsonCoef)

## Compute the t-statistical test comparing the two correlation vectors
cor.test(M1,M2,na.action=na.omit)
