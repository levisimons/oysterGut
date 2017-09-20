library(ggplot2)
library(reshape2)
library(phyloseq)
library(tidyr)
library(plyr)

setwd("~/Desktop/OysterMicrobiome")

## Read in OTU count data by sample file and the taxonomy file from MOTHUR output.
OTUCount <- import_mothur(mothur_shared_file = "stability.opti_mcc.shared", mothur_constaxonomy_file = "stability.cons.taxonomy")
## Scale OTU by relative sequence abundance.
OTUCount = transform_sample_counts(OTUCount,function(x) x / sum(x))
## Filter out OTUs with a relative abundance below a given cutoff.
OTUCount = filter_taxa(OTUCount, function(x) mean(x) > 1e-2, TRUE)
## Read in experimental design variables, such as the week the sample was taken.
factors <- read.table("MicrobiomeFactors.tsv", header=TRUE, sep="\t",as.is=T)
## Convert experimental design variables into a dataframe.
factors <- sample_data(data.frame(factors, row.names=sample_names(OTUCount)))
## Merge relative OTU abundance data with design variables into a Phyloseq object.
microbiome <- merge_phyloseq(OTUCount,factors)
## Update the taxonomic order names.
colnames(tax_table(microbiome)) <- c("kingdom", "phylum", "class", "order", "family",  "genus")
## Convert the Phyloseq object back to a dataframe.
microbiomeDF <- psmelt(microbiome)
## Create a subset dataframe with the OTU abundance data, and the genus name
## merged into a single row variable.  The resulting matrix will then have this
## row variable by week of sampling as the column variables.  The entries
## will be the relative OTU abundance.
microbiomeAbundance <- spread(microbiomeDF,WeekFromStart,Abundance,fill=NA)
microbiomeAbundance$OTUInfo <- paste(microbiomeAbundance$OTU,microbiomeAbundance$genus,microbiomeAbundance$SampleID,sep="_")
microbiomeAbundance <- microbiomeAbundance[,c(29,1:28)]
## Select and output experimental variables for each OTU node for later
## Cytoscape analysis.
LSANodes <-subset(microbiomeAbundance,select=c(1,5:8,11:17))
write.table(LSANodes,file="eLSANodes.csv",quote=FALSE,sep=",",row.names=FALSE)
## Finish aggregating data so that the first column contains:
## OTU number_genus_SampleID
## and the remaining columns are the relative OTU abundance by week.
microbiomeAbundance <- microbiomeAbundance[,-c(2:17)]
colnames(microbiomeAbundance) <- c("OTUInfo","Week1","Week2","Week3","Week4","Week5","Week6","Week7","Week8","Week9","Week10","Week11","Week12")
microbiomeAbundance <- aggregate(microbiomeAbundance,by=list(microbiomeAbundance$OTUInfo),FUN=function(x) na.omit(x)[1])[,-1]
LSAInput <- microbiomeAbundance
colnames(LSAInput) <- c("#","Week1","Week2","Week3","Week4","Week5","Week6","Week7","Week8","Week9","Week10","Week11","Week12")
write.table(LSAInput,file="LSAInput.tsv",quote=FALSE,sep="\t",row.names=FALSE)
