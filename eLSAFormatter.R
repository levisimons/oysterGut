library(ggplot2)
library(reshape2)
library(phyloseq)

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
eLSAMicrobiome <- subset(microbiomeDF, select=c("OTU","genus","SampleID"))
eLSAMicrobiome$OTUAndTaxa <- as.character((interaction(eLSAMicrobiome,sep="-")))
OTUAndTaxa <- subset(eLSAMicrobiome, select=c("OTUAndTaxa"))
AbundanceAndTime <- subset(microbiomeDF,select=c("Abundance","WeekFromStart"))
eLSAMicrobiome <- cbind(OTUAndTaxa,AbundanceAndTime)
eLSAMicrobiomeDF <- reshape(eLSAMicrobiome,idvar="OTUAndTaxa",timevar="WeekFromStart",direction="wide")
names(eLSAMicrobiomeDF) <- gsub(x = names(eLSAMicrobiomeDF), pattern = "Abundance.", replacement = "Week ")
## Re-order the week columns so they are in chronological order.
eLSAMicrobiomeDF <- eLSAMicrobiomeDF[c(1,7,10,13,8,4,2,6,11,9,3,5,12)]
## Output data as a .tsv file for extended Local Similarity Analysis (eLSA)
write.table(eLSAMicrobiomeDF,file="eLSAMicrobiome.tsv",quote=FALSE,sep="\t")
