library(ggplot2)
library(reshape2)
library(phyloseq)

setwd("~/Desktop/OysterMicrobiome")

OTUCount <- import_mothur(mothur_shared_file = "stability.opti_mcc.shared", mothur_constaxonomy_file = "stability.cons.taxonomy")
OTUCount = transform_sample_counts(OTUCount,function(x) x / sum(x))
OTUCount = filter_taxa(OTUCount, function(x) mean(x) > 1e-2, TRUE)

factors <- read.table("MicrobiomeFactors.tsv", header=TRUE, sep="\t",as.is=T)

factors <- sample_data(data.frame(factors, row.names=sample_names(OTUCount)))

microbiome <- merge_phyloseq(OTUCount,factors)
colnames(tax_table(microbiome)) <- c("kingdom", "phylum", "class", "order", "family",  "genus")
microbiomeDF <- psmelt(microbiome)
eLSAMicrobiome <- subset(microbiomeDF, select=c("OTU","genus","SampleID"))
eLSAMicrobiome$OTUAndTaxa <- as.character((interaction(eLSAMicrobiome,sep="-")))
OTUAndTaxa <- subset(eLSAMicrobiome, select=c("OTUAndTaxa"))
AbundanceAndTime <- subset(microbiomeDF,select=c("Abundance","WeekFromStart"))
eLSAMicrobiome <- cbind(OTUAndTaxa,AbundanceAndTime)
reshape(eLSAMicrobiome,idvar="OTUAndTaxa",timevar="WeekFromStart",direction="wide")
