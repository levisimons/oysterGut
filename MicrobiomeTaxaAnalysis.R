library("plyr")
library("ggplot2")
library(phyloseq)
library("ape")
library("vegan")
library("devtools")
library("microbiome")

setwd("~/Desktop/OysterMicrobiome")

## Read in OTU count data by sample file and the taxonomy file from MOTHUR output.
OTUCount <- import_mothur(mothur_shared_file = "stability.opti_mcc.shared",
                          mothur_constaxonomy_file = "stability.cons.taxonomy")
## Import phylogenetic tree into microbiome data set.
phy_tree(OTUCount) <- rtree(ntaxa(OTUCount), rooted=TRUE, tip.label=taxa_names(OTUCount))
## Determine a sequence count by sample threshold for subsampling.
OTUThreshold = as.integer(1e-3*sum(sample_sums(OTUCount)))

## Read in experimental design variables, such as the week the sample was taken.
factors <- read.table("MicrobiomeFactors.csv", header=TRUE, sep=",",as.is=T)
## Convert experimental design variables into a dataframe.
## The design variables for the oyster microbiome study are:
## FeedType (TET, ISO, or CHAE), SampleType (FEED, GUT, or FECAL)
## WeekFromStart (1-12), SampleID (1-5, TETFeed, ISOFeed, and CHAEFeed)
## Group (All sample group names), ExperimentalStatus (EXP or CON)
factors <- sample_data(data.frame(factors, row.names=sample_names(OTUCount)))
## Merge relative OTU abundance data with design variables into a Phyloseq object.
microbiomeRaw <- merge_phyloseq(OTUCount,factors)
## Update the taxonomic order names.
colnames(tax_table(microbiomeRaw)) <- c("kingdom", "phylum", "class", "order", "family",  "genus")

## Filter out samples which contain less than a certain fraction of the total
## sequence count.  This is used to help make sampling uniform across samples.
microbiomeRaw = prune_samples(sample_sums(microbiomeRaw) > OTUThreshold, microbiomeRaw )
set.seed(42)
microbiomeRaw<-rarefy_even_depth(microbiomeRaw, sample.size = OTUThreshold,rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

# Subset samples to look at core microbiomes by groups
#microbiome=subset_samples(microbiomeRaw,SampleType=="FECAL")
#microbiome=subset_samples(microbiome,FeedType=="TET")
# Get relative abundance for each sample in each group.
#microbiome <- microbiome::transform(microbiome,"compositional")
# Detection frequency for OTUs which are observed with at least a non-zero sequence count.
#head(prevalence(microbiome,detection=0,sort=TRUE),n=30)
#microbiome.core <- core(microbiome, detection=0,prevalence=0.99)
#otu_table(microbiome.core)
# Taxonomic labels of core microbiome.
#tax_table(microbiome.core)
# Total core abundance in each sample (sum of abundances of the core members)
#core.abundance <-sample_sums(core(microbiome,detection=0,prevalence=0.9))
# With compositional (relative) abundances
#det <- c(0, 0.1, 0.5, 2, 5, 20, 50)/100
#prevalences <- seq(.05, 1, .05)
plot_core(microbiome, prevalences = prevalences, detections = det, plot.type = "lineplot") + xlab("Relative Abundance (%)")

# Select feed type to subsamples on
feed="TET"
# Subset samples to look at core microbiomes by groups
microbiome=subset_samples(microbiomeRaw,FeedType==feed)
microbiome=subset_samples(microbiome,SampleType=="FECAL")
# Get relative abundance for each sample in each group.
microbiome <- microbiome::transform(microbiome,"compositional")
# What relative abundance threshold counts as an observation.
abundanceThreshold=0.001
# What fraction of selected samples counts as an observation
prevalenceThreshold=0.90
# Select a core microbiome given these thresholds
microbiome.core <- core(microbiome, detection=abundanceThreshold,prevalence=prevalenceThreshold)
# Calculate the number of core taxa and the number of samples they're found in.
numCore = dim(otu_table(microbiome.core))[1]
numSamples = dim(otu_table(microbiome.core))[2]
# Subset this core microbiome from the large sample set.
MBCoreSubset = prune_taxa(names(sort(taxa_sums(microbiome.core)))[1:numCore],microbiome)
# Re-select original sample set to track the change in relative OTU abundance
# between sample sources for the same feed type.
microbiome=subset_samples(microbiomeRaw,FeedType==feed)
# Get relative abundance for each sample in each group.
microbiome <- microbiome::transform(microbiome,"compositional")
MBSubset = prune_taxa(names(sort(taxa_sums(MBCoreSubset))),microbiome)

## Plot the most abundant OTUs
b = plot_bar(MBSubset,fill="family", facet_grid = ~SampleType~FeedType)
b = b + geom_bar(aes(color=genus, fill=family), stat="identity", position="stack")
b = b + ggtitle(paste(numCore,"OTUs observed in",numSamples,"samples","\n Relative abundance cutoff:",abundanceThreshold,"\n Fraction of samples cutoff:",prevalenceThreshold,sep=" "))
b = b + labs(y = "Relative sequence abundance")
b = b + scale_x_discrete(labels=NULL)
b

