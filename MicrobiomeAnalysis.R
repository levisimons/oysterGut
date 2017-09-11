library("plyr")
library("ggplot2")
library(phyloseq)
library("ape")
library("vegan")

setwd("~/Desktop/OysterMicrobiome")

## Read in OTU count data by sample file and the taxonomy file from MOTHUR output.
OTUCount <- import_mothur(mothur_shared_file = "stability.opti_mcc.shared",
                          mothur_constaxonomy_file = "stability.cons.taxonomy")
## Import phylogenetic tree into microbiome data set.
phy_tree(OTUCount) <- rtree(ntaxa(OTUCount), rooted=TRUE, tip.label=taxa_names(OTUCount))
## Determine a sequence count by sample threshold for subsampling.
OTUThreshold = as.integer(1e-4*sum(sample_sums(OTUCount)))

## Read in experimental design variables, such as the week the sample was taken.
factors <- read.table("MicrobiomeFactors.csv", header=TRUE, sep=",",as.is=T)
## Convert experimental design variables into a dataframe.
factors <- sample_data(data.frame(factors, row.names=sample_names(OTUCount)))
## Merge relative OTU abundance data with design variables into a Phyloseq object.
microbiomeRaw <- merge_phyloseq(OTUCount,factors)
## Update the taxonomic order names.
colnames(tax_table(microbiomeRaw)) <- c("kingdom", "phylum", "class", "order", "family",  "genus")

## Filter out samples which contain less than a certain fraction of the total
## sequence count.  This is used to help make sampling uniform across samples.
microbiomeRaw = prune_samples(sample_sums(microbiomeRaw) > OTUThreshold, microbiomeRaw )
set.seed(42)
microbiomeRaw<-rarefy_even_depth(microbiomeRaw, sample.size = 1e-4*sum(sample_sums(microbiomeRaw)),rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

## If you want to scale OTU by relative sequence abundance.
## This is needed for any beta diversity analysis.
## OTUs with a rarity below a given cutoff are removed from this analysis.
microbiome = transform_sample_counts(microbiomeRaw,function(x) x / sum(x))
## If you want to filter out OTUs with a relative abundance below a given cutoff.
microbiome = filter_taxa(microbiomeRaw, function(x) mean(x) > 1e-3, TRUE)

## Load distance methods.
dist_methods <- unlist(distanceMethodList)
print(dist_methods)
## Remove the user-defined distances.
## Choose 2 for weight unifrac and 8 for bray-curtis.
dist_methods <- dist_methods[c(2)]
print(dist_methods)

## Loop through each distance method, save each plot to a list, called plist.
plist <- vector("list", length(dist_methods))
names(plist) = dist_methods
for( i in dist_methods ){
  # Calculate distance matrix
  iDist <- distance(microbiome, method=i)
  # Calculate ordination
  iMDS  <- ordinate(microbiome, "MDS", distance=iDist)
  ## Make plot
  # Don't carry over previous plot (if error, p will be blank)
  p <- NULL
  # Create plot, store as temp variable, p
  p <- plot_ordination(microbiome, iMDS, color="SampleType")
  # Add title to each plot
  p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
  # Save the graphic to file.
  plist[[i]] = p
}
p
## If you want to plot the distance metrics
## in a grid of plots using available data.
df = ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
p = ggplot(df, aes(Axis.1, Axis.2, color=SampleType))
p = p + geom_point(size=1.5, alpha=0.5)
p = p + facet_wrap(~distance, scales="free")
p = p + ggtitle("Multidimensional Distance Scalings \n various distance metrics for oyster microbiome data")
p

## Plot Shannon index of samples.  Color by a design variable.
plot_richness(microbiomeRaw, measures=c("Shannon","Chao1"),color="SampleType",title="Alpha diversity metrics \n oyster microbiome data")

## Store alpha diversity metrics, concatenated with experimental variables,
## in a single dataframe for subsequent significance tests.
alphaDiversity <- estimate_richness(microbiomeRaw,measures=c("Shannon","Chao1"))
alphaDiversity <- as.data.frame(alphaDiversity)
alphaDiversity$Group <- rownames(alphaDiversity)
alphaDiversity <- alphaDiversity[,c(which(colnames(alphaDiversity)=="Group"),which(colnames(alphaDiversity)!="Group"))]
factors <- as.data.frame(factors)
alphaDiversity <- merge(factors,alphaDiversity,by=0)
alphaDiversity <- subset(alphaDiversity, select=-c(1,8))
colnames(alphaDiversity)[colnames(alphaDiversity)=="Group.x"] <- "Group"
alphaDiversity <- as.data.frame(alphaDiversity)

## Store alpha diversity metrics and select them by value for a particular design variable.
a <- alphaDiversity$Shannon[alphaDiversity$SampleType=="FEED"]
b <- alphaDiversity$Shannon[alphaDiversity$SampleType!="FEED"]

## Perform a Wilcoxon statistical test of significance on
## alpha diversity metrics separated by a design variable.
wilcox.test(a,b)

## Perform a PERMANOVA using a set number of permutations on a particular
## beta diversity metric and the significance of a particular design variable.
microbiomeDF = as(sample_data(microbiome), "data.frame")
microbiomeDist = distance(microbiome,method="wunifrac")
microbiomeAdonis = adonis(microbiomeDist ~ FeedType, microbiomeDF,permutations = 10000)
microbiomeAdonis
