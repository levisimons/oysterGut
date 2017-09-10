library("plyr")
library("ggplot2")
library(phyloseq)
library("ape")

setwd("~/Desktop/OysterMicrobiome")

## Read in OTU count data by sample file and the taxonomy file from MOTHUR output.
OTUCount <- import_mothur(mothur_shared_file = "stability.opti_mcc.shared",
                          mothur_constaxonomy_file = "stability.cons.taxonomy")

phy_tree(OTUCount) <- rtree(ntaxa(OTUCount), rooted=TRUE, tip.label=taxa_names(OTUCount))

## If you want to scale OTU by relative sequence abundance.
OTUCount = transform_sample_counts(OTUCount,function(x) x / sum(x))
## If you want to filter out OTUs with a relative abundance below a given cutoff.
OTUCount = filter_taxa(OTUCount, function(x) mean(x) > 1e-2, TRUE)
## Read in experimental design variables, such as the week the sample was taken.
factors <- read.table("MicrobiomeFactors.csv", header=TRUE, sep=",",as.is=T)
## Convert experimental design variables into a dataframe.
factors <- sample_data(data.frame(factors, row.names=sample_names(OTUCount)))
## Merge relative OTU abundance data with design variables into a Phyloseq object.
microbiome <- merge_phyloseq(OTUCount,factors)
## Update the taxonomic order names.
colnames(tax_table(microbiome)) <- c("kingdom", "phylum", "class", "order", "family",  "genus")

## Load distance methods.
dist_methods <- unlist(distanceMethodList)
## Remove the user-defined distances
dist_methods <- dist_methods[(2:8)]
dist_methods <- dist_methods[-(2:6)]
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
  p <- plot_ordination(microbiome, iMDS, color="FeedType", shape="SampleType")
  # Add title to each plot
  p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
  # Save the graphic to file.
  plist[[i]] = p
}

## Plot the distance metrics in a grid of plots using available data.
df = ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
p = ggplot(df, aes(Axis.1, Axis.2, color=FeedType, shape=SampleType))
p = p + geom_point(size=3, alpha=0.5)
p = p + facet_wrap(~distance, scales="free")
p = p + ggtitle("MDS on various distance metrics for oyster microbiome data")
p