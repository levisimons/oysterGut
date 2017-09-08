library("plyr")
library("ggplot2")
library(phyloseq)
library("ape")

setwd("~/Desktop/OysterMicrobiome")

## Read in OTU count data by sample file and the taxonomy file from MOTHUR output.
OTUCount <- import_mothur(mothur_shared_file = "stability.opti_mcc.shared",
                          mothur_constaxonomy_file = "stability.cons.taxonomy")

phy_tree(OTUCount) <- rtree(ntaxa(OTUCount), rooted=TRUE, tip.label=taxa_names(OTUCount))

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

MB.ord <- ordinate(microbiome,"NMDS","wunifrac")
p = plot_ordination(microbiome,MB.ord,type="samples", color="FeedType", title="Samples by feed type")
print(p)
