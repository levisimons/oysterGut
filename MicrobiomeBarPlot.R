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

microbiome = subset_taxa(microbiome,phylum=="Proteobacteria")
p = plot_bar(microbiome, fill="genus", facet_grid = ~WeekFromStart)
p + geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")

plot_heatmap(microbiome, taxa.label="class")
