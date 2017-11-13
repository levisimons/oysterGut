library("plyr")
library(dplyr)
library("ggplot2")
library(phyloseq)
library("ape")
library("vegan")
library("microbiome")
library(data.table)

setwd("~/Desktop/OysterMicrobiome")

## Read in OTU count data by sample file and the taxonomy file from MOTHUR output.
OTUCount <- import_mothur(mothur_shared_file = "stability.opti_mcc.shared",
                          mothur_constaxonomy_file = "stability.cons.taxonomy")
## Import phylogenetic tree into microbiome data set.
phy_tree(OTUCount) <- rtree(ntaxa(OTUCount), rooted=TRUE, tip.label=taxa_names(OTUCount))

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

## Plot log of read depth per sample.
n <- sample_sums(microbiomeRaw)
bn <- barplot(n,xaxt="n",log="y",main="Log of the total sequence count per sample",ylab="Log sequence total per sample",xlab="Samples")

## If you want to scale OTU by relative sequence abundance.
## This is needed for any beta diversity analysis.
## OTUs with a rarity below a given cutoff are removed from this analysis.
#microbiome <- subset_samples(microbiomeRaw,SampleType=="FECAL")
microbiome <- subset_samples(microbiomeRaw,PhaseAndStatus=="EXP_FECAL_MONTH_1"|PhaseAndStatus=="EXP_FECAL_MONTH_2"|PhaseAndStatus=="EXP_FECAL_MONTH_3"|PhaseAndStatus=="CON_FECAL_MONTH_1"|PhaseAndStatus=="CON_FECAL_MONTH_2"|PhaseAndStatus=="CON_FECAL_MONTH_3")
microbiome <- microbiome::transform(microbiome,"compositional")

## Load distance methods.
dist_methods <- unlist(distanceMethodList)
print(dist_methods)
## Remove the user-defined distances.
## Choose 2 for weighted unifrac and 8 for bray-curtis.
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
  ordinate(microbiome,"MDS",distance=iDist)
  ## Make plot
  # Don't carry over previous plot (if error, p will be blank)
  p <- NULL
  # Create plot, store as temp variable, p
  p <- plot_ordination(microbiome, iMDS, color="PhaseAndStatus")
  # Increase font size
  p <- p + theme(text = element_text(size = 10))+ stat_ellipse(aes(group=PhaseAndStatus))
  # Save the graphic to file.
  plist[[i]] = p
}
p + geom_point(size=1, alpha=1)
pdf()
pdf('TreatmentAndControlMicrobiomeWunifrac.pdf',7,7)
plot(p)
dev.off()


## Plot Shannon index of samples.  Color by a design variable.
plot_richness(microbiomeRaw, measures=c("Shannon","Observed"),color="SampleType",title="Alpha diversity metrics \n oyster microbiome data")

## Store alpha diversity metrics, concatenated with experimental variables,
## in a single dataframe for subsequent significance tests.
alphaDiversity <- estimate_richness(microbiomeRaw,measures=c("Shannon","Observed"))
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

## Plot alpha diversity metrics versus time
a <- ggplot(alphaDiversity, aes(WeekFromStart, Shannon, color = ExperimentalStatus))
a = a + labs(x="Week", y="Sample Shannon index")
a = a + geom_point(size=3.5, alpha=1)
a = a + ggtitle("Sample Shannon index over time")
a
b <- ggplot(alphaDiversity, aes(WeekFromStart, Observed, color = ExperimentalStatus))
b = b + labs(x="Week", y="Sample OTU richness")
b = b + geom_point(size=3.5, alpha=1)
b = b + ggtitle("Sample OTU richness index over time")
b

## Plot evenness of samples by OTU abundance.
TaxaEvenness <- evenness(microbiomeRaw,index="all",zeroes=TRUE)
TaxaEvenness <- merge(TaxaEvenness,factors,by="row.names")
t <- ggplot(data=TaxaEvenness,aes(x=Group,y=pielou,color=FeedType))+geom_point(size=3.5, alpha=1)
t <- t + labs(x="Sample") + ggtitle(paste("Sample evenness\nRarefied sampling depth:",OTUThreshold,"reads",sep=" "))+theme(axis.text.x=element_blank())
t 

## Perform a Wilcoxon statistical test of significance on
## alpha diversity metrics separated by a design variable.
wtest <- wilcox.test(a,b)
wtest
## Compute Z score of test.
qnorm(wtest$p.value)

## Perform a PERMANOVA using a set number of permutations on a particular
## beta diversity metric and the significance of a particular design variable.
MBSubset = subset_samples(microbiome, Phase=="3")
microbiomeDF = as(sample_data(MBSubset), "data.frame")
microbiomeAdonis = adonis(distance(MBSubset,method="bray")~ExperimentalStatus,data=microbiomeDF,permutations = 10000)
microbiomeAdonis

## If you want to plot the beta diversity distance metrics
## in a grid of plots using available data.
df = ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
p = ggplot(df, aes(Axis.1, Axis.2, color=SampleType))
p = p + geom_point(size=1.5, alpha=0.5)
p = p + facet_wrap(~distance, scales="free")
p = p + ggtitle("Multidimensional Distance Scalings \n various distance metrics for oyster microbiome data")
p

## Bar plot of relative OTU abundances.
g = plot_bar(microbiome, fill="phylum", facet_grid = ~PhaseAndStatus)
g = g + geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack")
g = g + ggtitle(paste("Relative OTU abundance by feed \n Even sampling depth of ", OTUThreshold, "sequences",sep=" "))
g

## Merge relative OTU abundance data with design variables into a Phyloseq object.
## This is done to select the most abundant taxa from the dataset for
## more downstream analysis.
microbiomeRaw <- merge_phyloseq(OTUCount,factors)
colnames(tax_table(microbiomeRaw)) <- c("kingdom", "phylum", "class", "order", "family",  "genus")

## Find the most abundant OTUs as selected by an experiment variable
## such as feed type.
aTaxa = 10
microbiome <- subset_samples(microbiomeRaw, Phase=="3")
microbiome <- subset_samples(microbiome, SampleType!="GUT")
microbiome <- merge_samples(microbiome,"PhaseAndStatus")
microbiome <- microbiome::transform(microbiome,"compositional")
MBAbundant = sort(taxa_sums(microbiome), TRUE)[1:aTaxa]
MBSubset = prune_taxa(names(MBAbundant), microbiome)

## Plot the most abundant OTUs
b = plot_bar(MBSubset, fill="genus")
b = b + theme(legend.text = element_text(size = 10))
b = b + theme(plot.title = element_text(size = rel(2))) + theme(legend.title = element_text(size=10))
b = b + theme(axis.text.x = element_text(size=10))+ theme(axis.text.y = element_text(size=10))
b = b + theme(axis.title.y = element_text(size = rel(1.5), angle = 90)) + theme(axis.title.x = element_text(size = rel(1.5)))
b = b + theme(axis.ticks = element_line(size = 2))
b = b + geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")
#b = b + ggtitle(paste("Most abundant",aTaxa,"OTUs by sample origin.  Weeks 9 to 12.",sep=" "))
b = b + labs(y = "Relative sequence\nabundance", x= "Sample groups")
b
pdf()
pdf('TopTaxaPhase3.pdf',7,7)
plot(b)
dev.off()
