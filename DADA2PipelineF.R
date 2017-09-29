## This script assumes that the data you are starting with meets certain criteria:
## 1. Non-biological nucleotides have been removed (primers/adapters/barcodes…)
## 2. Samples are demultiplexed (split into individual per-sample fastqs)
## 3. If paired-end sequencing, the forward and reverse fastqs contain reads in matched order
library(dada2); packageVersion("dada2")
library("ape")
library("vegan")
library(ggplot2); packageVersion("ggplot2")
library(phyloseq); packageVersion("phyloseq")
library("plyr")

## Set working directory and path variable
setwd("~/Desktop/DADA2Analysis")
path <- "/Users/levisimons/Desktop/DADA2Analysis"
list.files(path)

## First we read in the names of the fastq files, and perform some string manipulation to get lists of the forward and reverse fastq files in matched order.
# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files(path, pattern="_R1_001.fastq"))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
# Specify the full path to the fnFs
fnFs <- file.path(path, fnFs)

## Examine quality profile of forward reads to check where read quality drops.
##fcut is where the quality profile drops on the graph.
plotQualityProfile(fnFs[1:2])
fcut=210

## Define the filenames for the filtered fastq.gz files.
filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))

## Filter the forward reads
out <- filterAndTrim(fnFs, filtFs, truncLen=fcut,
                     maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

## Learn the error rates of each sample.
## Determine the error rate for forward reads:
errF <- learnErrors(filtFs, multithread=TRUE)

## Visualize the error rates
plotErrors(errF, nominalQ=TRUE)

## Dereplicate the filtered fastq files to get the counts of unique reads.
derepFs <- derepFastq(filtFs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names

## Infer the sequence variants in each sample
## First with the forward reads:
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)

## construct a “sequence table” of our samples, a higher-resolution version of the “OTU table”
seqtab <- makeSequenceTable(dadaFs)

## Remove chimeric sequences
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
## Get dimensions for the seqtab.nochim matrix.
dim(seqtab.nochim)
## Get the percentage of reads which are non-chimeric sequences
sum(seqtab.nochim)/sum(seqtab)

## Check the number of reads that made it through each step in the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)

## Assign taxonomy down to the genus level.
set.seed(100) # Initialize random number generator for reproducibility
taxa <- assignTaxonomy(seqtab.nochim, "/Users/levisimons/Desktop/DADA2Analysis/silva_nr_v128_train_set.fa.gz", minBoot=80)
microbiomeRaw <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), tax_table(taxa))

## Prune away samples.
## F2R15 is poorly sampled.  F15R16 is the even mock community.  F16R16 is the staggered mock community.
microbiomeRaw <- prune_samples(sample_names(microbiomeRaw) != "F2R15"
                               & sample_names(microbiomeRaw) != "F15R16"
                               & sample_names(microbiomeRaw) != "F16R16", microbiomeRaw)

## Load the phylogenetic tree structure and experimental design variables
## into the Phyloseq object.
phy_tree(microbiomeRaw) <- rtree(ntaxa(microbiomeRaw), rooted=TRUE, tip.label=taxa_names(microbiomeRaw))
factors <- read.table("MicrobiomeFactors.csv", header=TRUE, sep=",",as.is=T)
factors <- sample_data(data.frame(factors, row.names=sample_names(microbiomeRaw)))
microbiomeRaw <- merge_phyloseq(microbiomeRaw,factors)
## Check taxonomic assignments to each sequence down to the genus level.
unname(taxa)

## Save output for downstream Phyloseq analysis.
saveRDS(microbiomeRaw,"/Users/levisimons/Desktop/DADA2Analysis/microbiomeRaw.rds")

## Assign species level taxonomy
taxa.plus <- addSpecies(taxa, "/Users/levisimons/Desktop/DADA2Analysis/silva_species_assignment_v128.fa.gz", verbose=TRUE)

## Save output again for downstream Phyloseq analysis.
saveRDS(microbiomeRaw,"/Users/levisimons/Desktop/DADA2Analysis/microbiomeRaw.rds")

## If the DADA2 output needs to be loaded for subsequent Phyloseq analysis
## load file back as a Phyloseq object
microbiomeRaw <- readRDS("/Users/levisimons/Desktop/DADA2Analysis/microbiomeRaw.rds")

## Subset the microbiome data by various experimental variables.
microbiome <- subset_samples(microbiomeRaw,PhaseAndStatus!="FEED1" & PhaseAndStatus!="FEED2" & PhaseAndStatus!="FEED3")
## Scale the reads to relative abundance within each sample.
microbiome <- transform_sample_counts(microbiome, function(x) x/sum(x))

## Perform a PERMANOVA using a set number of permutations on a particular
## beta diversity metric and the significance of a particular design variable.
MBSubset = microbiome
microbiomeDF = as(sample_data(MBSubset), "data.frame")
microbiomeAdonis = adonis(distance(MBSubset,method="wunifrac")~PhaseAndStatus,data=microbiomeDF,permutations = 10000)
microbiomeAdonis

## Load beta diversity distance methods.
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
  p <- p + theme(text = element_text(size = 20))+ stat_ellipse(aes(group=PhaseAndStatus))
  # Add title to each plot
  p <- p + ggtitle(paste("MDS using distance method ", i,sep=""))
  # Save the graphic to file.
  plist[[i]] = p
}
p + geom_point(size=1, alpha=1)
