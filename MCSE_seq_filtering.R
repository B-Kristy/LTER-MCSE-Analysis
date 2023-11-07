### Post-processing sequence filtering for [Insert Manuscript Title Here]
### DOI:
### Please contact Brandon Kristy (kristybr@msu.edu) for any questions
 
rm(list=ls()) # Remove all objects from Environment
# MAC OS
# Note -- Change this to the GitHub Repository Directory once downloading onto your PC
setwd("~/Dropbox/Projects/LTER_MCSE/data/seq_data/")
# Required Libraries
library(data.table)
library(phyloseq)
library(plyr)
library(microbiome)
library(tidyverse)
library(vegan)



## Raw phyloseq object -------------------------------------------------------------------
# Load OTU table 
otu_table <- read.csv("otu_table.csv", header=T, row.names = 1)
otu <- otu_table(otu_table, taxa_are_rows = TRUE)

# Load metadata
metadata <- read.csv("sample_metadata.csv", row.names = 1)
samp <- sample_data(metadata)

# Reformat taxonomy; Remove Unassigned at the Kingdom level
# 60559 taxa; 55628 were classified 
tax_table <- read.table("nifH_taxonomy.out", sep = ";")

tax_table <- tax_table %>%
  subset(select = c(V1, V3)) %>%
  filter(!grepl('k__unassigned', V3)) %>%
  mutate(V3=gsub("[a-z]__", "", V3)) %>%
  separate_wider_delim(V3, delim = "|", names = c("Kingdom","Phylum","Class","Order","Family","Genus","Species")) %>%
  remove_rownames %>%
  column_to_rownames(var='V1')

tax_table <- as.matrix(tax_table)
tax <- tax_table(tax_table)

# Prepare phyloseq object
phylo <- merge_phyloseq(otu, samp, tax)
# Output raw phyloseq object 
save(phylo, file = "raw_phyloseq.RData")

## Filtered phyloseq object  -------------------------------------------------------------------
# Remove singletons
otus.mine <- prune_taxa(taxa_sums(phylo) > 1, phylo)
tax.remove <- ntaxa(phylo) - ntaxa(otus.mine) #0 OTUs removed

# Remove plant contaminants
bact.p <- subset_taxa(otus.mine, Kingdom != "\tEukaryota")

n.filtered <- ntaxa(otus.mine)- ntaxa(bact.p) # 248 OTUs removed

# remove low sequence coverage samples
## all samples 
sorted <- sort(sample_sums(bact.p))  
min <- min(sample_sums(bact.p)) #49
max <- max(sample_sums(bact.p)) # 83933
mean <- mean(sample_sums(bact.p)) #26,402.8
median <- median(sample_sums(bact.p)) #26601

bact.p_filtered <- prune_samples(sample_sums(bact.p)>=1000, bact.p)
save(bact.p_filtered, file = "phylo_filtered.RData")

## Rarefied phyloseq object -------------------------------------------------------------------
set.seed(37920)
# Rarefaction depth determined from rarefaction curve -- see Figure S1
bact.rfy <- rarefy_even_depth(bact.p_filtered, sample.size = 1289, replace = FALSE, rngseed = FALSE, verbose = FALSE)
save(bact.rfy, file= "phylo_rarefied.RData")

