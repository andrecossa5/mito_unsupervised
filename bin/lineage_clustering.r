# LINEAGE clustering method

library(tidyverse)
library(data.table)
library(LINEAGE)

# Paths
args <- commandArgs(trailingOnly=TRUE)

sample_name <- args[1]
path_afm <- args[2]
path_meta <- args[3]
n_repeats <- args[4]
n_cores <- args[5]

# Read data
afm <- fread(path_afm) %>% as.data.frame()
meta <- fread(path_meta) %>% as.data.frame()

# Add as row.names
row.names(afm) <- afm[, 1]
row.names(meta) <- meta[, 1]
afm <- afm %>% select(-index) %>% t() # Transpose
meta <- meta %>% select(-index)


##

extract_alt_wt <- function(x) {

    l1 <- strsplit(x, "_")[[1]][2]
    wt <- strsplit(l1, ">")[[1]][1]
    alt <- strsplit(l1, ">")[[1]][2]

    return(c(wt, alt))
}
alleles <- sapply(afm %>% row.names(), function(x) { extract_alt_wt(x) })

##

# Add alleles as columns
afm <- afm %>% as.data.frame() %>% mutate(altAllele=alleles[1,], refAllele=alleles[2,])

# LINEAGE workflow
afm[is.na(afm)] <- 0 # NA as 0
results <- lineage(data=afm, repeats=n_repeats, thread=n_cores)

# Format and write results
labels <- results$label
colnames(labels) <- c('label', 'cell')
row.names(labels) <- labels$cell
labels <- labels %>% select(-cell)

# Write
gt <- meta %>% select(GBC)
write.csv(gt, 'ground_truth.csv')
write.csv(labels, 'lineage_labels.csv')
