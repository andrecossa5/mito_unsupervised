# Phylo script for first trees with ape and ggtree

# Code
library(Rfast)
library(ape)
library(phangorn)
library(tidyverse)
library(ggtree)
library(TreeTools)
library(treeio)
library(viridis)
library(phytools)


##


################################################################ Utils
annotate_tree <- function(tree, df, cov='GBC') {

    if (is.character(df[[cov]])) {
        labels <- df[[cov]] %>% unique()
        L <- lapply( labels, function(x) { x <- df %>% filter(.data[[cov]] == x) %>% row.names() } )
        names(L) <- labels
        tree <- groupOTU(tree, L)

    } else if (is.numeric(df[[cov]])) {
        tree <- full_join(
            tree %>% as_tibble(), 
            df %>% select(cov) %>% rownames_to_column('label'), 
            by = 'label'
        ) %>% as.treedata()
    }
    return(tree)
}


##


plot_tree_cov <- function(tree, cov='group', colors=NULL) {
    p <- ggtree(treeNJ, aes(color=.data[[cov]]), layout='circular') +
    theme(legend.position="none") + 
    scale_color_manual(values=colors, name=NULL)
    return(p)
}
################################################################


##


# Set path
path_main <- '/Users/IEO5505/Desktop/MI_TO/images_DIPA/images'

# Load data
X_char <- read.csv(paste0(path_main, '/AML_clones_afm.csv'), row.names=1)
meta <- read.csv(paste0(path_main, '/AML_clones_meta.csv'), row.names=1)
cells <- row.names(X_char)
df <- cbind(X_char, meta)

# Compute NJ tree
treeNJ  <- NJ(dist(X_char))

# Annotate with clone type and top important muts for clone ...
clone <- 'CCGCAATAAATGCGACTT'
df <- df %>% mutate(
    clone=case_when(
        GBC == clone ~ clone, 
        GBC != clone ~ 'other'
    )
)

# Plot 
pdf(paste0(path_main, '/NJ_', clone, '.pdf'), width=15, height=15)
treeNJ <- annotate_tree(treeNJ, df, clone) 
p <- ggtree(treeNJ, aes(color='clone'), layout='circular') +
    #theme(legend.position="none") +
    scale_color_manual(values=c('black', 'red'), name='Clone type')
print(p)
dev.off()
