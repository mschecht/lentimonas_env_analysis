#!/usr/bin/env Rscript

# This script will plot the mitags placed on the 
# reference tree

#######################################
# Check if basic packages are installed 
# if not here then install them
#######################################

# Packages needed for script
packages <- c("tidyverse", "ggtree")

# Select packages that are NOT already installed
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]

# Install them
if(length(new.packages)) install.packages(new.packages)

###############
# Load packages
###############
suppressMessages(lapply(packages, library, character.only = TRUE))

#############
# Upload data
#############

# Reference tree
NWK_ref <- read.tree("data/external/lenti_16s_masked_phylib_phyml_tree.txt")

# Guppy to_csv + edpl
info <- read_csv("data/processed/mitags/query_info.csv", col_names = TRUE)

###########################
# Prepare data for plotting
###########################

# Calculate the amount of query sequences 
# placed on each edge of reference tree
ABUND_TABLE <- info %>% 
  group_by(edge_num) %>% 
  count() %>% 
  rename("edge" = edge_num, "mitags_placed" = n) %>% 
  mutate(color = "indianred")

# Tidy format the reference tree
tree_tidy <- fortify(NWK_ref)

# create metadata 
ABUND_TABLE$edge2map <- ABUND_TABLE$edge
ALL_EDGES <- data.frame(node = NWK_ref$edge[,2], 
                        parent = NWK_ref$edge[,1],
                        edge= 1:nrow(NWK_ref$edge)) %>% as_tibble()

ALL_EDGES_extended <- ALL_EDGES %>% 
  left_join(ABUND_TABLE, by = "edge") %>% 
  left_join(tree_tidy) %>% 
  mutate(bkg = case_when(grepl("rRNA", label)  ~ 'rRNA'))

######
# Plot
######

# Height and width
w <- 15
h <- 15

# Plot
p <- ggtree(NWK_ref, layout='circular') %<+% ALL_EDGES_extended +
  geom_tippoint(aes(x=branch, size = mitags_placed, color = color),       # Make bubbles on edges
                vjust=0, 
                alpha = 0.7) +
  geom_tiplab(aes(angle = angle, color = color, fill = bkg),              # Label tips
              size = 2, 
              align = TRUE, 
              linesize = 0.2, 
              linetype = "dotted") +
  theme(legend.position = "bottom", 
        legend.key = element_blank()) # no keys

p + guides(colour = "none") +                                             # remove color legend
  geom_cladelabel(node = 215,                                             # Label Lentimonas clones
                  "Lentimonas\n  Clones", 
                  offset=0.09, 
                  barsize=2, 
                  offset.text=0.005, 
                  fontsize=3)


pdf("~/Desktop/lentimonas_placed.pdf", height = h, width = w)

dev.off()