# Scripts to color phylogenetic tree by substrate specificity (Figure 2)
# Serina Robinson
# March 5, 2020

# Install packages
pacman::p_load("treeio", "ggtree", "stringr", "readxl", "RColorBrewer", "ggplot2")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/adenylpred_analysis/")

# Read in the RaxML output
# RaxML version 8.2.9_sse3_pthread was used
# raxmlHPC-PTHREADS-SSE3 -f a -x 1234 -m PROTGAMMAJTT -s 671_seqs_HMMaligned_converted_uppercased.fasta -n 671_seqs_ANL_20190708_bootstrapped.nwk -# 100 -p 1234
phylo <- treeio::read.newick("data/RAxML_bipartitions.671_seqs_ANL_20190708_bootstrapped.nwk")
  
# Make a phylogenetic tree
pl <- ggtree(phylo, layout = "circular")

# Append metadata
dd <- data.frame(label = pl$data$label,
                 grp = ifelse(pl$data$isTip, word(pl$data$label, -1, sep = "_"), "NONE"), 
                 size = ifelse(as.numeric(pl$data$label) > 75, 0.5, -1)) # NAs introduced by coercion error is ok
p2 <- pl %<+% dd

# Link substrate to tree tips
rawdat <- read_excel("data/combined_adenylate_forming_training_set_for_db_20191404.xlsx")
rawdat <- rawdat[-grep(paste0(c("HOLDOUT", "OTHER", "CAR", "NRPS"), collapse = "|"), rawdat$functional_class),] # remove extraneous classes

# Set the color palette
pal2 <- c("#E41A1C", "#92D050", "#377EB8", "#984EA3",
          "#FF7F00", "goldenrod", "gray68", "#A65628",
          "#F781BF", "blue1", "darkorchid1", "navy", 
          "gray68", "plum1", "blue1",
          "deepskyblue", "gold", "darkorchid1", 
          "deeppink2", "lightslateblue",
          "lightblue2", "darkseagreen1")

# Plot the tree
pdf("output/fig2.pdf", width = 12, height = 12)
par(mar=c(0.001,0.001,0.001,0.001))
ptree <- p2 +
  aes(color = grp) +
  scale_color_manual(values = pal2) +
  geom_nodepoint(size = p2$data$size[!p2$data$isTip], color = "gray40") +
  theme(legend.position = "right") +
  geom_treescale(offset = 1, fontsize = 4, x=3, y=1)
ptree
dev.off()

