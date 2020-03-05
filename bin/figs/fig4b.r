# Scripts to plot phylogenetic tree of nocardiolactone biosynthetic gene clusters (Figure 4b)
# Serina Robinson
# March 5, 2020

# Load packages
pacman::p_load("genoPlotR", "scales",
               "RColorBrewer", "ggthemes",
               "ggplot2","Biostrings",
               "ape","phangorn", "ade4")

#Set working directory
setwd("~/Documents/Wageningen_UR/github/adenylpred_analysis/")

# Read in NltC phylogenetic tree and convert to phylogeny
my_tre <- read.table("data/Nocardia_NltC_phylogeny_25seqs.nwk", header = F)
my_nwk <- ade4::newick2phylog(x.tre = my_tre[,1])

# Read in the genome neighborhoods (output from RODEO web tool http://rodeo.scs.illinois.edu/)
co <- read.csv("data/Nocardia_genome_neighborhoods.csv", stringsAsFactors=F)
acc <- substr(names(my_nwk$leaves), start = 1, stop = 12) # extract accession numbers
co <- co[grep(paste0(acc,collapse="|"), co$query.accession),] # find accession numbers in genome neighborhoods

# Create a list of genome neighborhoods
nb <- list()
cacc <- unique(co$query.accession)
for(i in 1:length(cacc)){
  nb[[i]]<-co[co$query.accession==cacc[i],]
}

# Create DNA segment objects
source("lib/get_dna_segs.r")
raw_segs <- get_dna_segs(nb)

# Trim the first gene in each neighborhood
for(i in 1:length(raw_segs)){
  raw_segs[[i]]<-raw_segs[[i]][2:(nrow(raw_segs[[i]])-1),]
}

# Get enriched PFAMs across all genome neighborhoods
source("lib/get_enriched_pfams.r")
pfam_desc <- get_enriched_pfams(co)
pfams <- pfam_desc$pfam
pfam_desc$desc[2:length(pfam_desc$desc)] <- c("NltA/NltB", "NltC", "NltD",
                                              "MMPL transporter", "TetR family regulator",
                                              "Short chain dehydrogenase", "Cytochrome p450",
                                              "4HBT Thioesterase superfamily", "Galactokinase galactose-binding signature",
                                              "Galactose-1-phosphate uridyl transferase", "Zinc-binding dehydrogenase")
# Color segs
source("lib/color_segs.r")
segs <- color_segs(raw_segs, pfams)
segs


# Fix names
fixnams <- lapply(1:length(segs), function(x){ gsub("\\.","_",segs[[x]]$name) })
tonam <- unlist(lapply(1:length(fixnams), function(x){fixnams[[x]][1]}))

for(i in 1:length(tonam)){
  tmp<-grep(tonam[i],names(my_nwk$leaves))
  names(segs)[i]<-names(my_nwk$leaves)[tmp]
}

# Reorder the gene clusters to match tree tip labels
source("lib/reorder_list.r")
segs2 <- reorder_list(segs, names(my_nwk$leaves))
names(segs2) <- names(my_nwk$leaves)

# Plot phylogenetic tree with gene clusters
pdf(paste0("output/fig4.pdf"), width=8, height=6)
plot_gene_map(segs2, tree = my_nwk, tree_branch_labels_cex = 30)
dev.off()

# Make PFAM legend
pdf(paste0("output/fig4_legend.pdf"), width=10, height=10)
plot.new()
pal <- c("gray", "#FFB341", "#92D050", "#F37FCC", "dodgerblue", "#E31A1C", "#CAB2D6",
         "#6A3D9A", "#FFFF99", "#B15928", "#34A02C", "black")
l <- legend("bottom", legend=pfam_desc$desc, fill=pal,
            col = "black", border = "black",
            bty="n", ncol = 2)
dev.off()
