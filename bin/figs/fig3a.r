# Scripts to plot and color clustering diagram (Figure 2A)
# Serina Robinson
# March 5, 2020

# Install packages
pacman::p_load("data.table", "scales", "tidyverse", "stringr", 
               "Biostrings","DECIPHER","igraph","RColorBrewer")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/adenylpred_analysis/")

# Read in the all-vs-all BLAST results
# cd /home/wackett/robi0916/wageningen/data/antismashv2/
# module load ncbi_blast+
# makeblastdb -in 3319_comb_cdhit40_for_ssn.fasta -dbtype prot -out db/3319_comb_cdhit_blastdb
# blastp -db db/3319_comb_cdhit_blastdb -query 3319_comb_cdhit40_for_ssn.fasta -outfmt 6 -out 3319_comb_cdhit40_all_v_all.tsv -num_threads 8
allvall <- fread("data/3319_comb_cdhit40_all_v_all.tsv", stringsAsFactors = F, data.table = F)

# Read in the random forest predictions
clusassn <- read_csv("data/2344_antismashv2_full_rf_prediction_results.csv") 
colnames(clusassn)[1] <- "V1"

# Join with predictions by the name
join_assnV1 <- allvall %>%
  left_join(., clusassn, by = "V1") %>%
  mutate(prob0.6V1 = coalesce(prob0.6, V1)) %>%
  mutate(prob0.5V1 = coalesce(prob0.5, V1)) %>%
  mutate(prob0.4V1 = coalesce(prob0.4, V1)) %>%
  mutate(prob0.3V1 = coalesce(prob0.3, V1)) %>%
  mutate(prob0.2V1 = coalesce(prob0.2, V1)) %>%
  mutate(prob0.1V1 = coalesce(prob0.1, V1))

colnames(clusassn)[1] <- "V2"
join_assnV2 <- allvall %>%
  left_join(., clusassn, by = "V2") %>%
  mutate(prob0.6V2 = coalesce(prob0.6, V2)) %>%
  mutate(prob0.5V2 = coalesce(prob0.5, V2)) %>%
  mutate(prob0.4V2 = coalesce(prob0.4, V2)) %>%
  mutate(prob0.3V2 = coalesce(prob0.3, V2)) %>%
  mutate(prob0.2V2 = coalesce(prob0.2, V2)) %>%
  mutate(prob0.1V2 = coalesce(prob0.1, V2))


# Color by functional class
namvec <- c("ARYL", "BLS", "FAAL", "LACS", "LUCIFERASE", "MACS", "NRPS", "SACS", "VLACSBILE")

# Color AMP clustering diagram
pal <- c("#E41A1C", "#92D050", "#377EB8", "#984EA3", "#FF7F00",   
          "goldenrod", "#A65628", "#F781BF", "blue1", 
          "gray68", "darkorchid1", "navy", "plum1",
          "deepskyblue", "gold", "deeppink2", "lightslateblue",
          "lightblue2", "darkseagreen1")

# Color by probability of 0.6
source("lib/make_cluster_diagram.r")
shrt <- allvall[,(c(1:2,11))]
colnames(shrt) <- c("prot1", "prot2", "eval")
gr <- make_cluster(shrt, thresh = 36, bynum = 1, namvec, pal)
net <- gr[[1]]
g <- gr[[2]]

# Manually color some that were incorrectly identified
V(g)$color[grep("Schizosaccharomyces_pombe_C6.through.C11_", V(g)$name)] <- "#A65628"
V(g)$color[grep("cerevisiae_C6.through.C11_", V(g)$name)] <- "#A65628"
V(g)$color[grep("X5IJ97_9ACTN_Streptomyces_sp._C6.through.C11_mediumchain_MACS", V(g)$name)] <- "#A65628"
V(g)$color[grep("ALA09345.1_Micromonospora_sp._C2.through.C5_shortchain_NRPS", V(g)$name)] <- "#F781BF"
thresh <- 36
l <- layout_components(g)

jpeg(file = paste0("output/fig3a.jpg"), width = 2000, height = 2000)
par(mar=c(1,1,1,1))

pl<-plot(g, vertex.label = NA, 
         vertex.size = 1.5,
         layout=l, edge.color = "gray40", edge.width=0.3)
title(main = paste0("BLAST e-value cut-off: 1e-",thresh))
dev.off()

cleg <- gr[[3]]
cleg_nams <- c("Aryl-CoA ligases \n (ARYL)", "Medium chain acyl-CoA \n synthetases (MACS) \n C6 - C12", 
               "Short chain acyl-CoA \n synthetases (SACS) \n C2 - C5", "Nonribosomal peptide \n synthases (NRPS)",
               "Long chain acyl-CoA \n synthetases (LACS) \n C13 - C17", "Fatty acyl-AMP \n ligases (FAAL)", 
               "Very-long chain and bile \n acyl-CoA synthetases \n (VLACSBILE) C18+", "Luciferases (LUC)",
               "Beta-lactone \n synthetases (BLS)", "None")

pdf(file=paste0("fig3_legend.pdf"), width = 40, height = 20)
plot.new()
legend("center", legend=cleg_nams, fill = as.character(cleg$colors), 
       bty="n", horiz = T)
dev.off() 


