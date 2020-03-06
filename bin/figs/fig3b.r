# Install packages
pacman::p_load("ggplot2", "readxl", "tidyverse", "RColorBrewer", "ggforce", "ggmosaic")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/adenylpred_analysis")

# Read in the antismash distribution table
raw <- read_excel("data/rf_predictions_antismash.xlsx") 

# Take only the rows of interest
tbl <- raw[1:8,]
tbl$Enzyme_class[tbl$Enzyme_class == "VLACSBILE"] <- "VLACS" # fix names

# Set color palette
pal1 <- c("#8B4513", "#E41A1C", "#377EB8", "#F781BF", "goldenrod", "#984EA3", "blue1", "#92d050")

tbl$Enzyme_class <- factor(tbl$Enzyme_class, levels = tbl$Enzyme_class)
tbl$log_Num_predicted <- log10(tbl$Num_predicted)

pdf("output/fig3b.pdf", width = 5, height = 3)
ggplot(data = tbl, aes(y = Num_predicted, x = Enzyme_class))+
  geom_bar(stat="identity", fill = pal1) +
  theme_bw() +
  scale_x_discrete(name = "Enzyme function") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_log10(expand = c(0, 0), name = "Log10 # of genes predicted") +
  annotation_logticks(sides = "l")  
dev.off()


