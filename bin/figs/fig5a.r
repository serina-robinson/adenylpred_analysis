# Install packages
pacman::p_load("readxl", "ggplot2", "tidyverse", "dplyr")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/adenylpred_analysis/")

# Read in the data
rawdat <- read_excel("data/NltC_activity_C20_C14.xlsx")
tm <- rawdat %>%
  dplyr::pull(Time)
rawdat

rawdat_long <- rawdat %>%
  dplyr::filter(Time != 1) %>%
  dplyr::filter(!Time %in% c(7:23)) %>%
  reshape2::melt(., id.vars = c("Time"), value.name = "C", variable.name = "id") %>%
  dplyr::mutate(substrate = word(id, sep = "_", 1),
                cofactor = word(id, sep = "_", 2),
                value = word(id, sep = "_", 3)) %>%
  dplyr::filter(cofactor != "ADP")

davg <- rawdat_long %>%
  dplyr::filter(value == "avg")  %>%
  dplyr::filter(cofactor == "AMP")

dstd <- rawdat_long %>%
  dplyr::filter(value == "stdev") %>%
  dplyr::filter(cofactor == "AMP")

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

pdf("output/fig5a.pdf", width = 4, height = 4)
ggplot(davg, aes(x=Time, y=C, color=substrate, group = substrate)) +
  geom_point(aes(shape = substrate)) +
  geom_line() +
  labs( y="AMP produced (µM)", x="Time (hours)")+
  geom_errorbar(data = dstd, aes(ymax=davg$C + dstd$C, ymin = davg$C - dstd$C), width=0.2,size=0.5) +
  theme(legend.title=element_blank(), axis.line=element_line(color="black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        panel.background=element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        text = element_text(size = 12, color = "black"),
        legend.key= element_rect(fill=NA, color=NA),
        legend.position=c(0.9, 0.2)) +
  scale_color_manual(values = c("orange", "navy",  "gray")) +
  ylim(0, 100)
dev.off()



