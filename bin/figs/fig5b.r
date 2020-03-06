# Install packages
pacman::p_load("ggplot2", "data.table", "tidyverse", "RColorBrewer")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/adenylpred_analysis/")

# Read in the dataset
filpath <- "data/NltAB_OleA_condensation/"

fils <- list.files(filpath, full.names = T)
filall <- fils[grep("_2\\.csv", fils)]

filshrt <- gsub("\\.csv", "", list.files(filpath)[grep("_2\\.csv", list.files(filpath))])

dat <- lapply(paste0(filall), function(x) {
  fread(x, data.table = F)}) 
head(dat)
names(dat) <- filshrt

dff <- dat %>%
  map2_df(names(dat), ~ mutate(.x, ID = .y)) %>%
  bind_rows() 
head(dff)


df2 <- dff %>%
  dplyr::filter(!grepl("NltAplusB", ID)) %>%
  dplyr::filter(dplyr::between(V1, left = 21.9,
                       right = 22.2)) %>%
  mutate(grp = case_when(str_detect(ID, "OleA") ~ "OleA",
                         str_detect(ID, "NltAB") ~ "NltA + NltB",
                         str_detect(ID, "NltB") ~ "NltB",
                         str_detect(ID, "_NltA_") ~ "NltA"))


# Set color palette
pal2 <- c("#E41A1C", "#92D050", "dodgerblue", "#984EA3") 
colnames(df2)[2] <- "V2"

df2$grp <- factor(df2$grp, levels = c("OleA",
                                      "NltA + NltB",
                                      "NltB",
                                      "NltA"))

pl2 <- ggplot(data = df2,  aes(x = V1, y = V2, colour = grp, group = ID)) +
  theme_bw() +
  xlab("Retention time (min)") +
  ylab("GC-FID Peak Area") +
  facet_grid(rows = vars(grp)) +
  geom_point(size = 0.0000000001, alpha = 0.01) +
  geom_path(size = 0.75) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 14, color = "black"),
        legend.text = element_text(size = 12, color = "black"), legend.title = element_blank(),
        legend.position = "none",
        axis.ticks.x = element_line(), axis.ticks.y = element_line(),
        plot.margin = unit(c(0, 1, 0, 1), "cm"),
        panel.border = element_blank()) +
  scale_color_manual(values = pal2)
pl2

pdf("output/fig5b.pdf", width = 3.5, height = 3.5)
pl2
dev.off()


