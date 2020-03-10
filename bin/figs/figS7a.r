# Install packages
pacman::p_load("ggplot2", "data.table", "tidyverse", "RColorBrewer", 
               "ggpubr")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/adenylpred_analysis/")

# Read in the dataset
filpath <- "data/1-24-18 Nlt Ole Mix/"
fils <- list.files(filpath, full.names = T)
filall <- fils[!grepl("MS 1", fils)]

dat <- lapply(paste0(filall), function(x) {
  fread(paste0(x, "/tic_front.csv"), data.table = F)}) 

names(dat) <- word(sep = "1-23-18 ", filall, 2)
names(dat)

dff <- dat %>%
  map2_df(names(dat), ~ mutate(.x, ID = .y)) %>%
  bind_rows() 
head(dff)

shift.df <- dff
head(shift.df)

df2 <- shift.df %>% 
   dplyr::filter(dplyr::between(V1, left = 10.5, 
                         right = 10.8)) %>%
  mutate(grp = case_when(str_detect(ID, "MS 2") ~ "OleA + NltD",
                         str_detect(ID, "MS 3") ~ "OleA + NltD + NltC",
                         str_detect(ID, "MS 4") ~ "OleA + NltD + OleC"))

pal2 <- c("firebrick", "dodgerblue", "#98D277")
colnames(df2)[2] <- "V2"

df2$grp <- factor(df2$grp, levels = c("OleA + NltD + NltC",
  "OleA + NltD + OleC",
  "OleA + NltD"))

pl2 <- ggplot(data = df2,  aes(x = V1, y = V2, colour = grp, group = ID)) +
  theme_bw() +
  xlab("Retention time (min)") +
  ylab("GC-FID Peak Area") +
  facet_grid(rows = vars(grp)) +
  geom_point(size = 0.0000000001, alpha = 0.01) +
  geom_path(size = 0.75) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size=10), axis.title = element_text(size = 14),
        legend.text = element_text(size = 14), legend.title = element_blank(),
        legend.position = "none",
        axis.ticks.x = element_blank(), axis.ticks.y = element_line(),
        plot.margin = unit(c(0, 1, 0, 1), "cm"),
        panel.border = element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_blank()) +
  scale_color_manual(values = pal2) +
  theme_pubr()
pl2

pdf("output/figS7A_blactone.pdf", width = 3.5, height = 4)
pl2
dev.off()

df2 <- shift.df %>% 
  dplyr::filter(dplyr::between(V1, left = 11.9, 
                               right = 12.3)) %>%
  mutate(grp = case_when(
    str_detect(ID, "MS 2") ~ "OleA + NltD",
    str_detect(ID, "MS 3") ~ "OleA + NltD + NltC",
    str_detect(ID, "MS 4") ~ "OleA + NltD + OleC"))
colnames(df2)[2] <- "V2"

df2$grp <- factor(df2$grp, levels = c("OleA + NltD + NltC",
  "OleA + NltD + OleC",
  "OleA + NltD"))

pl2 <- ggplot(data = df2,  aes(x = V1, y = V2, colour = grp, group = ID)) +
  theme_bw() +
  xlab("Retention time (min)") +
  ylab("GC-FID Peak Area") +
  facet_grid(rows = vars(grp)) +
  geom_point(size = 0.0000000001, alpha = 0.01) +
  geom_path(size = 0.75) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size=10), axis.title = element_text(size = 14),
        legend.text = element_text(size = 14), legend.title = element_blank(),
        legend.position = "none",
        axis.ticks.x = element_blank(), axis.ticks.y = element_line(),
        plot.margin = unit(c(0, 1, 0, 1), "cm"),
        panel.border = element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        strip.text.x = element_blank()) +
  scale_color_manual(values = pal2) +
  theme_pubr()
pl2
dev.off()

pdf("output/figS7A_bhydroxy_acid.pdf", width = 3.5, height = 4)
pl2 
dev.off()
