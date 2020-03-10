# Install packages
pacman::p_load("ggplot2", "data.table", "tidyverse", 
               "RColorBrewer", "baseline")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/adenylpred_analysis/")

# Read in the dataset
filpath <- "data/20180626_NltAB/"
fils <- list.files(filpath, full.names = T)
fils
mut <- fils[!grepl("old|BPER|WT", fils)]
wt <- fils[grepl("redo", fils)]
filall <- c(wt, mut)[grep("_2", c(wt, mut))]

filshrt <- gsub("\\.csv", "", filall)
dat <- lapply(paste0(filall), function(x) {
  fread(paste0(x, "/tic_front.csv"), data.table = F)}) 
warnings()

names(dat) <- filshrt

dff <- dat %>%
  map2_df(names(dat), ~ mutate(.x, ID = .y)) %>%
  bind_rows() 
head(dff)

df2 <- dff %>%
  dplyr::filter(!grepl("NltAplusB", ID)) %>%
  dplyr::filter(dplyr::between(V1, left = 17, 
                        right = 18)) %>%
  mutate(grp = case_when(str_detect(ID, "WT_NltAB") ~ "NltA + NltB",
                         str_detect(ID, "E57A") ~ "E57A",
                         str_detect(ID, "E57Q") ~ "E57Q"))

colnames(df2)[2] <- "V2"
df2$grp <- factor(df2$grp, levels = c("NltA + NltB",
                                      "E57A",
                                      "E57Q"))


olea_test <- df2 %>%
  dplyr::filter(ID == unique(ID)[1]) %>%
  dplyr::select(V2, V1) %>%
  t() 

# Baseline correction
olea.irls <- baseline(olea_test)
olea_new_base <- olea.irls@corrected %>%
  t() %>%
  as.data.frame()
olea_new_time <- olea.irls@spectra %>%
  t() %>%
  as.data.frame()


un_ids <- unique(df2$ID)
ll <- list()
time <- list()

for(i in 1:length(un_ids)) {
  orig <- df2 %>%
    dplyr::filter(ID == un_ids[i]) %>%
    dplyr::select(V2, V1) %>%
    t() 
  
  orig.irls <- baseline(orig)
  orig_new_base <- orig.irls@corrected %>%
    t() %>%
    as.data.frame() %>% 
    dplyr::select(V2) 
  
    
  ll[[i]] <- orig_new_base
  time[[i]] <- orig.irls@spectra %>%
      t() %>%
      as.data.frame() %>%
      dplyr::select(V1)
}


tm_df <- bind_rows(time)

names(ll) <- un_ids
dff2 <- ll %>%
  map2_df(names(ll), ~ mutate(.x, ID = .y)) %>%
  bind_rows() %>%
  mutate(grp = case_when(str_detect(ID, "WT_NltAB") ~ "NltA + NltB",
                         str_detect(ID, "E57A") ~ "E57A",
                         str_detect(ID, "E57Q") ~ "E57Q")) %>%
  dplyr::mutate(V1 = tm_df$V1)
  
  
dff2$grp <- factor(dff2$grp, levels = c(
  "NltA + NltB",
  "E57A",
  "E57Q"))

pal3 <- c("#4B0082", "forestgreen", "orange")

pdf("output/figS5b.pdf", height = 3, width = 4)
pl2 <- ggplot(data = dff2,  aes(x = V1, y = V2, colour = grp, group = ID)) +
  theme_bw() +
  xlab("Retention time (min)") +
  ylab("GC-FID Peak Area") +
  facet_grid(rows = vars(grp)) +
  geom_point(size = 0.0000000001, alpha = 0.01) +
  geom_path(size = 0.75) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 8), axis.title = element_text(size = 11),
        legend.text = element_text(size = 10), legend.title = element_blank(),
        legend.position = "none",
        axis.ticks.x = element_line(), axis.ticks.y = element_line(),
        plot.margin = unit(c(0, 1, 0, 1), "cm"),
        panel.border = element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_blank())+
  scale_color_manual(values = pal3) +
  ylim(NA, 20000)
pl2
dev.off()

