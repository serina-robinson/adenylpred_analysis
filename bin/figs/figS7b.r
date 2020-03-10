# Install packages
pacman::p_load("MSnbase", "data.table", "InterpretMSSpectrum", "cowplot")

# Read in files
setwd("~/Documents/Wageningen_UR/github/adenylpred_analysis/")

# Try reading in the mass spectrum
ms1 <- fread("data/1-24-18 Nlt Ole Mix/1-23-18 MS 4B.D/SPECTAB.CSV", sep = ",", data.table = F,
             fill=TRUE, skip = 3)
PlotSpec(ms1, masslab=NULL, neutral_losses=NA, xlim = c(50,300))

enz <- fread("data/20180125 Lactone with BME and hydroxylamine/B 20min.D/SPECTAB.CSV", sep = ",", data.table = F,
             fill=TRUE, skip = 3)

pdf("output/figS7b.pdf", width = 9, height = 6)
  par(mfrow = c(2,1))
  par(mar = c(0,0,0,0))
  PlotSpec(ms1, masslab=NULL, neutral_losses=NA, xlim = c(50,275), cols = rep("red", nrow(ms1)))
  PlotSpec(enz, masslab=NULL, neutral_losses=NA, xlim = c(50,275), mz_prec = 1, ylim = c(16000,0), cols = rep("black", nrow(enz)))
dev.off()

