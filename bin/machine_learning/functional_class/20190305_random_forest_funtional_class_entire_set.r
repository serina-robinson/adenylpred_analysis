## Install packages
pacman::p_load("caret", "data.table", "Biostrings", "seqinr", "DECIPHER", "cowplot", "tidymodels", "ranger", "tree", "rsample", "tidyverse", "randomForest","gbm","nnet","e1071","svmpath","lars","glmnet","svmpath")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Set seed 
set.seed(20190403)

# Read in the data
rawdat <- read_csv("data/703_training_sqs_with_loop_extracted.csv")
table(duplicated(rawdat[,2:ncol(rawdat)])) # check no duplicates
colnames(rawdat)[1] <- "nms"
rawdat$clf <- word(rawdat$nms, 2, sep = "_")
table(rawdat$clf)

# Remove rejects
dat <- rawdat[-grep(paste0(c("HOLDOUT", "OTHER", "amino.acid", "reject"), collapse = "|"), rawdat$nms),]

# Train a random forest model with optimal parameters, no max depth
x_train <- dat[,!colnames(dat) %in% c("nms", "clf")]
colnames(x_train)
y_train <- dat$clf
form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat$nms)

tunegrid <- expand.grid(.mtry = as.integer(sqrt(ncol(x_train))), .splitrule = 'gini', .min.node.size = 1)

rf_full_fc <- caret::train(
  x = x_train,
  y = y_train,
  method = "ranger",
  trControl = trainControl(method = "repeatedcv", number = 10,
                           repeats = 5,
                           verboseIter = T, classProbs = T,
                           savePredictions = "final"),
  num.trees = 1000,
  tuneGrid = tunegrid,
  verbose = TRUE,
  importance = "permutation")
rf_full_fc$results
# saveRDS(rf_full_fc, "data/20190305_rf_models/rf_fullset_10xvalidated_funct_class.rds")

# Now train on the entire dataset
rf_full_fc_noxval <- ranger(y_train ~., data = form_train, num.trees = 1000, splitrule = "gini",
                mtry = as.integer(sqrt(ncol(x_train))), min.node.size = 1,
                importance = "permutation", probability = TRUE)
rf_full_fc_noxval # OOB prediction error is 17.2%
saveRDS(rf_full_fc_noxval, "data/20190305_rf_fullset_fc_noxval.rds")
