## Install packages
pacman::p_load("caret", "Biostrings", "phangorn", "ape", "seqinr",
"DECIPHER", "cowplot", "tidymodels", "ranger", "tree", "rsample", 
"tidyverse", "randomForest","gbm","nnet","e1071","svmpath","lars","glmnet","svmpath")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/adenylpred_analysis/")

# Set random seed 
set.seed(20191304)

# Read in the data
rawdat <- read_csv("data/797_seqs_510_feats_for_supervised_centered_scaled_20191104.csv")
colnames(rawdat)[1] <- "nms"
rawdat$clf <- word(rawdat$nms, -1, sep = "_")
table(rawdat$clf)

# Remove the holdout test predictions
dat <- rawdat[-grep(paste0(c("HOLDOUT", "OTHER", "CAR", "amino.acid"), collapse = "|"), rawdat$nms),] # 658 observations
table(dat$clf)

dat_split <- initial_split(dat, strata = "clf")
dat_train <- training(dat_split)
dat_test  <- testing(dat_split)
nrow(dat_train)/nrow(dat) # 75 %

# Define our response
x_train <- dat_train[,!colnames(dat_train) %in% c("nms", "clf")]
x_test <- dat_test[,!colnames(dat_test) %in% c("nms", "clf")]
y_train <- as.factor(dat_train$clf)
y_test <- as.factor(dat_test$clf)
table(y_test)

# Make a data frame for prediction
df_train <- data.frame(x_train, stringsAsFactors = F, row.names = dat_train$nms)

# Complete dataset for training and testing
form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat_train$nms)
form_test <- data.frame(cbind(x_test, y_test), stringsAsFactors = F, row.names = dat_test$nms)

# Try without weights
tunegrid <- expand.grid(.splitrule = "gini",
.mtry = as.integer(sqrt(ncol(x_train))), .min.node.size = 1)
tunegrid

rf <- caret::train(
  x = x_train,
  y = y_train,
  method = "ranger",
  trControl = trainControl(method = "repeatedcv", number = 10,
                           repeats = 3,
                           verboseIter = T, classProbs = T,
                           savePredictions = "final"),
  num.trees = 500,
  tuneGrid = tunegrid,
  verbose = TRUE,
  max.depth = 9,
  importance = "permutation")
getTrainPerf(rf)

# Try with weights
classwts <- nrow(dat)/(length(unique(dat$clf)) * table(dat$clf))
classwts

rf_weighted <- caret::train(
  x = x_train,
  y = y_train,
  method = "ranger",
  trControl = trainControl(method = "repeatedcv", number = 10,
                           repeats = 3,
                           verboseIter = T, classProbs = T,
                           savePredictions = "final"),
  num.trees = 500,
  tuneGrid = tunegrid,
  verbose = TRUE,
  max.depth = 9,
  class.weights = classwts,
  importance = "permutation")
getTrainPerf(rf_weighted)

modellist <- list()
modellist[[1]] <- rf
modellist[[2]] <- rf_weighted
names(modellist) <- c("Unweighted", "Weighted")
results <- resamples(modellist)
summary(results)
dotplot(results)


# Now compare performance on the test set

# Rf unweighted
rf_nowts <- ranger(y_train ~., data = form_train, num.trees = 1000, splitrule = "gini",
                        mtry = as.integer(sqrt(ncol(x_train))), min.node.size = 1,
                        importance = "permutation", max.depth = 9, probability = TRUE)
saveRDS(rf_nowts, "output/rf_functional_class_1000trees_probability_centered_scaled_unweighted.rds")
rf_pred <- predict(rf_nowts, data = form_test)
length(rf_pred)
cm_rf_nowts <- confusionMatrix(rf_pred$predictions, as.factor(dat_test$clf))

cm_rf_nowts$overall

# Rf weighted
classwts <- nrow(dat)/(length(unique(dat$clf)) * table(dat$clf))
rf_weighted <- ranger(y_train ~., data = form_train, num.trees = 1000, splitrule = "gini",
                   mtry = as.integer(sqrt(ncol(x_train))), min.node.size = 1,
                   importance = "permutation", max.depth = 9, class.weights = classwts)
rf_pred_wts <- predict(rf_weighted, data = form_test)
cm_rf_wtd <- confusionMatrix(rf_pred_wts$predictions, as.factor(dat_test$clf))


tow <- data.frame(rownames(cm_rf_nowts$byClass),
            cm_rf_nowts$byClass[,"F1"], # F1 score for by class
           cm_rf_nowts$byClass[,"Balanced Accuracy"],
           cm_rf_wtd$byClass[,"F1"],
           cm_rf_wtd$byClass[,"Balanced Accuracy"])

colnames(tow) <- c("Classname", "No_weights_F1", "No_weights_accuracy", "Weighted_F1", "Weighted_accuracy")
write_csv(tow, "output/effects_of_class_balancing.csv")

# RF weighted probability
rf_weighted <- ranger(y_train ~., data = form_train, num.trees = 1000, splitrule = "gini",
                      mtry = as.integer(sqrt(ncol(x_train))), min.node.size = 1,
                      importance = "permutation", max.depth = 9, class.weights = classwts,
                      probability = TRUE)
saveRDS(rf_weighted, "data/rf_functional_class_1000trees_probability_centered_scaled_weighted.rds")
