# Install packages
pacman::p_load("caret", "Biostrings", "DECIPHER", "tidymodels", "ranger", 
               "tree", "rsample", "tidyverse", "randomForest", "nnet")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/adenylpred_analysis/")

# Read in the data
rawdat <- read_csv("data/703_training_sqs_with_loop_extracted.csv")
table(duplicated(rawdat[,2:ncol(rawdat)])) # check no duplicates
colnames(rawdat)[1] <- "nms"
rawdat$clf <- word(rawdat$nms, 2, sep = "_")

# Remove the holdout test predictions
dat <- rawdat[!grepl(paste0(c("HOLDOUT", "OTHER", "CAR", "amino.acid", "reject"), collapse = "|"), rawdat$clf),] # 658 observations
table(dat$clf)

nb_models <- list()
nb_pred <- list()
dat_test <- list()
nb_cm <- list()

for(i in 1:10) {
  
  set.seed(i)
  
  # Split into test and train
  dat_split <- initial_split(dat, prop = 3/4, strata = "clf")
  dat_train <- training(dat_split)
  dat_test  <- testing(dat_split)
  nrow(dat_train)/nrow(dat) # 75 %
  
  # Define our response
  x_train <- dat_train[,!colnames(dat_train) %in% c("nms", "clf")]
  x_test <- dat_test[,!colnames(dat_test) %in% c("nms", "clf")]
  y_train <- as.factor(dat_train$clf)
  y_test <- as.factor(dat_test$clf)
  print(table(y_test))
  
  form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat_train$nms)
  form_test <- data.frame(cbind(x_test, y_test), stringsAsFactors = F, row.names = dat_test$nms)
  
  # Train models
  fitControl <- trainControl(method = "none", classProbs = FALSE)
  nb_grid <- expand.grid(usekernel = TRUE, fL = 0, adjust = 1)
  
  nb_models[[i]] <-  nb <- train(
    x = x_train,
    y = y_train,
    method = "nb",
    tuneGrid = nb_grid,
    trControl = fitControl)
  
  nb_pred[[i]] <- predict(nb_models[[i]], newdata = x_test)
  nb_cm[[i]] <- confusionMatrix(nb_pred[[i]], as.factor(dat_test$clf))
  dat_test[[i]] <- dat_test
}

# Calculate the overall test statistics
length(nb_cm)
xx <- list()
for (i in 1:length(nb_cm)) {
  xx[[i]] <- nb_cm[[i]]$overall[1]
}
mean(unlist(xx)) # 0.820 classification accuracy
sd(unlist(xx)) # 0.019 standard deviation

