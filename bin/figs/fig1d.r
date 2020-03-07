# Install packages
pacman::p_load("ggplot2", "multiROC", "ranger", "dummies", "tidymodels", "caret", "tidyverse", "ggplot2")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/adenylpred_analysis/")

# Read in the data
rawdat <- read_csv("data/1553_training_sqs_with_loop_extracted.csv")
colnames(rawdat)[1] <- "nms"
rawdat$clf <- word(rawdat$nms, -1, sep = "_")

# Remove the holdout test predictions
dat <- rawdat[!grepl(paste0(c("HOLDOUT", "OTHER", "CAR", "amino.acid", "reject"), collapse = "|"), rawdat$clf),] # 658 observations

rf_models <- list()
rf_pred <- list()
dat_test <- list()
rf_cm <- list()

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
  
  # Complete dataset for training and testing
  form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat_train$nms)
  form_test <- data.frame(cbind(x_test, y_test), stringsAsFactors = F, row.names = dat_test$nms)
  
  classwts <- nrow(dat)/(length(unique(dat$clf)) * table(dat$clf))
  dtf <- data.frame(cbind(table(dat$clf), classwts))
  colnames(dtf) <- c("class_size", "class_weights")
  
  rf_models[[i]] <- ranger(y_train ~., data = form_train, num.trees = 1000, splitrule = "gini",
                           mtry = as.integer(sqrt(ncol(x_train))), min.node.size = 1,
                           importance = "permutation", class.weights = classwts)
  
  rf_pred[[i]] <- predict(rf_models[[i]], data = x_test)
  rf_cm[[i]] <- confusionMatrix(rf_pred[[i]]$predictions, as.factor(dat_test$clf))
  dat_test[[i]] <- dat_test
}
as.factor(dat_test$clf)
colnames(rf_pred[[i]]$predictions)[which.max(rf_pred[[i]]$predictions)]

# Calculate average test set error
xx <- list()
for (i in 1:length(rf_cm)) {
  xx[[i]] <- rf_cm[[i]]$overall[1]
}
mean(unlist(xx))  # 0.825 classification accuracy
sd(unlist(xx))  # 0.023 standard deviation
ind <- which.max(unlist(xx)) # ind = 4
dtf <- data.frame(rf_cm[[ind]]$table)

# Plot a confusion matrix
pdf("output/fig1d.pdf", width = 8, height = 8)
ggplot(dtf, aes(Prediction, Reference)) +
  geom_tile(aes(fill = Freq)) + 
  geom_text(aes(label = round(Freq, 1))) +
  scale_fill_gradient(low = "white", high = "red") +
  theme(axis.text.x = element_text(angle = 90), 
        axis.title.x = element_blank(),
        legend.position = "none")
dev.off()

