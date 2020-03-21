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
                           importance = "permutation", class.weights = classwts, probability = T)
  
  rf_pred[[i]] <- predict(rf_models[[i]], data = x_test)
  dat_test[[i]] <- dat_test
}

# Calculate the overall test statistics
ind <- 10
rf_df <- rf_pred[[ind]]$predictions
colnames(rf_df) <- paste0(colnames(rf_df), "_pred_RF")

true_label <- dummies::dummy(dat_test$clf, sep = ".")
true_label <- data.frame(true_label)
colnames(true_label) <- paste0(colnames(true_label), "_true")
true_label <- data.frame(true_label)

final_df <- data.frame(cbind(true_label, rf_df), stringsAsFactors = F)
colnames(final_df) <- gsub("clf\\.", "", colnames(final_df))

roc_res <- multi_roc(final_df, force_diag=T)
pr_res <- multi_pr(final_df, force_diag=T)
plot_roc_df <- plot_roc_data(roc_res)
plot_pr_df <- plot_pr_data(pr_res)

pal2 <- c("#E41A1C", "#92D050", "#377EB8", "#984EA3",
          "#FF7F00", "goldenrod", "#A65628",   "#F781BF",
          "blue1", "gray68", "darkorchid1", "navy", 
          "plum1","deepskyblue", "gold", 
          "deeppink2", "lightslateblue",
          "lightblue2", "darkseagreen1")

# Make AUROC plot (Figure 1C)
pdf("output/fig1c.pdf", width = 4, height = 4)
ggplot(plot_roc_df, aes(x = 1-Specificity, y=Sensitivity)) +
  geom_path(aes(color = Group), size = 1) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
               colour='grey', linetype = 'dotdash') +
  theme_bw() + 
  scale_color_manual(values = pal2) +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "none",
        legend.title = NULL)
dev.off()

# Calculate confidence intervals
roc_ci_res_micro <- roc_ci(final_df, conf= 0.95, type='basic', R = 1000, index = 17)
roc_ci_res_micro$t0
roc_ci_res_micro

roc_ci_res_macro <- roc_ci(final_df, conf= 0.95, type='basic', R = 1000, index = 16)
roc_ci_res_macro$t0


