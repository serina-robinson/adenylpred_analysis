# Install pacakges
pacman::p_load("caret", "data.table", "RColorBrewer", "Biostrings", 
               "multiROC", "tidymodels", "ggpubr", 
               "ranger", "tree", "rsample", "tidyverse")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/adenylpred_analysis/")

# Set seed 
set.seed(20193105)

# Read in the data
rawdat <- read_csv("data/703_training_sqs_with_loop_extracted.csv")
table(duplicated(rawdat[,2:ncol(rawdat)])) # check no duplicates
colnames(rawdat)[1] <- "nms"
rawdat$clf <- word(rawdat$nms, 2, sep = "_")

# Remove the holdout test predictions
dat <- rawdat[!grepl(paste0(c("HOLDOUT", "OTHER", "CAR", "amino.acid", "reject"), collapse = "|"), rawdat$clf),] # 658 observations
table(dat$clf)

# Split into test and train
dat_split <- initial_split(dat, prop = 3/4, strata = "clf")
dat_train <- training(dat_split)
dat_test  <- testing(dat_split)
nrow(dat_train)/nrow(dat) # 75 %
table(word(dat$nms, sep = "_", 2))

# Define our response
x_train <- dat_train[,!colnames(dat_train) %in% c("nms", "clf")]
x_test <- dat_test[,!colnames(dat_test) %in% c("nms", "clf")]
y_train <- as.factor(dat_train$clf)
y_test <- as.factor(dat_test$clf)

# Complete dataset for training and testing
form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat_train$nms)
form_test <- data.frame(cbind(x_test, y_test), stringsAsFactors = F, row.names = dat_test$nms)

classwts <- nrow(dat)/(length(unique(dat$clf)) * table(dat$clf))
dtf <- data.frame(cbind(table(dat$clf), classwts))
colnames(dtf) <- c("class_size", "class_weights")

# Train a model
rf <- ranger(y_train ~., data = form_train, num.trees = 1000, splitrule = "gini",
                           mtry = as.integer(sqrt(ncol(x_train))), min.node.size = 1,
                           class.weights = classwts, importance = "permutation", probability = T)

rf_pred1 <- predict(rf, data = form_test)

# Train model for confusion matrix
rf2 <- ranger(y_train ~., data = form_train, num.trees = 1000, splitrule = "gini",
             mtry = as.integer(sqrt(ncol(x_train))), min.node.size = 1,
             class.weights = classwts, importance = "permutation")

rf_pred2 <- predict(rf2, data = form_test)

# Make a confusion matrix
rf_cm <- confusionMatrix(rf_pred2$predictions, as.factor(dat_test$clf))
dtf <- data.frame(rf_cm$table, stringsAsFactors = F)

dtf$Prediction <- as.character(dtf$Prediction)
dtf$Reference <- as.character(dtf$Reference)

dtf$Prediction[dtf$Prediction == "VLACSBILE"] <- "VLACS"
dtf$Reference[dtf$Reference == "VLACSBILE"] <- "VLACS"
dtf$Prediction[dtf$Prediction == "LUCIFERASE"] <- "LUC"
dtf$Reference[dtf$Reference == "LUCIFERASE"] <- "LUC"

# Plot a confusion matrix for the best result
pdf("output/figS2a.pdf", width = 4, height = 4)
ggplot(dtf, aes(Prediction, Reference)) +
  geom_tile(aes(fill = Freq)) + 
  geom_text(aes(label = round(Freq, 1))) +
  scale_fill_gradient(low = "white", high = "red") +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 90, color = "black"), 
        axis.text.y = element_text(color = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") 
dev.off()

rf_pred <- data.frame(rf_pred1$predictions)
rf_pred
colnames(rf_pred) <- paste(colnames(rf_pred), "_pred_RF")
colnames(rf_pred)
true_label <- dummies::dummy(dat_test$clf, sep = ".")
true_label <- data.frame(true_label)
colnames(true_label) <- gsub(".*?\\.", "", colnames(true_label))
colnames(true_label) <- paste(colnames(true_label), "_true")
true_label <- data.frame(true_label)

final_df <- cbind(true_label, rf_pred)
final_df
colnames(final_df) <- gsub("\\.", " ", colnames(final_df))
colnames(final_df)

roc_res <- multi_roc(final_df, force_diag = T)
pr_res <- multi_pr(final_df, force_diag = T)
plot_roc_df <- plot_roc_data(roc_res)

pal2 <- c("#E41A1C", "#92D050", "#377EB8", "#984EA3",
                  "#FF7F00", "goldenrod", "#A65628",   "#F781BF",
                  "blue1", "gray68", "black", "navy", 
                  "plum1","deepskyblue", "gold", 
                  "deeppink2", "lightslateblue",
                  "lightblue2", "darkseagreen1", "darkorchid1")


pdf("output/figS2b.pdf", width = 4, height = 4)
ggplot(plot_roc_df, aes(x = 1-Specificity, y=Sensitivity)) +
  geom_path(aes(color = Group), size = 1) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
               colour='grey', linetype = 'dotdash') +
  theme_pubr() + 
  scale_color_manual(values = pal2) +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "none",
        legend.title = NULL)
dev.off()

pdf("output/figS2b_with_legend.pdf", width = 5, height = 5)
ggplot(plot_roc_df, aes(x = 1-Specificity, y=Sensitivity)) +
  geom_path(aes(color = Group), size = 1) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
               colour='grey', linetype = 'dotdash') +
  theme_pubr() + 
  scale_color_manual(values = pal2) +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "right",
        legend.title = NULL)
dev.off()

# Micro
roc_ci_res_micro <- roc_ci(final_df, conf= 0.95, type='basic', R = 1000, index = 11)
roc_ci_res_micro$t0
roc_ci_res_micro

# Macro
roc_ci_res_macro <- roc_ci(final_df, conf= 0.95, type='basic', R = 1000, index = 10)
roc_ci_res_macro$t0
