# Install packages
pacman::p_load("caret", "Biostrings", "RColorBrewer", 
                "tidymodels", "ranger", "DECIPHER",
                "rsample", "tidyverse", "multiROC")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/adenylpred_analysis/")

# Read in the data
rawdat <- read_csv("data/703_training_sqs_with_loop_extracted.csv")
table(duplicated(rawdat[,2:ncol(rawdat)])) # check no duplicates
colnames(rawdat)[1] <- "nms"
rawdat$clf <- word(rawdat$nms, 2, sep = "_")
dat <- rawdat[!grepl(paste0(c("HOLDOUT", "OTHER", "CAR"), collapse = "|"), rawdat$clf),] 

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

# xx <- list()
# for (i in 1:length(rf_cm)) {
#   xx[[i]] <- rf_cm[[i]]$overall[1]
# }
# mean(unlist(xx))
# sd(unlist(xx))
# 
# ind <- which.max(unlist(xx))

# rf_pred <- data.frame(rf_pred[[4]]$predictions)
rf_df <- rf_pred[[ind]]$predictions
colnames(rf_df) <- paste0(colnames(rf_df), "_pred_RF")
head(rf_df)

true_label <- dummies::dummy(dat_test$clf, sep = ".")
true_label <- data.frame(true_label)
# colnames(true_label) <- gsub(".*?\\.", "", colnames(true_label))
colnames(true_label) <- paste0(colnames(true_label), "_true")
head(true_label)
true_label <- data.frame(true_label)

final_df <- cbind(true_label, rf_df)
colnames(final_df) <- gsub("clf\\.", "", colnames(final_df))

roc_res <- multi_roc(final_df, force_diag=T)
pr_res <- multi_pr(final_df, force_diag=T)
# colnames(rf_pred) <- paste(colnames(rf_pred), "_pred_RF")

plot_roc_df <- plot_roc_data(roc_res) %>%
  dplyr::filter(!Group %in% c("Micro", "Macro"))

table(plot_roc_df$Group)
plot_pr_df <- plot_pr_data(pr_res)

pal1 <- colorRampPalette(colors=brewer.pal(8, "Set1"))(8)
pal1[pal1 == "#377EB8"] <- "#92D050"
# pal1[7] <- pal1[8]
# pal1[8] <- "#A65628"
# pal1[pal1 == "#A65628"] <- "gray68"
# pal1[pal1 == "#F781BF"] <- "#A65628"
pal1[pal1 == "#4DAF4A"] <- "#377EB8"
pal1[pal1 == "#FFFF33"] <- "goldenrod"
pal2 <- c(pal1, "blue1", "gray68", "darkorchid1", "navy", #"black", 
          "plum1",
          "deepskyblue", "gold", 
          "deeppink2", "lightslateblue",
          "lightblue2", "darkseagreen1")
palette(pal2)
length(pal2)

pdf("output/figS2.pdf", width = 4, height = 4)
ggplot(plot_roc_df, aes(x = 1-Specificity, y=Sensitivity)) +
  geom_path(aes(color = Group), size = 1) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
               colour='grey', linetype = 'dotdash') +
  theme_bw() + 
  scale_color_manual(values = pal2) +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.justification=c(1, 0), legend.position=c(.95, .05),
        legend.title = NULL)
dev.off()

unlist(roc_res$AUC)
roc_ci_res_micro <- roc_ci(final_df, conf= 0.95, type='basic', R = 100, index = 17)
roc_ci_res_micro

roc_ci_res_macro <- roc_ci(final_df, conf= 0.95, type='basic', R = 100, index = 1)

roc_auc_with_ci_res <- roc_auc_with_ci(final_df, conf= 0.95, type='basic', R = 100)

ggplot(plot_pr_df, aes(x=Recall, y=Precision)) + 
  geom_path(aes(color = Group, linetype=Method), size=1.5) + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        legend.justification=c(1, 0), legend.position=c(.95, .05),
        legend.title=element_blank(), 
        legend.background = element_rect(fill=NULL, size=0.5, 
                                         linetype="solid", colour ="black"))

