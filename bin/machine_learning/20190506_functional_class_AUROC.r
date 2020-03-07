# Install packages
pacman::p_load("caret", "data.table", "RColorBrewer", "tidyverse",
               "rsample", "ranger", "multiROC", "stringr")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/adenylpred_analysis/")

# Set random seed 
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
table(y_test)

# Complete dataset for training and testing
form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat_train$nms)
form_test <- data.frame(cbind(x_test, y_test), stringsAsFactors = F, row.names = dat_test$nms)

classwts <- nrow(dat)/(length(unique(dat$clf)) * table(dat$clf))
classwts
dtf <- data.frame(cbind(table(dat$clf), classwts))
colnames(dtf) <- c("class_size", "class_weights")

# writev(dtf, "output/small_sub_grp_class_sizes_and_weights.csv", quote = F, row.names = T)

### Try it with max depth 15 and tune grid determined
tunegrid <- expand.grid(.splitrule = "gini", .mtry = as.integer(sqrt(ncol(x_train))), .min.node.size = 1)

rf_weighted <- caret::train(
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
  max.depth = 20,
  class.weights = classwts,
  importance = "permutation")

getTrainPerf(rf_weighted)
rf_weighted$results
saveRDS(rf_weighted, "data/20193105_rf_weighted_with_loop_fc.rds")


### Try it with max depth 20 and tune grid determined
tunegrid <- expand.grid(.splitrule = "gini", .mtry = as.integer(sqrt(ncol(x_train))), .min.node.size = 1)

rf_unweighted <- caret::train(
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
  max.depth = 20,
  importance = "permutation")
getTrainPerf(rf_unweighted)
saveRDS(rf_unweighted, "data/20193105_rf_unweighted_with_loop_fc.rds")

##### Now try the exact same but modify to remove the extra amino acids from the loop sequence and compare performance
x_train <- x_train[,-grep("^X", colnames(x_train)),]
y_train <- y_train[,-grep("^X", colnames(y_train)),]

# Weighted
rf_weighted_no_loop <- caret::train(
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
  max.depth = 20,
  class.weights = classwts,
  importance = "permutation")
getTrainPerf(rf_weighted_no_loop)
# saveRDS(rf_weighted_no_loop, "data/20193105_rf_weighted_no_loop_fc.rds")



# Unweighted
rf_unweighted_no_loop <- caret::train(
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
  max.depth = 20,
  importance = "permutation")
getTrainPerf(rf_unweighted_no_loop)
# saveRDS(rf_unweighted_no_loop, "data/20193105_rf_unweighted_no_loop_fc.rds")



rf_fc_weighted_no_loop <- readRDS("data/20193105_rf_weighted_no_loop_fc.rds")
rf_fc_weighted_loop <- readRDS("data/20193105_rf_weighted_with_loop_fc.rds")
rf_fc_weighted_loop$results

dim(form_train)
rf_weighted_loop <- ranger(y_train ~., data = form_train, num.trees = 1000, splitrule = "gini",
                              mtry = as.integer(sqrt(ncol(x_train))), min.node.size = 1,
                              importance = "permutation", probability = T)
rf <- rf_weighted_loop


# ROC curve
approx_roc_curve <- function(x, label) {
  x %>%
    pluck("pred") %>%
    roc_curve(obs, y_train) %>%
    mutate(model = label)
}
# ROC curve is not looking good...
approx_roc_curve(rf, "Random Forest") %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_path()  +
  geom_abline(col = "red", alpha = .5)

dim(form_test)
rf

rf_pred <- predict(rf, data = form_test)
rf_pred$predictions
# rf_pred
rf_cm <- confusionMatrix(rf_pred, as.factor(dat_test$clf))
# rf_cm$overall
# rf_cm$byClass

rf_pred <- data.frame(rf_pred$predictions)
colnames(rf_pred) <- paste(colnames(rf_pred), "_pred_RF")
head(rf_pred)
true_label <- dummies::dummy(dat_test$clf, sep = ".")
true_label <- data.frame(true_label)
colnames(true_label) <- gsub(".*?\\.", "", colnames(true_label))
colnames(true_label) <- paste(colnames(true_label), "_true")
head(true_label)
true_label <- data.frame(true_label)

final_df <- cbind(true_label, rf_pred)
head(final_df)
colnames(final_df) <- gsub("\\.", " ", colnames(final_df))
head(final_df)
# roc_res <- multiclass.roc(final_df)
head(final_df)
write.csv(final_df, quote = F, "output/predictions_for_multiROC_analysis.csv")
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
pal2

pdf("output/functional_class_multiclass_AUROC_with_legend.pdf", width = 3, height = 3)
ggplot(plot_roc_df, aes(x = 1-Specificity, y=Sensitivity)) +
  geom_path(aes(color = Group), size = 1) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
               colour='grey', linetype = 'dotdash') +
  theme_bw() + 
  scale_color_manual(values = pal2[1:9]) +
  theme(plot.title = element_text(hjust = 0.5), 
        # legend.position=c(.95, .05),
        legend.position = "right",
        legend.title = NULL)
        # axis.text = element_text(size = 12),
       # legend.background = element_rect(fill=NULL, size=0.5, 
        #                                 linetype="solid", colour ="black"))

dev.off()

unlist(roc_res$AUC)
head(final_df)
roc_ci_res_micro <- roc_ci(final_df, conf= 0.95, type='basic', R = 100, index = 11)
roc_ci_res_micro$t0

roc_ci_res_macro <- roc_ci(final_df, conf= 0.95, type='basic', R = 100, index = 1)
roc_ci_res_macro$t0


roc_auc_with_ci_res <- roc_auc_with_ci(final_df, conf= 0.95, type='basic', R = 100)

ggplot(plot_pr_df, aes(x=Recall, y=Precision)) + 
  geom_path(aes(color = Group, linetype=Method), size=1.5) + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        legend.justification=c(1, 0), legend.position=c(.95, .05),
        legend.title=element_blank(), 
        legend.background = element_rect(fill=NULL, size=0.5, 
                                         linetype="solid", colour ="black"))

