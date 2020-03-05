# Install packages
pacman::p_load("data.table", "ranger", "muscle", "tidyverse", 
               "Biostrings","DECIPHER","igraph","RColorBrewer")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/adenylpred_analysis/")

# Now try with the 63,395 sequences (unclustered)
antismash_unclustered <- readAAStringSet("data/antismashv2/63395_antismashDBv2_AMP_binding_hits.fa")
length(antismash_unclustered)
# Remove duplicates
antismash_dedup <- antismash_unclustered[!duplicated(antismash_unclustered)]
length(antismash_dedup)
ex34 <- readAAStringSet("output/antismashv2_unclustered_34aa_extracted.faa")
sqs <- antismash_dedup[!names(antismash_dedup) %in% names(ex34)]
length(sqs)

# Extract 34 active site aa
source("src/extract_34_aa.r")
query_fils <- sapply(1:length(sqs), function(x) {tempfile(pattern = "", fileext = ".fasta")})
sapply(1:length(sqs), function(x) {writeXStringSet(sqs[x], query_fils[x])})
extract_34_list <- lapply(1:length(sqs), function(x) { extract_34_aa(query_fils[x]) })

nlist <- list()
grp1 <- readAAStringSet('output/last6251_antismashv2_unclustered_34aa_extracted.faa')
grp2 <- readAAStringSet('output/antismashv2_unclustered_34aa_extracted.faa')
grp_comb <- c(grp1, grp2)
names(grp_comb)

res <- strsplit(as.character(grp_comb), "")
nlist <- lapply(1:length(grp_comb), function(x) convert_seq_15aap(as.character(res[[x]])))
head(nlist)
head(nlist)
extract_34_df <- data.frame(matrix(unlist(nlist), nrow = length(res), byrow=T), stringsAsFactors=FALSE)

# ex34 <- readAAStringSet("output/antismashv2_unclustered_34aa_extracted.faa")

colnames(extract_34_df) <- as.character(fread("data/feature_names.txt", data.table = F)[,1])
dim(extract_34_df) # 2344 predictions by 510
rownames(extract_34_df) <- names(grp_comb)
# write_csv(data.frame(names(grp_comb)), "output/48250_antismashv2names.csv")
# write_csv(extract_34_df, "output/48250_antismashv2seqs_510_feats.csv")

# Read in the antismash training set
extract_34_df <- read_csv("output/48250_antismashv2seqs_510_feats.csv") 

# rownames(extract_34_df) <- names(sqs)
head(extract_34_df)
dim(extract_34_df)
colnames(extract_34_df)
# Predict functional class
rf_funct_class <- readRDS("data/rf_functional_class_1000trees_probability_centered_scaled.rds")
rf_funct_class_pred <- predict(rf_funct_class, data = extract_34_df, predict.all = F)
str(rf_funct_class_pred)
rf_fc_pred <- rf_funct_class_pred$predictions

# rf_fc_pred

res_fc_prob <- apply(rf_fc_pred, 1, max)
res_fc <- colnames(rf_fc_pred)[apply(rf_fc_pred, 1, which.max)]
full_prediction_results <- data.frame(names(grp_comb), res_fc, res_fc_prob)
colnames(full_prediction_results) <- c("names", "prediction", "probability")
fpr <- full_prediction_results
fpr$prob0.6 <- case_when(fpr$probability > 0.6 ~ paste0(fpr$names, "_", fpr$prediction),
                         TRUE ~ paste0(fpr$names, "_NOPRED"))
fpr$prob0.5 <- case_when(fpr$probability > 0.5 ~ paste0(fpr$names, "_", fpr$prediction),
                         TRUE ~ paste0(fpr$names, "_NOPRED"))
fpr$prob0.4 <- case_when(fpr$probability > 0.4 ~ paste0(fpr$names, "_", fpr$prediction),
                         TRUE ~ paste0(fpr$names, "_NOPRED"))
fpr$prob0.3 <- case_when(fpr$probability > 0.3 ~ paste0(fpr$names, "_", fpr$prediction),
                         TRUE ~ paste0(fpr$names, "_NOPRED"))
fpr$prob0.2 <- case_when(fpr$probability > 0.2 ~ paste0(fpr$names, "_", fpr$prediction),
                         TRUE ~ paste0(fpr$names, "_NOPRED"))
fpr$prob0.1 <- case_when(fpr$probability > 0.1 ~ paste0(fpr$names, "_", fpr$prediction),
                         TRUE ~ paste0(fpr$names, "_NOPRED"))

# write_csv(fpr, "output/48250_antismashv2_full_rf_prediction_results.csv")

# write_csv(full_prediction_results, "output/full_fc_prediction_results_old_train.csv")
preddf_nowts <- data.frame(rbind(table(full_prediction_results$prediction[full_prediction_results$probability > 0.6]),
                     table(full_prediction_results$prediction[full_prediction_results$probability > 0.5]),
                     table(full_prediction_results$prediction[full_prediction_results$probability > 0.4]),
                     table(full_prediction_results$prediction[full_prediction_results$probability > 0.3])))
table(full_prediction_results$prediction[full_prediction_results$probability > 0.5])
preddf_to_write <- cbind(c(0.6, 0.5, 0.4, 0.3), preddf_nowts)
# write_csv(preddf_to_write, "output/48250_antismashv2_predictions.csv")
full_prediction_results$prediction[full_prediction_results$probability > 0.5,]

nrow(fpr[fpr$prediction == "BLS",])
# write.csv(fpr[fpr$prediction == "BLS",], "output/166_BLS_predictions_from_antiSMASH.csv")
bls <- read_csv("output/166_BLS_predictions_from_antiSMASH.csv")
bls$names[grep("mibig", bls$names)]

925/2344 # about 37% confidence
table(high_conf$prediction)

high_conf <- full_prediction_results[full_prediction_results$probability > 0.4,]
dim(high_conf)
1386/2344 # 52%


high_conf <- full_prediction_results[full_prediction_results$probability > 0.3,]

dim(high_conf)


high_conf <- full_prediction_results[full_prediction_results$probability > 0.6,]
dim(high_conf)
675/2344

all_conf <- full_prediction_results
table(all_conf$prediction)
length(all_conf$prediction)
all_conf$prediction[all_conf$probability < 0.5] <- "NOPRED"
# write_csv(all_conf[all_conf$prediction == "BLS",], "output/BLS_predictions.csv")

sqnams_for_coloring <- paste0(names(grp_cmb), " ", as.character(all_conf$prediction))
write_csv(data.frame(sqnams_for_coloring), "output/antismashv2_predictions_for_clustering.csv")

# Try it with the weighted classes
rf_funct_class <- readRDS("data/rf_functional_class_1000trees_probability_centered_scaled_weighted.rds")
rf_funct_class_pred <- predict(rf_funct_class, data = extract_34_df, predict.all = F)
rf_fc_pred <- rf_funct_class_pred$predictions

res_fc_prob <- apply(rf_fc_pred, 1, max)
res_fc <- colnames(rf_fc_pred)[apply(rf_fc_pred, 1, which.max)]
full_prediction_results <- data.frame(names(grp_comb), res_fc, res_fc_prob)
colnames(full_prediction_results) <- c("names", "prediction", "probability")


summary(full_prediction_results$probability)
# write_csv(full_prediction_results, "output/full_fc_prediction_results_old_train.csv")

high_conf <- full_prediction_results[full_prediction_results$probability > 0.5,]
dim(high_conf)
966/2344 # about 37% confidence
table(high_conf$prediction)

high_conf <- full_prediction_results[full_prediction_results$probability > 0.4,]
dim(high_conf)
1559/2344 # 52%
table(high_conf$prediction)

high_conf <- full_prediction_results[full_prediction_results$probability > 0.3,]
table(high_conf$prediction)

preddf_wtd <- data.frame(rbind(table(full_prediction_results$prediction[full_prediction_results$probability > 0.6]),
                                     table(full_prediction_results$prediction[full_prediction_results$probability > 0.5]),
                                     table(full_prediction_results$prediction[full_prediction_results$probability > 0.4]),
                                     table(full_prediction_results$prediction[full_prediction_results$probability > 0.3])))

tbpred <- data.frame(cbind(preddf_nowts, preddf_wtd))
head(tbpred)



# 804 minutes or...13 hours... whoa
source("src/extract_34_aa.r")
query_fils <- sapply(1:length(sqs), function(x) {tempfile(pattern = "", fileext = ".fasta")})
sapply(1:length(sqs), function(x) {writeXStringSet(sqs[x], query_fils[x])})
extract_34_list <- lapply(1:length(sqs), function(x) { extract_34_aa(query_fils[x]) })
extract_34_df <- data.frame(matrix(unlist(extract_34_list), nrow = length(sqs), byrow=T),
                            stringsAsFactors=FALSE)
colnames(extract_34_df) <- as.character(fread("data/feature_names.txt", data.table = F)[,1])
dim(extract_34_df)
rf_full_pred <- predict(antismash_unclustered, data = extract_34_df, predict.all = F)