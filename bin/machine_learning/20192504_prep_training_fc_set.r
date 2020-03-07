# Install packages
pacman::p_load("readxl", "tidyverse", "DECIPHER", "data.table")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/adenylpred_analysis/")

# Set seed 
set.seed(20190304)

# Read in the full length sequences
rawdat <- read_excel("data/combined_adenylate_forming_training_set_for_db_20191404.xlsx") %>%
  dplyr::filter(!functional_class %in% c("OTHER", "CAR")) %>%
  dplyr::filter(small_substrate_group != "unknown.other") %>%
  dplyr::mutate(org_clned = paste0(word(organism, 1, sep = " "), "_", word(organism, 2, sep = " "))) %>%
  dplyr::mutate(acc_clned = gsub("_", "", acc))

sqs <- AAStringSet(rawdat$aa_seq)
sqnams_fc <- gsub(" ", "_", paste0(rawdat$acc_clned, "_", rawdat$functional_class, "_", rawdat$org_clned, "_", rawdat$small_substrate_group, "_", rawdat$small_substrate_group))
sqnams_ss <- gsub(" ", "_", paste0(rawdat$acc_clned, "_",  rawdat$small_substrate_group, "_", rawdat$org_clned, "_", rawdat$small_substrate_group, "_", rawdat$functional_class))
names(sqs) <- sqnams_fc # fix sqnams

# Extract the 34 amino acids AND loop 
source("lib/extract_34_aa_loop.r")
query_fils <- sapply(1:length(sqs), function(x) {tempfile(pattern = "", fileext = ".fasta")})
sapply(1:length(sqs), function(x) {writeXStringSet(sqs[x], query_fils[x])})
extract_34_list <- lapply(1:length(sqs), function(x) { extract_34_aa_loop(query_fils[x]) })
extract_34_df <- data.frame(matrix(unlist(extract_34_list), nrow = length(extract_34_list), byrow=T), stringsAsFactors=FALSE)
colnames(extract_34_df)[1:510] <- as.character(fread("data/feature_names.txt", data.table = F)[,1])
rownames(extract_34_df) <- names(sqs)
write.csv(extract_34_df, "data/655_uniprot_train_loop_extracted.csv", row.names = T, quote = F)

# Combine with the NRPS sequences
aa_dat <- readAAStringSet("data/sp2.adomains.faa")
aa_grps <- readAAStringSet("data/sp2_34extract_names_fixed_large_grps.faa")
names(aa_dat) <- names(aa_grps)

# Sample 150 sequences at random to prevent proportional bias from NRPS A-domains
ran_nums <- sample(1:length(aa_dat), 150, replace = F)
aa_samp <- aa_dat[ran_nums]
sqs <- aa_samp
names(sqs) <- paste0("NRPS_NRPS_", names(sqs))

# Extract the 34 amino acids AND loop 
source("lib/extract_34_aa_loop.r")
query_fils <- sapply(1:length(sqs), function(x) {tempfile(pattern = "", fileext = ".fasta")})
sapply(1:length(sqs), function(x) {writeXStringSet(sqs[x], query_fils[x])})
extract_34_list <- lapply(1:length(sqs), function(x) { extract_34_aa_loop(query_fils[x]) })
extract_34_df_NRPS <- data.frame(matrix(unlist(extract_34_list), nrow = length(extract_34_list), byrow=T), stringsAsFactors=FALSE)
colnames(extract_34_df_NRPS)
head(extract_34_df)
colnames(extract_34_df_NRPS)[1:510] <- as.character(fread("data/feature_names.txt", data.table = F)[,1])
rownames(extract_34_df_NRPS) <- names(sqs)
write.csv(extract_34_df_NRPS, "data/150_NRPS_train_loop_extracted_new.csv", quote = F, row.names = T)

# Combine the NRPS and the Uniprot
comb_df <- read_csv("data/655_uniprot_train_loop_extracted.csv") %>%
  column_to_rownames("X1")

comb_df2 <- read_csv("data/150_NRPS_train_loop_extracted_new.csv") %>%
  column_to_rownames("X1")

total_dat <- data.frame(rbind(comb_df, comb_df2), stringsAsFactors = F)

table(duplicated(total_dat)) # 105 duplicates
total_dedup <- total_dat[!duplicated(total_dat),]
dim(total_dedup)
table(word(rownames(total_dedup), 2, sep = "_"))
write.csv(total_dedup, "data/703_training_sqs_with_loop_extracted.csv", quote = F, row.names = T) 

# 