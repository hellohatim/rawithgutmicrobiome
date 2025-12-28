setwd("C:/Users/superadmin/OneDrive/Documents/Visa/EB1A/PUBLICATIONS/IEEE/RA_Study")

library(dplyr)
library(tibble)

# Load RDS files
asv_profile <- readRDS("1.ASV.profile.rds")
tax_info    <- readRDS("1.taxonomy.info.rds")

# Convert rownames (ASVs) to explicit column
asv_profile_df <- asv_profile %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV")

tax_info_df <- tax_info %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV")

# Identify sample columns
sample_cols <- setdiff(colnames(asv_profile_df), "ASV")

# Adjust patterns if needed
hc_cols <- grep("^HC", sample_cols, value = TRUE)
ra_cols <- grep("^RA", sample_cols, value = TRUE)

cat("HC samples:", length(hc_cols), "\n")
cat("RA samples:", length(ra_cols), "\n")

# Identify ASVs present in HC (>0 in any HC sample)
hc_asv_ids <- asv_profile_df %>%
  select(ASV, all_of(hc_cols)) %>%
  mutate(any_HC = rowSums(across(all_of(hc_cols)) > 0, na.rm = TRUE) > 0) %>%
  filter(any_HC) %>%
  pull(ASV)

# Identify ASVs present in RA (>0 in any RA sample)
ra_asv_ids <- asv_profile_df %>%
  select(ASV, all_of(ra_cols)) %>%
  mutate(any_RA = rowSums(across(all_of(ra_cols)) > 0, na.rm = TRUE) > 0) %>%
  filter(any_RA) %>%
  pull(ASV)

# Create HC and RA taxonomy-only tables
HC_tax_info <- tax_info_df %>%
  filter(ASV %in% hc_asv_ids)

RA_tax_info <- tax_info_df %>%
  filter(ASV %in% ra_asv_ids)

# Optional: also create taxonomy + abundance tables
HC_tax_abund <- HC_tax_info %>%
  left_join(asv_profile_df %>% select(ASV, all_of(hc_cols)), by = "ASV")

RA_tax_abund <- RA_tax_info %>%
  left_join(asv_profile_df %>% select(ASV, all_of(ra_cols)), by = "ASV")

library(stringr)

# Extract genus from taxonomy string
HC_genera <- HC_tax_info %>%
  mutate(Genus = str_extract(Taxon, "g__[^;]+")) %>%
  distinct(Genus)

RA_genera <- RA_tax_info %>%
  mutate(Genus = str_extract(Taxon, "g__[^;]+")) %>%
  distinct(Genus)

# Shared between HC and RA
shared_genera <- intersect(HC_genera$Genus, RA_genera$Genus)

# Only in Healthy Controls
hc_only_genera <- setdiff(HC_genera$Genus, RA_genera$Genus)

# Only in RA Patients
ra_only_genera <- setdiff(RA_genera$Genus, HC_genera$Genus)
################################################################################

#calculate 5% thresholds automatically
hc_thresh <- ceiling(0.05 * length(hc_cols))  # ~61
ra_thresh <- ceiling(0.05 * length(ra_cols))  # ~52

hc_thresh
ra_thresh

#find ASV with >=5% prevalence in HC
hc_asv_ids_5pct <- asv_profile_df %>%
  select(ASV, all_of(hc_cols)) %>%
  mutate(
    hc_prevalence = rowSums(across(all_of(hc_cols)) > 0, na.rm = TRUE)
  ) %>%
  filter(hc_prevalence >= hc_thresh) %>%
  pull(ASV)
#find ASV with >=5% prevalence in RA
ra_asv_ids_5pct <- asv_profile_df %>%
  select(ASV, all_of(ra_cols)) %>%
  mutate(
    ra_prevalence = rowSums(across(all_of(ra_cols)) > 0, na.rm = TRUE)
  ) %>%
  filter(ra_prevalence >= ra_thresh) %>%
  pull(ASV)
#create strictly filtered taxonomy tables
HC_tax_info_5pct <- tax_info_df %>%
  filter(ASV %in% hc_asv_ids_5pct)

RA_tax_info_5pct <- tax_info_df %>%
  filter(ASV %in% ra_asv_ids_5pct)

# genus level list
library(stringr)

HC_genera_5pct <- HC_tax_info_5pct %>%
  mutate(Genus = str_extract(Taxon, "g__[^;]+")) %>%
  distinct(Genus)

RA_genera_5pct <- RA_tax_info_5pct %>%
  mutate(Genus = str_extract(Taxon, "g__[^;]+")) %>%
  distinct(Genus)

#shared versus unique bacteria
shared_genera_5pct <- intersect(HC_genera_5pct$Genus, RA_genera_5pct$Genus)
hc_only_genera_5pct <- setdiff(HC_genera_5pct$Genus, RA_genera_5pct$Genus)
ra_only_genera_5pct <- setdiff(RA_genera_5pct$Genus, HC_genera_5pct$Genus)

#sanity check outputs
nrow(HC_tax_info_5pct)
nrow(RA_tax_info_5pct)

length(HC_genera_5pct$Genus)
length(RA_genera_5pct$Genus)
################################################################################
#heat map to show top 15 shared and only RA and HC only bacterial genus
library(dplyr)
library(tidyr)
library(stringr)
library(pheatmap)

# 1) Map ASV -> Genus from taxonomy
genus_map <- tax_info_df %>%
  mutate(Genus = str_extract(Taxon, "g__[^;]+")) %>%
  filter(!is.na(Genus)) %>%
  select(ASV, Genus)

# 2) Attach genus to ASV abundance table
asv_genus_profile <- asv_profile_df %>%
  left_join(genus_map, by = "ASV") %>%
  filter(!is.na(Genus))

library(tidyr)
# 3) Long format: one row per ASV–Sample
asv_long <- asv_genus_profile %>%
  pivot_longer(
    cols = -c(ASV, Genus),
    names_to = "Sample",
    values_to = "Abundance"
  ) %>%
  mutate(
    Group = case_when(
      Sample %in% hc_cols ~ "HC",
      Sample %in% ra_cols ~ "RA",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Group))

# 4) Aggregate to genus × group mean abundance
genus_group_means <- asv_long %>%
  group_by(Genus, Group) %>%
  summarise(
    mean_abundance = mean(Abundance, na.rm = TRUE),
    .groups = "drop"
  )

# 5) Wide format: rows = Genus, columns = HC / RA
genus_mat_full <- genus_group_means %>%
  pivot_wider(
    names_from = Group,
    values_from = mean_abundance,
    values_fill = 0
  ) %>%
  as.data.frame()

rownames(genus_mat_full) <- genus_mat_full$Genus
genus_mat_full$Genus <- NULL

# 6) Pick top 15 shared genera by overall mean abundance
shared_genus_mat <- genus_mat_full[shared_genera_5pct, , drop = FALSE]

# compute overall mean (HC & RA)
shared_genus_mat$overall_mean <- rowMeans(shared_genus_mat[, c("HC", "RA")], na.rm = TRUE)

top15_shared <- shared_genus_mat %>%
  as.data.frame() %>%
  arrange(desc(overall_mean)) %>%
  head(15) %>%
  rownames()

# 7) Build final genus set: all HC-only, all RA-only, top 15 shared
all_focus_genera <- unique(c(hc_only_genera_5pct, ra_only_genera_5pct, top15_shared))

genus_mat <- genus_mat_full[all_focus_genera, c("HC", "RA"), drop = FALSE]

# 8) Log-transform for nicer heatmap (add small pseudocount)
genus_mat_log <- log10(genus_mat + 1e-6)

# 9) Row annotation: Shared / HC_only / RA_only
genus_category <- data.frame(
  Genus = rownames(genus_mat_log),
  Category = case_when(
    rownames(genus_mat_log) %in% top15_shared          ~ "Shared_top15",
    rownames(genus_mat_log) %in% hc_only_genera_5pct   ~ "HC_only",
    rownames(genus_mat_log) %in% ra_only_genera_5pct   ~ "RA_only",
    TRUE ~ "Other"
  ),
  stringsAsFactors = FALSE
)

rownames(genus_category) <- genus_category$Genus
genus_category$Genus <- NULL

# 10) Plot heatmap
pheatmap(
  genus_mat_log,
  annotation_row = genus_category,
  cluster_cols = FALSE,          # keep HC vs RA as-is
  show_rownames = TRUE,
  fontsize_row = 7,
  main = "Genus-level abundance (HC vs RA)\nHC-only, RA-only, and top 15 shared genera"
)



################################################################################
library(dplyr)
library(tidyr)
library(stringr)
install.packages("tidymodels")
library(tidymodels)

tidymodels::tidymodels_prefer()
# 1) ASV → Genus mapping from taxonomy
genus_map <- tax_info_df %>%
  mutate(Genus = str_extract(Taxon, "g__[^;]+")) %>%
  filter(!is.na(Genus)) %>%
  select(ASV, Genus)

# 2) Attach genus to ASV abundances
asv_genus_profile <- asv_profile_df %>%
  left_join(genus_map, by = "ASV") %>%
  filter(!is.na(Genus))

# 3) Long format: ASV–Sample, then aggregate to Genus–Sample
genus_sample <- asv_genus_profile %>%
  select(-ASV) %>%
  pivot_longer(
    cols = -Genus,
    names_to = "Sample",
    values_to = "Abundance"
  ) %>%
  group_by(Sample, Genus) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(
    names_from = Genus,
    values_from = Abundance,
    values_fill = 0
  )
#add coutcome label(HC Vs RA) and normalize within samples
# Label samples from their IDs
genus_sample <- genus_sample %>%
  mutate(
    label = case_when(
      Sample %in% hc_cols ~ "HC",
      Sample %in% ra_cols ~ "RA",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(label))

# Convert to relative abundance per sample
genus_cols <- setdiff(colnames(genus_sample), c("Sample", "label"))

genus_sample <- genus_sample %>%
  rowwise() %>%
  mutate(total = sum(c_across(all_of(genus_cols)))) %>%
  ungroup() %>%
  mutate(across(all_of(genus_cols), ~ .x / ifelse(total == 0, 1, total))) %>%
  select(-total)
#train/test split
set.seed(123)

data_split <- initial_split(genus_sample, prop = 0.8, strata = label)
train_data <- training(data_split)
test_data  <- testing(data_split)
#recipe(pre-processing for ML) 
rec <- recipe(label ~ ., data = train_data) %>%
  update_role(Sample, new_role = "id") %>%    # don't use Sample as predictor
  step_zv(all_predictors()) %>%              # drop zero-variance genera
  step_normalize(all_predictors())           # standardize features
#model spec: lasso logistic regression(RA probability) 
log_spec <- logistic_reg(
  penalty = tune(),   # λ will be tuned
  mixture = 1         # 1 = pure lasso
) %>%
  set_engine("glmnet")
#workflow + cross-validation tuning
set.seed(123)
folds <- vfold_cv(train_data, v = 5, strata = label)

wf <- workflow() %>%
  add_model(log_spec) %>%
  add_recipe(rec)

lambda_grid <- grid_regular(penalty(), levels = 30)

tuned <- tune_grid(
  wf,
  resamples = folds,
  grid = lambda_grid,
  metrics = metric_set(roc_auc, accuracy)
)

collect_metrics(tuned)
best_lambda <- select_best(tuned, metric = "roc_auc")
best_lambda

#final model fit and evaluation
final_wf <- finalize_workflow(wf, best_lambda)

final_fit <- fit(final_wf, data = train_data)

# Predict on held-out test set
test_pred <- predict(final_fit, test_data, type = "prob") %>%
  bind_cols(test_data %>% select(label, Sample))

# RA probability is column ".pred_RA"
head(test_pred)

# Classification (RA vs HC)
class_pred <- predict(final_fit, test_data, type = "class") %>%
  bind_cols(test_data %>% select(label))

library(dplyr)
library(yardstick)

# Combine class and probability predictions + truth
all_pred <- test_data %>%
  select(Sample, label) %>%
  bind_cols(
    predict(final_fit, test_data, type = "prob"),
    predict(final_fit, test_data, type = "class")
  )

head(all_pred)
all_pred <- test_data %>%
  select(Sample, label) %>%
  bind_cols(
    predict(final_fit, test_data, type = "prob"),
    predict(final_fit, test_data, type = "class")
  ) %>%
  mutate(label = factor(label, levels = c("HC", "RA")))

class_metrics <- metric_set(accuracy, sens, spec)

class_metrics(all_pred, truth = label, estimate = .pred_class)
roc_auc(all_pred, truth = label, .pred_HC)

###############################################################################
library(dplyr)
library(tidyr)

# Transpose: samples in rows, ASVs in columns
asv_mat <- t(as.matrix(asv_profile))   # now: rows = samples, cols = ASVs

# Build metadata
meta <- data.frame(
  Sample = rownames(asv_mat),
  Group  = ifelse(rownames(asv_mat) %in% hc_cols, "HC",
                  ifelse(rownames(asv_mat) %in% ra_cols, "RA", NA))
) %>%
  filter(!is.na(Group)) %>%
  mutate(Group = factor(Group, levels = c("HC", "RA")))

# Keep only samples with labels in the matrix
asv_mat <- asv_mat[meta$Sample, , drop = FALSE]

n_samples <- nrow(asv_mat)
prev <- colSums(asv_mat > 0)
keep_asvs <- names(prev[prev >= 0.05 * n_samples])

asv_mat_filt <- asv_mat[, keep_asvs, drop = FALSE]


library(vegan)
library(ggplot2)

# Relative abundance (row-wise)
asv_rel <- sweep(asv_mat_filt, 1, rowSums(asv_mat_filt), "/")
asv_rel[is.na(asv_rel)] <- 0

# Bray-Curtis distance
dist_bc <- vegdist(asv_rel, method = "bray")

# PERMANOVA
set.seed(123)
perm <- adonis2(dist_bc ~ Group, data = meta, permutations = 999)
perm  # report R2 and p-value in paper

# PCoA
pcoa_res <- cmdscale(dist_bc, k = 2, eig = TRUE)
pcoa_df <- data.frame(
  Sample = rownames(asv_rel),
  PC1    = pcoa_res$points[,1],
  PC2    = pcoa_res$points[,2]
) %>%
  left_join(meta, by = "Sample")

ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(alpha = 0.7, size = 2) +
  theme_minimal() +
  labs(
    x = "PCoA1",
    y = "PCoA2",
    color = "Group",
    title = "PCoA of ASV-level gut microbiome (Bray–Curtis)"
  )

BiocManager::install("DESeq2")

library(DESeq2)

# DESeq2 expects counts: ASVs in rows, samples in columns
count_mat <- as.matrix(asv_profile)
count_mat <- count_mat[, meta$Sample, drop = FALSE]  # match meta order

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = meta,
  design    = ~ Group
)
dds <- estimateSizeFactors(dds, type = "poscounts")
dds <- DESeq(dds)



res <- results(dds, contrast = c("Group", "RA", "HC"))
res_df <- as.data.frame(res) %>%
  mutate(ASV = rownames(.)) %>%
  arrange(padj)

# Join taxonomy
tax_df <- tax_info %>%
  rename(ASV = 1)  # if first column is ASV ID; adjust if needed

res_annot <- res_df %>%
  left_join(tax_df, by = "ASV")

# Significant ASVs
sig_asv <- res_annot %>%
  filter(!is.na(padj), padj < 0.05)

head(sig_asv)


library(tidymodels)
tidymodels::tidymodels_prefer()

asv_df <- as.data.frame(asv_mat_filt) %>%
  mutate(Sample = rownames(.)) %>%
  relocate(Sample)

ml_df <- asv_df %>%
  left_join(meta, by = "Sample") %>%
  filter(!is.na(Group))

ml_df$Group <- factor(ml_df$Group, levels = c("HC", "RA"))

set.seed(123)
split <- initial_split(ml_df, prop = 0.8, strata = Group)
train_data <- training(split)
test_data  <- testing(split)

asv_cols <- base::setdiff(colnames(train_data), c("Sample", "Group"))

rec_asv <- recipe(Group ~ ., data = train_data) %>%
  update_role(Sample, new_role = "id") %>%
  step_zv(all_predictors()) %>%
  step_log(all_predictors(), offset = 1e-6) %>%   # log-transform counts
  step_normalize(all_predictors())

# 1) Lasso logistic
log_spec <- logistic_reg(
  penalty = tune(),
  mixture = 1
) %>% set_engine("glmnet")

# 2) Random Forest (ranger)
rf_spec <- rand_forest(
  mtry  = tune(),
  trees = 500,
  min_n = tune()
) %>% 
  set_engine("ranger") %>%
  set_mode("classification")

# 3) XGBoost
xgb_spec <- boost_tree(
  trees = 1000,
  learn_rate = tune(),
  mtry = tune(),
  tree_depth = tune(),
  min_n = tune(),
  loss_reduction = tune()
) %>%
  set_engine("xgboost") %>%
  set_mode("classification")

# 4) SVM (radial)
svm_spec <- svm_rbf(
  cost = tune(),
  rbf_sigma = tune()
) %>%
  set_engine("kernlab") %>%
  set_mode("classification")

set.seed(123)
folds <- vfold_cv(train_data, v = 5, strata = Group)

ctrl <- control_grid(save_pred = TRUE)

log_wf <- workflow() %>% add_model(log_spec) %>% add_recipe(rec_asv)
rf_wf  <- workflow() %>% add_model(rf_spec)  %>% add_recipe(rec_asv)
xgb_wf <- workflow() %>% add_model(xgb_spec) %>% add_recipe(rec_asv)
svm_wf <- workflow() %>% add_model(svm_spec) %>% add_recipe(rec_asv)

# Example: tune logistic first
log_grid <- grid_regular(penalty(range = c(-4, 0)), levels = 20)

log_tuned <- tune_grid(
  log_wf,
  resamples = folds,
  grid = log_grid,
  metrics = yardstick::metric_set(yardstick::roc_auc, yardstick::accuracy),
  control = ctrl
)

collect_metrics(log_tuned)
best_log <- select_best(log_tuned, metric = "roc_auc")
best_log


final_log_wf <- finalize_workflow(log_wf, best_log)
final_log_fit <- fit(final_log_wf, data = train_data)

log_pred <- test_data %>%
  select(Sample, Group) %>%
  bind_cols(
    predict(final_log_fit, test_data, type = "prob"),
    predict(final_log_fit, test_data, type = "class")
  )

log_pred <- log_pred %>%
  mutate(Group = factor(Group, levels = c("HC", "RA")))

# Metrics
class_metrics <- yardstick::metric_set(yardstick::accuracy,
                                       yardstick::sens, 
                                       yardstick::spec)
class_metrics(log_pred, truth = Group, estimate = .pred_class)

yardstick::roc_auc(log_pred, truth = Group, .pred_RA, event_level = "second")







