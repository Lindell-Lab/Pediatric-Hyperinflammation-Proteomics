##============================================================================##
# Olinking_FINAL.R
# Comparative plasma and serum proteomics in pediatric sepsis, MIS-C, and CRS
# Critical Care Explorations
##============================================================================##

# Load packages
library(tidyverse)
library(dplyr)
library(here)

## Load and process data
set.seed(3158)
olinking_df <- read_csv(here::here("data", "final_df_meta_FINAL.csv"))
olinking_df <- olinking_df %>% mutate(case_control = case_when(
  Condition=="COVID_healthy" ~ "Control",
  Condition=="COVID_minimal" ~ "Case",
  Condition=="COVID_MISC" ~ "Case",
  Condition=="COVID_severe" ~ "Case",
  Condition=="CRS_minimal" ~ "Case",
  Condition=="CRS_severe" ~ "Case",
  Condition=="Sepsis_Healthy Control" ~ "Control",
  Condition=="Sepsis_MODS" ~ "Case"))
olinking_df$case_control <- factor(olinking_df$case_control, levels=c("Case", "Control"))

# Impute IL6 and IL6R a priori
olinking_df <- olinking_df %>% select(-"IL6R", -"IL6")

# Select the peak timepoint for each patient
# Peak timepoint is Timepoint=="1" for most patients (time of diagnosis and enrollment)
# CRS patients with pre-CAR-T (Timepoint=="1") and post-CAR-T (Timepoint=="2") samples per protocol
# Timepoint=="3" is onset of CRS, some severe CRS patients reach peak at Timepoint=="4" per Penn Scale
olinking_df <- olinking_df %>% mutate(peak_timepoint = case_when(
  Condition=="COVID_healthy" & Timepoint=="1" ~ "Peak",
  Condition=="COVID_minimal" & Timepoint=="1" ~ "Peak",
  Condition=="COVID_MISC" & Timepoint=="1" ~ "Peak",
  Condition=="COVID_severe" & Timepoint=="1" ~ "Peak",
  Condition=="CRS_minimal" & Timepoint=="3" ~ "Peak",
  Condition=="CRS_severe" & Timepoint=="3" ~ "Peak",
  Condition=="CRS_severe" & Timepoint=="4" & analysis_sample == TRUE ~ "Peak",
  Condition=="Sepsis_Healthy Control" & Timepoint=="0" ~ "Peak",
  Condition=="Sepsis_MODS" & Timepoint=="1" ~ "Peak"))
olinking_df %>% group_by(peak_timepoint) %>% count()
olinking_df %>% group_by(peak_timepoint, Condition) %>% count()
olinking_df_peak <- filter(olinking_df, peak_timepoint=="Peak")

olinking_data_labels <- read_csv(here::here("data", "olink_data_labels.csv"))
olinking_data_labels <- olinking_data_labels %>% filter(Assay != "IL6")
olinking_df_peak_long <- olinking_df_peak %>%
  pivot_longer(cols = where(is.numeric), names_to = "Assay", values_to = "NPX")
olinking_df_peak_long <- left_join(olinking_df_peak_long, olinking_data_labels)
rm(olinking_data_labels)

olinking_df_peak_wide <- olinking_df_peak_long %>%
  select("SubjectID", "Condition", "case_control", "Assay", "NPX") %>%
  distinct() %>%                                       
  pivot_wider(names_from = Assay, values_from = NPX)

table(olinking_df_peak_wide$Condition)

##----------------------------------------------------------------------------##

## PCA and UMAP for dimensionality reduction and visualization
library("umap")
library("FactoMineR")
library("factoextra")

# Define custom colors
custom_colors <- c(
  "CRS_severe" = "#CC33FF",
  "Sepsis_MODS" = "#E6004C",
  "COVID_MISC" = "#FF9900",
  "COVID_severe" = "#FFCC00",
  "COVID_minimal" = "#33CC33",
  "CRS_minimal" = "#33CCCC",
  "Sepsis_Healthy Control" = "#3399FF")

olinking_df_peak$Condition <- factor(olinking_df_peak$Condition,
                                     levels = names(custom_colors))

# PCA for overall data architecture
protein_matrix <- olinking_df_peak %>%
  select(-SubjectID, -Condition, -case_control, -Timepoint, -peak_timepoint) %>%
  select(where(is.numeric))

olinking_pca <- PCA(protein_matrix, ncp = 20, scale.unit = TRUE, graph = FALSE)

# Extract PC scores for first two components
pca_scores <- as.data.frame(olinking_pca$ind$coord[, 1:2])
pca_scores$Condition <- olinking_df_peak$Condition

colnames(pca_scores)[1:2] <- c("PC1", "PC2")

# Extract PC % variance
pc1_var <- round(olinking_pca$eig[1, 2], 1)
pc2_var <- round(olinking_pca$eig[2, 2], 1)

# Plot with ggplot
ggplot(pca_scores, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(
    shape = 21, 
    size = 3,
    stroke = 0.7,
    fill = custom_colors[pca_scores$Condition]
  ) +
  scale_color_manual(values = custom_colors) +
  labs(
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)"),
    color = "Condition"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold")
  )

# Extract PCA loadings
loadings <- as.data.frame(olinking_pca$var$coord)

# Add protein names
loadings$Protein <- rownames(loadings)

# Top 10 for PC1
top10_PC1 <- loadings %>%
  dplyr::select(Protein, Dim.1) %>%
  dplyr::arrange(desc(abs(Dim.1))) %>%
  dplyr::slice(1:10)

# Top 10 for PC2
top10_PC2 <- loadings %>%
  dplyr::select(Protein, Dim.2) %>%
  dplyr::arrange(desc(abs(Dim.2))) %>%
  dplyr::slice(1:10)

top10_PC1
top10_PC2

# Eigenvalues
fviz_eig(olinking_pca, addlabels = TRUE, ylim = c(0, 50), ncp=15) #at PC15 Eigenvalue is <1%
summary(olinking_pca) #cum var >70% at Dim17, cum var >80% at Dim23

# Set up PCA data as matrix
pca_data <- as.matrix(olinking_pca$ind$coord)

# Run UMAP
olinking_UMAP <- umap::umap(pca_data, n_neighbors = 15, min_dist = 0.1, spread = 1.0)

# Create a data frame for UMAP results
umap_df <- as.data.frame(olinking_UMAP$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")

# Add condition labels
umap_df$Condition <- olinking_df_peak$Condition

# Ensure the levels of the Condition factor are ordered and renamed
umap_df$Condition <- factor(umap_df$Condition, levels = c("Sepsis_Healthy Control", "COVID_minimal", "COVID_severe", "COVID_MISC","CRS_minimal", "CRS_severe", "Sepsis_MODS"))

# Plot UMAP results
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Condition)) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_color_manual(values = custom_colors) +
  theme_minimal() +
  labs(title = "UMAP of PCA Results", x = "UMAP1", y = "UMAP2") +
  theme(
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),   
        axis.title = element_text(size = 16), 
        axis.text = element_text(size = 16),
        plot.title = element_text(size = 16, hjust = 0.5)
       ) +
  guides(color = guide_legend(order = 1))

##----------------------------------------------------------------------------##

## Extract top features in UMAP1 and UMAP2
library(xgboost)
library(ggplot2)

X_prot <- olinking_df_peak %>%
  select(-SubjectID, -Condition, -case_control, -Timepoint, -peak_timepoint) %>%
  select(where(is.numeric))
stopifnot(nrow(X_prot) == nrow(umap_df))

X_mat <- scale(as.matrix(X_prot))
feature_names <- colnames(X_mat)

# Spearman correlation across UMAP axes
spearman_rank_axis <- function(axis_vec, X_mat, feature_names) {
  rho <- apply(X_mat, 2, function(v) suppressWarnings(cor(v, axis_vec, method = "spearman")))
  p   <- apply(X_mat, 2, function(v) suppressWarnings(cor.test(v, axis_vec, method = "spearman")$p.value))
  out <- data.frame(
    protein = feature_names,
    rho = as.numeric(rho),
    p = as.numeric(p),
    fdr = p.adjust(p, method = "fdr")
  ) %>% arrange(desc(abs(rho)))
  out
}

spearman_umap1 <- spearman_rank_axis(umap_df$UMAP1, X_mat, feature_names)
spearman_umap2 <- spearman_rank_axis(umap_df$UMAP2, X_mat, feature_names)

head(spearman_umap1, 20)
head(spearman_umap2, 20)

##----------------------------------------------------------------------------##

## Run ssGSVA and prepare matrix for heatmap
library("GSVA")
library("msigdbr")

# Define gene set
hallmark_immune_pathways <- c(
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING"
)

hallmark_msigdb <- msigdbr(species = "Homo sapiens", category = "H")
GSEA_gene_list <- split(hallmark_msigdb$gene_symbol, hallmark_msigdb$gs_name)
GSEA_gene_list <- lapply(GSEA_gene_list, unique)
GSEA_gene_list <- GSEA_gene_list[hallmark_immune_pathways]

# Create matrix
proteomic_matrix <- olinking_df_peak_long %>%
  select(SubjectID, Assay, NPX) %>%
  distinct() %>%
  pivot_wider(names_from = Assay, values_from = NPX) %>%
  column_to_rownames("SubjectID") %>%
  as.matrix()

common_genes <- intersect(rownames(proteomic_matrix), unlist(GSEA_gene_list))
length(common_genes)

GSEA_matrix <- t(proteomic_matrix)

# Create the parameter object
param <- gsvaParam(
  exprData = GSEA_matrix,
  geneSets = GSEA_gene_list,
  kcdf = "Gaussian"
)

# Run GSVA
gsva_scores <- gsva(param)

# Convert gsva_scores to dataframe
gsva_df <- as.data.frame(t(gsva_scores))

# Align metadata to GSVA sample order
sample_ids <- colnames(gsva_scores)

meta <- olinking_df_peak_wide %>%
  select(SubjectID, Condition)

meta_ordered <- meta[match(sample_ids, meta$SubjectID), ]

# Add condition vector to GSVA df
condition_vector <- dplyr::recode(
  meta_ordered$Condition,
  
  # COVID groups
  "COVID_minimal" = "COVID Minimal (n=26)",
  "COVID_severe" = "COVID Severe (n=15)",
  "COVID_MISC" = "MISC (n=22)",
  
  # CRS groups
  "CRS_minimal" = "CRS Minimal (n=13)",
  "CRS_severe" = "CRS Severe (n=13)",
  
  # Sepsis groups
  "Sepsis_Healthy Control" = "Healthy Controls (n=26)",
  "Sepsis_MODS" = "Sepsis (n=35)"
)

gsva_df$Condition <- factor(
  condition_vector,
  levels = c(
    "Healthy Controls (n=26)",
    "COVID Minimal (n=26)",
    "CRS Minimal (n=13)",
    "Sepsis (n=35)",
    "COVID Severe (n=15)",
    "MISC (n=22)",
    "CRS Severe (n=13)"
  )
)

# Label pathway columns
hallmark_cols <- grep("^HALLMARK_", names(gsva_df), value = TRUE)

# Kruskal-Wallis omnibus test
library(rstatix)
library(tidyr)

gsva_long <- gsva_df %>%
  select(all_of(hallmark_cols), Condition) %>%
  pivot_longer(cols = all_of(hallmark_cols), names_to = "pathway", values_to = "score")

kw_results <- gsva_long %>%
  group_by(pathway) %>%
  kruskal_test(score ~ Condition) %>%
  ungroup()

fdr_method <- "BH"

kw_results$p_adj <- p.adjust(kw_results$p, method = fdr_method)
kw_results <- kw_results %>% arrange(p_adj)
kw_results

library(writexl)
write_xlsx(kw_results, here::here("output", "Supplemental Digital content 2.xlsx"))

# Dunn's post-hoc cross-condition pairwise comparisons w/ BH correction
dunn_raw <- gsva_long %>%
  group_by(pathway) %>%
  dunn_test(score ~ Condition, p.adjust.method = "none") %>%
  ungroup()

dunn_raw$p_adj_global <- p.adjust(dunn_raw$p, method = fdr_method)

dunn_results <- dunn_raw %>% arrange(p_adj_global)
dunn_results

dunn_results %>% filter(p_adj_global < 0.05)

write_xlsx(dunn_results, here::here("output", "Supplemental Digital content 3.xlsx"))

# Row normalize and then create heatmap
# Define row-wise min-max normalization function
min_max_scale <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep(0, length(x)))  # Avoid division by zero
  return(2 * (x - rng[1]) / (rng[2] - rng[1]) - 1)
}

# Apply to each row of the matrix
gsva_scores_rescaled <- apply(gsva_scores, 1, min_max_scale)

# Create the heatmap
library(ComplexHeatmap)
library(circlize)
library(grid)

# Subset and rename matrix rows
gsva_scores_rescaled <- t(gsva_scores_rescaled)

# Define condition factor for column splitting and relabel
condition_factor <- factor(
  olinking_df_peak_wide$Condition,
  levels = c(
    "Sepsis_Healthy Control", "COVID_minimal", "CRS_minimal", "Sepsis_MODS",
    "COVID_severe", "COVID_MISC", "CRS_severe"
  )
)

levels(condition_factor) <- c(
  "Healthy Controls (n=26)",
  "COVID Minimal (n=26)",
  "CRS Minimal (n=13)",
  "Sepsis (n=35)",
  "COVID Severe (n=15)",
  "MISC (n=22)",
  "CRS Severe (n=13)"
)

# Define group colors
group_colors <- c(
  "CRS Severe (n=13)" = "#CC33FF",
  "Sepsis (n=35)" = "#E6004C",
  "MISC (n=22)" = "#FF9900",
  "COVID Severe (n=15)" = "#FFCC00",
  "COVID Minimal (n=26)" = "#33CC33",
  "CRS Minimal (n=13)" = "#33CCCC",
  "Healthy Controls (n=26)" = "#3399FF"
)

# Create column & row annotation
col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

ha_col <- HeatmapAnnotation(
  Condition = condition_factor,
  col = list(Condition = group_colors),
  annotation_name_side = "left",
  show_legend = TRUE
)

# Define condition levels explicitly
condition_levels <- c(
  "Healthy Controls (n=26)",
  "COVID Minimal (n=26)",
  "CRS Minimal (n=13)",
  "Sepsis (n=35)",
  "COVID Severe (n=15)",
  "MISC (n=22)",
  "CRS Severe (n=13)"
)

condition_factor <- factor(condition_factor, levels = condition_levels)

# Sort columns within each group by average NES (from blue to red)
column_order_by_nes <- unlist(lapply(levels(condition_factor), function(group) {
  group_cols <- which(condition_factor == group)
  avg_nes <- colMeans(gsva_scores_rescaled[, group_cols, drop = FALSE])
  group_cols[order(avg_nes)]  # ascending: blue to red
}))

# Reorder all inputs
gsva_scores_ordered <- gsva_scores_rescaled[, column_order_by_nes]
condition_factor_ordered <- condition_factor[column_order_by_nes]
ha_col_ordered <- ha_col[column_order_by_nes]
      
# Plot with manual order (no additional clustering to maintain in-group sorting)
ComplexHeatmap::Heatmap(gsva_scores_ordered,
    name = "NES",
    col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
    top_annotation = ha_col_ordered,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    column_split = condition_factor_ordered,
    cluster_column_slices = FALSE,
    show_column_names = FALSE,
    show_row_names = TRUE,
    column_names_gp = gpar(fontsize = 12),
    row_names_gp = gpar(fontsize = 12),
    row_dend_gp = gpar(lwd = 3),
    heatmap_legend_param = list(title = "Normalized NES", legend_height = unit(4, "cm")))

## Radar graph to see overlap between Sepsis, MIS-C, and CRS Severe
library(fmsb)

# Convert GSVA matrix to long format
GSVA_df_long <- gsva_scores_ordered %>%
  as.data.frame() %>%
  t() %>% 
  as.data.frame() %>%
  mutate(Condition = condition_factor_ordered) %>%
  group_by(Condition) %>%
  summarise(across(everything(), mean)) %>%
  as.data.frame()

df_radar <- GSVA_df_long %>% 
  filter(Condition %in% c("Healthy Controls (n=26)", "Sepsis (n=35)", "MISC (n=22)", "CRS Severe (n=13)"))                        

# Extract only numeric pathway columns
pathway_values <- df_radar %>% select(-Condition)

# Create max/min rows
max_row <- apply(pathway_values, 2, max)
min_row <- apply(pathway_values, 2, min)

# Build final radar df
df_radar_final <- rbind(max_row, min_row, pathway_values)
rownames(df_radar_final) <- c("max", "min", df_radar$Condition)

# Create the radar graph
radarchart(
  df_radar_final,
  axistype = 1,
  pcol = c("#3399FF", "#E6004C", "#FF9900", "#CC33FF"),
  pfcol = c("#3399FF40", "#E6004C40", "#FF990040", "#CC33FF40"),
  plwd = 3,
  plty = 1,
  cglcol = "grey80",
  cglty = 1,
  cglwd = 0.6,
  axislabcol = "grey20",
  vlcex = 0.5,
  cex = 1.2
)

legend(
  x = 1.35, y = 1,   # move legend far to the right
  legend = rownames(df_radar_final)[3:5],
  col = c("#3399FF", "#E6004C", "#FF9900", "#CC33FF"),
  lwd = 3,
  bty = "n",
  cex = 1.4
)

##----------------------------------------------------------------------------##

## Use XGBoost to identify minimal set of proteins that differentiate MIS-C from Sepsis
# Remove rows & columns not being used for feature ranking
binary_df_MISC_sepsis <- olinking_df_peak_wide %>%
  filter(Condition %in% c("COVID_MISC", "Sepsis_MODS")) %>%
  select(-SubjectID, -case_control, -Timepoint)

binary_df_MISC_sepsis$Condition <- factor(binary_df_MISC_sepsis$Condition)
levels(binary_df_MISC_sepsis$Condition)

binary_df_MISC_sepsis <- binary_df_MISC_sepsis %>%
  mutate(
    Condition_binary = dplyr::recode(
      Condition,
      "Sepsis_MODS" = 0,
      "COVID_MISC" = 1
    )
  )

# Prepare features & labels
features <- binary_df_MISC_sepsis %>% select(-Condition, -Condition_binary)
labels <- binary_df_MISC_sepsis$Condition

# Convert to matrix
features_matrix <- as.matrix(features)
labels_vector <- as.numeric(binary_df_MISC_sepsis[["Condition_binary"]])

print(dim(features_matrix))
print(length(labels_vector))

library(caret)

set.seed(3158)
feature_names <- colnames(features)

# Define model parameters
params <- list(
  objective = "binary:logistic",
  eval_metric = "aucpr",
  max_depth = 4,
  eta = 0.1,
  subsample = 1,
  colsample_bytree = 1,
  lambda = 5,
  alpha = 5,
  gamma = 1
)

# Train XGBoost using repeated stratified 5-fold CV (no bootstrap) for feature-importance
set.seed(3158)
feature_names <- colnames(features)

n_reps <- 40   # 40 repeats x nfold = 200 total fold-fits (matches old n_boot = 200)
top_n <- 10
nfold <- 5
nrounds_max <- 500
early_stopping_rounds <- 20

feature_counts <- setNames(numeric(length(feature_names)), feature_names)
best_nrounds_used <- c()
cv_aucpr <- c()

for (rep in seq_len(n_reps)) {
  fold_ids <- caret::createFolds(factor(labels_vector), k = nfold, list = FALSE)
  
  for (k in seq_len(nfold)) {
    test_idx  <- which(fold_ids == k)
    train_idx <- setdiff(seq_len(length(labels_vector)), test_idx)
    
    dtrain <- xgb.DMatrix(data = features_matrix[train_idx, , drop = FALSE],
                          label = labels_vector[train_idx])
    dtest  <- xgb.DMatrix(data = features_matrix[test_idx, , drop = FALSE],
                          label = labels_vector[test_idx])
    
    fold_fit <- xgb.train(
      params = params,
      data = dtrain,
      nrounds = nrounds_max,
      watchlist = list(val = dtest),
      early_stopping_rounds = early_stopping_rounds,
      verbose = 0
    )
    best_nrounds_used <- c(best_nrounds_used, fold_fit$best_iteration)
    
    metric_col <- grep("val.*aucpr$", colnames(fold_fit$evaluation_log), value = TRUE)[1]
    cv_aucpr <- c(cv_aucpr, fold_fit$evaluation_log[[metric_col]][fold_fit$best_iteration])
    
    importance <- xgb.importance(feature_names = feature_names, model = fold_fit)
    top_feats <- head(importance$Feature, top_n)
    feature_counts[top_feats] <- feature_counts[top_feats] + 1
  }
}

n_boot <- n_reps * nfold

cat("Mean CV-selected nrounds:", mean(best_nrounds_used), "\n")
cat("Median CV AUCPR:", median(cv_aucpr, na.rm = TRUE), "\n")

# Most frequently selected features
top10_MISC <- data.frame(
  Feature = names(feature_counts),
  times_in_topN = as.integer(feature_counts),
  pct_of_resamples = round(100 * feature_counts / n_boot, 1)
) %>%
  as_tibble() %>%
  arrange(desc(pct_of_resamples)) %>%
  dplyr::slice(1:10)

top10_MISC

##----------------------------------------------------------------------------##

## Use XGBoost to identify minimal set of proteins that differentiate CRS Severe from Sepsis
# Remove rows & columns not being used for feature ranking
binary_df_CRS_sepsis <- olinking_df_peak_wide %>%
  filter(Condition %in% c("CRS_severe", "Sepsis_MODS")) %>%
  select(-SubjectID, -case_control, -Timepoint)

binary_df_CRS_sepsis$Condition <- factor(binary_df_CRS_sepsis$Condition)
levels(binary_df_CRS_sepsis$Condition)

binary_df_CRS_sepsis <- binary_df_CRS_sepsis %>%
  mutate(
    Condition_binary = dplyr::recode(
      Condition,
      "Sepsis_MODS" = 0,
      "CRS_severe" = 1
    )
  )

# Prepare features & labels
features1 <- binary_df_CRS_sepsis %>% select(-Condition, -Condition_binary)
labels1 <- binary_df_CRS_sepsis$Condition

# Convert to matrix
features_matrix1 <- as.matrix(features1)
labels_vector1 <- as.numeric(binary_df_CRS_sepsis[["Condition_binary"]])

print(dim(features_matrix1))
print(length(labels_vector1))

set.seed(3158)
features_names1 <- colnames(features1)

# Define model parameters
params <- list(
  objective = "binary:logistic",
  eval_metric = "aucpr",
  max_depth = 4,
  eta = 0.1,
  subsample = 1,
  colsample_bytree = 1,
  lambda = 5,
  alpha = 5,
  gamma = 1
)

# Train XGBoost using repeated stratified 5-fold CV (no bootstrap) for feature-importance
set.seed(3158)
features_names1 <- colnames(features1)

n_reps1 <- 40
top_n1 <- 10
nfold1 <- 5
nrounds_max1 <- 500
early_stopping_rounds1 <- 20

feature_counts1 <- setNames(numeric(length(features_names1)), features_names1)
best_nrounds_used1 <- c()
cv_aucpr1 <- c()

for (rep in seq_len(n_reps1)) {
  fold_ids1 <- caret::createFolds(factor(labels_vector1), k = nfold1, list = FALSE)
  
  for (k in seq_len(nfold1)) {
    test_idx1  <- which(fold_ids1 == k)
    train_idx1 <- setdiff(seq_len(length(labels_vector1)), test_idx1)
    
    dtrain1 <- xgb.DMatrix(data = features_matrix1[train_idx1, , drop = FALSE],
                           label = labels_vector1[train_idx1])
    dtest1  <- xgb.DMatrix(data = features_matrix1[test_idx1, , drop = FALSE],
                           label = labels_vector1[test_idx1])
    
    fold_fit1 <- xgb.train(
      params = params,
      data = dtrain1,
      nrounds = nrounds_max1,
      watchlist = list(val = dtest1),
      early_stopping_rounds = early_stopping_rounds1,
      verbose = 0
    )
    best_nrounds_used1 <- c(best_nrounds_used1, fold_fit1$best_iteration)
    
    metric_col1 <- grep("val.*aucpr$", colnames(fold_fit1$evaluation_log), value = TRUE)[1]
    cv_aucpr1 <- c(cv_aucpr1, fold_fit1$evaluation_log[[metric_col1]][fold_fit1$best_iteration])
    
    importance1 <- xgb.importance(feature_names = features_names1, model = fold_fit1)
    top_feats1 <- head(importance1$Feature, top_n1)
    feature_counts1[top_feats1] <- feature_counts1[top_feats1] + 1
  }
}

n_boot <- n_reps1 * nfold1

cat("Mean CV-selected nrounds:", mean(best_nrounds_used1), "\n")
cat("Median CV AUCPR:", median(cv_aucpr1, na.rm = TRUE), "\n")

# Most frequently selected features
top10_CRS <- data.frame(
  Feature = names(feature_counts1),
  times_in_topN = as.integer(feature_counts1),
  pct_of_resamples = round(100 * feature_counts1 / n_boot, 1)
) %>%
  as_tibble() %>%
  arrange(desc(pct_of_resamples)) %>%
  dplyr::slice(1:10)

top10_CRS
