# Load packages
library(tidyverse)
library(dplyr)

## Load and process data ##

set.seed(3158)
setwd("C:/Users/USER/Documents/Olinking")
olinking_df <- read_csv("C:/Users/USER/Documents/Olinking/final_df_meta_012826.csv")
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
olinking_df <- olinking_df %>% select(-"IL6R")

# Select the peak timepoint for each patient
olinking_df <- olinking_df %>% mutate(peak_timepoint = case_when(
  Condition=="COVID_healthy" & Timepoint=="1" ~ "Peak",
  Condition=="COVID_minimal" & Timepoint=="1" ~ "Peak",
  Condition=="COVID_MISC" & Timepoint=="1" ~ "Peak",
  Condition=="COVID_severe" & Timepoint=="1" ~ "Peak",
  Condition=="CRS_minimal" & Timepoint=="3" ~ "Peak",
  Condition=="CRS_severe" & Timepoint=="3" ~ "Peak",
  Condition=="CRS_severe" & Timepoint=="4" & (SubjectID=="CHP959-120" | SubjectID=="CHP959-131" | SubjectID=="CHP959-139") ~ "Peak",
  Condition=="Sepsis_Healthy Control" & Timepoint=="0" ~ "Peak",
  Condition=="Sepsis_MODS" & Timepoint=="1" ~ "Peak"))
olinking_df %>% group_by(peak_timepoint) %>% count()
olinking_df %>% group_by(peak_timepoint, Condition) %>% count()
olinking_df_peak <- filter(olinking_df, peak_timepoint=="Peak")
olinking_df_peak <- olinking_df_peak %>% filter(Condition != "COVID_healthy")

olinking_data_labels <- read_csv("C:/Users/USER/Documents/Olinking/olink_data_labels.csv")
olinking_data_labels <- olinking_data_labels %>% filter(Assay != "IL6")
olinking_df_peak_long <- olinking_df_peak %>%
  pivot_longer(cols = where(is.numeric), names_to = "Assay", values_to = "NPX")
olinking_df_peak_long <- left_join(olinking_df_peak_long, olinking_data_labels)
rm(olinking_data_labels)
olinking_df_peak_long <- olinking_df_peak_long %>% filter(Condition != "COVID_healthy")

olinking_df_peak_wide <- olinking_df_peak_long %>%
  select("SubjectID", "Condition", "case_control", "Assay", "NPX") %>%
  distinct() %>%                                       
  pivot_wider(names_from = Assay, values_from = NPX)

table(olinking_df_peak_wide$Condition)

##----------------------------------------------------------------------------##

## PCA and UMAP for dimensionality reduction and visualization
library("OlinkAnalyze")
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
olinking_pca <- PCA(olinking_df_peak[,2:545], ncp = 20, scale.unit = TRUE, graph = FALSE)

# Extract PC scores for first two components
pca_scores <- as.data.frame(olinking_pca$ind$coord[, 1:2])
pca_scores$Condition <- olinking_df_peak$Condition

colnames(pca_scores)[1:2] <- c("PC1", "PC2")

# Extract PC % variance
pc1_var <- round(olinking_pca$eig[1, 2], 1)
pc2_var <- round(olinking_pca$eig[2, 2], 1)

# Plot with ggplot: all solid circles, custom colors
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

# Extract loadings (coordinates of variables on PCs)
loadings <- olinking_pca$var$coord

# Top 10 positive for PC1
top_pc1 <- sort(abs(loadings[, "Dim.1"]), decreasing = TRUE)[1:10]
names_pc1 <- names(top_pc1)

top_loadings_df_pc1 <- as_tibble(loadings, rownames = "Protein") %>%
  mutate(abs_PC1 = abs(Dim.1), abs_PC2 = abs(Dim.2)) %>%
  arrange(desc(abs_PC1)) %>%
  slice(1:10) %>%
  select(Protein, Dim.1)

# Top 10 most negative for PC1
neg_pc1 <- sort(loadings[, "Dim.1"])[1:10]
names_neg_pc1 <- names(neg_pc1)

neg_loadings_pc1 <- as_tibble(loadings, rownames = "Protein") %>%
  arrange(Dim.1) %>%
  slice(1:10) %>%
  select(Protein, Dim.1)

# Top 10 positive for PC2
top_pc2 <- sort(abs(loadings[, "Dim.2"]), decreasing = TRUE)[1:10]
names_pc2 <- names(top_pc2)

top_loadings_df_pc2 <- as_tibble(loadings, rownames = "Protein") %>%
  mutate(abs_PC2 = abs(Dim.2)) %>%
  arrange(desc(abs_PC2)) %>%
  slice(1:10) %>%
  select(Protein, Dim.2)

# Top 10 most negative for PC2
neg_pc2 <- sort(loadings[, "Dim.2"])[1:10]
names_neg_pc2 <- names(neg_pc2)

neg_loadings_pc2 <- as_tibble(loadings, rownames = "Protein") %>%
  arrange(Dim.2) %>%
  slice(1:10) %>%
  select(Protein, Dim.2)

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


## PCA and UMAP for dimensionality reduction and visualization with Healthy omitted

olinking_df_filtered <- olinking_df_peak %>%
  filter(Condition != "Sepsis_Healthy Control")

# Define custom colors
custom_colors <- c(
  "CRS_severe" = "#CC33FF",
  "Sepsis_MODS" = "#E6004C",
  "COVID_MISC" = "#FF9900",
  "COVID_severe" = "#FFCC00",
  "COVID_minimal" = "#33CC33",
  "CRS_minimal" = "#3399FF")

olinking_df_filtered$Condition <- factor(olinking_df_filtered$Condition,
                                     levels = names(custom_colors))

# PCA for overall data architecture
olinking_pca <- PCA(olinking_df_filtered[,2:545], ncp = 20, scale.unit = TRUE, graph = FALSE)

fviz_pca_ind(olinking_pca,
             geom.ind = "point",
             pointsize = 2,
             alpha = 0.8,
             col.ind = olinking_df_filtered$Condition,
             addEllipses = TRUE,
             ellipse.level = 0.95,
             ellipse.type = "confidence",
             repel = TRUE,
             title = "",
             legend.title = "Condition") +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  scale_shape_manual(values = c(8, 12, 13, 15, 2, 17, 16)) +
  theme(
    text = element_text(size = 16),               # Global text size
    axis.title = element_text(size = 16),         # Axis titles
    axis.text = element_text(size = 16),          # Axis tick labels
    legend.title = element_text(size = 16),       # Legend title
    legend.text = element_text(size = 16),        # Legend labels
    plot.title = element_text(size = 16, face = "bold")  # Plot title
  )

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
umap_df$Condition <- olinking_df_filtered$Condition

# Ensure the levels of the Condition factor are ordered and renamed
umap_df$Condition <- factor(umap_df$Condition, levels = c("COVID_minimal", "COVID_severe", "COVID_MISC","CRS_minimal", "CRS_severe", "Sepsis_MODS"))

# Plot UMAP results
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Condition)) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_color_manual(values = custom_colors) +
  theme_minimal() +
  labs(title = "UMAP of PCA Results with Healthy Patients Omitted", x = "UMAP1", y = "UMAP2") +
  theme(
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),   
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 16, hjust = 0.5)
  ) +
  guides(color = guide_legend(order = 1))


##----------------------------------------------------------------------------##


## Extract top features in UMAP1 and UMAP2 (with healthies included)
library(dplyr)
library(xgboost)
library(ggplot2)

X_prot <- olinking_df_peak %>%  select(,2:545) %>% select(where(is.numeric))
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

# XGBoost to predict UMAP1/UMAP2 from proteins

xgb_cv_train_axis <- function(y, X_mat, seed = 3158) {
  set.seed(seed)
  d <- xgb.DMatrix(data = X_mat, label = y)
  
  params_reg <- list(
    objective = "reg:squarederror",
    eval_metric = "rmse",
    max_depth = 4,
    eta = 0.1,
    subsample = 1,
    colsample_bytree = 1,
    lambda = 5,
    alpha = 5,
    gamma = 1
  )
  
  cv <- xgb.cv(
    params = params_reg,
    data = d,
    nrounds = 2000,
    nfold = 5,
    early_stopping_rounds = 25,
    verbose = 0
  )
  
  best_nrounds <- cv$best_iteration
  
  model <- xgb.train(
    params = params_reg,
    data = d,
    nrounds = best_nrounds,
    verbose = 0
  )
  
  imp <- xgb.importance(feature_names = colnames(X_mat), model = model)
  
  list(cv = cv, best_nrounds = best_nrounds, model = model, importance = imp)
}

fit_umap1 <- xgb_cv_train_axis(umap_df$UMAP1, X_mat, seed = 3158)
fit_umap2 <- xgb_cv_train_axis(umap_df$UMAP2, X_mat, seed = 3158)

# Top proteins by XGBoost gain
head(fit_umap1$importance, 20)
head(fit_umap2$importance, 20)


##----------------------------------------------------------------------------##


## Run ssGSVA and prepare matrix for heatmap
library("GSVA")

# Define gene set
GSEA_gene_list <- list(
  HALLMARK_TNFA_SIGNALING_VIA_NFKB = c("ABCA1", "ACKR3", "AREG", "ATF3", "ATP2B1", "B4GALT1", "B4GALT5", "BCL2A1", "BCL3", "BCL6", "BHLHE40", "BIRC2", "BIRC3", "BMP2", "BTG1", "BTG2", "BTG3", "CCL2", "CCL20", "CCL4", "CCL5", "CCN1", "CCND1", "CCNL1", "CCRL2", "CD44", "CD69", "CD80", "CD83", "CDKN1A", "CEBPB", "CEBPD", "CFLAR", "CLCF1", "CSF1", "CSF2", "CXCL1", "CXCL10", "CXCL11", "CXCL2", "CXCL3", "CXCL6", "DENND5A", "DNAJB4", "DRAM1", "DUSP1", "DUSP2", "DUSP4", "DUSP5", "EDN1", "EFNA1", "EGR1", "EGR2", "EGR3", "EHD1", "EIF1", "ETS2", "F2RL1", "F3", "FJX1", "FOS", "FOSB", "FOSL1", "FOSL2", "FUT4", "G0S2", "GADD45A", "GADD45B", "GCH1", "GEM", "GFPT2", "GPR183", "HBEGF", "HES1", "ICAM1", "ICOSLG", "ID2", "IER2", "IER3", "IER5", "IFIH1", "IFIT2", "IFNGR2", "IL12B", "IL15RA", "IL18", "IL1A", "IL1B", "IL23A", "IL6", "IL6ST", "IL7R", "INHBA", "IRF1", "IRS2", "JAG1", "JUN", "JUNB", "KDM6B", "KLF10", "KLF2", "KLF4", "KLF6", "KLF9", "KYNU", "LAMB3", "LDLR", "LIF", "LITAF", "MAFF", "MAP2K3", "MAP3K8", "MARCKS", "MCL1", "MSC", "MXD1", "MYC", "NAMPT", "NFAT5", "NFE2L2", "NFIL3", "NFKB1", "NFKB2", "NFKBIA", "NFKBIE", "NINJ1", "NR4A1", "NR4A2", "NR4A3", "OLR1", "PANX1", "PDE4B", "PDLIM5", "PER1", "PFKFB3", "PHLDA1", "PHLDA2", "PLAU", "PLAUR", "PLEK", "PLK2", "PLPP3", "PMEPA1", "PNRC1", "PPP1R15A", "PTGER4", "PTGS2", "PTPRE", "PTX3", "RCAN1", "REL", "RELA", "RELB", "RHOB", "RIGI", "RIPK2", "RNF19B", "SAT1", "SDC4", "SERPINB2", "SERPINB8", "SERPINE1", "SGK1", "SIK1", "SLC16A6", "SLC2A3", "SLC2A6", "SMAD3", "SNN", "SOCS3", "SOD2", "SPHK1", "SPSB1", "SQSTM1", "STAT5A", "TANK", "TAP1", "TGIF1", "TIPARP", "TLR2", "TNC", "TNF", "TNFAIP2", "TNFAIP3", "TNFAIP6", "TNFAIP8", "TNFRSF9", "TNFSF9", "TNIP1", "TNIP2", "TRAF1", "TRIB1", "TRIP10", "TSC22D1", "TUBB2A", "VEGFA", "YRDC", "ZBTB10", "ZC3H12A", "ZFP36"),
  HALLMARK_IL6_JAK_STAT3_SIGNALING = c("A2M", "ACVR1B", "ACVRL1", "BAK1", "CBL", "CCL7", "CCR1", "CD14", "CD36", "CD38", "CD44", "CD9", "CNTFR", "CRLF2", "CSF1", "CSF2", "CSF2RA", "CSF2RB", "CSF3R", "CXCL1", "CXCL10", "CXCL11", "CXCL13", "CXCL3", "CXCL9", "DNTT", "EBI3", "FAS", "GRB2", "HAX1", "HMOX1", "IFNAR1", "IFNGR1", "IFNGR2", "IL10RB", "IL12RB1", "IL13RA1", "IL15RA", "IL17RA", "IL17RB", "IL18R1", "IL1B", "IL1R1", "IL1R2", "IL2RA", "IL2RG", "IL3RA", "IL4R", "IL6", "IL6ST", "IL7", "IL9R", "INHBE", "IRF1", "IRF9", "ITGA4", "ITGB3", "JUN", "LEPR", "LTB", "LTBR", "MAP3K8", "MYD88", "OSMR", "PDGFC", "PF4", "PIK3R5", "PIM1", "PLA2G2A", "PTPN1", "PTPN11", "PTPN2", "REG1A", "SOCS1", "SOCS3", "STAM2", "STAT1", "STAT2", "STAT3", "TGFB1", "TLR2", "TNF", "TNFRSF12A", "TNFRSF1A", "TNFRSF1B", "TNFRSF21", "TYK2"),
  HALLMARK_IL2_STAT5_SIGNALING = c("ABCB1", "ADAM19", "AGER", "AHCY", "AHNAK", "AHR", "ALCAM", "AMACR", "ANXA4", "APLP1", "ARL4A", "BATF", "BATF3", "BCL2", "BCL2L1", "BHLHE40", "BMP2", "BMPR2", "CA2", "CAPG", "CAPN3", "CASP3", "CCND2", "CCND3", "CCNE1", "CCR4", "CD44", "CD48", "CD79B", "CD81", "CD83", "CD86", "CDC42SE2", "CDC6", "CDCP1", "CDKN1C", "CISH", "CKAP4", "COCH", "COL6A1", "CSF1", "CSF2", "CST7", "CTLA4", "CTSZ", "CXCL10", "CYFIP1", "DCPS", "DENND5A", "DHRS3", "DRC1", "ECM1", "EEF1AKMT1", "EMP1", "ENO3", "ENPP1", "EOMES", "ETFBKMT", "ETV4", "F2RL2", "FAH", "FGL2", "FLT3LG", "FURIN", "GABARAPL1", "GADD45B", "GALM", "GATA1", "GBP4", "GLIPR2", "GPR65", "GPR83", "GPX4", "GSTO1", "GUCY1B1", "HIPK2", "HK2", "HOPX", "HUWE1", "HYCC2", "ICOS", "IFITM3", "IFNGR1", "IGF1R", "IGF2R", "IKZF2", "IKZF4", "IL10", "IL10RA", "IL13", "IL18R1", "IL1R2", "IL1RL1", "IL2RA", "IL2RB", "IL3RA", "IL4R", "IRF4", "IRF6", "IRF8", "ITGA6", "ITGAE", "ITGAV", "ITIH5", "KLF6", "LCLAT1", "LIF", "LRIG1", "LRRC8C", "LTB", "MAFF", "MAP3K8", "MAP6", "MAPKAPK2", "MUC1", "MXD1", "MYC", "MYO1C", "MYO1E", "NCOA3", "NCS1", "NDRG1", "NFIL3", "NFKBIZ", "NOP2", "NRP1", "NT5E", "ODC1", "P2RX4", "P4HA1", "PDCD2L", "PENK", "PHLDA1", "PHTF2", "PIM1", "PLAGL1", "PLEC", "PLIN2", "PLPP1", "PLSCR1", "PNP", "POU2F1", "PRAF2", "PRKCH", "PRNP", "PTCH1", "PTGER2", "PTH1R", "PTRH2", "PUS1", "RABGAP1L", "RGS16", "RHOB", "RHOH", "RNH1", "RORA", "RRAGD", "S100A1", "SCN9A", "SELL", "SELP", "SERPINB6", "SERPINC1", "SH3BGRL2", "SHE", "SLC1A5", "SLC29A2", "SLC2A3", "SLC39A8", "SMPDL3A", "SNX14", "SNX9", "SOCS1", "SOCS2", "SPP1", "SPRED2", "SPRY4", "ST3GAL4", "SWAP70", "SYNGR2", "SYT11", "TGM2", "TIAM1", "TLR7", "TNFRSF18", "TNFRSF1B", "TNFRSF21", "TNFRSF4", "TNFRSF8", "TNFRSF9", "TNFSF10", "TNFSF11", "TRAF1", "TTC39B", "TWSG1", "UCK2", "UMPS", "WLS", "XBP1"),
  HALLMARK_INTERFERON_GAMMA_RESPONSE = c("ADAR", "APOL6", "ARID5B", "ARL4A", "AUTS2", "B2M", "BANK1", "BATF2", "BPGM", "BST2", "BTG1", "C1R", "C1S", "CASP1", "CASP3", "CASP4", "CASP7", "CASP8", "CCL2", "CCL5", "CCL7", "CD274", "CD38", "CD40", "CD69", "CD74", "CD86", "CDKN1A", "CFB", "CFH", "CIITA", "CMKLR1", "CMPK2", "CMTR1", "CSF2RB", "CXCL10", "CXCL11", "CXCL9", "DDX60", "DHX58", "EIF2AK2", "EIF4E3", "EPSTI1", "FAS", "FCGR1A", "FGL2", "FPR1", "GBP4", "GBP6", "GCH1", "GPR18", "GZMA", "HELZ2", "HERC6", "HIF1A", "HLA-A", "HLA-B", "HLA-DMA", "HLA-DQA1", "HLA-DRB1", "HLA-G", "ICAM1", "IDO1", "IFI27", "IFI30", "IFI35", "IFI44", "IFI44L", "IFIH1", "IFIT1", "IFIT2", "IFIT3", "IFITM2", "IFITM3", "IFNAR2", "IL10RA", "IL15", "IL15RA", "IL18BP", "IL2RB", "IL4R", "IL6", "IL7", "IRF1", "IRF2", "IRF4", "IRF5", "IRF7", "IRF8", "IRF9", "ISG15", "ISG20", "ISOC1", "ITGB7", "JAK2", "KLRK1", "LAP3", "LATS2", "LCP2", "LGALS3BP", "LY6E", "LYSMD2", "MARCHF1", "METTL7B", "MT2A", "MTHFD2", "MVP", "MX1", "MX2", "MYD88", "NAMPT", "NCOA3", "NFKB1", "NFKBIA", "NLRC5", "NMI", "NOD1", "NUP93", "OAS2", "OAS3", "OASL", "OGFR", "P2RY14", "PARP12", "PARP14", "PDE4B", "PELI1", "PFKP", "PIM1", "PLA2G4A", "PLSCR1", "PML", "PNP", "PNPT1", "PSMA2", "PSMA3", "PSMB10", "PSMB2", "PSMB8", "PSMB9", "PSME1", "PSME2", "PTGS2", "PTPN1", "PTPN2", "PTPN6", "RAPGEF6", "RBCK1", "RIGI", "RIPK1", "RIPK2", "RNF213", "RNF31", "RSAD2", "RTP4", "SAMD9L", "SAMHD1", "SECTM1", "SELP", "SERPING1", "SLAMF7", "SLC25A28", "SOCS1", "SOCS3", "SOD2", "SP110", "SPPL2A", "SRI", "SSPN", "ST3GAL5", "ST8SIA4", "STAT1", "STAT2", "STAT3", "STAT4", "TAP1", "TAPBP", "TDRD7", "TNFAIP2", "TNFAIP3", "TNFAIP6", "TNFSF10", "TOR1B", "TRAFD1", "TRIM14", "TRIM21", "TRIM25", "TRIM26", "TXNIP", "UBE2L6", "UPP1", "USP18", "VAMP5", "VAMP8", "VCAM1", "WARS1", "XAF1", "XCL1", "ZBP1", "ZNFX1"),
  HALLMARK_PI3K_AKT_MTOR_SIGNALING = c("ACACA", "ACTR2", "ACTR3", "ADCY2", "AKT1", "AKT1S1", "AP2M1", "ARF1", "ARHGDIA", "ARPC3", "ATF1", "CAB39", "CAB39L", "CALR", "CAMK4", "CDK1", "CDK2", "CDK4", "CDKN1A", "CDKN1B", "CFL1", "CLTC", "CSNK2B", "CXCR4", "DAPP1", "DDIT3", "DUSP3", "E2F1", "ECSIT", "EGFR", "EIF4E", "FASLG", "FGF17", "FGF22", "FGF6", "GNA14", "GNGT1", "GRB2", "GRK2", "GSK3B", "HRAS", "HSP90B1", "IL2RG", "IL4", "IRAK4", "ITPR2", "LCK", "MAP2K3", "MAP2K6", "MAP3K7", "MAPK1", "MAPK10", "MAPK8", "MAPK9", "MAPKAP1", "MKNK1", "MKNK2", "MYD88", "NCK1", "NFKBIB", "NGF", "NOD1", "PAK4", "PDK1", "PFN1", "PIK3R3", "PIKFYVE", "PIN1", "PITX2", "PLA2G12A", "PLCB1", "PLCG1", "PPP1CA", "PPP2R1B", "PRKAA2", "PRKAG1", "PRKAR2A", "PRKCB", "PTEN", "PTPN11", "RAC1", "RAF1", "RALB", "RIPK1", "RIT1", "RPS6KA1", "RPS6KA3", "RPTOR", "SFN", "SLA", "SLC2A1", "SMAD2", "SQSTM1", "STAT2", "TBK1", "THEM4", "TIAM1", "TNFRSF1A", "TRAF2", "TRIB3", "TSC2", "UBE2D3", "UBE2N", "VAV3", "YWHAB")
)

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


## Statistical analyses
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

gsva_df <- as.data.frame(t(gsva_scores))

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


## Row normalize and then create heatmap
# Define row-wise min-max normalization function
min_max_scale <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep(0, length(x)))  # Avoid division by zero
  return(2 * (x - rng[1]) / (rng[2] - rng[1]) - 1)
}

# Apply to each row of the matrix
gsva_scores_rescaled <- apply(gsva_scores, 1, min_max_scale)


## Create the heatmap
library(ComplexHeatmap)
library(circlize)
library(grid)

# Select pathways of interest then rename
selected_pathways <- c(
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB", 
  "HALLMARK_IL6_JAK_STAT3_SIGNALING", 
  "HALLMARK_IL2_STAT5_SIGNALING", 
  "HALLMARK_INTERFERON_GAMMA_RESPONSE", 
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING"
)

geneset_names <- c(
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB" = "Hallmark TNFa Signaling via NFKB",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING" = "Hallmark IL6-JAK-STAT3 Signaling",
  "HALLMARK_IL2_STAT5_SIGNALING" = "Hallmark IL2-STAT5 Signaling",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE" = "Hallmark IFNg Response",
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING" = "Hallmark PI3K-AKT-mTOR Signaling"
)

# Subset and rename matrix rows
gsva_scores_rescaled <- t(gsva_scores_rescaled)

gsva_scores_subset <- gsva_scores_rescaled[selected_pathways, ]

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
column_annotation1 <- HeatmapAnnotation(
  Group = condition_factor,
  col = list(Group = group_colors),
  show_legend = TRUE,
  annotation_legend_param = list(
    title = "Condition",
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontsize = 10)
  )
)

col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

ha_col <- HeatmapAnnotation(
  Condition = condition_factor,
  col = list(Condition = group_colors),
  annotation_name_side = "left",
  show_legend = TRUE
)

row_ha <- rowAnnotation(Pathway = geneset_names)

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
library(dplyr)
library(tidyr)
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


par(mar = c(4, 4, 4, 4), xpd = TRUE)
par(mar = c(10, 10, 10, 10))   # huge margins push labels off the plot

par(mar = c(12, 12, 12, 12), xpd = TRUE)

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
  vlcex = 0.0001,     # shrinks pathway labels to invisible
  cex = 1.2
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

## Train XGBoost on all patients
library(xgboost)
library(dplyr)
library(ggplot2)

set.seed(3158)
dall <- xgb.DMatrix(data = features_matrix, label = labels_vector)
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

# Train model 
xgb_model <- xgb.train(
params = params,
data = dall,
nrounds = 500,
verbose = 1
)

## Top Proteins by Gain
importance <- xgb.importance(feature_names = feature_names, model = xgb_model)

library(dplyr)

top10_MISC <- importance %>%
  as_tibble() %>%
  arrange(desc(Gain)) %>%
  dplyr::slice(1:10)


library(ggplot2)
library(Ckmeans.1d.dp)

xgb.ggplot.importance(importance_matrix = importance, top_n = 10, n_clusters = 1) +
  geom_bar(stat = "identity", fill = "#3399FF") +  # Override default bar color
  labs(title = "Top 10 Feature Importances", x = "Feature", y = "Importance (Gain)") +
  theme_minimal(base_size = 18) +  # Increase base font size
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    plot.title = element_text(size = 18, face = "bold")
  )


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

## Train XGBoost on all patients
library(xgboost)
library(dplyr)
library(ggplot2)
library(caret)

set.seed(3158)
dall1 <- xgb.DMatrix(data=features_matrix1, label=labels_vector1)
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

# Train model
xgb_model1 <- xgb.train(
  params = params,
  data = dall1,
  nrounds = 500,
  verbose = 1
)


## Top Proteins by Gain
importance1 <- xgb.importance(feature_names = features_names1, model = xgb_model1)

library(dplyr)

top10_CRS <- importance1 %>%
  as_tibble() %>%
  arrange(desc(Gain)) %>%
  dplyr::slice(1:10)


library(ggplot2)
library(Ckmeans.1d.dp)

xgb.ggplot.importance(importance_matrix = importance1, top_n = 10, n_clusters = 1) +
  geom_bar(stat = "identity", fill = "#3399FF") +  # Override default bar color
  labs(title = "Top 10 Feature Importances", x = "Feature", y = "Importance (Gain)") +
  theme_minimal(base_size = 18) +  # Increase base font size
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    plot.title = element_text(size = 18, face = "bold")
  )