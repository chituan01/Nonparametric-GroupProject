# Load libraries and data -------------------------------------------------
library(vegan)
library(ggplot2)
library(ggrepel)
library(dunn.test) 
library(dplyr)
library(tidyr)
library(stringr)

#Load DF, remove sex, age, weight
cleanMRITotal <- read.csv("cleanMRITotal_final.csv")
cleanMRITotal <- cleanMRITotal[,c(-1,-3,-4,-5)] # drop SubjectID, Sex, Age, Weight
head(cleanMRITotal)
cleanMRITotal <- cleanMRITotal %>%
  rename_with(~ str_remove_all(.x, "(?i)total.|\\.Volume_mm3"))

## Not merge group ---------------------------------------------------------
features <- scale(cleanMRITotal[, -1])
mrpp <- mrpp(d = features, 
             grouping = cleanMRITotal$Group, 
             distance = "euclidean",  # can also use other distance metrics
             permutations = 9999)

# View global significance.
print(mrpp)
# We conclude there is a difference: Significance of delta: 1e-04 
# But Chance corrected within-group agreement A: 0.0131 
# => The features, collectively, cannot effectively separate any of the five groups well

## Now we do it for each feature
## 
features <- colnames(features)

# Initialize a vector to store p-values
pvals <- numeric(length(features))

# Loop through features
for (i in seq_along(features)) {
  # Run Kruskal-Wallis test for each feature
  test <- kruskal.test(cleanMRITotal[[features[i]]] ~ cleanMRITotal$Group)
  pvals[i] <- test$p.value
}

# Combine results
results <- data.frame(feature = features, p_value = pvals)
print(results)

# False discovery rate correction
results$adj_pval <- p.adjust(results$p_value, method = "fdr") 
significant_features <- results[results$adj_pval < 0.05, ]

cleanMRITotal$Group <- factor(
  cleanMRITotal$Group,
  levels = c("CN", "EMCI", "MCI", "LMCI", "AD")
)


for (feat in significant_features$feature) {
  # 1. Assign the plot to a variable (e.g., 'p')
  p <- ggplot(cleanMRITotal, aes(x = Group, y = .data[[feat]])) +
    geom_boxplot() +
    labs(title = paste("Feature:", feat)) +
    theme_minimal()
  
  print(p)
}
# 3 main ideas to support groupping
# - From boxplot: The visual difference between these three specific classes is minimal.
# - MRPP Effect Size: Chance corrected within-group agreement A: 0.0131 
# - Clinical Rationale
#CN vs MCI and CN vs AD----------
## Merge group ------------------------------------------------------------
cleanMRITotal$Group[cleanMRITotal$Group %in% c("EMCI", "MCI", "LMCI")] <- "MCI"
cleanMRITotal$Group <- factor(cleanMRITotal$Group)
head(cleanMRITotal)

# Re-run mrpp to confirm
mrpp <- mrpp(d = scale(cleanMRITotal[, -1]), 
             grouping = cleanMRITotal$Group, 
             distance = "euclidean",  # can also use other distance metrics
             permutations = 9999)

print(mrpp)
# p-value is the same
# A is smaller!!! 0.01004 < 0.0131!!!
# => - Still group the classes because the features couldn't effectively separate them 
# in the first place (as shown by the boxplot overlap)
#    - Confirm the fact that the features, as a whole, are weak separators for these categories.
# => Find significant features to distinguish

# First scenario (Dunn's Test) ----------------------------------------------------------
reference_group <- "CN"
results_list <- list()
num_MRITotal <- cleanMRITotal[,-1]
features <-  significant_features$feature

### Functions -------------------------------------------------------------------
calc_correlation <- function(df){
  # Remove Group column if exists
  roi_data <- df 
  cor_matrix <- cor(roi_data, use = "pairwise.complete.obs") # method: spearman
  return(cor_matrix)
}

calc_meff <- function(cor_matrix){
  eigenvalues <- eigen(cor_matrix, symmetric = TRUE)$values
  # Li & Ji formula
  m_eff <- sum(ifelse(eigenvalues >= 1, 1, eigenvalues))
  return(m_eff) 
}

# Li and Ji's method
cor_mat <- calc_correlation(num_MRITotal)
m_eff <- calc_meff(cor_mat)
alpha_eff <- 0.05/round(m_eff)

# Loop through features to perform post-hoc tests and calculate effect sizes
for (feat in features) {
  # 1. Perform Dunn's Test (Post-hoc for Kruskal-Wallis)
  #    method = "none" gives raw p-values for us to compare against alpha_eff
  dunn_res <- dunn.test(x = cleanMRITotal[[feat]], 
                        g = cleanMRITotal$Group, 
                        method = "none", # Use raw p-value
                        table = FALSE,
                        altp = FALSE)
  
  # 2. Extract P-values and match with comparison pairs
  p_vals_df <- data.frame(
    comparison = dunn_res$comparisons,
    p_value = dunn_res$P
  )
  
  # 3. Compute mean per group
  group_means <- tapply(cleanMRITotal[[feat]], cleanMRITotal$Group, mean)
  
  # 4. Process all comparisons against the reference group
  for (i in 1:nrow(p_vals_df)) {
    comp_str <- p_vals_df[i, "comparison"] # e.g., "AD - CN"
    p_val <- p_vals_df[i, "p_value"]
    
    groups <- unlist(strsplit(comp_str, " - "))
    grp1 <- groups[1]
    grp2 <- groups[2]
    
    # We only care about comparisons against the reference group ("CN")
    if (grp1 == reference_group | grp2 == reference_group) {
      test_group <- ifelse(grp1 == reference_group, grp2, grp1)
      
      # Effect size = Test Group Mean - Reference Group Mean
      effect_size <- group_means[test_group] - group_means[reference_group]
      
      # Store results using the P-value from the specific pairwise comparison
      results_list <- append(results_list, list(
        data.frame(
          feature = feat, 
          group = test_group, 
          p_value = p_val,
          effect_size = effect_size
        )
      ))
    }
  }
}

# Final data frame
temp_df <- bind_rows(results_list)
temp_df$is_significant <- temp_df$p_value < alpha_eff

# Extract final significant features
sig_AD <- subset(temp_df, group == "AD" & is_significant == TRUE)
sig_MCI <- subset(temp_df, group == "MCI" & is_significant == TRUE)

print("--- Final Significant Features (Corrected) ---")
print("Features Significant for AD vs. CN:")
print(unique(sig_AD$feature))

print("Features Significant for MCI vs. CN:")
print(unique(sig_MCI$feature))

# Second scenario: Wilcoxon Signed-Rank Test (Mean-centering)  -------------------------------------------------------------
reference_group <- "CN"
results_list <- list()
wilcox_results_list <- list() 
num_MRITotal <- cleanMRITotal[, -1] 
features <- significant_features$feature

mean_vec <- colMeans(cleanMRITotal[cleanMRITotal$Group == reference_group, features])
centered_MRITotal <- cleanMRITotal 
centered_MRITotal[, features] <- sweep(
  x = centered_MRITotal[, features], # Only select the numeric feature columns (X)
  MARGIN = 2,                    # Apply operation column-wise
  STATS = mean_vec,              # The vector of CN means
  FUN = "-"                      # Subtraction function
)

# 3. Wilcoxon Signed-Rank Test Loop 
for (feat in features) {
  
  # Loop through groups to be tested against 0 (AD and MCI)
  for (test_group in setdiff(levels(centered_MRITotal$Group), reference_group)) {
    
    # Isolate the centered data for the current test group
    test_data_centered <- centered_MRITotal[[feat]][centered_MRITotal$Group == test_group]
    
    # Remove NAs for accurate median and test calculation
    test_data_centered <- na.omit(test_data_centered)
    
    # Skip if insufficient data remains
    if (length(test_data_centered) < 10) next 
    
    # Perform Wilcoxon Signed-Rank Test (against median mu = 0)
    wilcox_res <- wilcox.test(test_data_centered, mu = 0, exact = FALSE) 
    
    # Calculate the Effect Size (Median Deviation from 0)
    median_deviation <- median(test_data_centered)
    
    # Store results
    wilcox_results_list <- append(wilcox_results_list, list(
      data.frame(
        feature = feat, 
        group = test_group, 
        p_value = wilcox_res$p.value, 
        effect_size = median_deviation # Median deviation from CN mean (0)
      )
    ))
  }
}

# Final data frame for analysis 
temp_wilcox <- bind_rows(wilcox_results_list)
temp_wilcox$is_significant <- temp_wilcox$p_value < alpha_eff

# --- Final Filtering and Output ---
sig_AD_wilcox <- subset(temp_wilcox, group == "AD" & is_significant == TRUE)
sig_MCI_wilcox <- subset(temp_wilcox, group == "MCI" & is_significant == TRUE)

print("--- Final Significant Features (Wilcoxon Signed-Rank Test) ---")
print("Features Significant for AD (Median different from CN Baseline):")
print(sig_AD_wilcox$feature)

print("Features Significant for MCI (Median different from CN Baseline):")
#print(sig_MCI_wilcox[order(sig_MCI_wilcox$p_value), c("feature", "p_value", "effect_size")])
print(sig_MCI_wilcox$feature)

#MCI vs AD----------
# Subset: Only MCI and AD patients
data_stage <- subset(cleanMRITotal, Group %in% c("MCI", "AD"))
data_stage$Group <- factor(data_stage$Group, levels = c("MCI", "AD")) 

# 2. REMOVAL OF CONFOUNDING VARIABLES (Sex, Age, Weight)
# Define columns to remove (ensure names match your dataset)
vars_to_remove <- c("sex", "age", "Weight") 

# Remove them from the main dataset for this step
data_stage <- data_stage[, !(names(data_stage) %in% vars_to_remove)]

cat("Removed variables:", paste(vars_to_remove, collapse=", "), "\n")
cat("Remaining columns (first 5):", head(colnames(data_stage), 5), "\n")

# 3. Numeric matrix preparation for correlations
# Remove 'Group' column to keep only anatomical features
num_Stage <- data_stage[, !(names(data_stage) %in% c("Group"))]
features <- colnames(num_Stage)


# ==============================================================================
# M_eff CALCULATION (Specific for MCI\AD)
# ==============================================================================

calc_correlation <- function(df){
  cor_matrix <- cor(df, use = "pairwise.complete.obs") 
  return(cor_matrix)
}

calc_meff <- function(cor_matrix){
  eigenvalues <- eigen(cor_matrix, symmetric = TRUE)$values
  m_eff <- sum(ifelse(eigenvalues >= 1, 1, eigenvalues))
  return(m_eff) 
}

cor_mat_stage <- calc_correlation(num_Stage)
m_eff_stage <- calc_meff(cor_mat_stage)
alpha_eff_stage <- 0.05 / round(m_eff_stage)

cat("\n--- Statistical Parameters Step 2 ---\n")
cat("Number of Anatomical Features:", length(features), "\n")
cat("M_eff (Independent Tests):", round(m_eff_stage, 2), "\n")
cat("Corrected Alpha (P-value threshold):", format(alpha_eff_stage, scientific=TRUE), "\n")


# ==============================================================================
# SCENARIO 1: DUNN'S TEST (MCI vs AD)
# ==============================================================================
cat("\n--- Execution Scenario 1: Dunn (MCI vs AD) ---\n")
reference_group <- "MCI"
results_list <- list()

for (feat in features) {
  # Dunn's Test
  dunn_res <- dunn.test(x = data_stage[[feat]], 
                        g = data_stage$Group, 
                        method = "none", 
                        kw = FALSE,
                        table = FALSE,
                        altp = FALSE)
  
  p_vals_df <- data.frame(comparison = dunn_res$comparisons, p_value = dunn_res$P)
  group_means <- tapply(data_stage[[feat]], data_stage$Group, mean)
  
  for (i in 1:nrow(p_vals_df)) {
    comp_str <- p_vals_df[i, "comparison"]
    p_val <- p_vals_df[i, "p_value"]
    
    groups <- unlist(strsplit(comp_str, " - "))
    grp1 <- groups[1]
    grp2 <- groups[2]
    
    if (grp1 == reference_group | grp2 == reference_group) {
      test_group <- ifelse(grp1 == reference_group, grp2, grp1) # This will be AD
      effect_size <- group_means[test_group] - group_means[reference_group]
      
      results_list <- append(results_list, list(
        data.frame(feature = feat, group = test_group, p_value = p_val, effect_size = effect_size)
      ))
    }
  }
}

volcano_df_dunn <- bind_rows(results_list)
volcano_df_dunn$is_significant <- volcano_df_dunn$p_value < alpha_eff_stage

# Extract Dunn Features
sig_AD_vs_MCI_Dunn <- subset(volcano_df_dunn, is_significant == TRUE)
features_Dunn_Step2 <- unique(sig_AD_vs_MCI_Dunn$feature)

cat("Significant Features (DUNN):", length(features_Dunn_Step2), "\n")
print(features_Dunn_Step2)


cat("\n--- Execution Scenario 2: Centered Wilcoxon (Baseline MCI) ---\n")

wilcox_results_list <- list() 

# 1. Center on MCI group means
mean_vec <- colMeans(data_stage[data_stage$Group == "MCI", features])
centered_Stage <- data_stage 
centered_Stage[, features] <- sweep(x = centered_Stage[, features], MARGIN = 2, STATS = mean_vec, FUN = "-")

# 2. Wilcoxon Loop
for (feat in features) {
  test_group <- "AD"
  test_data_centered <- centered_Stage[[feat]][centered_Stage$Group == test_group]
  test_data_centered <- na.omit(test_data_centered)
  
  if (length(test_data_centered) < 10) next 
  
  # Test against 0 (which represents the MCI mean)
  wilcox_res <- wilcox.test(test_data_centered, mu = 0, exact = FALSE) 
  median_deviation <- median(test_data_centered)
  
  wilcox_results_list <- append(wilcox_results_list, list(
    data.frame(feature = feat, group = test_group, p_value = wilcox_res$p.value, effect_size = median_deviation)
  ))
}

volcano_df_wilcox <- bind_rows(wilcox_results_list)
volcano_df_wilcox$is_significant <- volcano_df_wilcox$p_value < alpha_eff_stage


# Extract Wilcoxon Features
sig_AD_vs_MCI_Wilcox <- subset(volcano_df_wilcox, is_significant == TRUE)
features_Wilcox_Step2 <- unique(sig_AD_vs_MCI_Wilcox$feature)

cat("Significant Features (WILCOXON):", length(features_Wilcox_Step2), "\n")
print(features_Wilcox_Step2)

