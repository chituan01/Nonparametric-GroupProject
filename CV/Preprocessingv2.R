#Preprocessing 2026 - Nov 
library(dplyr) #dataframe playtime
library(vegan) # mrpp
source("utils.R")

#Load data mri + metadata 2026 v3
mrifull <- read_csv("mrifullv3.csv")

# Volume adjustment/correction ----
volume_indices <- c(2:10, 12:46, 52:167)
target_vols <- colnames(mrifull)[volume_indices]

mrifull_adjusted <- residualize_volumes(
  train_df = mrifull,
  test_df = mrifull,
  vol_cols = target_vols
)

#Clear worspace
rm(list = setdiff(ls(), "mrifull_adjusted"))

# MRPP ----
features <- scale(mrifull_adjusted[, 197:ncol(mrifull_adjusted)])
mrpp <- mrpp(d = features, 
             grouping = mrifull_adjusted$`Research Group`, 
             distance = "euclidean",  # can also use other distance metrics
             permutations = 9999)

# View global significance.
print(mrpp)
# We conclude there is a difference: Significance of delta: 1e-04 
# But Chance corrected within-group agreement A: 0.022
# => The features, collectively, cannot effectively separate any of the five groups well

# Kruskal Wallis ----
features <- colnames(features)

# Initialize a vector to store p-values
pvals <- numeric(length(features))

# Loop through features
for (i in seq_along(features)) {
  # Run Kruskal-Wallis test for each feature
  test <- kruskal.test(mrifull_adjusted[[features[i]]] ~  mrifull_adjusted$`Research Group`)
  pvals[i] <- test$p.value
}

# Combine results
results <- data.frame(feature = features, p_value = pvals)
print(results)

# False discovery rate correction
results$adj_pval <- p.adjust(results$p_value, method = "fdr") 
significant_features <- results[results$adj_pval < 0.05, ]

mrifull_adjusted$`Research Group` <- factor(
  mrifull_adjusted$`Research Group`,
  levels = c("CN", "EMCI", "MCI", "LMCI", "AD")
)


for (feat in significant_features$feature) {
  # 1. Assign the plot to a variable (e.g., 'p')
  p <- ggplot(mrifull_adjusted, aes(x = Group, y = .data[[feat]])) +
    geom_boxplot() +
    labs(title = paste("Feature:", feat)) +
    theme_minimal()
  
  print(p)
}


#Li and Ji's method ----
# Li and Ji's method
cor_mat <- calc_correlation(mrifull_adjusted[, 197:ncol(mrifull_adjusted)])
m_eff <- calc_meff(cor_mat)
alpha_eff <- 0.05/round(m_eff)

