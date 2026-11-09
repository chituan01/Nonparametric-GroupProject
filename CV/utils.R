# ==============================================================================
# File: utils.R
# Description: Helper function for FastSurfer volume residualization
# ==============================================================================

residualize_volumes <- function(train_df, test_df = train_df, vol_cols, 
                                icv_col = "eTIV", 
                                age_col = "Subject Age", group_col = "Research Group",
                                cn_label = "CN") {
  
  # Filter training controls
  cn_train <- train_df %>% filter(.data[[group_col]] == cn_label)
  
  if (nrow(cn_train) == 0) {
    stop(paste("ERROR: No subjects found where", group_col, "==", cn_label))
  }
  
  # Reference baseline means from training CN cohort
  mean_icv <- mean(cn_train[[icv_col]], na.rm = TRUE)
  mean_age <- mean(cn_train[[age_col]], na.rm = TRUE)
  
  # Create a copy of the target dataset to hold adjusted values
  adjusted_df <- test_df
  
  for (col in vol_cols) {
    # Fit regression model on CN train cohort only
    formula_str <- as.formula(paste0("`", col, "` ~ `", icv_col, "` + `", age_col, "`"))
    #print(formula_str)
    fit <- lm(formula_str, data = cn_train, na.action = na.exclude)
    
    # Extract coefficients
    alpha    <- coef(fit)[1]
    beta_icv <- coef(fit)[2]
    beta_age <- coef(fit)[3]
    
    #print(beta_icv)
    #print(beta_age)
    
    # Mean-Centered Version V_adj = V_raw - B1 (icv-mean_icv) + B2 (age-mean_age)
    adj_col_name <- paste0(col, "_adj")
    adjusted_df[[adj_col_name]] <- test_df[[col]] - 
      (beta_icv * (test_df[[icv_col]] - mean_icv)) - 
      (beta_age * (test_df[[age_col]] - mean_age))
    # Residual Version V_adj = V_raw - (alpha  + B1 icv + B2 age)
    adjusted_df[[adj_col_name]] <- test_df[[col]] - 
      (alpha + beta_icv * test_df[[icv_col]] + 
         beta_age * test_df[[age_col]])
  }
  
  return(adjusted_df)
} 
  
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
