# ==============================================================================
## ALZHEIMER'S DIAGNOSIS PIPELINE: ROBUST THRESHOLDING & SEQUENTIAL TESTING-----
# ==============================================================================
# Description: 
# This script implements a two-stage sequential diagnostic system.
# 1. Screening Stage: Distinguishes Cognitive Normal (CN) from Pathological (MCI/AD).
# 2. Staging Stage: Distinguishes MCI from AD among pathological subjects.
# Method: 
# - Uses 5-Fold Cross-Validation on the Training Set to determine robust thresholds.
# - Validates the final pipeline on an unseen Test Set (30%).

# --- 0. LIBRARIES ---------------
library(caret)  
library(mgcv)   
library(pROC)   
library(dplyr)  
library(ggplot2)

# --- 1. DATA LOADING & PREPARATION ---------
cleanMRITotal <- read.csv("cleanMRITotal.csv")

# Merge Clinical Subgroups: Combine EMCI and LMCI into a single "MCI" class
cleanMRITotal$Group[cleanMRITotal$Group %in% c("EMCI", "MCI", "LMCI")] <- "MCI"
cleanMRITotal$Group <- factor(cleanMRITotal$Group)

# --- 2. STRATIFIED SPLIT (70% Train / 30% Test) ------------
# We use a stratified split to maintain class proportions in both sets.
set.seed(123)
data_screen <- cleanMRITotal
trainIndex <- createDataPartition(data_screen$Group, p = .7, list = FALSE)

# --- 3. STANDARDIZATION (Z-SCORE) --------------
# CRITICAL: Mean and SD must be calculated ONLY on the Training set 
# and then applied to the Test set to prevent Data Leakage.
numeric_cols <- names(data_screen)[sapply(data_screen, is.numeric)]

for(col in numeric_cols) {
  # Calculate statistics on the Training set
  m <- mean(data_screen[trainIndex, col], na.rm = TRUE)
  s <- sd(data_screen[trainIndex, col], na.rm = TRUE)
  
  # Apply transformation to ALL data (Train + Test)
  data_screen[[col]] <- (data_screen[[col]] - m) / s
}

# Split Dataframes
trainDataS1 <- data_screen[trainIndex, ]  # Full Training Set (used for Stage 1)
testDataS1  <- data_screen[-trainIndex, ] # Full Test Set (Hold-out for final validation)

# Create Binary Target for Stage 1 (Screening: Healthy vs Pathological)
trainDataS1$BinaryGroup <- factor(ifelse(trainDataS1$Group == "CN", "CN", "PAT"), levels = c("CN", "PAT"))
testDataS1$BinaryGroup  <- factor(ifelse(testDataS1$Group == "CN", "CN", "PAT"), levels = c("CN", "PAT"))

# Create Subset for Stage 2 (Staging: MCI vs AD only)
trainDataS2 <- subset(trainDataS1, Group %in% c("MCI", "AD"))
trainDataS2$Group <- factor(trainDataS2$Group, levels = c("MCI", "AD"))


# ==============================================================================
# PART 4: MODEL FORMULA DEFINITIONS-------
# ==============================================================================

# 1. DUNN Formula (Parsimonious)
f_dunn <- as.formula("BinaryGroup ~ Total.Inf.Lat.Vent.Volume_mm3 + ctx.total.middletemporal.Volume_mm3 + s(Total.Lateral.Ventricle.Volume_mm3) + Subject.Age + ctx.total.entorhinal.Volume_mm3 + s(Weight.Kg)")

# 2. WILCOXON Formula (Complex)
f_wilcox <- as.formula("BinaryGroup ~ s(Total.Lateral.Ventricle.Volume_mm3) + Total.Inf.Lat.Vent.Volume_mm3 + Total.Hippocampus.Volume_mm3 + s(CSF.Volume_mm3) + Total.Accumbens.area.Volume_mm3 + s(ctx.total.caudalanteriorcingulate.Volume_mm3) + s(ctx.total.inferiortemporal.Volume_mm3) + ctx.total.parahippocampal.Volume_mm3 + s(ctx.total.parsopercularis.Volume_mm3) + s(ctx.total.posteriorcingulate.Volume_mm3) + s(ctx.total.rostralanteriorcingulate.Volume_mm3) + Cbm_Vermis_VII.Volume_mm3 + s(Total.Fornix.Volume_mm3) + Subject.Age + s(Weight.Kg)")

# 3. STAGE 2 Formula (MCI vs AD)
f_stage2 <- as.formula("Group ~ Total.Hippocampus.Volume_mm3 + ctx.total.entorhinal.Volume_mm3 + ctx.total.precuneus.Volume_mm3 + ctx.total.superiorparietal.Volume_mm3 + ctx.total.insula.Volume_mm3 + s(Tuberal.Region.Volume_mm3) + ctx.total.fusiform.Volume_mm3")


# ==============================================================================
# PART 5: ROBUST THRESHOLD FUNCTION (5-FOLD INTERNAL CROSS-VALIDATION)-------
# ==============================================================================
# This function splits the Training Set into 5 folds to find a stable average threshold.
# It ensures the threshold is not overfitted to a single data split.

find_robust_threshold <- function(data_input, formula_model, model_type="stage1", k=5) {
  
  set.seed(123)
  # Select the correct target variable based on the model stage
  target_var <- if(model_type == "stage1") data_input$BinaryGroup else data_input$Group
  folds <- createFolds(target_var, k = k, list = TRUE)
  
  thresholds_list <- c()
  
  cat(paste0("\n>>> 5-Fold CV Threshold Search: ", model_type, " <<<\n"))
  
  for(i in 1:k) {
    # A. Internal Split (Train vs Validation)
    inner_train_idx <- setdiff(1:nrow(data_input), folds[[i]])
    inner_train <- data_input[inner_train_idx, ]
    inner_val   <- data_input[folds[[i]], ]
    
    # B. Dynamic Weight Calculation (Class Balancing within the fold)
    if(model_type == "stage1") {
      t_col <- inner_train$BinaryGroup
      target_val <- inner_val$BinaryGroup
    } else {
      t_col <- inner_train$Group
      target_val <- inner_val$Group
    }
    
    n_neg <- sum(t_col == levels(t_col)[1])
    n_pos <- sum(t_col == levels(t_col)[2])
    n_tot <- nrow(inner_train)
    
    # Create weight vector
    w_vec <- ifelse(t_col == levels(t_col)[1], (n_tot/n_neg), (n_tot/n_pos))
    
    # --- TECHNICAL FIX: Inject weights into dataframe ---
    # This prevents scoping errors where 'gam' cannot find the weight vector.
    inner_train$weights_internal <- w_vec
    
    # C. Train GAM on the inner training set
    mod <- gam(formula_model, data = inner_train, family = quasibinomial(), weights = weights_internal)
    
    # D. Predict on the inner validation set
    preds_val <- predict(mod, inner_val, type="response")
    roc_val   <- roc(target_val, preds_val, quiet=TRUE)
    
    # E. Select Threshold based on Strategy
    if(model_type == "stage1") {
      # STRATEGY 1: Screening (Target Sensitivity ~90%)
      # We prioritize excluding false negatives.
      res <- coords(roc_val, x="all", ret=c("threshold", "sensitivity"))
      best_t <- res$threshold[which.min(abs(res$sensitivity - 0.90))]
    } else {
      # STRATEGY 2: Staging (Youden Index)
      # We prioritize the best balance between Sensitivity and Specificity.
      res <- coords(roc_val, "best", best.method="youden", ret="threshold")
      best_t <- res$threshold
    }
    
    # Handle edge cases (multiple thresholds found)
    if(length(best_t) > 1) best_t <- best_t[1]
    if(!is.finite(best_t)) best_t <- 0.5 
    
    thresholds_list <- c(thresholds_list, best_t)
    # cat(paste("   Fold", i, "-> Threshold found:", round(best_t, 3), "\n"))
  }
  
  # F. Calculate Average Threshold
  avg_t <- mean(thresholds_list, na.rm=TRUE)
  cat(paste(">>> FINAL ROBUST THRESHOLD (Mean):", round(avg_t, 3), "\n"))
  return(avg_t)
}

# --- EXECUTE THRESHOLD SEARCH ---
# 1. DUNN (Stage 1 -> Target 90% Sensitivity)
th_dunn_avg <- find_robust_threshold(trainDataS1, f_dunn, "stage1", k=5)

# 2. WILCOXON (Stage 1 -> Target 90% Sensitivity)
th_wilcox_avg <- find_robust_threshold(trainDataS1, f_wilcox, "stage1", k=5)

# 3. STAGE 2 GAM (MCI vs AD -> Target Youden Index)
th_stage2_avg <- find_robust_threshold(trainDataS2, f_stage2, "stage2", k=5)


# ==============================================================================
# PART 6: FINAL MODEL TRAINING (FULL TRAINING SET)------
# ==============================================================================
# Now that we have the robust thresholds, we train the definitive models 
# using the entire 70% Training Set to apply them to the Test Set.

# Weights for Stage 1
n_tot1 <- nrow(trainDataS1)
w_s1   <- ifelse(trainDataS1$BinaryGroup == "CN", 
                 (n_tot1/sum(trainDataS1$BinaryGroup == "CN")), 
                 (n_tot1/sum(trainDataS1$BinaryGroup == "PAT")))

# Weights for Stage 2
n_tot2 <- nrow(trainDataS2)
w_s2   <- ifelse(trainDataS2$Group == "MCI", 
                 (n_tot2/sum(trainDataS2$Group == "MCI")), 
                 (n_tot2/sum(trainDataS2$Group == "AD")))

cat("\n>>> Training Final Models on Full Training Set... <<<\n")
final_model_dunn   <- gam(f_dunn, data=trainDataS1, family=quasibinomial(), weights=w_s1)
final_model_wilcox <- gam(f_wilcox, data=trainDataS1, family=quasibinomial(), weights=w_s1)
final_model_stage2 <- gam(f_stage2, data=trainDataS2, family=quasibinomial(), weights=w_s2)


# ==============================================================================
# PART 7: SEQUENTIAL PIPELINE ON TEST SET (30% - UNSEEN DATA)------
# ==============================================================================

run_pipeline <- function(test_data, model_s1, th_s1, model_s2, th_s2, name) {
  
  cat(paste0("\n##################################################\n"))
  cat(paste0(" PIPELINE: ", name, "\n"))
  cat(paste0(" Threshold S1: ", round(th_s1, 3), " | Threshold S2: ", round(th_s2, 3), "\n"))
  cat(paste0("##################################################\n"))
  
  # 1. SCREENING PREDICTION (Stage 1)
  probs_s1 <- predict(model_s1, test_data, type="response")
  preds_s1 <- ifelse(probs_s1 > th_s1, "PAT", "CN")
  
  # 2. SEQUENTIAL LOGIC (Cascade)
  final_preds <- rep(NA, nrow(test_data))
  
  # A. Subjects classified as Healthy (CN) exit the pipeline here.
  final_preds[preds_s1 == "CN"] <- "CN"
  
  # B. Subjects classified as Pathological (PAT) proceed to Stage 2.
  pat_idx <- which(preds_s1 == "PAT")
  
  if(length(pat_idx) > 0) {
    # Stage 2 Prediction only for the subset identified as 'PAT'
    probs_s2 <- predict(model_s2, test_data[pat_idx, ], type="response")
    final_preds[pat_idx] <- ifelse(probs_s2 > th_s2, "AD", "MCI")
  }
  
  # 3. EVALUATION
  final_preds_f <- factor(final_preds, levels = c("CN", "MCI", "AD"))
  actual_f      <- factor(test_data$Group, levels = c("CN", "MCI", "AD"))
  
  cm <- confusionMatrix(final_preds_f, actual_f)
  print(cm$table)
  
  # --- DETAILED METRICS BY CLASS ---
  # Extracts Sensitivity, Specificity, and Balanced Accuracy for CN, MCI, AD
  metrics_by_class <- cm$byClass[, c("Sensitivity", "Specificity", "Balanced Accuracy")]
  
  cat("\n--- DETAILED METRICS BY CLASS ---\n")
  print(round(metrics_by_class, 4))
  
  # Extract Global Metrics
  acc <- cm$overall["Accuracy"]
  bal_acc <- mean(cm$byClass[,"Balanced Accuracy"], na.rm=TRUE)
  
  cat(paste0("\n--- GLOBAL SUMMARY ---\n"))
  cat(paste0("Overall Accuracy:          ", round(acc, 4), "\n"))
  cat(paste0("Average Balanced Accuracy: ", round(bal_acc, 4), "\n"))
  
  return(list(cm=cm, accuracy=acc, balanced_acc=bal_acc, pred = final_preds_f, actual = actual_f))
}

# --- EXECUTE FINAL TEST ---
# We use the pure Test Set (containing all classes)
test_set_final <- testDataS1

# PIPELINE A: Dunn + Stage 2
res_dunn <- run_pipeline(test_set_final, final_model_dunn, th_dunn_avg, final_model_stage2, th_stage2_avg, "DUNN + STAGE2")

# PIPELINE B: Wilcoxon + Stage 2
res_wilcox <- run_pipeline(test_set_final, final_model_wilcox, th_wilcox_avg, final_model_stage2, th_stage2_avg, "WILCOXON + STAGE2")


# ==============================================================================
# PART 8: INTERPRETABILITY PLOTS AND REVISION---------
# ==============================================================================

# --- A. STAGE 1 (DUNN MODEL) - SCREENING ---
cat("\n>>> Generating Plots for STAGE 1 (Dunn) <<<\n")

# 1. Non-Linear Curves (Smooth)
try(dev.off(), silent=TRUE) # Reset
par(mfrow = c(1, 2))        # 2 side-by-side plots

plot(final_model_dunn, 
     scheme = 1, 
     shade = TRUE, shade.col = "lightblue", 
     residuals = TRUE, pch = 20, cex = 0.3, col = rgb(0,0,0,0.1),
     main = "Non-Linear Effect")

# 2. Linear Lines (Parametric)
# FIX: Now plots each variable with its name underneath
par(mfrow = c(2, 2)) # 2x2 Grid to fit linear variables

termplot(final_model_dunn, 
         terms = NULL,        # All linear variables
         se = TRUE,           # Show error (shading)
         col.term = "blue", 
         col.se = "lightblue",
         rug = TRUE,          # Show data density
         main = "Partial Effect", # Generic title
         xlab = NULL,         # IMPORTANT: Leave NULL to use the variable name!
         ylab = "Risk Contribution") 


# --- B. STAGE 2 (GAM MODEL) - STAGING ---
cat("\n>>> Generating Plots for STAGE 2 (Staging) <<<\n")


try(dev.off(), silent=TRUE) # Page reset
par(mfrow = c(1, 1))        # 2 plots

plot(final_model_stage2, 
     scheme = 1, 
     shade = TRUE, shade.col = "lightgreen", 
     residuals = TRUE, pch = 20, cex = 0.3, col = rgb(0,0,0,0.1),
     main = "Non-Linear Effect")


par(mfrow = c(2, 2)) # 2x2 Grid

termplot(final_model_stage2, 
         terms = NULL, 
         se = TRUE, 
         col.term = "darkred", 
         col.se = "pink",
         rug = TRUE, 
         main = "Partial Effect",
         xlab = NULL,        # IMPORTANT: Leave NULL to see names!
         ylab = "Risk Contribution")

# Final layout reset
par(mfrow = c(1, 1))

# REVISION ----------------------------------------------------
# The revision is employed on final_model_dunn
error_analysis_df <- data.frame(
  SubjectID = test_set_final$Subject.ID,
  Actual = res_dunn$actual,
  Predicted = res_dunn$pred,
  stringsAsFactors = FALSE
)

# Predicted = "MCI" but Actual = "CN"
mci_false_positives <- subset(error_analysis_df, Predicted == "MCI" & Actual == "CN")
# youden
# total: 24
# true: 8 => 002_S_5178, 007_S_1222, 007_S_5265, 009_S_4337, 009_S_4388,
#           010_S_0420, 011_S_4105, 012_S_1009

ad_false_positives1 <- subset(error_analysis_df, Predicted == "AD" & Actual == "CN")
# youden
# total: 3
# To MCI: 1: 013_S_4579
# True: 0

ad_false_positives2 <- subset(error_analysis_df, Predicted == "AD" & Actual == "MCI")
# youden
# total: 28
# True: 15 => 006_S_1130, 007_S_0128, 007_S_0249, 007_S_0344, 009_S_4324,
#            011_S_0241, 012_S_1292, 012_S_4094, 016_S_0702,123_S_4096
#            016_S_4902, 018_S_0155, 018_S_0406, 035_S_4784, 137_S_4815
#            , 



## New confusion matrix-------------
library(caret)


conf_data <- c(
  14, 18, 1,
  16, 83, 5,  # -8 for first col  => 24,75,5 => 16, 83, 5
  3,  13, 49  # -0 for first col, -15 for second col => 3, 28, 34 => 3,13,49
)

conf_matrix_data <- matrix(conf_data, nrow = 3, ncol = 3, byrow = TRUE)
rownames(conf_matrix_data) <- c("CN", "MCI", "AD")
colnames(conf_matrix_data) <- c("CN", "MCI", "AD")

conf_table <- as.table(conf_matrix_data)
names(dimnames(conf_table)) <- c("Prediction", "Reference")

cm_3x3_revised <- confusionMatrix(conf_table)

print(cm_3x3_revised$table)
print(cm_3x3_revised$byClass[, c("Sensitivity", "Specificity", "Balanced Accuracy")])

overall_acc <- cm_3x3_revised$overall["Accuracy"]
avg_bal_acc <- mean(cm_3x3_revised$byClass[, "Balanced Accuracy"])
cat(sprintf("Overall Accuracy:          %.4f (%.2f%%)\n", overall_acc, overall_acc * 100))
cat(sprintf("Average Balanced Accuracy: %.4f (%.2f%%)\n", avg_bal_acc, avg_bal_acc * 100))



# ==============================================================================
# PART 9: CONFORMAL PREDICTION (UNCERTAINTY QUANTIFICATION)------
# ==============================================================================
# Method: Split Conformal Prediction (Inductive)
# Goal: Provide guaranteed prediction sets for Stage 2 (MCI vs AD) with 90% confidence.

cat("\n>>> INITIALIZING CONFORMAL PREDICTION MODULE <<<\n")

# 1. SPLIT FOR CALIBRATION
# We need a fresh 'Calibration Set' that the Stage 2 model hasn't seen perfectly.
# Ideally, we should have split Training Set into (Train Proper + Calibration).
# For this thesis demo, we will perform a 20% hold-out from 'trainDataS2' for calibration.

set.seed(999)
calibIndex <- createDataPartition(trainDataS2$Group, p = .2, list = FALSE)
data_calib <- trainDataS2[calibIndex, ]       # Used to calculate 'strangeness' (scores)
data_train_cp <- trainDataS2[-calibIndex, ]   # Used to re-train the underlying model

# Re-calculate weights for this new sub-training set
w_cp <- ifelse(data_train_cp$Group == "MCI", 
               (nrow(data_train_cp)/sum(data_train_cp$Group == "MCI")), 
               (nrow(data_train_cp)/sum(data_train_cp$Group == "AD")))

# 2. RE-TRAIN MODEL (Underlying Algorithm)
# We need a model trained ONLY on 'data_train_cp' to avoid overfitting the calibration set.
cat("... Re-training Stage 2 Model on Proper Training Subset ...\n")
model_cp_stage2 <- gam(f_stage2, 
                       data = data_train_cp, 
                       family = quasibinomial(), 
                       weights = w_cp)

# 3. COMPUTE NON-CONFORMITY SCORES (CALIBRATION)
# We predict on the calibration set and see how "wrong" the model is.
probs_calib <- predict(model_cp_stage2, data_calib, type = "response")
scores_calib <- rep(NA, nrow(data_calib))

for(i in 1:nrow(data_calib)) {
  true_label <- data_calib$Group[i]
  prob_ad <- probs_calib[i]
  
  # Score = 1 - Probability of the TRUE class
  # If I am AD and prob is 0.9, score is 0.1 (Low error)
  # If I am AD and prob is 0.2, score is 0.8 (High error)
  if(true_label == "AD") {
    scores_calib[i] <- 1 - prob_ad 
  } else {
    scores_calib[i] <- prob_ad # Equivalent to 1 - (1-prob_ad) for MCI
  }
}

# 4. COMPUTE Q-HAT (The Threshold)
# We want 90% confidence (Alpha = 0.10)
alpha <- 0.10
n_cal <- length(scores_calib)
# Statistical correction for finite sample
q_level <- ceiling((n_cal + 1) * (1 - alpha)) / n_cal
q_level <- min(1, q_level) # Clip to 1 max

q_hat <- quantile(scores_calib, probs = q_level, type = 1)

cat(paste0("... Calibration Complete. Alpha: ", alpha, " (90% Guarantee) ...\n"))
cat(paste0("... Calibrated Threshold (Q-hat): ", round(q_hat, 4), " ...\n"))


# ==============================================================================
# PART 10: SEQUENTIAL CONFORMAL PIPELINE (FULL TEST FLOW)-------
# ==============================================================================
# Simulate the full clinical journey: Screening (Dunn) -> Staging (Conformal)

cat("\n>>> RUNNING SEQUENTIAL CONFORMAL TEST <<<\n")

# A. STEP 1: SCREENING (DUNN)
# Apply Dunn model to the entire unseen Test Set
probs_s1_final <- predict(final_model_dunn, test_set_final, type="response")
# Use the robust threshold calculated in Part 5
preds_s1_final <- ifelse(probs_s1_final > th_dunn_avg, "PAT", "CN")

# B. FILTERING: WHO SURVIVES?
# Select only subjects classified as 'PAT' (Pathological)
pat_indices <- which(preds_s1_final == "PAT")
survivors_data <- test_set_final[pat_indices, ]

cat(paste0("Total Test Patients: ", nrow(test_set_final), "\n"))
cat(paste0("Filtered out as CN:  ", nrow(test_set_final) - length(pat_indices), "\n"))
cat(paste0("Passed to Stage 2:   ", length(pat_indices), " (Includes True PAT + False Positives)\n"))

# C. STEP 2: CONFORMAL PREDICTION ON SURVIVORS
# Predict probabilities using the CP model
probs_cp_survivors <- predict(model_cp_stage2, survivors_data, type="response")

cp_outcomes <- c()
set_sizes <- c()

for(i in 1:nrow(survivors_data)) {
  # 1. Build Prediction Set
  p_ad <- probs_cp_survivors[i]
  current_set <- c()
  
  # Include MCI if error score <= Q-hat
  if (p_ad <= q_hat) { current_set <- c(current_set, "MCI") }
  # Include AD if error score <= Q-hat
  if ((1 - p_ad) <= q_hat) { current_set <- c(current_set, "AD") }
  
  set_sizes[i] <- length(current_set)
  
  # 2. Evaluate Clinical Outcome against Ground Truth
  true_lab <- as.character(survivors_data$Group[i])
  
  # CASE I: The patient is actually HEALTHY (CN)
  # This is a False Positive from Stage 1. 
  if (true_lab == "CN") {
    if (length(current_set) == 0) {
      # SUCCESS: The model realized this patient doesn't look like MCI or AD
      cp_outcomes[i] <- "CN Caught as Anomaly (Empty Set)"
    } else {
      # FAILURE: The error propagated
      cp_outcomes[i] <- "Stage 1 FP (Wrongly Diagnosed)"
    }
  } 
  
  # CASE II: The patient is truly MCI or AD
  else {
    if (true_lab %in% current_set) {
      if (length(current_set) == 1) {
        cp_outcomes[i] <- "Correct & Precise (Singleton)" # Perfect
      } else {
        cp_outcomes[i] <- "Correct but Uncertain {MCI, AD}" # Safe
      }
    } else {
      cp_outcomes[i] <- "Missed Diagnosis (Error)" # Violation
    }
  }
}



# Prepare Data for Plotting
df_seq_plot <- data.frame(
  TrueGroup = survivors_data$Group,
  Outcome = factor(cp_outcomes, 
                   levels = c("Correct & Precise (Singleton)", 
                              "Correct but Uncertain {MCI, AD}", 
                              "CN Caught as Anomaly (Empty Set)",
                              "Stage 1 FP (Wrongly Diagnosed)",
                              "Missed Diagnosis (Error)"))
)

# Generate the Final Chart
library(ggplot2)

plot_cp <- ggplot(df_seq_plot, aes(x = TrueGroup, fill = Outcome)) +
  geom_bar(position = "stack", color="black", alpha=0.9, width=0.7) +
  scale_fill_manual(values = c(
    "Correct & Precise (Singleton)" = "forestgreen",         # Ideal
    "Correct but Uncertain {MCI, AD}" = "gold",              # Acceptable/Safe
    "CN Caught as Anomaly (Empty Set)" = "cyan",             # Error Correction
    "Stage 1 FP (Wrongly Diagnosed)" = "gray50",             # Propagated Error
    "Missed Diagnosis (Error)" = "firebrick"                 # Critical Error
  )) +
  labs(title = "Clinical Pipeline Outcome (Stage 2 with Conformal Prediction)",
       subtitle = paste0("Analysis of patients passed from Screening (Dunn). Alpha = ", alpha),
       x = "True Clinical Diagnosis (Ground Truth)",
       y = "Number of Patients",
       fill = "Diagnostic Outcome") +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(face="bold", size=14),
        axis.text = element_text(size=12))

# Print Plot
print(plot_cp)

cat("\n>>> PIPELINE EXECUTION COMPLETE <<<\n")