#GAM models for ADNI
#18-Dec-2025
#
library(ISLR2)
#library(car)
library(mgcv)
library(rgl)
library(splines)
library(pbapply)
library(dplyr) #to clean the names of the vars
library(stringr) #for renaming the col
library(caret) #Partition test train
#Load the dataset without the IDs
cleanMRITotal <- read.csv("D:/POLIMI/2Y/AppStat/Project/Repository/ADNIAlzheimer/Rcode/Processing/cleanMRITotal.csv")
cleanMRITotal <- cleanMRITotal[,-c(2,3,4,76)] #age sex weight, Third ventricle which is duplicated ,52:67 if you want Cbm duplicated

cleanMRITotal <- cleanMRITotal %>%
  rename_with(~ str_remove_all(.x, "(?i)total.|\\.Volume_mm3"))#Let us try a simple GAM using 3 regressors with cubic splines
#As factor for the label otherwise it messes up the label
cleanMRITotal$Group <- as.factor(cleanMRITotal$Group)
#Merge MCIs
levels(cleanMRITotal$Group)[levels(cleanMRITotal$Group) %in% c("EMCI", "LMCI")] <- "MCI"
cleanMRITotal$Group <- factor(cleanMRITotal$Group, levels = c("CN", "MCI", "AD"))




# Partitioning ------------------------------------------------------------

#Partitioning
#For proportional train vs test 70% 30%
set.seed(2025)

# 1. Create the index for a 70% stratified split
# 'y' is your outcome variable (the merged Group column)
trainIndex <- createDataPartition(cleanMRITotal$Group, p = 0.7, list = FALSE)

# 2. Subset the data
trainData <- cleanMRITotal[trainIndex, ]
testData  <- cleanMRITotal[-trainIndex, ]

# 3. Verify the proportions are the same
prop.table(table(cleanMRITotal$Group)) # Original
prop.table(table(trainData$Group))     # Training
prop.table(table(testData$Group))      # Testing




# Simple GAM test ---------------------------------------------------------


#For gam classification
trainData$Group <- as.numeric(trainData$Group) - 1


#Gam solely with intercept
gamIntercept <- gam(list(Group ~ 1,~ 1),
                    data = trainData, 
                    family = multinom(K = 2))

#Gam with one volume + intercept
gamHippo <- gam(list(Group ~ s(Hippocampus, bs='cr'),
         ~ s(Hippocampus, bs='cr')),
    data = trainData, 
    family = multinom(K = 2))


# Gam with two smooth volumes + intercept
gamHippoCSF <- gam(list(Group ~ s(Hippocampus, bs='cr') + s(CSF, bs='cr'),
                     ~ s(Hippocampus, bs='cr') + s(CSF, bs='cr')),
                data = trainData, 
                family = multinom(K = 2))

#Anova to test for explanability of the variance
#How much variability of the data is explained by x vols

anova(gamIntercept, gamHippo, gamHippoCSF, test = "Chisq")


#AIC to test for significance of the features being added

AIC(gamIntercept, gamHippo, gamHippoCSF)

#GLM 
#
# 1. Identify predictors (all columns except 'Group')
target <- "Group"
predictors <- setdiff(names(trainData), target)

# 2. Create the smooth terms string: "s(Var1, bs='cr', k=5) + s(Var2, bs='cr', k=5) ..."
# We use a smaller k=5 because you have 300 observations and many features.
  #smooth_terms <- paste0("s(", predictors, ")", collapse = " + ")
  smooth_terms <- paste0("s(", predictors, ", bs='cr', k=2)", collapse = " + ")
  
# 3. Construct the list of formulas for multinom(K=2)
# We need two identical formulas (one for each non-reference group)
formula_1 <- as.formula(paste(target, "~", smooth_terms))
formula_2 <- as.formula(paste("~", smooth_terms))

formula_list <- list(formula_1, formula_2)

# 4. Fit the model
# We use select = TRUE to help the model drop useless features (Automatic Selection)
gamGLM <- gam(formula_list, 
              data = trainData, 
              family = multinom(K = 2),
              method = "REML",
              select = TRUE)

# 1. Identify all predictor names (everything except Group)
target <- "Group"
predictors <- setdiff(names(trainData), target)

# 2. Create a single string of all predictors: "Var1 + Var2 + Var3..."
# This is the standard 'linear' way without s()
rhs <- paste(predictors, collapse = " + ")

# 3. Create the list of formulas for 3 groups (K=2)
# Both formulas use the same 'rhs' string
formula_list <- list(
  as.formula(paste(target, "~", rhs)), # Predicts Group 2 vs 1
  as.formula(paste("~", rhs))          # Predicts Group 3 vs 1
)

# 4. Fit the parametric model
# Since there are no smooths, this is technically a Multinomial Logistic Regression
gam_glm_version <- gam(formula_list, 
                       data = trainData, 
                       family = multinom(K = 2))

summary(gam_glm_version)

numeric_cols <- sapply(trainData, is.numeric)

# 2. Identify the names of these columns, but remove "Group"
# (Even if Group is numeric, we don't want to scale the category labels)
cols_to_scale <- setdiff(names(trainData)[numeric_cols], "Group")

# 3. Apply the scale function
# We use as.vector() because scale() returns a matrix which can sometimes 
# cause issues in model formulas
trainData[cols_to_scale] <- lapply(trainData[cols_to_scale], function(x) as.vector(scale(x)))

## GLM GAM extract p values  -----------------------------------------------
# 1. Get the summary object
sumGam <- summary(gam_glm_version)

# 2. Extract the parametric coefficients table
# This contains Estimate, Std. Error, z value, and Pr(>|z|)
#paramCoeffs <- sumGam$p.table
paramCoeffs <- sumGam$p.table

# 3. Filter for p-values <= 0.01
significantTerms <- paramCoeffs[paramCoeffs[, "Pr(>|z|)"] < 0.05, , drop = FALSE]

# 4. View the results
print(significantTerms)


## GLS GAM significant only ------------------------------------------------
#From the significant of GLM GAM
# Create the right-hand side of the formula: "var1 + var2 + var3"
alpha <- 0.05
significantVars <- rownames(significantTerms[significantTerms[, "Pr(>|z|)"] < alpha, ])
significantVars <- significantVars[significantVars != "(Intercept).1"]
cleanNames <- gsub("\\.[0-9]+$", "", significantVars)
cleanNames <- unique(cleanNames)
rhs <- paste(cleanNames, collapse = " + ")

# Create the full formula (replace 'y' with your actual response variable)
formula_list <- list(
  as.formula(paste(target, "~", rhs)), # Predicts Group 2 vs 1
  as.formula(paste("~", rhs))          # Predicts Group 3 vs 1
)


# Refit the model
gamGLMsignificant <- gam(formula_list, 
                         data = trainData, 
                         family = multinom(K = 2))


# Calculate weights inversely proportional to class frequency
targetFreq <- table(trainData$Group)
weightsVec <- 1 / targetFreq[as.character(trainData$Group)]

# Refit with weights
gam_weighted <- gam(formula_list, 
                    data = trainData, 
                    family = multinom(K = 2),
                    weights = weightsVec)
  
# MCI ---------------------------------------------------------------------

#The GAM classification model is as
#y = b0 + b1x1 ... 
#x1,x2,... vols significant in the approach 1

## 1st Approach --------------------------------------------------------------


# This includes 7 vars
#Now, for approach 1, run the code then ... 
volumeVars <- sig_MCI$feature
trainData$Group <- as.numeric(as.factor(trainData$Group)) - 1

#trainData$Group <- as.factor(trainData$Group)

# 2. Initialize with the Null Model (Intercept only for both predictors)
nullModel <- gam(list(Group ~ 1, ~ 1), data = trainData, family = multinom(K = 2))
bestAic <- AIC(nullModel)
selectedVars <- c()

# 3. The Loop
repeat {
  bestCandidateAic <- Inf
  bestCandidate <- NULL
  
  for(var in volumeVars) {
    if(!(var %in% selectedVars)) {
      
      # Build the list of formulas dynamically
      # We add the new variable to BOTH linear predictors to match your example
      base_f <- paste(c("Group ~ 1", selectedVars), collapse = " + ")
      if(length(selectedVars) > 0) {
        # Add splines for existing vars
        f1 <- as.formula(paste("Group ~", paste0("s(", c(selectedVars, var), ", bs='cr')", collapse = " + ")))
        f2 <- as.formula(paste("~", paste0("s(", c(selectedVars, var), ", bs='cr')", collapse = " + ")))
      } else {
        # First iteration
        f1 <- as.formula(paste("Group ~ s(", var, ", bs='cr')"))
        f2 <- as.formula(paste("~ s(", var, ", bs='cr')"))
      }
      
      # Fit using the list format
      modelTry <- try(gam(list(f1, f2), data = trainData, family = multinom(K = 2)), silent = TRUE)
      
      if(!inherits(modelTry, "try-error")) {
        currentAic <- AIC(modelTry)
        if(currentAic < bestCandidateAic) {
          bestCandidateAic <- currentAic
          bestCandidate <- var
        }
      }
    }
  }
  
  # Decision Rule
  if(bestCandidateAic < bestAic) {
    bestAic <- bestCandidateAic
    selectedVars <- c(selectedVars, bestCandidate)
    print(paste("Added:", bestCandidate, "| New AIC:", round(bestAic, 2)))
  } else {
    print("No further improvement. Stopping.")
    break
  }
}


cleanVars <- paste(selectedVars, collapse = " + ")

# Combine with the dependent variable
formulaText <- paste("Group ~", cleanVars)
print(formulaText)
# "Group ~ Amygdala + Inf.Lat.Vent + ctx.middletemporal + ctx.entorhinal 
#         + Lateral.Ventricle""


## 2nd Approach ------------------------------------------------------------

#nullModel <- gam(list(Group ~ 1, ~ 1), data = trainData, family = multinom(K = 2))
#
sig_MCI_wilcox_sorted <- sig_MCI_wilcox %>%
  arrange(match(feature, sig_MCI$feature))

volumeVars <- sig_MCI_wilcox_sorted$feature
bestAic <- AIC(nullModel)
selectedVars <- c()

# 3. The Loop
repeat {
  bestCandidateAic <- Inf
  bestCandidate <- NULL
  
  for(var in volumeVars) {
    if(!(var %in% selectedVars)) {
      
      # Build the list of formulas dynamically
      # We add the new variable to BOTH linear predictors to match your example
      base_f <- paste(c("Group ~ 1", selectedVars), collapse = " + ")
      if(length(selectedVars) > 0) {
        # Add splines for existing vars
        f1 <- as.formula(paste("Group ~", paste0("s(", c(selectedVars, var), ", bs='cr')", collapse = " + ")))
        f2 <- as.formula(paste("~", paste0("s(", c(selectedVars, var), ", bs='cr')", collapse = " + ")))
      } else {
        # First iteration
        f1 <- as.formula(paste("Group ~ s(", var, ", bs='cr')"))
        f2 <- as.formula(paste("~ s(", var, ", bs='cr')"))
      }
      
      # Fit using the list format
      modelTry <- try(gam(list(f1, f2), data = trainData, family = multinom(K = 2)), silent = TRUE)
      
      if(!inherits(modelTry, "try-error")) {
        currentAic <- AIC(modelTry)
        if(currentAic < bestCandidateAic) {
          bestCandidateAic <- currentAic
          bestCandidate <- var
        }
      }
    }
  }
  
  # Decision Rule
  if(bestCandidateAic < bestAic) {
    bestAic <- bestCandidateAic
    selectedVars <- c(selectedVars, bestCandidate)
    print(paste("Added:", bestCandidate, "| New AIC:", round(bestAic, 2)))
  } else {
    print("No further improvement. Stopping.")
    break
  }
}


cleanVars <- paste(selectedVars, collapse = " + ")

# Combine with the dependete("Group ~", cleanVars)
print(formulaText)

# NOT SORTED
# "Group ~ Amygdala + Inf.Lat.Vent + ctx.inferiorparietal + Accumbens.area
#           + CSF + Cbm_VIIIb + Post.Hypothalamus + ctx.inferiortemporal 
#           + ctx.parsopercularis 

#SORTED 
#"Group ~ Amygdala + Inf.Lat.Vent + ctx.inferiorparietal + Accumbens.area
#           + CSF + Cbm_VIIIb + Post.Hypothalamus + ctx.inferiortemporal 
#           + ctx.parsopercularis"

final_f1 <- as.formula(paste("Group ~", paste0("s(", selectedVars, ", bs='cr')", collapse = " + ")))
final_f2 <- as.formula(paste("~", paste0("s(", selectedVars, ", bs='cr')", collapse = " + ")))
finalModel <- gam(list(final_f1, final_f2), data = trainData, family = multinom(K = 2))
summary(finalModel)



# CN vs AD ----------------------------------------------------------------

volumeVars <- sig_AD_wilcox$feature
bestAic <- AIC(nullModel)
selectedVars <- c()

# 3. The Loop
repeat {
  bestCandidateAic <- Inf
  bestCandidate <- NULL
  
  for(var in volumeVars) {
    if(!(var %in% selectedVars)) {
      
      # Build the list of formulas dynamically
      # We add the new variable to BOTH linear predictors to match your example
      base_f <- paste(c("Group ~ 1", selectedVars), collapse = " + ")
      if(length(selectedVars) > 0) {
        # Add splines for existing vars
        f1 <- as.formula(paste("Group ~", paste0("s(", c(selectedVars, var), ", bs='ts')", collapse = " + ")))
        f2 <- as.formula(paste("~", paste0("s(", c(selectedVars, var), ", bs='ts')", collapse = " + ")))
      } else {
        # First iteration
        f1 <- as.formula(paste("Group ~ s(", var, ", bs='ts')"))
        f2 <- as.formula(paste("~ s(", var, ", bs='ts')"))
      }
      
      # Fit using the list format
      modelTry <- try(gam(list(f1, f2), data = trainData, family = multinom(K = 2)), silent = TRUE)
      
      if(!inherits(modelTry, "try-error")) {
        currentAic <- AIC(modelTry)
        if(currentAic < bestCandidateAic) {
          bestCandidateAic <- currentAic
          bestCandidate <- var
        }
      }
    }
  }
  
  # Decision Rule
  if(bestCandidateAic < bestAic) {
    bestAic <- bestCandidateAic
    selectedVars <- c(selectedVars, bestCandidate)
    print(paste("Added:", bestCandidate, "| New AIC:", round(bestAic, 2)))
  } else {
    print("No further improvement. Stopping.")
    break
  }
}


cleanVars <- paste(selectedVars, collapse = " + ")

# Combine with the dependent variable
formulaText <- paste("Group ~", cleanVars)
print(formulaText)

# NOT SORTED
# "Group ~ Amygdala + Inf.Lat.Vent + ctx.inferiorparietal + Accumbens.area
#           + CSF + Cbm_VIIIb + Post.Hypothalamus + ctx.inferiortemporal 
#           + ctx.parsopercularis 

#SORTED 
#"Group ~ Amygdala + Inf.Lat.Vent + ctx.inferiorparietal + Accumbens.area
#           + CSF + Cbm_VIIIb + Post.Hypothalamus + ctx.inferiortemporal 
#           + ctx.parsopercularis"

           
selectedVars <- c("Amygdala")
selectedVars <- c("Amygdala", "Inf.Lat.Vent", "ctx.inferiorparietal", "Accumbens.area"
                  , "CSF", "Cbm_VIIIb", "Post.Hypothalamus", "ctx.inferiortemporal" 
                  , "ctx.parsopercularis" )


#These vars are taken from the AD vs CN 61 vars.
selectedVars <- c(
  "Amygdala", "Inf.Lat.Vent", "ctx.precuneus", "Putamen", 
  "ctx.precentral", "ctx.paracentral", "ctx.lateraloccipital", 
  "Post.Hypothalamus", "CSF", "ctx.lingual", "Accumbens.area", 
  "Third.Ventricle", "ctx.postcentral", "ctx.inferiortemporal", 
  "VentralDC", "Hippocampus", "Cbm_Vermis_VII", "Cbm_Vermis_IX", 
  "ctx.middletemporal", "Ant.Commisure", "ctx.parsorbitalis", 
  "N.opticus", "Optic.tract", "Cerebellum.Cortex"
)

final_f1 <- as.formula(paste("Group ~", paste0("s(", selectedVars, ", bs='tp')", collapse = " + ")))
final_f2 <- as.formula(paste("~", paste0("s(", selectedVars, ", bs='tp')", collapse = " + ")))
finalModel <- gam(list(final_f1, final_f2), data = trainData, family = multinom(K = 2), method="REML", select = TRUE)
summary(finalModel)

#finalModel <- gam(list(final_f1, final_f2), data = trainData, family = multinom(K = 2))

gam.check(finalModel)

# Generate probability predictions for the test set
prob_preds <- predict(finalModel, newdata = testData, type = "response")

labels <- c("CN", "MCI", "AD")
trainData$Group <- factor(trainData$Group, 
                          levels = 1:3, 
                          labels = labels)

colnames(prob_preds) <- levels(trainData$Group)

# Identify the class with the highest probability for each row
predicted_classes <- colnames(prob_preds)[apply(prob_preds, 1, which.max)]

# Convert to factor with the same levels as your original data
predicted_classes <- factor(predicted_classes, levels = levels(testData$Group))



# Create the confusion matrix
conf_matrix <- confusionMatrix(predicted_classes, testData$Group)

# Print the results
print(conf_matrix)




#THis is the base model 
cols_names <- setdiff(names(cleanMRITotal), "Group")
volumeVars <- cols_names

##
currentFormula <- "Group ~ 1"
selectedVars <- c() 
bestAic <- Inf 

# Iterate to find the best variables
repeat {
  bestCandidateAic <- Inf
  bestCandidate <- NULL
  
  # Try every variable not yet in the model
  for(var in volumeVars) {
    print(var)
    if(!(var %in% selectedVars)) {
      
      # Construct formula: Logit(P) = Current_Model + s(New_Var)
      # We use 'bs="cr"' (Cubic Regression Spline) as defined in Source 1 [6]
      formulaTry <- paste(currentFormula, "+ s(", var, ", bs=ts )")
      
      # Fit the GAM
      # family = binomial() is used for binary classification (AD vs CN) 
      # similar to the wage>250 example in Source 1 [9]
      modelTry <- try(gam(list(formulaTry_f1, formulaTry_f2), 
                          data = trainData, 
                          family = multinom(K = 2),
                          method = "REML"),    # Enable extra penalty for null space
                      silent = TRUE)
      
      # Check if model converged and get AIC (proxy for CV error [2])
      if(!inherits(modelTry, "try-error")) {
        currentAic <- AIC(modelTry)
        
        # Look for the lowest error
        if(currentAic < bestCandidateAic) {
          bestCandidateAic <- currentAic
          bestCandidate <- var
        }
      }
    }
  }
  
  # Decision Rule: Did we improve the model?
  if(bestCandidateAic < bestAic) {
    # Update the best score and add the variable to our "Building Blocks" [3]
    bestAic <- bestCandidateAic
    selectedVars <- c(selectedVars, bestCandidate)
    currentFormula <- paste(currentFormula, "+ s(", bestCandidate, ", bs=' r')")
    print(paste("Added:", bestCandidate, "| New AIC:", round(bestAic, 2)))
  } else {
    # If no variable improves the model, stop the process
    print("No further improvement. Stopping selection.")
    break
  }
}

