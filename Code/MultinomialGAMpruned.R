library(ISLR2)
library(car)
library(mgcv)
library(rgl)
library(splines)
library(pbapply)
library(dplyr) #to clean the names of the vars
library(stringr) #for renaming the col
library(caret) #Partition test train

#Post processing
library(ggplot2)
library(tidyr)


#Load the dataset without the IDs
cleanMRITotal <- read.csv("D:/POLIMI/2Y/AppStat/Project/Repository/ADNIAlzheimer/Rcode/Processing/cleanMRITotal.csv")
#cleanMRITotal <- cleanMRITotal[,-c(2,3,4,11)] #age sex weight, Third ventricle which is duplicated ,52:67 if you want Cbm duplicated

cleanMRITotal <- cleanMRITotal %>%
  rename_with(~ str_remove_all(.x, "(?i)total.|\\.Volume_mm3"))#Let us try a simple GAM using 3 regressors with cubic splines
#As factor for the label otherwise it messes up the label
cleanMRITotal$Group <- as.factor(cleanMRITotal$Group)
#Merge MCIs
levels(cleanMRITotal$Group)[levels(cleanMRITotal$Group) %in% c("EMCI", "LMCI")] <- "MCI"
#CN vs all
cleanMRITotal$Group <- ordered(cleanMRITotal$Group, levels = c("CN", "MCI", "AD"))
#MCI vs All
#
#
#For proportional train vs test 70% 30%
set.seed(123)

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



# Calculate weights inversely proportional to class frequency
targetFreq <- table(trainData$Group)
weightsVec <- 1 / targetFreq[as.character(trainData$Group)]
weightsVec <- weightsVec * (nrow(trainData) / sum(weightsVec))


#For gam classification
#
trainData$Group <- as.numeric(trainData$Group) - 1




#####Semi-Parametric significant 
final_f1 <- as.formula(paste("Group ~ Inf.Lat.Vent + ctx.middletemporal +
                              Subject.Age "))


final_f2 <- as.formula(paste(" ~ ctx.paracentral + Accumbens.area + 
                              Third.Ventricle + Hippocampus +
                              Subject.Age +
                              s(Amygdala) + s(ctx.precuneus) + s(CSF) + 
                              s(ctx.lingual) + s(ctx.inferiortemporal) + 
                              Weight.Kg"))

finalModelGrouped <- gam(list(final_f1, final_f2), data = trainData,
                         family = multinom(K = 2), method="REML", ,
                         weights = weightsVec,
                         select = TRUE)
summary(finalModelGrouped)

# Generate probability predictions for the test set
prob_preds <- predict(finalModelGrouped, newdata = testData, type = "response")

labels <- c("CN", "MCI", "AD")

trainData$Group <- factor(trainData$Group, 
                          levels = labels)

colnames(prob_preds) <- levels(trainData$Group)



# Identify the class with the highest probability for each row
predicted_classes <- colnames(prob_preds)[apply(prob_preds, 1, which.max)]

# Convert to factor with the same levels as your original data
predicted_classes <- factor(predicted_classes, levels = levels(testData$Group))

# Create the confusion matrix
conf_matrix <- confusionMatrix(predicted_classes, testData$Group)

# Print the results
print(conf_matrix)

plot(finalModelGrouped, pages = 1, shades = TRUE, seWithMean = TRUE)
plot(finalModelGrouped, pages = 1, shades = TRUE, seWithMean = FALSE)
