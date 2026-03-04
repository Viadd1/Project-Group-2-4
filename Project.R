# ==============================================================================
# Script Name: Host Classification using MCFS and Rough Set Theory (ROSETTA)
# Description: This pipeline performs data preprocessing, feature selection using 
#              Monte Carlo Feature Selection (MCFS), and rule-based classification 
#              using R.ROSETTA to distinguish between 'Host' classes (Avian vs. Human).
#              Results are visualized as a network using Cytoscape.
# ==============================================================================

# ------------------------------------------------------------------------------
# SECTION 1: Setup and Data Loading
# ------------------------------------------------------------------------------

# Set the working directory to the location of the dataset.
# (Update this path based on the machine running the script)
# setwd("/Users/dipc/Developer/uu/1MB416_KnowledgeBasedSystem/Project")
setwd("C:/Users/david/OneDrive - Uppsala universitet/Dokument/Uppsala Universitet/Biologi och bioteknik/Kunskapsbaserade system/Projekt")

# Load required libraries
library(rJava)      # Required for R.ROSETTA and rmcfs (Java-based packages)
library(rmcfs)      # For Monte Carlo Feature Selection
library(R.ROSETTA)  # For Rough Set-based rule generation and classification
library(devtools)   # For sourcing scripts or loading development packages
library(csVisuNet)  # For visualizing rule networks

# Read the dataset. 
# colClasses = "character" prevents automatic conversion of string features (like amino acids) to factors prematurely.
data <- read.csv("Project2.csv", sep="\t", colClasses = "character")

# ------------------------------------------------------------------------------
# SECTION 2: Data Preprocessing and Cleaning
# ------------------------------------------------------------------------------

# Remove the 'id' column as it is an identifier, not a predictive feature.
data$id <- NULL

# Filter out features with excessive missing values.
# In this dataset, missing values are denoted by "?". We retain only the columns 
# where less than 50% of the entries are missing to prevent noise.
cols_to_keep <- sapply(data, function(col) {
  sum(col == "?") / length(col) < 0.50 
})

# Apply the filter to create a clean subset
data_clean <- data[, cols_to_keep]

# Inspect the target variable ('Host') to understand the class distribution 
# (e.g., class imbalance between Avian and Human)
data_clean$Host
table(data_clean$Host)

# # ------------------------------------------------------------------------------
# # SECTION 3: MCFS Hyperparameter Optimization
# # ------------------------------------------------------------------------------
# 
# # Function to test combinations of 'projections' and 'projectionSize'
# test_mcfs_params <- function(clean_data) {
#   # Define the parameter grids to search over
#   proj_sizes <- c(0.05, 0.10, 0.20, 0.30) # Percentage of features evaluated per subset (5% to 30%)
#   projections_list <- c(500, 1000, 1500)  # Number of random subsets/trees to build
#   
#   # Initialize an empty dataframe to log performance metrics
#   results <- data.frame()
#   
#   for (p in projections_list) {
#     for (ps in proj_sizes) {
#       cat("Running MCFS -> Projections:", p, "| Projection Size:", ps, "\n")
#       
#       # Run MCFS with current grid parameters
#       temp_model <- mcfs(Host ~ ., 
#                          data = clean_data, 
#                          projections = p, 
#                          projectionSize = ps, 
#                          splits = 5,             # Defaulting to 5 splits for this test
#                          threadsNumber = 8,      # Parallel processing (adjust based on CPU)
#                          cutoffPermutations = 3) # Permutation tests to find significance cutoff
#       
#       # Extract key evaluation metrics
#       num_features <- temp_model$cutoff_value
#       top_feat <- as.character(temp_model$RI$attribute[1]) # The feature with the highest Relative Importance
#       
#       # Append current run's metrics to the results table
#       results <- rbind(results, data.frame(
#         Projections = p, 
#         ProjSize = ps, 
#         Num_Significant_Features = num_features, 
#         Top_Feature = top_feat
#       ))
#     }
#   }
#   return(results)
# }
# 
# # Run the parameter tuning and view results
# mcfs_experiment_results <- test_mcfs_params(data_clean)
# print(mcfs_experiment_results)
# 
# 
# # Function to test the effect of the 'splits' parameter on feature selection
# test_mcfs_splits <- function(clean_data) {
#   # Fix projection size and count based on optimal findings from previous test
#   proj_size <- 0.10 
#   projections <- 1500  
#   
#   # Test Shallow (2), Medium (5), and Deep (10) tree splits
#   splits_list <- c(2, 5, 10)  
#   
#   results <- data.frame()
#   
#   for (s in splits_list) {
#     cat("Running MCFS -> Splits:", s, "\n")
#     
#     temp_model <- mcfs(Host ~ ., data = clean_data, 
#                        projections = projections, 
#                        projectionSize = proj_size, 
#                        splits = s,           # Dynamic parameter
#                        threadsNumber = 8, 
#                        cutoffPermutations = 3) 
#     
#     num_features <- temp_model$cutoff_value
#     top_feat <- as.character(temp_model$RI$attribute[1])
#     
#     results <- rbind(results, data.frame(
#       Splits = s, 
#       Num_Significant_Features = num_features, 
#       Top_Feature = top_feat
#     ))
#   }
#   return(results)
# }
# 
# # Run the split experiment and evaluate
# split_experiment_results <- test_mcfs_splits(data_clean)
# print(split_experiment_results)

# ------------------------------------------------------------------------------
# SECTION 4: Final MCFS Execution and Feature Selection
# ------------------------------------------------------------------------------

# Run the final MCFS model utilizing the optimal parameters discovered above
mcfs_result <- mcfs(Host ~ ., 
                    data = data_clean, 
                    projections = 1500, 
                    projectionSize = 0.1, 
                    splits = 10,
                    threadsNumber = 8,
                    cutoffPermutations = 3)

# Display the top-ranked features based on Relative Importance (RI)
head(mcfs_result$RI)

# Plot convergence graph to verify if 1500 projections were sufficient to stabilize the model
plot(mcfs_result, type="distances")

# Plot Relative Importance graph to visualize the significance cutoff point
plot(mcfs_result, type = "ri", size = 47)

# Extract features that pass the permutation-based significance threshold (RI > cutoff_value)
mcfs_result2 <- mcfs_result$RI[1:mcfs_result$cutoff_value,]
top_features <- mcfs_result$RI[1:mcfs_result$cutoff_value, "attribute"]
print("Significant Features Selected:")
print(top_features)

# Subset the dataset to retain only the significant features and the target variable
selected_data <- data_clean[, c(as.character(top_features), "Host")]

# # ------------------------------------------------------------------------------
# # SECTION 5: ROSETTA Modeling (Rough Set Theory)
# # ------------------------------------------------------------------------------
# 
# # Helper function to evaluate different R.ROSETTA configurations
# run_rosetta_model <- function(data, reducer, classifier, cv_folds) {
#   model <- rosetta(
#     data,
#     discrete = TRUE,         # Features are treated as categorical/discrete (e.g., amino acids)
#     reducer = reducer,       # Algorithm for rule reduction (finding reducts)
#     classifier = classifier, # Algorithm for voting/classification
#     cv = cv_folds            # Cross-validation folds
#   )
#   
#   # Return combined metadata and quality metrics (Accuracy, AUC, etc.)
#   return(data.frame(
#     Reducer = reducer,
#     Classifier = classifier,
#     CV = cv_folds,
#     model$quality
#   ))
# }
# 
# # Define grid of ROSETTA parameters to test
# classifiers <- c("StandardVoter", "ObjectTrackingVoter", "NaiveBayesClassifier")
# cv_values <- c(5, 10)
# 
# # Evaluate using the "Johnson" reducer heuristic
# rosetta_summary_table_johnson <- do.call(rbind, lapply(classifiers, function(clf) {
#   do.call(rbind, lapply(cv_values, function(cv) {
#     run_rosetta_model(selected_data, "Johnson", clf, cv)
#   }))
# }))
# print("Johnson Reducer Results:")
# print(rosetta_summary_table_johnson)
# 
# # Evaluate using the "Genetic" reducer heuristic
# rosetta_summary_table_genetic <- do.call(rbind, lapply(classifiers, function(clf) {
#   do.call(rbind, lapply(cv_values, function(cv) {
#     run_rosetta_model(selected_data, "Genetic", clf, cv)
#   }))
# }))
# print("Genetic Reducer Results:")
# print(rosetta_summary_table_genetic)

# Train the final rule-based model based on optimal parameters 
# (Johnson reducer, Standard Voter, 5-fold CV chosen here)
rosetta_model <- rosetta(selected_data, 
                         discrete = TRUE,
                         reducer = "Johnson",
                         classifier = "StandardVoter",
                         cv = 5)

# Output overall model rules and cross-validation performance metrics
rules <- rosetta_model$main
qual <- rosetta_model$quality
print("Final Model Quality:")
print(qual)

# Show how the discrete values were interpreted by the model
slevels <- rosetta_model$main$levels

# ------------------------------------------------------------------------------
# SECTION 6: Rule Analysis and Clustering
# ------------------------------------------------------------------------------

# Filter for statistically significant rules (p-value < 0.05) to ensure 
# biological or scientific relevance
sig_rules <- rosetta_model$main[rosetta_model$main$pValue < 0.05, ]

# Load custom utility functions for rule handling
# (Note: These files must exist in the working directory)
source("recalculateRules.R")
source("clusterRules.R")

# Recalculate support and accuracy metrics for the extracted rules based on the subset data
recalc_rules <- recalculateRules(selected_data, rosetta_model$main, discrete=TRUE)

# Cluster the recalculated rules to group similar rules together.
# WARNING FIXED: Passing the full training df (selected_data) instead of sig_rules 
# avoids structural mismatch errors during clustering.
clustered_rules <- clusterRules(selected_data, recalc_rules)

# Sort rules by p-value to find the strongest, most significant predictors first
sorted_rules <- recalc_rules[order(recalc_rules$pValue), ]

# Print the top 5 distinguishing rules for the "Avian" host class
print("TOP 5 AVIAN RULES:")
head(sorted_rules[sorted_rules$decision == "Avian", 
                  c("features", "levels", "supportRHS", "accuracyRHS", "pValue")], 5)

# Print the top 5 distinguishing rules for the "Human" host class
print("TOP 5 HUMAN RULES:")
head(sorted_rules[sorted_rules$decision == "Human", 
                  c("features", "levels", "supportRHS", "accuracyRHS", "pValue")], 5)

# ------------------------------------------------------------------------------
# SECTION 7: Network Visualization (Cytoscape)
# ------------------------------------------------------------------------------

# Initiate connection to Cytoscape. 
# IMPORTANT: The Cytoscape desktop application MUST be running locally before executing this.
library(RCy3)
cytoscapePing()

# Generate the initial VisuNet network based on the ROSETTA rules

vis <- visunetcyto(rosetta_model$main)

# Create specific VisuArc diagrams to see how the highest-impact features 
# (e.g., specific positions P590, P337) map directly to the Human and Avian classes.
visuArc(vis, 'Human', 'P590')
visuArc(vis, 'Avian', 'P337')