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

# ------------------------------------------------------------------------------
# SECTION 3: MCFS Hyperparameter Optimization
# ------------------------------------------------------------------------------

# Function to test combinations of 'projections' and 'projectionSize'
test_mcfs_params <- function(clean_data) {
  # Define the parameter grids to search over
  proj_sizes <- c(0.05, 0.10, 0.20, 0.30) # Percentage of features evaluated per subset (5% to 30%)
  projections_list <- c(500, 1000, 1500)  # Number of random subsets/trees to build

  # Initialize an empty dataframe to log performance metrics
  results <- data.frame()

  for (p in projections_list) {
    for (ps in proj_sizes) {
      cat("Running MCFS -> Projections:", p, "| Projection Size:", ps, "\n")

      # Run MCFS with current grid parameters
      temp_model <- mcfs(Host ~ .,
                         data = clean_data,
                         projections = p,
                         projectionSize = ps,
                         splits = 5,             # Defaulting to 5 splits for this test
                         threadsNumber = 8,      # Parallel processing (adjust based on CPU)
                         cutoffPermutations = 3) # Permutation tests to find significance cutoff

      # Extract key evaluation metrics
      num_features <- temp_model$cutoff_value
      top_feat <- as.character(temp_model$RI$attribute[1]) # The feature with the highest Relative Importance

      # Append current run's metrics to the results table
      results <- rbind(results, data.frame(
        Projections = p,
        ProjSize = ps,
        Num_Significant_Features = num_features,
        Top_Feature = top_feat
      ))
    }
  }
  return(results)
}

# Run the parameter tuning and view results
mcfs_experiment_results <- test_mcfs_params(data_clean)
print(mcfs_experiment_results)


# Function to test the effect of the 'splits' parameter on feature selection
test_mcfs_splits <- function(clean_data) {
  # Fix projection size and count based on optimal findings from previous test
  proj_size <- 0.10
  projections <- 1500

  # Test Shallow (2), Medium (5), and Deep (10) tree splits
  splits_list <- c(2, 5, 10)

  results <- data.frame()

  for (s in splits_list) {
    cat("Running MCFS -> Splits:", s, "\n")

    temp_model <- mcfs(Host ~ ., data = clean_data,
                       projections = projections,
                       projectionSize = proj_size,
                       splits = s,           # Dynamic parameter
                       threadsNumber = 8,
                       cutoffPermutations = 3)

    num_features <- temp_model$cutoff_value
    top_feat <- as.character(temp_model$RI$attribute[1])

    results <- rbind(results, data.frame(
      Splits = s,
      Num_Significant_Features = num_features,
      Top_Feature = top_feat
    ))
  }
  return(results)
}

# Run the split experiment and evaluate
split_experiment_results <- test_mcfs_splits(data_clean)
print(split_experiment_results)

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

# ------------------------------------------------------------------------------
# SECTION 5: ROSETTA Modeling (Rough Set Theory)
# ------------------------------------------------------------------------------

# Helper function to evaluate different R.ROSETTA configurations
run_rosetta_model <- function(data, reducer, classifier, cv_folds) {
  model <- rosetta(
    data,
    discrete = TRUE,         # Features are treated as categorical/discrete (e.g., amino acids)
    reducer = reducer,       # Algorithm for rule reduction (finding reducts)
    classifier = classifier, # Algorithm for voting/classification
    cv = cv_folds            # Cross-validation folds
  )

  # Return combined metadata and quality metrics (Accuracy, AUC, etc.)
  return(data.frame(
    Reducer = reducer,
    Classifier = classifier,
    CV = cv_folds,
    model$quality
  ))
}

# Define grid of ROSETTA parameters to test
classifiers <- c("StandardVoter", "ObjectTrackingVoter", "NaiveBayesClassifier")
cv_values <- c(5, 10)

# Evaluate using the "Johnson" reducer heuristic
rosetta_summary_table_johnson <- do.call(rbind, lapply(classifiers, function(clf) {
  do.call(rbind, lapply(cv_values, function(cv) {
    run_rosetta_model(selected_data, "Johnson", clf, cv)
  }))
}))
print("Johnson Reducer Results:")
print(rosetta_summary_table_johnson)

# Evaluate using the "Genetic" reducer heuristic
rosetta_summary_table_genetic <- do.call(rbind, lapply(classifiers, function(clf) {
  do.call(rbind, lapply(cv_values, function(cv) {
    run_rosetta_model(selected_data, "Genetic", clf, cv)
  }))
}))
print("Genetic Reducer Results:")
print(rosetta_summary_table_genetic)

# Train the final rule-based model based on optimal parameters 
# (Johnson reducer, Standard Voter, 5-fold CV chosen here)
rosetta_model <- rosetta(selected_data, 
                         discrete = TRUE,
                         reducer = "Johnson",
                         classifier = "StandardVoter",
                         cv = 10)

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

# Taking all the rules from each class
avian_rules <- sorted_rules[sorted_rules$decision == "Avian", ]
human_rules <- sorted_rules[sorted_rules$decision == "Human", ]

# Computing how many rules each class has
cat("Avian rules:", nrow(avian_rules), "\n")
cat("Human rules:", nrow(human_rules), "\n")

# Generate VisuNet networks based on the rules
# Create specific VisuArc diagrams to see how the highest-impact features 
# (e.g., position P590 in both networks) map directly to the Human and Avian classes.
vis_avian <- visunetcyto(avian_rules)
visuArc(vis_avian, 'Avian', 'P590')

# Save progress in avian network before running!
vis_human <- visunetcyto(human_rules)
visuArc(vis_human, 'Human', 'P590')





# Load the plotting library
library(ggplot2)

library(ggplot2)

slide_bg <- "#F3F7FF"   # rgb(243, 247, 255)
text_col <- "#100090"   # rgb(16, 0, 144)

# ---------------------------------------------------------
# PLOT 1: The MCFS Projection Size Experiment (Line Graph)
# ---------------------------------------------------------
plot_mcfs <- ggplot(mcfs_experiment_results, 
                    aes(x = as.factor(ProjSize), 
                        y = Num_Significant_Features, 
                        color = as.factor(Projections), 
                        group = Projections)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4) +
  theme_minimal(base_size = 14) +
  labs(title = "Effect of Projection Size and Projection Size on Feature Selection",
       x = "Projection Size (% of amino acids per tree)",
       y = "Number of Significant Features",
       color = "Total Projections") +
  scale_color_brewer(palette = "Set1") +
  theme(
    legend.position = "bottom",
    
    plot.background = element_rect(fill = slide_bg, color = NA),
    panel.background = element_rect(fill = slide_bg, color = NA),
    
    # Legend background blending (prevents white boxes around the keys)
    legend.background = element_rect(fill = slide_bg, color = NA),
    legend.key = element_rect(fill = slide_bg, color = NA),
    
    # TEXT COLOR UPDATES
    plot.title = element_text(color = text_col, face = "bold"),
    plot.subtitle = element_text(color = text_col),
    axis.title.x = element_text(color = text_col),
    axis.title.y = element_text(color = text_col),
    axis.text.x = element_text(color = text_col),
    axis.text.y = element_text(color = text_col),
    legend.text = element_text(color = text_col),
    legend.title = element_text(color = text_col)
  )

# Show and save Plot 1
print(plot_mcfs)
ggsave("MCFS_ProjectionSize_Custom.png", plot = plot_mcfs, width = 10, height = 6, dpi = 300)


# ---------------------------------------------------------
# PLOT 2: The MCFS Splits Experiment (Bar Chart)
# ---------------------------------------------------------
plot_splits <- ggplot(split_experiment_results, 
                      aes(x = as.factor(Splits), 
                          y = Num_Significant_Features)) +
  
  geom_col(fill = "#5e81ac", width = 0.6) + 
  
  geom_text(aes(label = Num_Significant_Features), vjust = -0.5, size = 5, color = text_col) + 
  
  theme_minimal(base_size = 14) +
  labs(title = "Effect of Tree Depth (Splits)",
       x = "Maximum Tree Depth (Splits)",
       y = "Number of Significant Features") +
  theme(
    # Background blending
    plot.background = element_rect(fill = slide_bg, color = NA),
    panel.background = element_rect(fill = slide_bg, color = NA),
    
    # Clean up vertical grid lines 
    panel.grid.major.x = element_blank(),
    
    # TEXT COLOR UPDATES
    plot.title = element_text(color = text_col, face = "bold"),
    plot.subtitle = element_text(color = text_col),
    axis.title.x = element_text(color = text_col),
    axis.title.y = element_text(color = text_col),
    axis.text.x = element_text(color = text_col),
    axis.text.y = element_text(color = text_col)
  )

# Show and save Plot 2
print(plot_splits)
ggsave("MCFS_Splits_Custom.png", plot = plot_splits, width = 10, height = 6, dpi = 300)


library(ggplot2)

# Extract the top 10 features
top_features_data <- head(mcfs_result$RI, 10)

# Define  slide colors 
slide_bg <- "#F3F7FF"   # rgb(243, 247, 255)
text_col <- "#100090"   # rgb(16, 0, 144)

# Create the Horizontal Bar Chart
plot_importance <- ggplot(top_features_data, 
                          aes(x = reorder(attribute, RI), y = RI)) +
  
  # The bars
  geom_col(fill = "#5e81ac", width = 0.7) +   
  
  # The numbers at the end of the bars
  geom_text(aes(label = round(RI, 3)), 
            hjust = -0.2, size = 5, color = text_col) + 
  
  coord_flip() + 
  theme_minimal(base_size = 14) +
  labs(title = "Top 10 Most Significant Amino Acid Positions",
       subtitle = "Ranked by Monte Carlo Feature Selection (MCFS)",
       x = "PB1 Protein Position",
       y = "Relative Importance (RI) Score") +
  theme(
    # Background blending
    plot.background = element_rect(fill = slide_bg, color = NA),
    panel.background = element_rect(fill = slide_bg, color = NA),
    
    # Grid lines cleanup
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    
    # TEXT COLOR UPDATES
    plot.title = element_text(color = text_col, face = "bold"),
    plot.subtitle = element_text(color = text_col),
    axis.title.x = element_text(color = text_col),
    axis.title.y = element_text(color = text_col),
    axis.text.x = element_text(color = text_col),
    axis.text.y = element_text(color = text_col)
  )

# Show the plot!
print(plot_importance)

# Save the high-res image
ggsave("MCFS_Top10_Features_Custom.png", plot = plot_importance, width = 10, height = 6, dpi = 300)


library(ggplot2)
master_rosetta

# 1. Combine two summary tables into one master table
master_rosetta <- rbind(rosetta_summary_table_johnson, rosetta_summary_table_genetic)

# 2. Define slide colors 
slide_bg <- "#F3F7FF"   # rgb(243, 247, 255)
text_col <- "#100090"   # rgb(16, 0, 144)

# 3. Create the Presentation Plot
plot_rosetta <- ggplot(master_rosetta, 
                       aes(x = Classifier, y = accuracyMean, fill = Reducer)) +
  
  # Create the grouped bars
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  
  # Add the error bars to show the Standard Deviation
  geom_errorbar(aes(ymin = accuracyMean - accuracyStd, ymax = accuracyMean + accuracyStd),
                position = position_dodge(width = 0.8), width = 0.25, color = text_col, linewidth = 0.8) +
  
  # Split the chart into two panels based on the CV folds
  facet_wrap(~ paste(CV, "- Fold Cross Validation")) +
  
  # Zoom in the Y-axis so the differences are highly visible
  coord_cartesian(ylim = c(0.80, 1.0)) +
  
  # Custom colors for the Reducers 
  scale_fill_manual(values = c("Johnson" = "#5e81ac", "Genetic" = "#88c0d0")) +
  
  theme_minimal(base_size = 14) +
  labs(title = "R.ROSETTA Model Performance Comparison",
       x = "Voting Classifier",
       y = "Mean Accuracy (+/- Standard Deviation)",
       fill = "Algorithm") +
  
  theme(
    # Background blending
    plot.background = element_rect(fill = slide_bg, color = NA),
    panel.background = element_rect(fill = slide_bg, color = NA),
    legend.background = element_rect(fill = slide_bg, color = NA),
    legend.key = element_rect(fill = slide_bg, color = NA),
    
    # Clean up grid lines
    panel.grid.major.x = element_blank(),
    
    # TEXT COLOR UPDATES
    plot.title = element_text(color = text_col, face = "bold"),
    plot.subtitle = element_text(color = text_col),
    axis.title.x = element_text(color = text_col),
    axis.title.y = element_text(color = text_col),
    # Angle the x-axis text so the long classifier names don't overlap!
    axis.text.x = element_text(color = text_col, angle = 45, hjust = 1), 
    axis.text.y = element_text(color = text_col),
    legend.text = element_text(color = text_col),
    legend.title = element_text(color = text_col),
    
    # Style the facet headers (the grey boxes saying "5 - Fold" and "10 - Fold")
    strip.background = element_rect(fill = "#2c3e50", color = NA),
    strip.text = element_text(color = "white", face = "bold")
  )

# Show and save the plot
print(plot_rosetta)
ggsave("Rosetta_Performance_Custom.png", plot = plot_rosetta, width = 11, height = 7, dpi = 300)