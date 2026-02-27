# Set working directory to where the dataset is located
setwd("C:/Users/david/OneDrive - Uppsala universitet/Dokument/Uppsala Universitet/Biologi och bioteknik/Kunskapsbaserade system/Projekt")

# Load appropriate packages
library(rJava)
library(rmcfs)
library(R.ROSETTA)
library(devtools)
library(csVisuNet)

# Read the data
data <- read.csv("Project2.csv", sep="\t", colClasses = "character")

# delete the id column
data$id <- NULL

# Calculate the threshold: Find columns where LESS THAN 50% of the data is "?"
# (You can adjust the 0.50 to whatever threshold you prefer)
cols_to_keep <- sapply(data, function(col) {
  sum(col == "?") / length(col) < 0.50 
})

# Subset the data to keep only the good columns
data_clean <- data[, cols_to_keep]

data_clean$Host
table(data_clean$Host)

# MCFS
mcfs_result <- mcfs(Host ~ ., 
                    data = data_clean, 
                    projections = 1500, 
                    projectionSize = 0.1, 
                    splits = 5,
                    threadsNumber = 8,
                    cutoffPermutations = 5
                    )

# Show selected features and plot convergence graph
head(mcfs_result$RI)
plot(mcfs_result, type="distances")

# Plot relative importance graph and save features with RI>cutoff_value
plot(mcfs_result, type = "ri", size = 50)
mcfs_result2 <- mcfs_result$RI[1:mcfs_result$cutoff_value,]
mcfs_result2

# 5. Extract only the significant features
top_features <- mcfs_result$RI[1:mcfs_result$cutoff_value, "attribute"]
print(top_features)

# Create a subset of the data containing ONLY the top features and the target 'Host'
selected_data <- data_clean[, c(as.character(top_features), "Host")]

selected_data

# Train the rule-based model using 10-fold cross-validation
rosetta_model <- rosetta(selected_data, 
                         discrete = TRUE,
                         reducer = "Johnson"
                         #classifier = "StandardVoter",
                         #roc = TRUE,
                       #  cv = 10,
                         #clroc = "Human"
                       )

# Retrieve rules and model performance (Accuracy, Sensitivity, Specificity)
rules <- rosetta_model$main
qual <- rosetta_model$quality
print(rules)
print(qual)

# Show discretization levels (bins) used in the rules
slevels <- rosetta_model$main$levels
slevels

# 1. Filter rules with a significant p-value
sig_rules <- rosetta_model$main[rosetta_model$main$pValue < 0.05, ]
sig_rules

# If these functions give errors when running, download from github to wd and load them manually 
#source("recalculateRules.R")
#source("clusterRules.R")

# Recalculate rules
recalc_rules <- recalculateRules(selected_data, rosetta_model$main, discrete=TRUE)

# NOTE: the clusterRules function takes training df and recalculated rules as arguments
# Previously, the sig_rules table was used as an argument instead of training df, which caused issues in the figure
clustered_rules <- clusterRules(selected_data, recalc_rules)

# Link to CytoScape (open CytoScape before running)
library(RCy3)
cytoscapePing ()
vis <- visunetcyto(sig_rules)