#!/usr/bin/env Rscript

# GateMeClass wrapper for omnibenchmark
# Compatible with simple random train/test split preprocessing

library(argparse)

# Auto-install GateMeClass if not available
if (!require("GateMeClass", quietly = TRUE)) {
  message("GateMeClass not found, installing from GitHub...")
  if (!require("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org", quiet = TRUE)
  }
  remotes::install_github("simo1c/GateMeClass", quiet = TRUE)
}

library(GateMeClass)
library(data.table)

library(argparse)
library(data.table)
library(GateMeClass)

# Parse arguments
parser <- ArgumentParser(description = "Run GateMeClass annotation")
parser$add_argument("--train_matrix", required = TRUE, help = "Training matrix file (.gz)")
parser$add_argument("--train_labels", required = TRUE, help = "Training labels file (.gz)")
parser$add_argument("--test_matrix", required = TRUE, help = "Test matrix file (.gz)")
parser$add_argument("--output", required = TRUE, help = "Output file for predictions")
parser$add_argument("--GMM_parameterization", default = "V", help = "GMM parameterization")
parser$add_argument("--sampling", type = "double", default = 0.1, help = "Sampling fraction")
parser$add_argument("--k", type = "integer", default = 20, help = "Number of neighbors")
parser$add_argument("--seed", type = "integer", default = 42, help = "Random seed")

args <- parser$parse_args()

cat("=== GateMeClass Configuration ===\n")
cat("Train matrix:", args$train_matrix, "\n")
cat("Train labels:", args$train_labels, "\n")
cat("Test matrix:", args$test_matrix, "\n")
cat("GMM:", args$GMM_parameterization, ", sampling:", args$sampling, ", k:", args$k, "\n")
cat("=================================\n\n")

# Set seed
set.seed(args$seed)

# Load training data
cat("Loading training data...\n")
train_dt <- fread(args$train_matrix)
train_labels_dt <- fread(args$train_labels)

# Remove 'col' column if present (preprocessing artifact)
if ("col" %in% names(train_dt)) {
  train_dt[, col := NULL]
}

# Standardize column names (remove special characters)
setnames(train_dt, gsub("[^A-Za-z0-9_]", "_", names(train_dt)))

cat("  Train matrix:", nrow(train_dt), "cells ×", ncol(train_dt), "markers\n")
cat("  Train labels:", nrow(train_labels_dt), "cells\n")

# Convert numeric IDs to CellType names
train_labels <- as.character(train_labels_dt[[1]])
train_labels_celltype <- paste0("CellType_", train_labels)
cat("  Unique cell types in training:", paste(sort(unique(train_labels)), collapse = ", "), "\n\n")

# Filter out unlabeled cells (empty strings, NA, or '""')
train_labels_char <- as.character(train_labels)
valid_idx <- !is.na(train_labels_char) & 
             train_labels_char != "" & 
             train_labels_char != '""'

if (sum(!valid_idx) > 0) {
  cat("  Removing", sum(!valid_idx), "unlabeled cells from training\n")
  train_dt <- train_dt[valid_idx, ]
  train_labels_celltype <- train_labels_celltype[valid_idx]
  cat("  After filtering:", nrow(train_dt), "labeled cells\n\n")
}

# Load test data
cat("Loading test data...\n")
test_dt <- fread(args$test_matrix)

# Remove 'col' column if present
if ("col" %in% names(test_dt)) {
  test_dt[, col := NULL]
}

# Standardize column names
setnames(test_dt, gsub("[^A-Za-z0-9_]", "_", names(test_dt)))

cat("  Test matrix:", nrow(test_dt), "cells ×", ncol(test_dt), "markers\n\n")

# Prepare matrices for GateMeClass (transpose to markers × cells)
cat("Preparing data for GateMeClass...\n")
train_matrix <- t(as.matrix(train_dt))
test_matrix <- t(as.matrix(test_dt))

cat("  Train matrix transposed:", nrow(train_matrix), "markers ×", ncol(train_matrix), "cells\n")
cat("  Test matrix transposed:", nrow(test_matrix), "markers ×", ncol(test_matrix), "cells\n\n")

# METHOD 2: Train model first, then classify
cat("=== Running GateMeClass Annotation (Method 2) ===\n")
cat("Step 1: Training model...\n")

# Train the model
trained_model <- GateMeClass_train(
  reference = train_matrix,
  labels = train_labels_celltype,
  marker_table = NULL,
  GMM_parameterization = args$GMM_parameterization,
  sampling = args$sampling,
  verbose = TRUE,
  seed = args$seed
)

cat("\nStep 2: Classifying test cells...\n")

# Classify test cells
predictions <- GateMeClass_classify(
  exp_matrix = test_matrix,
  trained_model = trained_model,
  k = args$k,
  verbose = TRUE,
  seed = args$seed
)

cat("\n=== GateMeClass Annotation Complete ===\n")

# Convert predictions back to numeric IDs
# Extract the numeric part from "CellType_X"
predicted_ids <- as.numeric(sub("CellType_", "", predictions))

cat("Prediction summary:\n")
cat("  Unique predicted types:", length(unique(predicted_ids)), "\n")
cat("  Predicted type counts:\n")
print(table(predicted_ids))

# Save predictions
output_dt <- data.table(label = paste0(predicted_ids, ".0"))
fwrite(output_dt, args$output, sep = "\t", col.names = FALSE)

cat("\nPredictions saved to:", args$output, "\n")


