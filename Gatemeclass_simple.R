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

library(argparse)
library(data.table)
library(GateMeClass)

library(argparse)
library(data.table)
library(GateMeClass)

# Parse arguments
parser <- ArgumentParser(description = "Run GateMeClass annotation")
parser$add_argument("--train.data.matrix", type="character", dest="train_matrix", 
                    help="Training matrix file (.gz)")
parser$add_argument("--labels_train", type="character", dest="train_labels",
                    help="Training labels file (.gz)")
parser$add_argument("--test.data.matrix", type="character", dest="test_matrix",
                    help="Test matrix file (.gz)")
parser$add_argument("--labels_test", type="character", dest="test_labels",
                    help="Test labels file (.gz)")
parser$add_argument("--output_dir", "-o", dest="output_dir", type="character",
                    help="Output directory", default=getwd())
parser$add_argument("--name", "-n", dest="name", type="character", 
                    help="Dataset name")

# GateMeClass parameters
parser$add_argument("--GMM_parameterization", default = "V", 
                    help = "GMM parameterization (V or E)")
parser$add_argument("--sampling", type = "double", default = 0.1, 
                    help = "Sampling fraction")
parser$add_argument("--k", type = "integer", default = 20, 
                    help = "Number of neighbors")
parser$add_argument("--seed", type = "integer", default = 42, 
                    help = "Random seed")

args <- parser$parse_args()

cat("=== GateMeClass Configuration ===\n")
cat("Train matrix:", args$train_matrix, "\n")
cat("Train labels:", args$train_labels, "\n")
cat("Test matrix:", args$test_matrix, "\n")
cat("GMM:", args$GMM_parameterization, ", sampling:", args$sampling, ", k:", args$k, "\n")
cat("=================================\n\n")

set.seed(args$seed)

# Load training data
cat("Loading training data...\n")
train_dt <- fread(args$train_matrix)
train_labels_dt <- fread(args$train_labels)

if ("col" %in% names(train_dt)) train_dt[, col := NULL]
setnames(train_dt, gsub("[^A-Za-z0-9_]", "_", names(train_dt)))

cat("  Train matrix:", nrow(train_dt), "cells ×", ncol(train_dt), "markers\n")

# Convert labels and filter unlabeled
train_labels <- as.character(train_labels_dt[[1]])
train_labels_celltype <- paste0("CellType_", train_labels)

valid_idx <- !is.na(train_labels) & train_labels != "" & train_labels != '""'
if (sum(!valid_idx) > 0) {
  cat("  Removing", sum(!valid_idx), "unlabeled cells\n")
  train_dt <- train_dt[valid_idx, ]
  train_labels_celltype <- train_labels_celltype[valid_idx]
}

cat("  Training cells:", nrow(train_dt), "\n")
cat("  Unique cell types:", length(unique(train_labels_celltype)), "\n\n")

# Load test data
cat("Loading test data...\n")
test_dt <- fread(args$test_matrix)
if ("col" %in% names(test_dt)) test_dt[, col := NULL]
setnames(test_dt, gsub("[^A-Za-z0-9_]", "_", names(test_dt)))
cat("  Test cells:", nrow(test_dt), "\n\n")

# Transpose for GateMeClass
train_matrix <- t(as.matrix(train_dt))
test_matrix <- t(as.matrix(test_dt))

# METHOD 3: Train and classify in one step
cat("=== Running GateMeClass (Method 3) ===\n")
cat("Training and classifying in one step...\n\n")

result <- GateMeClass_annotate(
  exp_matrix = test_matrix,
  marker_table = NULL,
  train_parameters = list(
    reference = train_matrix,
    labels = train_labels_celltype
  ),
  GMM_parameterization = args$GMM_parameterization,
  sampling = args$sampling,
  k = args$k,
  verbose = TRUE,
  seed = args$seed
)

cat("\n=== GateMeClass Complete ===\n")

# Extract predictions
predictions <- result$labels
predicted_ids <- as.numeric(sub("CellType_", "", predictions))

cat("Prediction summary:\n")
cat("  Unique predicted types:", length(unique(predicted_ids)), "\n")
print(table(predicted_ids))

# Save predictions
output_file <- file.path(args$output_dir, paste0(args$name, "_predicted_labels.txt"))
output_dt <- data.table(label = paste0(predicted_ids, ".0"))
fwrite(output_dt, output_file, sep = "\t", col.names = FALSE)

cat("\nPredictions saved to:", output_file, "\n")
