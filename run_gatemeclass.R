
#!/usr/bin/env Rscript

# Auto-install GateMeClass with all dependencies
if (!require("GateMeClass", quietly = TRUE)) {
  message("Installing GateMeClass and dependencies...")
  
  # Install devtools if needed
  if (!require("devtools", quietly = TRUE)) {
    install.packages("devtools", repos = "https://cloud.r-project.org/")
  }
  
  # Install GateMeClass from GitHub
  devtools::install_github("simo1c/GateMeClass", quiet = TRUE)
}

library(argparse)
library(GateMeClass)
library(data.table)

# Create argument parser
parser <- ArgumentParser(description = "Run GateMeClass with train/test split")

# Required omnibenchmark arguments
parser$add_argument("--output_dir", "-o", 
                   type = "character",
                   required = TRUE,
                   help = "Output directory where results will be saved")

parser$add_argument("--name", "-n",
                   type = "character", 
                   required = TRUE,
                   help = "Name of the dataset")

# Input files - SEPARATE train and test
parser$add_argument("--data.train_matrix",
                   type = "character",
                   required = TRUE,
                   help = "Path to training data matrix (.matrix.gz)")

parser$add_argument("--data.train_labels",
                   type = "character",
                   required = TRUE,
                   help = "Path to training labels (.true_labels.gz)")

parser$add_argument("--data.test_matrix",
                   type = "character",
                   required = TRUE,
                   help = "Path to test data matrix (.matrix.gz)")

# GateMeClass-specific parameters
parser$add_argument("--GMM_parameterization",
                   type = "character",
                   default = "V",
                   choices = c("V", "E"),
                   help = "GMM variance parameter: 'V' (Variable) or 'E' (Equal)")

parser$add_argument("--reject_option",
                   action = "store_true",
                   default = FALSE,
                   help = "Try to detect cell types not defined in marker table")

parser$add_argument("--sampling",
                   type = "double",
                   default = 0.1,
                   help = "Percentage of cells used for annotation (0.0-1.0)")

parser$add_argument("--k",
                   type = "integer",
                   default = 20,
                   help = "k parameter for k-NN label refinement")

parser$add_argument("--seed",
                   type = "integer",
                   default = 1,
                   help = "Random seed")

# Parse arguments
args <- parser$parse_args()

# Create output directory
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

# Print configuration
cat(strrep("=", 70), "\n")
cat("GateMeClass Train/Test Configuration\n")
cat(strrep("=", 70), "\n")
cat("Dataset name:", args$name, "\n")
cat("Training matrix:", args$data.train_matrix, "\n")
cat("Training labels:", args$data.train_labels, "\n")
cat("Test matrix:", args$data.test_matrix, "\n")
cat("Output directory:", args$output_dir, "\n")
cat("GMM parameterization:", args$GMM_parameterization, "\n")
cat("Reject option:", args$reject_option, "\n")
cat("Sampling:", args$sampling, "\n")
cat("k:", args$k, "\n")
cat("Seed:", args$seed, "\n")
cat(strrep("=", 70), "\n\n")

# ============================================================================
# Helper function: Load and convert numeric labels to cell type names
# ============================================================================
load_and_convert_labels <- function(label_file) {
  cat("Loading labels from:", label_file, "\n")
  
  # Read labels
  if (grepl("\\.gz$", label_file)) {
    labels_raw <- readLines(gzfile(label_file))
  } else {
    labels_raw <- readLines(label_file)
  }
  
  # Try to convert to numeric
  labels_numeric <- suppressWarnings(as.numeric(labels_raw))
  
  # Remove NA values (from empty lines or headers)
  if (any(is.na(labels_numeric))) {
    na_count <- sum(is.na(labels_numeric))
    cat("  Removing", na_count, "non-numeric/empty lines\n")
    labels_raw <- labels_raw[!is.na(labels_numeric)]
    labels_numeric <- labels_numeric[!is.na(labels_numeric)]
  }
  
  # Find the mapping file
  # Path structure:
  # Label file: .../preprocessing/data_preprocessing/.hash/data_import_train.true_labels.gz
  # Mapping:    .../.hash_base/data_import.input_labels.gz
  # Need to go up 4 directory levels!
  
  base_dir <- dirname(dirname(dirname(dirname(label_file))))  # Up 4 levels
  dataset_name <- sub("_train\\.true_labels\\.gz$", "", basename(label_file))
  dataset_name <- sub("_test\\.true_labels\\.gz$", "", dataset_name)
  
  mapping_file <- file.path(base_dir, paste0(dataset_name, ".input_labels.gz"))
  
  cat("  Label file path:", label_file, "\n")
  cat("  Base directory (4 levels up):", base_dir, "\n")
  cat("  Dataset name:", dataset_name, "\n")
  cat("  Looking for mapping file at:", mapping_file, "\n")
  
  if (!file.exists(mapping_file)) {
    cat("  ERROR: Mapping file not found!\n")
    cat("  Searched at:", mapping_file, "\n")
    cat("  Please verify the path structure.\n")
    stop("Label mapping file is required but not found")
  }
  
  # Load the mapping
  cat("  ✓ Found mapping file\n")
  label_mapping <- fread(mapping_file)
  
  cat("  Mapping file contents:\n")
  print(head(label_mapping))
  
  # Convert numeric IDs to cell type names
  labels_string <- label_mapping$population[match(labels_numeric, label_mapping$label)]
  
  cat("  ✓ Converted", length(labels_numeric), "numeric labels to cell type names\n")
  cat("  Unique cell types:", length(unique(labels_string)), "\n")
  cat("  Example labels:", paste(head(unique(labels_string), 3), collapse=", "), "...\n\n")
  
  return(labels_string)
}

# ============================================================================
# Load TRAINING data
# ============================================================================
cat("Loading TRAINING data...\n")
tryCatch({
  if (!file.exists(args$data.train_matrix)) {
    stop("Training matrix not found: ", args$data.train_matrix)
  }
  
  # Read training matrix WITH EXPLICIT HEADER
  cat("  Reading training matrix with header...\n")
  train_dt <- fread(args$data.train_matrix, header = TRUE)
  
  # Convert to matrix while preserving column names
  train_matrix <- as.matrix(train_dt)
  
  # Verify column names were preserved
  cat("  Training matrix dimensions:", nrow(train_matrix), "rows x", ncol(train_matrix), "cols\n")
  cat("  Column names preserved?:", !is.null(colnames(train_matrix)), "\n")
  cat("  Training markers:", paste(colnames(train_matrix), collapse = ", "), "\n\n")
  
  if (is.null(colnames(train_matrix)) || all(colnames(train_matrix) == "")) {
    stop("ERROR: Column names not found in training matrix!")
  }
  
}, error = function(e) {
  stop("Error loading training matrix: ", e$message)
})

# Load training labels
train_labels <- load_and_convert_labels(args$data.train_labels)

if (length(train_labels) != nrow(train_matrix)) {
  stop("Number of training labels (", length(train_labels), 
       ") does not match number of training cells (", nrow(train_matrix), ")")
}

cat("Training set cell type distribution:\n")
print(table(train_labels))
cat("\n")

# ============================================================================
# Load TEST data
# ============================================================================
cat("Loading TEST data...\n")
tryCatch({
  if (!file.exists(args$data.test_matrix)) {
    stop("Test matrix not found: ", args$data.test_matrix)
  }
  
  # Read test matrix WITH EXPLICIT HEADER
  cat("  Reading test matrix with header...\n")
  test_dt <- fread(args$data.test_matrix, header = TRUE)
  
  # Convert to matrix while preserving column names
  test_matrix <- as.matrix(test_dt)
  
  # Verify column names were preserved
  cat("  Test matrix dimensions:", nrow(test_matrix), "rows x", ncol(test_matrix), "cols\n")
  cat("  Column names preserved?:", !is.null(colnames(test_matrix)), "\n")
  cat("  Test markers:", paste(colnames(test_matrix), collapse = ", "), "\n\n")
  
  if (is.null(colnames(test_matrix)) || all(colnames(test_matrix) == "")) {
    stop("ERROR: Column names not found in test matrix!")
  }
  
  # Verify markers match between train and test
  cat("Verifying train/test marker compatibility...\n")
  cat("  Train markers:", paste(colnames(train_matrix), collapse = ", "), "\n")
  cat("  Test markers:", paste(colnames(test_matrix), collapse = ", "), "\n")
  cat("  Markers identical?:", identical(colnames(train_matrix), colnames(test_matrix)), "\n\n")
  
  if (!identical(colnames(train_matrix), colnames(test_matrix))) {
    cat("ERROR: Training and test matrices have different markers!\n")
    cat("Train only:", setdiff(colnames(train_matrix), colnames(test_matrix)), "\n")
    cat("Test only:", setdiff(colnames(test_matrix), colnames(train_matrix)), "\n")
    stop("Marker mismatch between training and test data")
  }
  
}, error = function(e) {
  stop("Error loading test matrix: ", e$message)
})

# ============================================================================
# Train on TRAINING set, Predict on TEST set
# ============================================================================
cat(strrep("=", 70), "\n")
cat("TRAINING AND PREDICTION\n")
cat(strrep("=", 70), "\n")
cat("Training on", nrow(train_matrix), "cells...\n")
cat("Will predict on", nrow(test_matrix), "test cells...\n\n")

res <- GateMeClass_annotate(
  exp_matrix = test_matrix,  # ← Predict on TEST set only!
  marker_table = NULL,
  train_parameters = list(
    reference = train_matrix,
    labels = train_labels
  ),
  GMM_parameterization = args$GMM_parameterization,
  reject_option = args$reject_option,
  sampling = args$sampling,
  k = args$k,
  verbose = TRUE,
  seed = args$seed
)

# ============================================================================
# Save results
# ============================================================================
cat("\n", strrep("=", 70), "\n")
cat("RESULTS\n")
cat(strrep("=", 70), "\n")
cat("Test set predictions:\n")
print(table(res$labels))
cat("\n")

# Main output: predictions for test cells
prediction_file <- file.path(args$output_dir, 
                             paste0(args$name, "_predicted_labels.txt"))

writeLines(res$labels, prediction_file)
cat("✓ Predicted labels saved to:", prediction_file, "\n")
cat("  Format: Plain text, one label per line (", length(res$labels), "test cells)\n")

# Auxiliary outputs
if (!is.null(res$marker_table)) {
  marker_table_path <- file.path(args$output_dir, 
                                  paste0(args$name, ".marker_table.gz"))
  fwrite(res$marker_table, marker_table_path, sep = "\t", compress = "gzip")
  cat("  Marker table saved to:", marker_table_path, "\n")
}

signatures_path <- file.path(args$output_dir, 
                             paste0(args$name, ".cell_signatures.gz"))
fwrite(res$cell_signatures, signatures_path, sep = "\t", compress = "gzip")
cat("  Cell signatures saved to:", signatures_path, "\n")

cat("\n✓ GateMeClass completed successfully!\n")
