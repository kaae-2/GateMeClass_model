#!/usr/bin/env Rscript

# GateMeClass wrapper for omnibenchmark
# Compatible with simple random train/test split preprocessing
# WITH outlier handling to prevent GMM failures

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

# ============================================================================
# Argument Parser
# ============================================================================

parser <- ArgumentParser(description="GateMeClass caller with transformation and outlier handling")

parser$add_argument('--train.data.matrix', type="character", dest="train_matrix")
parser$add_argument('--labels_train', type="character", dest="train_labels")
parser$add_argument('--test.data.matrix', type="character", dest="test_matrix")
parser$add_argument('--labels_test', type="character", dest="test_labels")
parser$add_argument('--seed', type="integer", default=42)
parser$add_argument("--output_dir", "-o", dest="output_dir", type="character", default=getwd())
parser$add_argument("--name", "-n", dest="name", type="character")
parser$add_argument("--GMM_parameterization", type="character", default="V")
parser$add_argument("--sampling", type="double", default=0.1)
parser$add_argument("--k", type="integer", default=20)
parser$add_argument("--cofactor", type="double", default=5)

args <- parser$parse_args()
set.seed(args$seed)

message("=== GateMeClass Configuration ===")
message("Train matrix: ", args$train_matrix)
message("Train labels: ", args$train_labels)
message("Test matrix: ", args$test_matrix)
message("GMM: ", args$GMM_parameterization, ", sampling: ", args$sampling, ", k: ", args$k)
message("Cofactor: ", args$cofactor, " (arcsinh transformation)")
message("=================================\n")

# ============================================================================
# Load Data
# ============================================================================

message("Loading training data...")
train_dt <- fread(args$train_matrix)
train_labels_numeric <- fread(args$train_labels, header=FALSE)[[1]]

# Standardize column names (two-step process)
names(train_dt) <- gsub("[^A-Za-z0-9_]", "_", names(train_dt))
names(train_dt) <- gsub("_+", "_", names(train_dt))  # Remove duplicate underscores

message("  Original train matrix: ", nrow(train_dt), " cells × ", ncol(train_dt), " columns")
message("  Original train labels: ", length(train_labels_numeric), " cells")

# ============================================================================
# Clean Training Data
# ============================================================================

# Remove 'col' column if it exists
if ("col" %in% names(train_dt)) {
  train_dt <- train_dt[, !names(train_dt) %in% "col", with=FALSE]
  message("  Removed 'col' column from training data")
}

# Remove unlabeled cells
train_labels_char <- as.character(train_labels_numeric)
valid_train_idx <- !is.na(train_labels_char) & train_labels_char != "" & train_labels_char != '""'
n_unlabeled <- sum(!valid_train_idx)

if (n_unlabeled > 0) {
  message("  Removing ", n_unlabeled, " unlabeled cells from training set")
  train_dt <- train_dt[valid_train_idx, ]
  train_labels_numeric <- train_labels_numeric[valid_train_idx]
  train_labels_char <- train_labels_char[valid_train_idx]
}

message("  Clean train matrix: ", nrow(train_dt), " cells × ", ncol(train_dt), " markers")

# Convert to cell type names
unique_ids <- sort(unique(train_labels_char))
id_to_name <- setNames(paste0("CellType_", unique_ids), unique_ids)
name_to_id <- setNames(unique_ids, paste0("CellType_", unique_ids))

train_labels_celltype <- id_to_name[train_labels_char]
names(train_labels_celltype) <- NULL

message("  Cell types: ", length(unique(train_labels_celltype)), " unique types")

# ============================================================================
# Load and Clean Test Data
# ============================================================================

message("\nLoading test data...")
test_dt <- fread(args$test_matrix)

names(test_dt) <- gsub("[^A-Za-z0-9_]", "_", names(test_dt))
names(test_dt) <- gsub("_+", "_", names(test_dt))

if ("col" %in% names(test_dt)) {
  test_dt <- test_dt[, !names(test_dt) %in% "col", with=FALSE]
}

message("  Test matrix: ", nrow(test_dt), " cells × ", ncol(test_dt), " markers")

# ============================================================================
# Apply Transformation WITH OUTLIER HANDLING
# ============================================================================

message("\nPreparing data for GateMeClass...")

train_matrix_raw <- as.matrix(train_dt)
test_matrix_raw <- as.matrix(test_dt)

# Handle extreme outliers BEFORE transformation
message("Handling outliers...")
upper_bound <- quantile(as.vector(train_matrix_raw), 0.999, na.rm = TRUE)
message("  99.9th percentile: ", round(upper_bound, 2))
message("  Max before capping: ", round(max(train_matrix_raw, na.rm = TRUE), 2))

# Cap extreme values at 99.9th percentile
train_matrix_raw[train_matrix_raw > upper_bound] <- upper_bound
test_matrix_raw[test_matrix_raw > upper_bound] <- upper_bound

# Remove any negative values (shouldn't exist in CyTOF but just in case)
train_matrix_raw[train_matrix_raw < 0] <- 0
test_matrix_raw[test_matrix_raw < 0] <- 0

message("  Max after capping: ", round(max(train_matrix_raw, na.rm = TRUE), 2))

# Now apply arcsinh transformation
message("Applying arcsinh transformation (cofactor = ", args$cofactor, ")...")
train_matrix_transformed <- asinh(train_matrix_raw / args$cofactor)
test_matrix_transformed <- asinh(test_matrix_raw / args$cofactor)

message("  Final range: [", round(min(train_matrix_transformed, na.rm = TRUE), 2), ", ", 
        round(max(train_matrix_transformed, na.rm = TRUE), 2), "]")

# Transpose for GateMeClass (markers as rows, cells as columns)
train_matrix <- t(train_matrix_transformed)
test_matrix <- t(test_matrix_transformed)

message("  Train: ", nrow(train_matrix), " markers × ", ncol(train_matrix), " cells")
message("  Test: ", nrow(test_matrix), " markers × ", ncol(test_matrix), " cells")

# ============================================================================
# Run GateMeClass
# ============================================================================

message("\n=== Running GateMeClass (Method 3) ===\n")

res <- GateMeClass_annotate(
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

message("\n=== GateMeClass Complete ===")

# ============================================================================
# Save Results
# ============================================================================

predictions_celltype <- res$labels
predictions_numeric <- as.numeric(name_to_id[predictions_celltype])

res_char <- paste0(predictions_numeric, ".0")
res_char[is.na(predictions_numeric)] <- NA

outfile <- file.path(args$output_dir, paste0(args$name, "_predicted_labels.txt"))
write.table(file=outfile, res_char, col.names=FALSE, row.names=FALSE, quote=FALSE, na='""')

message("Predictions saved to: ", outfile)
message("\nPrediction summary: ", length(unique(predictions_numeric)), " unique types")

message("\nGateMeClass analysis complete!")
