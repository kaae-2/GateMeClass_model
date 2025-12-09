#!/usr/bin/env Rscript

# GateMeClass wrapper for omnibenchmark
# Handles unlabeled cells in training data

library(argparse)

# Auto-install GateMeClass if not available
if (!require("GateMeClass", quietly = TRUE)) {
  message("GateMeClass not found, installing from GitHub...")
  if (!require("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org")
  }
  remotes::install_github("simo1c/GateMeClass")
}

library(GateMeClass)
library(data.table)

# ============================================================================
# Argument Parser
# ============================================================================

parser <- ArgumentParser(description="GateMeClass caller")

parser$add_argument('--train.data.matrix', type="character")
parser$add_argument('--labels_train', type="character")
parser$add_argument('--test.data.matrix', type="character")
parser$add_argument('--labels_test', type="character")
parser$add_argument('--seed', type="integer", default=42)
parser$add_argument("--output_dir", "-o", dest="output_dir", type="character", default=getwd())
parser$add_argument("--name", "-n", dest="name", type="character")
parser$add_argument("--GMM_parameterization", type="character", default="V")
parser$add_argument("--sampling", type="double", default=0.1)
parser$add_argument("--k", type="integer", default=20)

args <- parser$parse_args()
set.seed(args$seed)

message("=== GateMeClass Configuration ===")
message("Train matrix: ", args$`train.data.matrix`)
message("Train labels: ", args$labels_train)
message("Test matrix: ", args$`test.data.matrix`)
message("GMM: ", args$GMM_parameterization, ", sampling: ", args$sampling, ", k: ", args$k)
message("=================================\n")

# ============================================================================
# Load Data
# ============================================================================

message("Loading training data...")
train_dt <- fread(args$`train.data.matrix`)
train_labels_numeric <- fread(args$labels_train, header=FALSE)[[1]]

message("  Original train matrix: ", nrow(train_dt), " cells × ", ncol(train_dt), " columns")
message("  Original train labels: ", length(train_labels_numeric), " cells")

# ============================================================================
# Clean Training Data
# ============================================================================

# Remove 'col' column if it exists (bug in preprocessing)
if ("col" %in% names(train_dt)) {
  train_dt <- train_dt[, !names(train_dt) %in% "col", with=FALSE]
  message("  Removed 'col' column from training data")
}

# Remove unlabeled cells from training (empty strings)
valid_train_idx <- train_labels_numeric != ""
n_unlabeled <- sum(!valid_train_idx)

if (n_unlabeled > 0) {
  message("  Removing ", n_unlabeled, " unlabeled cells from training set")
  train_dt <- train_dt[valid_train_idx, ]
  train_labels_numeric <- train_labels_numeric[valid_train_idx]
}

message("  Clean train matrix: ", nrow(train_dt), " cells × ", ncol(train_dt), " markers")
message("  Clean train labels: ", length(train_labels_numeric), " cells")

# Convert numeric IDs to cell type names
unique_ids <- sort(unique(train_labels_numeric))
id_to_name <- setNames(paste0("CellType_", unique_ids), unique_ids)
name_to_id <- setNames(unique_ids, paste0("CellType_", unique_ids))

train_labels_celltype <- sapply(train_labels_numeric, function(id) {
  id_to_name[as.character(id)]
})
names(train_labels_celltype) <- NULL

message("  Cell types in training: ", paste(sort(unique(train_labels_celltype)), collapse=", "))

# ============================================================================
# Load and Clean Test Data
# ============================================================================

message("\nLoading test data...")
test_dt <- fread(args$`test.data.matrix`)

# Remove 'col' column if it exists
if ("col" %in% names(test_dt)) {
  test_dt <- test_dt[, !names(test_dt) %in% "col", with=FALSE]
  message("  Removed 'col' column from test data")
}

message("  Test matrix: ", nrow(test_dt), " cells × ", ncol(test_dt), " markers")

# Verify marker names match
if (!all(names(train_dt) == names(test_dt))) {
  stop("Marker names differ between train and test!")
}

# ============================================================================
# Prepare Data for GateMeClass
# ============================================================================

message("\nPreparing data for GateMeClass...")

# Transpose matrices (GateMeClass expects markers as ROWS, cells as COLUMNS)
train_matrix <- t(as.matrix(train_dt))
test_matrix <- t(as.matrix(test_dt))

message("  Train matrix transposed: ", nrow(train_matrix), " markers × ", ncol(train_matrix), " cells")
message("  Test matrix transposed: ", nrow(test_matrix), " markers × ", ncol(test_matrix), " cells")

# ============================================================================
# Run GateMeClass
# ============================================================================

message("\n=== Running GateMeClass Annotation ===")
message("Using Method 3: Training and classification in one step")

# Use train_parameters with reference and labels
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

predictions_celltype <- res$labels
message("\nGateMeClass completed!")
message("  Predicted ", length(predictions_celltype), " labels")
message("  Unique predictions: ", length(unique(predictions_celltype)))

# ============================================================================
# Save Results
# ============================================================================

message("\nSaving results...")

# Convert back to numeric IDs
predictions_numeric <- sapply(predictions_celltype, function(name) {
  id <- name_to_id[name]
  if (is.na(id)) {
    warning(sprintf("Unknown cell type: %s", name))
    return(NA)
  }
  return(as.numeric(id))
})
names(predictions_numeric) <- NULL

# Format with .0 suffix (matching random baseline output format)
res_char <- as.character(predictions_numeric)
res_char <- paste0(res_char, ".0")
res_char[is.na(predictions_numeric)] <- NA

# Save predictions
outfile <- file.path(args$output_dir, paste0(args$name, "_predicted_labels.txt"))
write.table(file=outfile, res_char, 
            col.names=FALSE, row.names=FALSE, quote=FALSE, na='""')

message("Predictions saved to: ", outfile)

# Show prediction distribution
message("\n=== Prediction Summary ===")
pred_table <- table(predictions_celltype)
for (i in seq_along(pred_table)) {
  message(sprintf("  %s: %d cells (%.1f%%)", 
                  names(pred_table)[i], 
                  pred_table[i],
                  100 * pred_table[i] / length(predictions_celltype)))
}
message("==========================\n")

message("GateMeClass analysis complete!")
