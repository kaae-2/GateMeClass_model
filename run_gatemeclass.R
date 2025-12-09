#!/usr/bin/env Rscript

# GateMeClass wrapper for omnibenchmark - Simple R script
# Works like the random baseline - no Python wrapper needed

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
# Argument Parser (matching random baseline pattern)
# ============================================================================

parser <- ArgumentParser(description="GateMeClass caller")

parser$add_argument('--train.data.matrix',
                    type="character",
                    help='gz-compressed CSV file with training data')
parser$add_argument('--labels_train',
                    type="character",
                    help='gz-compressed file with training labels (numeric IDs)')
parser$add_argument('--test.data.matrix',
                    type="character",
                    help='gz-compressed CSV file with test data')
parser$add_argument('--labels_test',
                    type="character",
                    help='gz-compressed file with test labels (for reference only)')
parser$add_argument('--seed',
                    type="integer",
                    help='Random seed',
                    default=42)
parser$add_argument("--output_dir", "-o",
                    dest="output_dir",
                    type="character",
                    help="output directory where files will be saved",
                    default=getwd())
parser$add_argument("--name", "-n",
                    dest="name",
                    type="character",
                    help="name of the dataset")
parser$add_argument("--GMM_parameterization",
                    type="character",
                    default="V",
                    help="GMM parameterization (V, E, or VV)")
parser$add_argument("--sampling",
                    type="double",
                    default=0.1,
                    help="Sampling fraction for training")
parser$add_argument("--k",
                    type="integer",
                    default=20,
                    help="Number of nearest neighbors")

args <- parser$parse_args()

# Set random seed
set.seed(args$seed)

message("=== GateMeClass Configuration ===")
message("Train matrix: ", args$`train.data.matrix`)
message("Train labels: ", args$labels_train)
message("Test matrix: ", args$`test.data.matrix`)
message("GMM parameterization: ", args$GMM_parameterization)
message("Sampling: ", args$sampling)
message("k: ", args$k)
message("Seed: ", args$seed)
message("=================================\n")

# ============================================================================
# Load Data
# ============================================================================

message("Loading training data...")
train_dt <- fread(args$`train.data.matrix`)
message("  Train matrix: ", nrow(train_dt), " cells × ", ncol(train_dt), " markers")

# Load training labels (numeric IDs)
train_labels_numeric <- fread(args$labels_train, header=FALSE)[[1]]
message("  Train labels: ", length(train_labels_numeric), " cells")
message("  Unique label IDs: ", paste(sort(unique(train_labels_numeric)), collapse=", "))

if (nrow(train_dt) != length(train_labels_numeric)) {
  stop("Train matrix rows (", nrow(train_dt), ") != train labels length (", 
       length(train_labels_numeric), ")")
}

# Convert numeric IDs to cell type names for GateMeClass
# Use simple naming: CellType_1, CellType_2, etc.
unique_ids <- sort(unique(train_labels_numeric))
id_to_name <- setNames(paste0("CellType_", unique_ids), unique_ids)
name_to_id <- setNames(unique_ids, paste0("CellType_", unique_ids))

train_labels_celltype <- sapply(train_labels_numeric, function(id) {
  id_to_name[as.character(id)]
})
names(train_labels_celltype) <- NULL

message("  Converted numeric IDs to cell type names")

# Load test data
message("\nLoading test data...")
test_dt <- fread(args$`test.data.matrix`)
message("  Test matrix: ", nrow(test_dt), " cells × ", ncol(test_dt), " markers")

# Verify marker compatibility
if (ncol(train_dt) != ncol(test_dt)) {
  stop("Train markers (", ncol(train_dt), ") != test markers (", ncol(test_dt), ")")
}
if (!all(names(train_dt) == names(test_dt))) {
  stop("Marker names differ between train and test")
}

# ============================================================================
# Prepare Data for GateMeClass
# ============================================================================

message("\nPreparing data for GateMeClass...")

# Transpose matrices (GateMeClass expects markers as ROWS, cells as COLUMNS)
train_matrix <- t(as.matrix(train_dt))
test_matrix <- t(as.matrix(test_dt))

message("  Train matrix transposed: ", nrow(train_matrix), " markers × ", 
        ncol(train_matrix), " cells")
message("  Test matrix transposed: ", nrow(test_matrix), " markers × ", 
        ncol(test_matrix), " cells")

# Prepare reference (train) data with cell type labels
train_reference <- data.frame(
  population = train_labels_celltype,
  stringsAsFactors = FALSE
)

message("  Training reference prepared with ", nrow(train_reference), " cells")

# ============================================================================
# Run GateMeClass
# ============================================================================

message("\n=== Running GateMeClass Annotation ===")

# Set training parameters
train_parameters <- list(
  GMM_parameterization = args$GMM_parameterization,
  sampling = args$sampling,
  k = args$k
)

message("Training parameters:")
message("  GMM_parameterization: ", train_parameters$GMM_parameterization)
message("  sampling: ", train_parameters$sampling)
message("  k: ", train_parameters$k)

# Run GateMeClass annotation
message("\nAnnotating test data...")
result <- GateMeClass_annotate(
  exp_matrix = test_matrix,
  reference = train_matrix,
  reference_labels = train_reference,
  train_parameters = train_parameters
)

predictions_celltype <- result$labels$population
message("\nGateMeClass completed!")
message("  Predicted ", length(predictions_celltype), " labels")

# ============================================================================
# Save Results (matching random baseline output format)
# ============================================================================

message("\nSaving results...")

# Convert predictions back to numeric IDs
predictions_numeric <- sapply(predictions_celltype, function(name) {
  id <- name_to_id[name]
  if (is.na(id)) {
    warning(sprintf("Unknown cell type: %s", name))
    return(NA)
  }
  return(as.numeric(id))
})
names(predictions_numeric) <- NULL

# Format as character with .0 suffix (matching random baseline)
res_char <- as.character(predictions_numeric)
res_char <- paste0(res_char, ".0")
res_char[is.na(predictions_numeric)] <- NA

# Save predicted labels (matching random baseline output format)
outfile <- file.path(args$output_dir, paste0(args$name, "_predicted_labels.txt"))
write.table(file=outfile, res_char, 
            col.names=FALSE, row.names=FALSE, quote=FALSE, na='""')

message("Predictions saved to: ", outfile)
message("\nGateMeClass analysis complete!")
