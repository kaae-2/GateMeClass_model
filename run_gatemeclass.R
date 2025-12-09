#!/usr/bin/env Rscript

# GateMeClass wrapper for omnibenchmark with train/test split
# Compatible with preprocessing script that outputs:
#   - Numeric label IDs in train/test labels files
#   - Label mapping file (ID <-> cell type name) as {dataset}.label_mapping.gz

# Auto-install GateMeClass if not available
if (!require("GateMeClass", quietly = TRUE)) {
  message("GateMeClass not found, installing from GitHub...")
  if (!require("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org")
  }
  remotes::install_github("simo1c/GateMeClass")
}
#!/usr/bin/env Rscript

# GateMeClass wrapper for omnibenchmark - SIMPLIFIED VERSION
# Works with numeric labels directly (no label mapping file needed)

# Auto-install GateMeClass if not available
if (!require("GateMeClass", quietly = TRUE)) {
  message("GateMeClass not found, installing from GitHub...")
  if (!require("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org")
  }
  remotes::install_github("LeidenCBC/GateMeClass")
}

suppressPackageStartupMessages({
  library(GateMeClass)
  library(data.table)
  library(optparse)
})

# ============================================================================
# Argument Parser
# ============================================================================

option_list <- list(
  make_option("--output_dir", type = "character", default = ".",
              help = "Output directory [default: %default]"),
  make_option("--name", type = "character", default = "predictions",
              help = "Base name for output files [default: %default]"),
  make_option("--data.train_matrix", type = "character",
              help = "Path to training data matrix (gzipped CSV)"),
  make_option("--data.train_labels", type = "character",
              help = "Path to training labels (gzipped, numeric IDs)"),
  make_option("--data.test_matrix", type = "character",
              help = "Path to test data matrix (gzipped CSV)"),
  make_option("--GMM_parameterization", type = "character", default = "V",
              help = "GMM parameterization (V, E, or VV) [default: %default]"),
  make_option("--sampling", type = "double", default = 0.1,
              help = "Sampling fraction for training [default: %default]"),
  make_option("--k", type = "integer", default = 20,
              help = "Number of nearest neighbors [default: %default]"),
  make_option("--seed", type = "integer", default = 42,
              help = "Random seed [default: %default]")
)

parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

# Validate required arguments (label_mapping is now optional)
required_args <- c("data.train_matrix", "data.train_labels", "data.test_matrix")
missing_args <- required_args[!sapply(required_args, function(x) !is.null(args[[x]]))]
if (length(missing_args) > 0) {
  stop("Missing required arguments: ", paste(missing_args, collapse = ", "))
}

# Set random seed
set.seed(args$seed)

message("=== GateMeClass Configuration ===")
message("Train matrix: ", args$`data.train_matrix`)
message("Train labels: ", args$`data.train_labels`)
message("Test matrix: ", args$`data.test_matrix`)
message("GMM parameterization: ", args$GMM_parameterization)
message("Sampling: ", args$sampling)
message("k: ", args$k)
message("Seed: ", args$seed)
message("=================================\n")

# ============================================================================
# Load Data
# ============================================================================

# Load training data
message("Loading training data...")
train_dt <- fread(args$`data.train_matrix`)
message("  Train matrix: ", nrow(train_dt), " cells × ", ncol(train_dt), " markers")

# Load training labels (numeric IDs)
train_labels_numeric <- fread(args$`data.train_labels`, header = FALSE)[[1]]
message("  Train labels: ", length(train_labels_numeric), " cells")
message("  Unique label IDs: ", paste(sort(unique(train_labels_numeric)), collapse = ", "))

if (nrow(train_dt) != length(train_labels_numeric)) {
  stop("Train matrix rows (", nrow(train_dt), ") != train labels length (", 
       length(train_labels_numeric), ")")
}

# Convert numeric IDs to cell type names for GateMeClass
# Use simple naming: CellType_1, CellType_2, etc.
unique_ids <- sort(unique(train_labels_numeric))
id_to_name <- setNames(paste0("CellType_", unique_ids), unique_ids)
name_to_id <- setNames(unique_ids, paste0("CellType_", unique_ids))

train_labels_celltype <- sapply(train_labels_numeric, function(id) id_to_name[as.character(id)])
names(train_labels_celltype) <- NULL

message("  Converted numeric IDs to cell type names:")
for (i in seq_along(unique_ids)) {
  message(sprintf("    %s -> %s", unique_ids[i], id_to_name[as.character(unique_ids[i])]))
}

# Load test data
message("\nLoading test data...")
test_dt <- fread(args$`data.test_matrix`)
message("  Test matrix: ", nrow(test_dt), " cells × ", ncol(test_dt), " markers")

# Verify marker compatibility
if (ncol(train_dt) != ncol(test_dt)) {
  stop("Train markers (", ncol(train_dt), ") != test markers (", ncol(test_dt), ")")
}
if (!all(names(train_dt) == names(test_dt))) {
  stop("Marker names differ between train and test")
}

message("  Markers: ", paste(head(names(train_dt), 10), collapse = ", "), 
        if (ncol(train_dt) > 10) "..." else "")

# ============================================================================
# Prepare Data for GateMeClass
# ============================================================================

message("\nPreparing data for GateMeClass...")

# CRITICAL: Transpose matrices (GateMeClass expects markers as ROWS, cells as COLUMNS)
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
message("  Cell types in training: ", paste(sort(unique(train_reference$population)), collapse = ", "))

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
message("  Unique predictions: ", length(unique(predictions_celltype)))
message("  Predicted cell types: ", paste(sort(unique(predictions_celltype)), collapse = ", "))

# ============================================================================
# Save Results
# ============================================================================

dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

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

# Save predicted labels as numeric IDs
predictions_file <- file.path(args$output_dir, paste0(args$name, "_predicted_labels.txt"))
write.table(
  predictions_numeric, 
  file = predictions_file, 
  quote = FALSE, 
  row.names = FALSE, 
  col.names = FALSE
)
message("Predictions (numeric) saved to: ", predictions_file)

# Print summary statistics
message("\n=== Summary ===")
message("Test cells: ", length(predictions_numeric))
message("Predicted cell types:")
pred_table <- table(predictions_celltype)
for (i in seq_along(pred_table)) {
  message(sprintf("  %s: %d cells (%.1f%%)", 
                  names(pred_table)[i], 
                  pred_table[i],
                  100 * pred_table[i] / length(predictions_celltype)))
}
message("===============\n")

message("GateMeClass analysis complete!")
