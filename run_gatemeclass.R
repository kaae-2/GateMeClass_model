#!/usr/bin/env Rscript

# GateMeClass wrapper for omnibenchmark with train/test split
# Accepts train matrix, train labels, and test matrix as inputs
# Outputs predicted labels for test set

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
              help = "Path to training data matrix (gzipped)"),
  make_option("--data.train_labels", type = "character",
              help = "Path to training labels (gzipped)"),
  make_option("--data.test_matrix", type = "character",
              help = "Path to test data matrix (gzipped)"),
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

# Validate required arguments
required_args <- c("data.train_matrix", "data.train_labels", "data.test_matrix")
missing_args <- required_args[!sapply(required_args, function(x) !is.null(args[[x]]))]
if (length(missing_args) > 0) {
  stop("Missing required arguments: ", paste(missing_args, collapse = ", "))
}

# Set random seed
set.seed(args$seed)

message("GateMeClass Configuration:")
message("  Train matrix: ", args$`data.train_matrix`)
message("  Train labels: ", args$`data.train_labels`)
message("  Test matrix: ", args$`data.test_matrix`)
message("  GMM parameterization: ", args$GMM_parameterization)
message("  Sampling: ", args$sampling)
message("  k: ", args$k)
message("  Seed: ", args$seed)

# ============================================================================
# Label Loading and Conversion
# ============================================================================

load_and_convert_labels <- function(label_file, dataset_name) {
  message("Loading labels from: ", label_file)
  
  # Read numeric labels
  labels_numeric <- fread(label_file, header = FALSE)[[1]]
  message("  Found ", length(labels_numeric), " labels")
  message("  Unique values: ", paste(unique(labels_numeric), collapse = ", "))
  
  # Load mapping file (go up 4 directory levels from label file)
  base_dir <- dirname(dirname(dirname(dirname(label_file))))
  mapping_file <- file.path(base_dir, paste0(dataset_name, ".input_labels.gz"))
  
  message("  Looking for mapping at: ", mapping_file)
  
  if (!file.exists(mapping_file)) {
    stop("Mapping file not found: ", mapping_file)
  }
  
  # Read mapping (numeric ID -> cell type name)
  mapping <- fread(mapping_file, header = FALSE)
  label_map <- setNames(mapping[[2]], mapping[[1]])
  message("  Loaded mapping with ", length(label_map), " entries")
  
  # Convert numeric IDs to cell type names
  labels_converted <- sapply(labels_numeric, function(id) {
    name <- label_map[as.character(id)]
    if (is.na(name)) {
      warning(sprintf("No mapping found for ID: %s", id))
      return(NA)
    }
    return(name)
  })
  
  message("  Converted labels to cell types")
  message("  Unique cell types: ", length(unique(labels_converted)))
  
  return(list(labels = labels_converted, mapping = label_map))
}

# ============================================================================
# Load Data
# ============================================================================

# Extract dataset name from file path
dataset_name <- sub("_train\\.matrix\\.gz$", "", basename(args$`data.train_matrix`))
message("\nDataset: ", dataset_name)

# Load training data
message("\nLoading training data...")
train_dt <- fread(args$`data.train_matrix`)
message("  Train matrix: ", nrow(train_dt), " cells × ", ncol(train_dt), " markers")

# Load and convert training labels
label_result <- load_and_convert_labels(args$`data.train_labels`, dataset_name)
train_labels <- label_result$labels
label_map <- label_result$mapping

if (nrow(train_dt) != length(train_labels)) {
  stop("Train matrix rows (", nrow(train_dt), ") != train labels length (", length(train_labels), ")")
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

message("  Markers: ", paste(names(train_dt), collapse = ", "))

# ============================================================================
# Prepare Data for GateMeClass
# ============================================================================

message("\nPreparing data for GateMeClass...")

# CRITICAL: Transpose matrices (GateMeClass expects markers as ROWS, cells as COLUMNS)
train_matrix <- t(as.matrix(train_dt))
test_matrix <- t(as.matrix(test_dt))

message("  Train matrix: ", nrow(train_matrix), " markers × ", ncol(train_matrix), " cells")
message("  Test matrix: ", nrow(test_matrix), " markers × ", ncol(test_matrix), " cells")

# Prepare reference (train) data
train_reference <- data.frame(
  population = train_labels,
  stringsAsFactors = FALSE
)

# ============================================================================
# Run GateMeClass
# ============================================================================

message("\nRunning GateMeClass annotation...")

# Set training parameters
train_parameters <- list(
  GMM_parameterization = args$GMM_parameterization,
  sampling = args$sampling,
  k = args$k
)

message("  Training parameters:")
message("    GMM_parameterization: ", train_parameters$GMM_parameterization)
message("    sampling: ", train_parameters$sampling)
message("    k: ", train_parameters$k)

# Run GateMeClass annotation
result <- GateMeClass_annotate(
  exp_matrix = test_matrix,
  reference = train_matrix,
  reference_labels = train_reference,
  train_parameters = train_parameters
)

predictions <- result$labels$population
message("\nGateMeClass completed!")
message("  Predicted ", length(predictions), " labels")
message("  Unique predictions: ", length(unique(predictions)))

# ============================================================================
# Save Results
# ============================================================================

dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

# Save predicted labels (convert back to numeric IDs for metrics)
predictions_file <- file.path(args$output_dir, paste0(args$name, "_predicted_labels.txt"))

# Convert predictions back to numeric IDs for metrics compatibility
# Create reverse mapping: cell type name -> numeric ID
reverse_map <- setNames(names(label_map), label_map)
predictions_numeric <- sapply(predictions, function(name) {
  id <- reverse_map[name]
  if (is.na(id)) {
    warning(sprintf("Unknown cell type: %s", name))
    return(NA)
  }
  return(as.numeric(id))
})

write.table(
  predictions_numeric, 
  file = predictions_file, 
  quote = FALSE, 
  row.names = FALSE, 
  col.names = FALSE
)
message("Predictions saved to: ", predictions_file)

# Save marker table (optional, for debugging)
marker_table_file <- file.path(args$output_dir, paste0(args$name, "_marker_table.txt"))
write.table(
  result$marker_table,
  file = marker_table_file,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
message("Marker table saved to: ", marker_table_file)

# Save cell signatures (optional, for debugging)
signatures_file <- file.path(args$output_dir, paste0(args$name, "_signatures.txt"))
write.table(
  result$cell_signatures,
  file = signatures_file,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
message("Cell signatures saved to: ", signatures_file)

message("\nGateMeClass analysis complete!")
