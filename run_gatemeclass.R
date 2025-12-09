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
  make_option("--data.label_mapping", type = "character",
              help = "Path to label mapping file (ID,cell_type)"),
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
required_args <- c("data.train_matrix", "data.train_labels", "data.test_matrix", "data.label_mapping")
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
message("Label mapping: ", args$`data.label_mapping`)
message("GMM parameterization: ", args$GMM_parameterization)
message("Sampling: ", args$sampling)
message("k: ", args$k)
message("Seed: ", args$seed)
message("=================================\n")

# ============================================================================
# Load Label Mapping
# ============================================================================

load_label_mapping <- function(mapping_file) {
  message("Loading label mapping from: ", mapping_file)
  
  if (!file.exists(mapping_file)) {
    stop("Mapping file not found: ", mapping_file)
  }
  
  # Read mapping (numeric ID, cell type name)
  mapping <- fread(mapping_file, header = FALSE, col.names = c("id", "cell_type"))
  
  message("  Loaded ", nrow(mapping), " cell type mappings:")
  for (i in seq_len(nrow(mapping))) {
    message(sprintf("    %d: %s", mapping$id[i], mapping$cell_type[i]))
  }
  
  # Create named vector: ID -> cell_type
  label_map <- setNames(mapping$cell_type, mapping$id)
  
  # Create reverse mapping: cell_type -> ID
  reverse_map <- setNames(mapping$id, mapping$cell_type)
  
  return(list(
    forward = label_map,    # ID -> cell_type
    reverse = reverse_map,   # cell_type -> ID
    mapping_df = mapping
  ))
}

# ============================================================================
# Convert Numeric Labels to Cell Type Names
# ============================================================================

convert_numeric_to_celltype <- function(numeric_labels, label_map) {
  message("Converting ", length(numeric_labels), " numeric labels to cell types...")
  
  # Convert to character for lookup
  numeric_labels_char <- as.character(numeric_labels)
  
  # Convert using mapping
  cell_types <- sapply(numeric_labels_char, function(id) {
    cell_type <- label_map[id]
    if (is.na(cell_type)) {
      warning(sprintf("No mapping found for ID: %s", id))
      return(NA)
    }
    return(cell_type)
  })
  
  # Remove names attribute
  names(cell_types) <- NULL
  
  message("  Converted to ", length(unique(cell_types)), " unique cell types")
  return(cell_types)
}

# ============================================================================
# Convert Cell Type Names to Numeric IDs
# ============================================================================

convert_celltype_to_numeric <- function(cell_types, reverse_map) {
  message("Converting ", length(cell_types), " cell type predictions to numeric IDs...")
  
  numeric_ids <- sapply(cell_types, function(cell_type) {
    id <- reverse_map[cell_type]
    if (is.na(id)) {
      warning(sprintf("No mapping found for cell type: %s", cell_type))
      return(NA)
    }
    return(id)
  })
  
  # Remove names attribute
  names(numeric_ids) <- NULL
  
  message("  Converted to numeric IDs")
  return(numeric_ids)
}

# ============================================================================
# Load Data
# ============================================================================

# Load label mapping first
label_maps <- load_label_mapping(args$`data.label_mapping`)

# Load training data
message("\nLoading training data...")
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
train_labels_celltype <- convert_numeric_to_celltype(train_labels_numeric, label_maps$forward)

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

# Convert predictions back to numeric IDs for metrics compatibility
predictions_numeric <- convert_celltype_to_numeric(predictions_celltype, label_maps$reverse)

# Save predicted labels as numeric IDs
predictions_file <- file.path(args$output_dir, paste0(args$name, ".predicted_labels.gz"))
write.table(
  predictions_numeric, 
  file = gzfile(predictions_file), 
  quote = FALSE, 
  row.names = FALSE, 
  col.names = FALSE
)
message("Predictions (numeric) saved to: ", predictions_file)

# Also save predictions as cell type names for inspection
predictions_celltype_file <- file.path(args$output_dir, 
                                       paste0(args$name, ".predicted_celltypes.txt"))
write.table(
  predictions_celltype,
  file = predictions_celltype_file,
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
message("Predictions (cell types) saved to: ", predictions_celltype_file)

# Save marker table (optional, for debugging)
marker_table_file <- file.path(args$output_dir, paste0(args$name, ".marker_table.txt"))
write.table(
  result$marker_table,
  file = marker_table_file,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
message("Marker table saved to: ", marker_table_file)

# Save cell signatures (optional, for debugging)
signatures_file <- file.path(args$output_dir, paste0(args$name, ".signatures.txt"))
write.table(
  result$cell_signatures,
  file = signatures_file,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
message("Cell signatures saved to: ", signatures_file)

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
