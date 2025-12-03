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


# Automated cell type annotation for flow/CyTOF data
# Output format compatible with CyAnno and flow_metrics

library(argparse)
library(GateMeClass)
library(data.table)

# Create argument parser
parser <- ArgumentParser(description = "Run GateMeClass annotation on flow/CyTOF data")

# Required omnibenchmark arguments
parser$add_argument("--output_dir", "-o", 
                   type = "character",
                   required = TRUE,
                   help = "Output directory where results will be saved")

parser$add_argument("--name", "-n",
                   type = "character", 
                   required = TRUE,
                   help = "Name of the dataset")

# Input files (match your YAML stage inputs)
parser$add_argument("--data.matrix",
                   type = "character",
                   required = TRUE,
                   help = "Path to preprocessed data matrix (.matrix.gz)")

parser$add_argument("--data.true_labels",
                   type = "character",
                   required = FALSE,
                   default = NULL,
                   help = "Path to true labels file (.true_labels.gz) - optional, for training mode")

# GateMeClass-specific parameters
parser$add_argument("--marker_table",
                   type = "character",
                   required = FALSE,
                   default = NULL,
                   help = "Path to marker table file (optional)")

parser$add_argument("--mode",
                   type = "character",
                   default = "train_and_annotate",
                   choices = c("annotate", "train", "train_and_annotate"),
                   help = "GateMeClass mode: 'annotate' (with marker table), 'train' (extract marker table), or 'train_and_annotate'")

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

parser$add_argument("--narrow_marker_table",
                   action = "store_true",
                   default = FALSE,
                   help = "Use narrow marker table format")

parser$add_argument("--seed",
                   type = "integer",
                   default = 1,
                   help = "Random seed")

# Parse arguments
args <- parser$parse_args()

# Create output directory if it doesn't exist
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

# Print configuration
cat(strrep("=", 70), "\n")
cat("GateMeClass Configuration\n")
cat(strrep("=", 70), "\n")
cat("Dataset name:", args$name, "\n")
cat("Mode:", args$mode, "\n")
cat("Input matrix:", args$data.matrix, "\n")
cat("Output directory:", args$output_dir, "\n")
cat("GMM parameterization:", args$GMM_parameterization, "\n")
cat("Reject option:", args$reject_option, "\n")
cat("Sampling:", args$sampling, "\n")
cat("k:", args$k, "\n")
cat("Seed:", args$seed, "\n")
cat(strrep("=", 70), "\n\n")

# Load data matrix with error handling
cat("Loading data matrix...\n")
tryCatch({
  if (!file.exists(args$data.matrix)) {
    stop("Data matrix file not found: ", args$data.matrix)
  }
  
  # Read the matrix file (gzipped CSV)
  dt <- fread(args$data.matrix)
  
  # Check if first column contains cell IDs (non-numeric)
  if (ncol(dt) > 0 && !is.numeric(dt[[1]])) {
    cat("Detected cell IDs in first column\n")
    cell_ids <- dt[[1]]
    dt <- dt[, -1, with = FALSE]
    exp_matrix <- as.matrix(dt)
    rownames(exp_matrix) <- cell_ids
  } else {
    exp_matrix <- as.matrix(dt)
  }
  
  # Ensure column names are preserved (marker names)
  if (is.null(colnames(exp_matrix))) {
    warning("No column names detected - markers will be numbered")
    colnames(exp_matrix) <- paste0("Marker_", 1:ncol(exp_matrix))
  }
  
  cat("Matrix dimensions:", nrow(exp_matrix), "cells x", ncol(exp_matrix), "markers\n")
  cat("Markers:", paste(colnames(exp_matrix), collapse = ", "), "\n\n")
  
}, error = function(e) {
  stop("Error loading data matrix: ", e$message)
})

# Helper function to load and convert labels
load_and_convert_labels <- function(label_file) {
  # Read numeric labels
  if (grepl("\\.gz$", label_file)) {
    labels_numeric <- readLines(gzfile(label_file))
  } else {
    labels_numeric <- readLines(label_file)
  }
  labels_numeric <- as.numeric(labels_numeric)
  
  # Find the mapping file
  # Try multiple possible locations
  mapping_paths <- c(
    gsub("/preprocessing/data_preprocessing/default/.*\\.true_labels\\.gz", ".input_labels.gz", label_file),
    gsub("\\.true_labels\\.gz", ".input_labels.gz", label_file),
    gsub("/data_preprocessing/default/", "/", gsub("\\.true_labels\\.gz", ".input_labels.gz", label_file))
  )
  
  label_mapping <- NULL
  mapping_file_used <- NULL
  for (path in mapping_paths) {
    if (file.exists(path)) {
      cat("Found label mapping at:", path, "\n")
      label_mapping <- fread(path)
      mapping_file_used <- path
      break
    }
  }
  
  if (is.null(label_mapping)) {
    cat("WARNING: Could not find label mapping file.\n")
    cat("Tried paths:\n")
    for (p in mapping_paths) cat("  -", p, "\n")
    cat("Returning numeric labels as strings...\n")
    return(as.character(labels_numeric))
  }
  
  # Convert numeric IDs to cell type names
  labels_string <- label_mapping$population[match(labels_numeric, label_mapping$label)]
  
  cat("Converted", length(labels_numeric), "numeric labels to cell type names\n")
  cat("Unique cell types:", length(unique(labels_string)), "\n")
  cat("Example labels:", paste(head(unique(labels_string), 5), collapse=", "), "...\n")
  
  return(labels_string)
}

# Initialize results
res <- NULL
marker_table <- NULL

# Execute based on mode
if (args$mode == "train") {
  # Training mode: extract marker table from annotated data
  cat("Running in TRAINING mode - extracting marker table\n")
  
  if (is.null(args$data.true_labels)) {
    stop("Training mode requires --data.true_labels argument")
  }
  
  # Load true labels
  cat("Loading true labels...\n")
  tryCatch({
    if (!file.exists(args$data.true_labels)) {
      stop("True labels file not found: ", args$data.true_labels)
    }
    
    # Load and convert labels
    true_labels <- load_and_convert_labels(args$data.true_labels)
    
    # Check if number of labels matches number of cells
    if (length(true_labels) != nrow(exp_matrix)) {
      stop("Number of labels (", length(true_labels), 
           ") does not match number of cells (", nrow(exp_matrix), ")")
    }
    
  }, error = function(e) {
    stop("Error loading true labels: ", e$message)
  })
  
  # Remove any unassigned cells
  assigned_idx <- true_labels != "unassigned"
  cat("Found", sum(!assigned_idx), "unassigned cells (will be removed from training)\n")
  exp_matrix <- exp_matrix[assigned_idx, , drop = FALSE]
  true_labels <- true_labels[assigned_idx]
  
  cat("Training on", length(true_labels), "labeled cells\n")
  cat("Cell types:", length(unique(true_labels)), "\n")
  cat("Cell type distribution:\n")
  print(table(true_labels))
  cat("\n")
  
  # Train GateMeClass
  cat("Extracting marker table...\n")
  marker_table <- GateMeClass_train(
    reference = exp_matrix,
    labels = true_labels,
    GMM_parameterization = args$GMM_parameterization,
    verbose = TRUE,
    seed = args$seed
  )
  
  cat("\nMarker table extracted successfully!\n\n")
  print(marker_table)
  
  # Save marker table (auxiliary output)
  marker_table_path <- file.path(args$output_dir, 
                                  paste0(args$name, ".marker_table.gz"))
  fwrite(marker_table, marker_table_path, sep = "\t", compress = "gzip")
  cat("Marker table saved to:", marker_table_path, "\n")
  
} else if (args$mode == "annotate") {
  # Annotation mode: use provided marker table
  cat("Running in ANNOTATION mode - using provided marker table\n")
  
  if (is.null(args$marker_table)) {
    stop("Annotation mode requires --marker_table argument")
  }
  
  # Load marker table
  cat("Loading marker table...\n")
  tryCatch({
    if (!file.exists(args$marker_table)) {
      stop("Marker table file not found: ", args$marker_table)
    }
    marker_table <- as.data.frame(fread(args$marker_table))
    marker_table[is.na(marker_table)] <- "*"  # Replace NA with wildcards
    
  }, error = function(e) {
    stop("Error loading marker table: ", e$message)
  })
  
  cat("Marker table loaded with", nrow(marker_table), "cell types\n\n")
  print(marker_table)
  
  # Annotate
  cat("\nAnnotating cells...\n")
  res <- GateMeClass_annotate(
    exp_matrix = exp_matrix,
    marker_table = marker_table,
    GMM_parameterization = args$GMM_parameterization,
    reject_option = args$reject_option,
    sampling = args$sampling,
    k = args$k,
    verbose = TRUE,
    narrow_marker_table = args$narrow_marker_table,
    seed = args$seed
  )
  
} else if (args$mode == "train_and_annotate") {
  # Combined mode: train and annotate in one step
  cat("Running in TRAIN_AND_ANNOTATE mode\n")
  
  if (is.null(args$data.true_labels)) {
    stop("Train and annotate mode requires --data.true_labels argument")
  }
  
  # Load true labels for training
  cat("Loading true labels for training...\n")
  tryCatch({
    if (!file.exists(args$data.true_labels)) {
      stop("True labels file not found: ", args$data.true_labels)
    }
    
    # Load and convert labels
    true_labels <- load_and_convert_labels(args$data.true_labels)
    
    # Check if number of labels matches number of cells
    if (length(true_labels) != nrow(exp_matrix)) {
      stop("Number of labels (", length(true_labels), 
           ") does not match number of cells (", nrow(exp_matrix), ")")
    }
    
  }, error = function(e) {
    stop("Error loading true labels: ", e$message)
  })
  
  # Remove any unassigned cells from training data
  assigned_idx <- true_labels != "unassigned"
  cat("Found", sum(!assigned_idx), "unassigned cells\n")
  training_matrix <- exp_matrix[assigned_idx, , drop = FALSE]
  training_labels <- true_labels[assigned_idx]
  
  cat("Training on", length(training_labels), "labeled cells\n")
  cat("Cell types:", length(unique(training_labels)), "\n")
  cat("Cell type distribution:\n")
  print(table(training_labels))
  cat("\n")
  
  cat("Full matrix markers:", paste(colnames(exp_matrix), collapse=", "), "\n")
  cat("Training matrix markers:", paste(colnames(training_matrix), collapse=", "), "\n\n")
  
  # Train and annotate
  cat("Training and annotating...\n")
  res <- GateMeClass_annotate(
    exp_matrix = exp_matrix,
    marker_table = NULL,
    train_parameters = list(
      reference = training_matrix,
      labels = training_labels
    ),
    GMM_parameterization = args$GMM_parameterization,
    reject_option = args$reject_option,
    sampling = args$sampling,
    k = args$k,
    verbose = TRUE,
    seed = args$seed
  )
  
  marker_table <- res$marker_table
}

# Save results if annotation was performed
if (!is.null(res)) {
  cat("\n", strrep("=", 70), "\n")
  cat("Annotation Results\n")
  cat(strrep("=", 70), "\n")
  cat("Label distribution:\n")
  print(table(res$labels))
  cat("\n")
  
  # ============================================================================
  # MAIN OUTPUT: Plain text file with predicted labels (matches CyAnno format)
  # ============================================================================
  # This is the file that flow_metrics will use for evaluation
  prediction_file <- file.path(args$output_dir, 
                               paste0(args$name, "_predicted_labels.txt"))
  
  writeLines(res$labels, prediction_file)
  cat("✓ Predicted labels saved to:", prediction_file, "\n")
  cat("  Format: Plain text, one label per line (", length(res$labels), "cells)\n")
  
  # ============================================================================
  # AUXILIARY OUTPUTS: Additional files for reference/debugging
  # ============================================================================
  
  # Save marker table used (for interpretation)
  if (!is.null(marker_table)) {
    marker_table_path <- file.path(args$output_dir, 
                                    paste0(args$name, ".marker_table.gz"))
    fwrite(marker_table, marker_table_path, sep = "\t", compress = "gzip")
    cat("  Marker table saved to:", marker_table_path, "\n")
  }
  
  # Save cell signatures (for debugging/analysis)
  signatures_path <- file.path(args$output_dir, 
                               paste0(args$name, ".cell_signatures.gz"))
  fwrite(res$cell_signatures, signatures_path, sep = "\t", compress = "gzip")
  cat("  Cell signatures saved to:", signatures_path, "\n")
}

