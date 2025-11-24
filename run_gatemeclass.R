#!/usr/bin/env Rscript

# GateMeClass module for omnibenchmark
# Automated cell type annotation for flow/CyTOF data

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
                   default = "annotate",
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

# Load data matrix
cat("Loading data matrix...\n")
exp_matrix <- as.matrix(fread(args$data.matrix))
cat("Matrix dimensions:", nrow(exp_matrix), "cells x", ncol(exp_matrix), "markers\n\n")

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
  true_labels <- readLines(gzfile(args$data.true_labels))
  
  # Remove any unassigned cells
  assigned_idx <- true_labels != "unassigned"
  exp_matrix <- exp_matrix[assigned_idx, ]
  true_labels <- true_labels[assigned_idx]
  
  cat("Training on", length(true_labels), "labeled cells\n")
  cat("Cell types:", length(unique(true_labels)), "\n\n")
  
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
  
  # Save marker table
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
  marker_table <- as.data.frame(fread(args$marker_table))
  marker_table[is.na(marker_table)] <- "*"  # Replace NA with wildcards
  
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
  true_labels <- readLines(gzfile(args$data.true_labels))
  
  # Remove any unassigned cells from training data
  assigned_idx <- true_labels != "unassigned"
  training_matrix <- exp_matrix[assigned_idx, ]
  training_labels <- true_labels[assigned_idx]
  
  cat("Training on", length(training_labels), "labeled cells\n")
  cat("Cell types:", length(unique(training_labels)), "\n\n")
  
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
  
  # Save predicted labels (main output for omnibenchmark)
  labels_path <- file.path(args$output_dir, 
                           paste0(args$name, ".labels.gz"))
  writeLines(res$labels, gzfile(labels_path))
  cat("Predicted labels saved to:", labels_path, "\n")
  
  # Save marker table used
  if (!is.null(marker_table)) {
    marker_table_path <- file.path(args$output_dir, 
                                    paste0(args$name, ".marker_table.gz"))
    fwrite(marker_table, marker_table_path, sep = "\t", compress = "gzip")
    cat("Marker table saved to:", marker_table_path, "\n")
  }
  
  # Save cell signatures
  signatures_path <- file.path(args$output_dir, 
                               paste0(args$name, ".signatures.gz"))
  fwrite(res$cell_signatures, signatures_path, sep = "\t", compress = "gzip")
  cat("Cell signatures saved to:", signatures_path, "\n")
}

cat("\n", strrep("=", 70), "\n")
cat("GateMeClass completed successfully!\n")
cat(strrep("=", 70), "\n")
