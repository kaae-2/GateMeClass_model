#!/usr/bin/env Rscript

# GateMeClass wrapper for omnibenchmark
# - Uses random train/test split from preprocessing
# - Drops "unassigned" label (largest id)
# - Drops technical markers (Time, DNA, etc.)
# - Applies arcsinh transform
# - Trains via GateMeClass_annotate(train_parameters = ...)
# - Outputs numeric labels with ".0" and "" for NA as expected by metrics

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
})

# ---------------------------------------------------------------------
# Load / install GateMeClass
# ---------------------------------------------------------------------
if (!requireNamespace("GateMeClass", quietly = TRUE)) {
  message("GateMeClass not found, installing from GitHub...")
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org", quiet = TRUE)
  }
  remotes::install_github("simo1c/GateMeClass", quiet = TRUE)
}

suppressPackageStartupMessages(library(GateMeClass))

# ---------------------------------------------------------------------
# Patch: force `type` in set_marker_expression_GMM to be scalar logical
# ---------------------------------------------------------------------
fix_set_marker_expression_GMM <- function() {
  ns <- asNamespace("GateMeClass")
  orig_fun <- get("set_marker_expression_GMM", envir = ns)

  patched_fun <- function(X, GMM_parameterization, type, RSS) {
    # Defensive fix: make sure `type` is length-1 TRUE/FALSE
    if (length(type) == 0L || is.na(type[1])) {
      # Fall back to simple 2-component gating (+/-)
      type <- TRUE
    } else {
      type <- isTRUE(type[1])
    }
    orig_fun(X, GMM_parameterization, type, RSS)
  }

  assignInNamespace(
    x      = "set_marker_expression_GMM",
    value  = patched_fun,
    ns     = "GateMeClass",
    unlock = TRUE
  )
}

fix_set_marker_expression_GMM()

# ---------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------
parser <- ArgumentParser(description = "GateMeClass wrapper for omnibenchmark")

parser$add_argument('--train.data.matrix',
                    dest = "train_matrix",
                    type = "character",
                    help = 'gz-compressed CSV with training feature matrix (cells × markers)')
parser$add_argument('--labels_train',
                    dest = "train_labels",
                    type = "character",
                    help = 'gz-compressed file with training labels (one per cell)')
parser$add_argument('--test.data.matrix',
                    dest = "test_matrix",
                    type = "character",
                    help = 'gz-compressed CSV with test feature matrix (cells × markers)')
parser$add_argument('--labels_test',
                    dest = "test_labels",
                    type = "character",
                    help = 'gz-compressed file with test labels (optional, not used here)')
parser$add_argument('--seed',
                    type = "integer",
                    default = 42,
                    help = 'Random seed for GateMeClass')
parser$add_argument("--output_dir", "-o",
                    dest = "output_dir",
                    type = "character",
                    default = getwd(),
                    help = "Output directory for predictions")
parser$add_argument("--name", "-n",
                    dest = "name",
                    type = "character",
                    help = "Dataset name for output files")
parser$add_argument("--GMM_parameterization",
                    dest = "GMM_parameterization",
                    type = "character",
                    default = "V",
                    help = "GMM parameterization: 'V' (variable) or 'E' (equal)")
parser$add_argument("--sampling",
                    type = "double",
                    default = 0.3,
                    help = "Fraction of cells to use for annotation (0-1)")
parser$add_argument("--k",
                    type = "integer",
                    default = 20,
                    help = "k-NN parameter for label refinement")
parser$add_argument("--cofactor",
                    type = "double",
                    default = 5,
                    help = "Cofactor for arcsinh transformation")

args <- parser$parse_args()
set.seed(args$seed)

message("=== GateMeClass Configuration ===")
message("Train matrix: ", args$train_matrix)
message("Train labels: ", args$train_labels)
message("Test matrix:  ", args$test_matrix)
message("GMM: ", args$GMM_parameterization,
        ", sampling: ", args$sampling,
        ", k: ", args$k)
message("Cofactor: ", args$cofactor, " (arcsinh transformation)")
message("=================================\n")

# ---------------------------------------------------------------------
# Load training data
# ---------------------------------------------------------------------
message("Loading training data...")
train_dt <- fread(args$train_matrix)
train_labels_numeric <- fread(args$train_labels, header = FALSE)[[1]]

message("  Original train matrix: ", nrow(train_dt), " cells × ", ncol(train_dt), " columns")
message("  Original train labels: ", length(train_labels_numeric), " cells")

# Sanity check
if (nrow(train_dt) != length(train_labels_numeric)) {
  stop("Training matrix rows don't match label count!")
}

# Clean column names a bit (keep technical marker names intact)
names(train_dt) <- gsub("[^A-Za-z0-9_]", "_", names(train_dt))
names(train_dt) <- gsub("_+", "_", names(train_dt))

# Remove possible 'col' column
if ("col" %in% names(train_dt)) {
  train_dt <- train_dt[, ! "col", with = FALSE]
  message("  Removed 'col' column from training data")
}

# Drop NA / empty labels if any
train_labels_char <- as.character(train_labels_numeric)
valid_idx <- !is.na(train_labels_char) & train_labels_char != "" & train_labels_char != '""'
n_dropped_empty <- sum(!valid_idx)

if (n_dropped_empty > 0) {
  message("  Dropping ", n_dropped_empty, " cells with empty/NA labels from TRAINING")
  train_dt <- train_dt[valid_idx, ]
  train_labels_numeric <- train_labels_numeric[valid_idx]
  train_labels_char <- train_labels_char[valid_idx]
}

# Drop "unassigned" population: assume it's the largest label id
label_max <- max(train_labels_numeric, na.rm = TRUE)
message("  Largest label id in training: ", label_max)

idx_unassigned <- which(train_labels_numeric == label_max)
if (length(idx_unassigned) > 0) {
  message("  Dropping ", length(idx_unassigned),
          " cells with label == ", label_max,
          " from TRAINING (debug: mimic population != 'unassigned').")
  train_dt <- train_dt[-idx_unassigned, ]
  train_labels_numeric <- train_labels_numeric[-idx_unassigned]
  train_labels_char <- train_labels_char[-idx_unassigned]
}

message("  Clean train matrix: ", nrow(train_dt), " cells × ", ncol(train_dt), " markers")

# Map numeric labels -> generic cell type names
unique_ids <- sort(unique(train_labels_char))
id_to_name <- setNames(paste0("CellType_", unique_ids), unique_ids)
name_to_id <- setNames(unique_ids, paste0("CellType_", unique_ids))

train_labels_celltype <- id_to_name[train_labels_char]
names(train_labels_celltype) <- NULL

message("  Cell types used for training: ",
        paste(sort(unique(train_labels_celltype)), collapse = ", "))

# ---------------------------------------------------------------------
# Load test data
# ---------------------------------------------------------------------
message("\nLoading test data...")
test_dt <- fread(args$test_matrix)
message("  Original test matrix: ", nrow(test_dt), " cells × ", ncol(test_dt), " columns")

names(test_dt) <- gsub("[^A-Za-z0-9_]", "_", names(test_dt))
names(test_dt) <- gsub("_+", "_", names(test_dt))

if ("col" %in% names(test_dt)) {
  test_dt <- test_dt[, ! "col", with = FALSE]
  message("  Removed 'col' column from test data")
}

# ---------------------------------------------------------------------
# Drop technical markers (Time, DNA, etc.)
# ---------------------------------------------------------------------
technical_markers <- c("Time", "Cell_length", "DNA1", "DNA2", "Viability", "event_number")

train_drop <- intersect(technical_markers, names(train_dt))
test_drop  <- intersect(technical_markers, names(test_dt))

if (length(train_drop) > 0) {
  message("  Dropping technical markers from training: ",
          paste(train_drop, collapse = ", "))
  train_dt <- train_dt[, ! names(train_dt) %in% train_drop, with = FALSE]
}

if (length(test_drop) > 0) {
  message("  Dropping technical markers from test: ",
          paste(test_drop, collapse = ", "))
  test_dt <- test_dt[, ! names(test_dt) %in% test_drop, with = FALSE]
}

message("  Final train matrix: ", nrow(train_dt), " cells × ", ncol(train_dt), " markers")
message("  Final test matrix : ", nrow(test_dt),  " cells × ", ncol(test_dt),  " markers")

# ---------------------------------------------------------------------
# Transform and transpose
# ---------------------------------------------------------------------
message("\nPreparing data for GateMeClass...")

train_matrix_raw <- as.matrix(train_dt)
test_matrix_raw  <- as.matrix(test_dt)

message("Applying arcsinh transformation (cofactor = ", args$cofactor, ")...")
train_matrix_transformed <- asinh(train_matrix_raw / args$cofactor)
test_matrix_transformed  <- asinh(test_matrix_raw  / args$cofactor)

message("  Range(train_matrix): [",
        round(min(train_matrix_transformed, na.rm = TRUE), 2), ", ",
        round(max(train_matrix_transformed, na.rm = TRUE), 2), "]")

# GateMeClass expects markers in rows, cells in columns
train_matrix <- t(train_matrix_transformed)
test_matrix  <- t(test_matrix_transformed)

message("  GateMeClass train: ", nrow(train_matrix), " markers × ", ncol(train_matrix), " cells")
message("  GateMeClass test : ", nrow(test_matrix),  " markers × ", ncol(test_matrix),  " cells")

# ---------------------------------------------------------------------
# Run GateMeClass (Method 3: training + annotation in one step)
# ---------------------------------------------------------------------
message("\n=== Running GateMeClass_annotate (train_parameters) ===\n")

res <- GateMeClass_annotate(
  exp_matrix       = test_matrix,
  marker_table     = NULL,
  train_parameters = list(
    reference = train_matrix,
    labels    = train_labels_celltype
  ),
  GMM_parameterization = args$GMM_parameterization,
  reject_option       = FALSE,
  sampling            = args$sampling,
  k                   = args$k,
  verbose             = TRUE,
  narrow_marker_table = TRUE,
  seed                = args$seed
)

message("\n=== GateMeClass complete ===")

# ---------------------------------------------------------------------
# Save results
# ---------------------------------------------------------------------
predictions_celltype <- res$labels
predictions_numeric  <- as.numeric(name_to_id[predictions_celltype])

# Omnibenchmark expects numeric labels with ".0" and "" for NA
res_char <- paste0(predictions_numeric, ".0")
res_char[is.na(predictions_numeric)] <- NA

outfile <- file.path(args$output_dir, paste0(args$name, "_predicted_labels.txt"))

write.table(
  file      = outfile,
  x         = res_char,
  col.names = FALSE,
  row.names = FALSE,
  quote     = FALSE,
  na        = '""'
)

message("Predictions saved to: ", outfile)
message("Prediction summary (unique numeric labels): ",
        paste(sort(unique(predictions_numeric)), collapse = ", "))

message("\nGateMeClass analysis complete!")
