#!/usr/bin/env Rscript

## GateMeClass wrapper for OmniBenchmark / ob-pipeline-cytof
## - Uses train/test splits from preprocessing
## - asinh transform with cofactor (default 5), no capping
## - Trains on train set, predicts labels for test set
## - Writes {name}_predicted_labels.txt in output_dir

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
})

## ----------------------------------------------------------------------
## Ensure GateMeClass is available
## ----------------------------------------------------------------------
if (!requireNamespace("GateMeClass", quietly = TRUE)) {
  message("GateMeClass not found, installing from GitHub (simo1c/GateMeClass)...")
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org", quiet = TRUE)
  }
  remotes::install_github("simo1c/GateMeClass", quiet = TRUE)
}
suppressPackageStartupMessages(library(GateMeClass))

## ----------------------------------------------------------------------
## Argument parser – match OmniBenchmark conventions
## ----------------------------------------------------------------------
parser <- ArgumentParser(description = "GateMeClass wrapper for OmniBenchmark")

parser$add_argument("--train.data.matrix",
                    dest  = "train_matrix",
                    type  = "character",
                    help  = "gz-compressed CSV with training cells (rows) × markers (cols).",
                    required = TRUE)

parser$add_argument("--labels_train",
                    dest  = "labels_train",
                    type  = "character",
                    help  = "gz-compressed labels file for training data (one label per row).",
                    required = TRUE)

parser$add_argument("--test.data.matrix",
                    dest  = "test_matrix",
                    type  = "character",
                    help  = "gz-compressed CSV with test cells (rows) × markers (cols).",
                    required = TRUE)

parser$add_argument("--labels_test",
                    dest  = "labels_test",
                    type  = "character",
                    help  = "gz-compressed labels file for test data (unused by wrapper).",
                    required = TRUE)

parser$add_argument("--seed",
                    dest  = "seed",
                    type  = "integer",
                    default = 42,
                    help = "Random seed.")

parser$add_argument("--output_dir", "-o",
                    dest  = "output_dir",
                    type  = "character",
                    default = getwd(),
                    help = "Directory to write prediction file.")

parser$add_argument("--name", "-n",
                    dest  = "name",
                    type  = "character",
                    help = "Dataset name used in output filename.",
                    required = TRUE)

## GateMeClass hyperparameters
parser$add_argument("--GMM_parameterization",
                    dest  = "GMM_parameterization",
                    type  = "character",
                    default = "V",
                    help = "GMM variance parameterization ('V' or 'E').")

parser$add_argument("--sampling",
                    dest  = "sampling",
                    type  = "double",
                    default = 0.1,
                    help = "Fraction of cells used in GateMeClass sampling step.")

parser$add_argument("--k",
                    dest  = "k",
                    type  = "integer",
                    default = 20,
                    help = "k for k-NN / MNN label refinement.")

parser$add_argument("--cofactor",
                    dest  = "cofactor",
                    type  = "double",
                    default = 5,
                    help = "Cofactor for arcsinh transformation (tutorial uses 5).")

args <- parser$parse_args()
set.seed(args$seed)

message("=== GateMeClass wrapper ===")
message(" train.data.matrix : ", args$train_matrix)
message(" labels_train      : ", args$labels_train)
message(" test.data.matrix  : ", args$test_matrix)
message(" labels_test       : ", args$labels_test, " (unused)")
message(" GMM_parameterization = ", args$GMM_parameterization)
message(" sampling             = ", args$sampling)
message(" k                    = ", args$k)
message(" cofactor             = ", args$cofactor)
message(" seed                 = ", args$seed)
message(" output_dir           = ", args$output_dir)
message(" name                 = ", args$name)
message("===========================\n")

dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

## ----------------------------------------------------------------------
## Helper I/O
## ----------------------------------------------------------------------
read_matrix <- function(path) {
  message("Reading matrix: ", path)
  dt <- fread(path, header = TRUE, data.table = FALSE)
  # Sanitize column names (GateMeClass expects "marker-like" names)
  cn <- colnames(dt)
  cn <- gsub("[^A-Za-z0-9_]", "_", cn)
  cn <- gsub("_+", "_", cn)
  colnames(dt) <- cn

  # Occasionally pipelines generate a stray 'col' column – drop if present
  if ("col" %in% colnames(dt)) {
    dt <- dt[, setdiff(colnames(dt), "col"), drop = FALSE]
    message("  Dropped spurious 'col' column.")
  }

  as.matrix(dt)
}

read_labels <- function(path) {
  message("Reading labels: ", path)
  x <- fread(path, header = FALSE, data.table = FALSE)[[1]]
  as.character(x)
}

## ----------------------------------------------------------------------
## Load matrices
## ----------------------------------------------------------------------
train_x_raw <- read_matrix(args$train_matrix)
test_x_raw  <- read_matrix(args$test_matrix)

message("Train matrix: ", nrow(train_x_raw), " cells × ", ncol(train_x_raw), " markers")
message("Test  matrix: ", nrow(test_x_raw),  " cells × ", ncol(test_x_raw),  " markers")

if (!identical(colnames(train_x_raw), colnames(test_x_raw))) {
  stop("Marker names differ between train and test matrices.")
}

## ----------------------------------------------------------------------
## Labels: numeric + \"\" → character labels for GateMeClass
## ----------------------------------------------------------------------
train_y_raw <- read_labels(args$labels_train)

if (nrow(train_x_raw) != length(train_y_raw)) {
  stop("Training matrix rows (", nrow(train_x_raw),
       ") != training labels length (", length(train_y_raw), ")")
}

# Treat NA, "" and '""' as unlabeled
is_empty <- is.na(train_y_raw) |
            train_y_raw == "" |
            train_y_raw == '""'

n_empty <- sum(is_empty)
message("Training labels: ", length(train_y_raw),
        " total, ", n_empty, " empty/NA")

if (n_empty > 0) {
  message("Dropping ", n_empty, " unlabeled cells from training set.")
  train_x_raw <- train_x_raw[!is_empty, , drop = FALSE]
  train_y_raw <- train_y_raw[!is_empty]
}

# Normalize numeric labels:
#   - strip trailing ".0" if present (e.g. "3.0" -> "3")
#   - require they are integers-as-strings
train_ids_chr <- sub("\\.0$", "", train_y_raw)

if (!all(grepl("^[0-9]+$", train_ids_chr))) {
  bad <- unique(train_ids_chr[!grepl("^[0-9]+$", train_ids_chr)])
  stop("Non-numeric training labels found after cleaning: ",
       paste(bad, collapse = ", "))
}

training_set_lab <- train_ids_chr

message("Unique training label IDs after cleaning: ",
        paste(sort(unique(training_set_lab)), collapse = ", "))

message("Cleaned train matrix: ", nrow(train_x_raw), " cells × ",
        ncol(train_x_raw), " markers")

## ----------------------------------------------------------------------
## Transformation: arcsinh with cofactor (no capping)
## ----------------------------------------------------------------------
message("\nApplying arcsinh transform with cofactor = ", args$cofactor, " (no capping)")

train_x <- asinh(train_x_raw / args$cofactor)
test_x  <- asinh(test_x_raw  / args$cofactor)

message("  Range(train_x): [",
        round(min(train_x, na.rm = TRUE), 3), ", ",
        round(max(train_x, na.rm = TRUE), 3), "]")

## GateMeClass expects markers as rows, cells as columns
train_mat_gc <- t(train_x)
test_mat_gc  <- t(test_x)

message("GateMeClass train matrix: ", nrow(train_mat_gc), " markers × ",
        ncol(train_mat_gc), " cells")
message("GateMeClass test matrix : ", nrow(test_mat_gc),  " markers × ",
        ncol(test_mat_gc),  " cells")

## ----------------------------------------------------------------------
## Run GateMeClass: training + annotation in one step (Method 3)
## ----------------------------------------------------------------------
set.seed(args$seed)
message("\n=== Running GateMeClass_annotate (train_parameters) ===")

res <- GateMeClass_annotate(
  exp_matrix       = test_mat_gc,
  marker_table     = NULL,
  train_parameters = list(
    reference = train_mat_gc,
    labels    = training_set_lab
  ),
  GMM_parameterization = args$GMM_parameterization,
  sampling             = args$sampling,
  k                    = args$k,
  verbose              = TRUE,
  seed                 = args$seed
)

pred_gate <- res$labels

if (length(pred_gate) != nrow(test_x_raw)) {
  stop("GateMeClass returned ", length(pred_gate),
       " labels for ", nrow(test_x_raw), " test cells.")
}

message("\nGateMeClass done. Unique predicted labels: ",
        paste(sort(unique(pred_gate)), collapse = ", "))

## ----------------------------------------------------------------------
## Map predictions back to numeric IDs with .0 and write file
## ----------------------------------------------------------------------
# Convert predictions (e.g. "1","2") -> integer; non-convertible -> NA
pred_ids <- suppressWarnings(as.integer(pred_gate))

# Metrics expect something like "1.0", "2.0", "" for missing
pred_str <- ifelse(is.na(pred_ids), '""', paste0(pred_ids, ".0"))

outfile <- file.path(args$output_dir,
                     paste0(args$name, "_predicted_labels.txt"))
message("Writing predictions to: ", outfile)

write.table(
  pred_str,
  file      = outfile,
  col.names = FALSE,
  row.names = FALSE,
  quote     = FALSE
)

message("Prediction file written.")
message("Prediction ID distribution (numeric):")
print(table(pred_ids, useNA = "ifany"))


