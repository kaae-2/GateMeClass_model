#!/usr/bin/env Rscript

# GateMeClass wrapper for OmniBenchmark / ob-pipeline-cytof
# - Uses train/test splits from the preprocessing stage
# - Does arcsinh(cofactor) transform + light outlier capping
# - Trains on train set, predicts labels for test set
# - Writes {name}_predicted_labels.txt in output_dir

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
})

## ----------------------------------------------------------------------
## Ensure GateMeClass is installed
## ----------------------------------------------------------------------
if (!requireNamespace("GateMeClass", quietly = TRUE)) {
  message("GateMeClass not found, installing from GitHub...")
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org", quiet = TRUE)
  }
  remotes::install_github("simo1c/GateMeClass", quiet = TRUE)
}
library(GateMeClass)

## ----------------------------------------------------------------------
## Argument parser – must match OmniBenchmark conventions
## ----------------------------------------------------------------------
parser <- ArgumentParser(description = "GateMeClass wrapper for OmniBenchmark")

parser$add_argument("--train.data.matrix", type = "character",
                    help = "gz-compressed CSV with training cells (rows) × markers (cols).")
parser$add_argument("--labels_train", type = "character",
                    help = "gz-compressed text file with one training label per row.")
parser$add_argument("--test.data.matrix", type = "character",
                    help = "gz-compressed CSV with test cells (rows) × markers (cols).")
parser$add_argument("--labels_test", type = "character",
                    help = "gz-compressed text file with one test label per row (unused).")

parser$add_argument("--seed", type = "integer", default = 42)
parser$add_argument("--output_dir", "-o", dest = "output_dir", type = "character",
                    default = getwd())
parser$add_argument("--name", "-n", dest = "name", type = "character")

# GateMeClass hyperparameters (fed from YAML `parameters`)
parser$add_argument("--GMM_parameterization", type = "character", default = "V")
parser$add_argument("--sampling", type = "double", default = 0.1)
parser$add_argument("--k", type = "integer", default = 20)
parser$add_argument("--cofactor", type = "double", default = 5,
                    help = "cofactor for arcsinh transformation")

args <- parser$parse_args()
set.seed(args$seed)

message("=== GateMeClass wrapper ===")
message(" train.data.matrix : ", args$`train.data.matrix`)
message(" labels_train      : ", args$labels_train)
message(" test.data.matrix  : ", args$`test.data.matrix`)
message(" labels_test       : ", args$labels_test, " (unused by wrapper)")
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
  # Standardize column names (avoid spaces, punctuation)
  colnames(dt) <- gsub("[^A-Za-z0-9_]", "_", colnames(dt))
  colnames(dt) <- gsub("_+", "_", colnames(dt))
  # Drop accidental 'col' column if present
  if ("col" %in% colnames(dt)) {
    dt <- dt[, setdiff(colnames(dt), "col"), drop = FALSE]
    message("  Dropped spurious 'col' column.")
  }
  as.matrix(dt)
}

read_labels <- function(path) {
  message("Reading labels: ", path)
  x <- fread(path, header = FALSE, data.table = FALSE)[[1]]
  x
}

## ----------------------------------------------------------------------
## Load data
## ----------------------------------------------------------------------
train_x_raw <- read_matrix(args$`train.data.matrix`)
train_y_raw <- read_labels(args$labels_train)

if (nrow(train_x_raw) != length(train_y_raw)) {
  stop("Training matrix rows (", nrow(train_x_raw),
       ") != training labels length (", length(train_y_raw), ")")
}

test_x_raw  <- read_matrix(args$`test.data.matrix`)
message("Train: ", nrow(train_x_raw), " cells × ", ncol(train_x_raw), " markers")
message("Test : ", nrow(test_x_raw),  " cells × ", ncol(test_x_raw),  " markers")

if (!identical(colnames(train_x_raw), colnames(test_x_raw))) {
  stop("Marker names differ between train and test matrices.")
}

## ----------------------------------------------------------------------
## Clean labels – drop unlabeled cells from training
## ----------------------------------------------------------------------
train_y_char <- as.character(train_y_raw)
valid_idx <- !is.na(train_y_char) & train_y_char != "" & train_y_char != '""'
n_dropped <- sum(!valid_idx)

if (n_dropped > 0) {
  message("Dropping ", n_dropped, " unlabeled cells from training set.")
  train_x_raw <- train_x_raw[valid_idx, , drop = FALSE]
  train_y_char <- train_y_char[valid_idx]
} else {
  message("No unlabeled cells in training set.")
}

message("Clean train: ", nrow(train_x_raw), " cells × ", ncol(train_x_raw), " markers")

## Map numeric labels -> pseudo cell type names (GateMeClass expects character labels)
unique_ids <- sort(unique(train_y_char))
id_to_name <- setNames(paste0("CellType_", unique_ids), unique_ids)
name_to_id <- setNames(unique_ids, paste0("CellType_", unique_ids))

train_y_celltype <- unname(id_to_name[train_y_char])

message("Training cell types: ", paste(sort(unique(train_y_celltype)), collapse = ", "))

## ----------------------------------------------------------------------
## Outlier handling + arcsinh transform
## ----------------------------------------------------------------------
message("\nTransforming data (arcsinh, cofactor = ", args$cofactor, ") with light outlier capping...")

train_vals <- as.vector(train_x_raw)
upper_bound <- stats::quantile(train_vals, 0.999, na.rm = TRUE)

message("  99.9% quantile of train intensities: ", round(upper_bound, 3))
message("  Max before capping: ", round(max(train_vals, na.rm = TRUE), 3))

train_x_raw[train_x_raw > upper_bound] <- upper_bound
test_x_raw[test_x_raw > upper_bound]   <- upper_bound

train_x_raw[train_x_raw < 0] <- 0
test_x_raw[test_x_raw < 0]   <- 0

train_x <- asinh(train_x_raw / args$cofactor)
test_x  <- asinh(test_x_raw  / args$cofactor)

message("  Range(train_x) after transform: [",
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
## Run GateMeClass (Method 3: training + annotation in one step)
## ----------------------------------------------------------------------
set.seed(args$seed)
message("\n=== Running GateMeClass_annotate (train_parameters) ===")

res <- GateMeClass_annotate(
  exp_matrix      = test_mat_gc,
  marker_table    = NULL,
  train_parameters = list(
    reference = train_mat_gc,
    labels    = train_y_celltype
  ),
  GMM_parameterization = args$GMM_parameterization,
  sampling            = args$sampling,
  k                   = args$k,
  verbose             = TRUE,
  seed                = args$seed
)

pred_celltypes <- res$labels
if (length(pred_celltypes) != nrow(test_x_raw)) {
  stop("GateMeClass returned ", length(pred_celltypes),
       " labels for ", nrow(test_x_raw), " test cells.")
}

message("\nGateMeClass done. Unique predicted types: ",
        length(unique(pred_celltypes)))

## ----------------------------------------------------------------------
## Map predictions back to numeric IDs and save
## ----------------------------------------------------------------------
pred_ids <- as.numeric(name_to_id[pred_celltypes])

# Metrics code generally expects something like "1.0", "2.0", ...
pred_str <- paste0(pred_ids, ".0")
pred_str[is.na(pred_ids)] <- NA_character_

outfile <- file.path(args$output_dir, paste0(args$name, "_predicted_labels.txt"))
message("Writing predictions to: ", outfile)

write.table(
  pred_str,
  file       = outfile,
  col.names  = FALSE,
  row.names  = FALSE,
  quote      = FALSE,
  na         = '""'   # consistent with how NAs are stored elsewhere in the pipeline
)

message("Prediction file written.")
message("Prediction distribution:")
print(table(pred_celltypes))

