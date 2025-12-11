#!/usr/bin/env Rscript

# GateMeClass wrapper for omnibenchmark
# Final version: Technical marker filtering + transformation
#  - Drops largest numeric label from TRAINING (treat as "unassigned")
#  - Drops technical markers: Time, Cell_length, DNA1, DNA2, Viability, event_number
#  - Drops markers with zero variance in training
#  - Uses arcsinh transformation on BIOLOGICAL markers only

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
})

# ---------------------------------------------------------------------------
# Ensure GateMeClass is installed
# ---------------------------------------------------------------------------

if (!require("GateMeClass", quietly = TRUE)) {
  message("GateMeClass not found, installing from GitHub...")
  if (!require("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org", quiet = TRUE)
  }
  remotes::install_github("simo1c/GateMeClass", quiet = TRUE)
}

suppressPackageStartupMessages(library(GateMeClass))

# ============================================================================
# Argument Parser
# ============================================================================

parser <- ArgumentParser(description = "GateMeClass - Technical marker filtering + transformation")

parser$add_argument('--train.data.matrix', type = "character", dest = "train_matrix")
parser$add_argument('--labels_train',      type = "character", dest = "train_labels")
parser$add_argument('--test.data.matrix',  type = "character", dest = "test_matrix")
parser$add_argument('--labels_test',       type = "character", dest = "test_labels")

parser$add_argument('--seed', type = "integer", default = 42)
parser$add_argument("--output_dir", "-o", dest = "output_dir", type = "character", default = getwd())
parser$add_argument("--name", "-n",        dest = "name",       type = "character")

parser$add_argument("--GMM_parameterization", type = "character", default = "V")
parser$add_argument("--sampling", type = "double", default = 0.3)
parser$add_argument("--k",        type = "integer", default = 20)
parser$add_argument("--cofactor", type = "double", default = 5)

args <- parser$parse_args()
set.seed(args$seed)

message("=== GateMeClass Configuration ===")
message("Train matrix: ", args$train_matrix)
message("Train labels: ", args$train_labels)
message("Test matrix:  ", args$test_matrix)
message("GMM: ", args$GMM_parameterization,
        ", sampling: ", args$sampling, ", k: ", args$k)
message("Cofactor: ", args$cofactor, " (arcsinh transformation)")
message("=================================\n")

# ============================================================================
# Load Training Data
# ============================================================================

message("Loading training data...")
train_dt <- fread(args$train_matrix)
train_labels_numeric <- fread(args$train_labels, header = FALSE)[[1]]

# Standardize column names
names(train_dt) <- gsub("[^A-Za-z0-9_]", "_", names(train_dt))
names(train_dt) <- gsub("_+", "_", names(train_dt))

message("  Original train matrix: ", nrow(train_dt), " cells × ", ncol(train_dt), " columns")

# Remove 'col' column if present
if ("col" %in% names(train_dt)) {
  train_dt <- train_dt[, !names(train_dt) %in% "col", with = FALSE]
}

# Convert labels to character
train_labels_char <- as.character(train_labels_numeric)

# Remove truly unlabeled cells
valid_train_idx <- !is.na(train_labels_char) &
                   train_labels_char != "" &
                   train_labels_char != '""'
n_unlabeled <- sum(!valid_train_idx)

if (n_unlabeled > 0) {
  message("  Removing ", n_unlabeled, " unlabeled cells from training set")
  train_dt             <- train_dt[valid_train_idx, ]
  train_labels_numeric <- train_labels_numeric[valid_train_idx]
  train_labels_char    <- train_labels_char[valid_train_idx]
}

# Drop largest numeric label from TRAINING (treat as "unassigned")
train_ids_num <- suppressWarnings(as.numeric(train_labels_char))
unassigned_id <- max(train_ids_num, na.rm = TRUE)

keep_idx <- train_ids_num != unassigned_id
n_dropped_unassigned <- sum(!keep_idx)

message("  Largest label id: ", unassigned_id)
message("  Dropping ", n_dropped_unassigned, " cells with label == ", unassigned_id)

train_dt             <- train_dt[keep_idx, ]
train_labels_numeric <- train_labels_numeric[keep_idx]
train_labels_char    <- train_labels_char[keep_idx]

if (nrow(train_dt) == 0L) {
  stop("After dropping unassigned cells, training matrix has 0 rows.")
}

message("  Clean train matrix: ", nrow(train_dt), " cells × ", ncol(train_dt), " markers")

# Map numeric IDs -> CellType_<id>
unique_ids <- sort(unique(train_labels_char))
id_to_name <- setNames(paste0("CellType_", unique_ids), unique_ids)
name_to_id <- setNames(unique_ids, paste0("CellType_", unique_ids))

train_labels_celltype <- id_to_name[train_labels_char]
names(train_labels_celltype) <- NULL

message("  Cell types: ", paste(sort(unique(train_labels_celltype)), collapse = ", "))

# ============================================================================
# Load Test Data
# ============================================================================

message("\nLoading test data...")
test_dt <- fread(args$test_matrix)

names(test_dt) <- gsub("[^A-Za-z0-9_]", "_", names(test_dt))
names(test_dt) <- gsub("_+", "_", names(test_dt))

if ("col" %in% names(test_dt)) {
  test_dt <- test_dt[, !names(test_dt) %in% "col", with = FALSE]
}

message("  Test matrix: ", nrow(test_dt), " cells × ", ncol(test_dt), " columns")

# ============================================================================
# CRITICAL: Drop technical markers BEFORE transformation
# ============================================================================

technical_markers <- c("Time", "Cell_length", "DNA1", "DNA2", "Viability", "event_number")

drop_train <- intersect(names(train_dt), technical_markers)
drop_test  <- intersect(names(test_dt),  technical_markers)

if (length(drop_train) > 0) {
  message("\n*** DROPPING TECHNICAL MARKERS from training: ",
          paste(drop_train, collapse = ", "))
  train_dt <- train_dt[, !names(train_dt) %in% drop_train, with = FALSE]
}

if (length(drop_test) > 0) {
  message("*** DROPPING TECHNICAL MARKERS from test: ",
          paste(drop_test, collapse = ", "))
  test_dt <- test_dt[, !names(test_dt) %in% drop_test, with = FALSE]
}

# Drop flat markers
marker_var <- sapply(train_dt, function(x) var(as.numeric(x), na.rm = TRUE))
flat_markers <- names(marker_var)[is.na(marker_var) | marker_var == 0]

if (length(flat_markers) > 0) {
  message("*** DROPPING flat markers: ", paste(flat_markers, collapse = ", "))
  train_dt <- train_dt[, !names(train_dt) %in% flat_markers, with = FALSE]
  test_dt  <- test_dt[, !names(test_dt)  %in% flat_markers, with = FALSE]
}

# Sanity check
if (!identical(names(train_dt), names(test_dt))) {
  stop("Train and test have different markers after cleaning.")
}

message("\nFinal matrices after filtering:")
message("  Train: ", nrow(train_dt), " cells × ", ncol(train_dt), " BIOLOGICAL markers")
message("  Test:  ", nrow(test_dt),  " cells × ", ncol(test_dt),  " BIOLOGICAL markers")

# ============================================================================
# Transform biological markers only
# ============================================================================

message("\nApplying arcsinh transformation (cofactor = ", args$cofactor, ")")
message("  This transformation is applied to BIOLOGICAL markers only")

train_matrix_raw <- as.matrix(train_dt)
test_matrix_raw  <- as.matrix(test_dt)

train_matrix_transformed <- asinh(train_matrix_raw / args$cofactor)
test_matrix_transformed  <- asinh(test_matrix_raw  / args$cofactor)

message("  Transformed range: [",
        round(min(train_matrix_transformed, na.rm = TRUE), 2), ", ",
        round(max(train_matrix_transformed, na.rm = TRUE), 2), "]")

# Transpose for GateMeClass (markers × cells)
train_matrix <- t(train_matrix_transformed)
test_matrix  <- t(test_matrix_transformed)

message("  GateMeClass format:")
message("    Train: ", nrow(train_matrix), " markers × ", ncol(train_matrix), " cells")
message("    Test:  ", nrow(test_matrix),  " markers × ", ncol(test_matrix),  " cells")

# ============================================================================
# Run GateMeClass
# ============================================================================

message("\n=== Running GateMeClass ===\n")

res <- GateMeClass_annotate(
  exp_matrix       = test_matrix,
  marker_table     = NULL,
  train_parameters = list(
    reference = train_matrix,
    labels    = train_labels_celltype
  ),
  GMM_parameterization = args$GMM_parameterization,
  sampling             = args$sampling,
  k                    = args$k,
  verbose              = TRUE,
  seed                 = args$seed
)

message("\n=== GateMeClass Complete ===")

# ============================================================================
# Save Results
# ============================================================================

predictions_celltype <- res$labels
predictions_numeric  <- as.numeric(name_to_id[predictions_celltype])

# Map to "k.0" format
res_char <- paste0(predictions_numeric, ".0")
res_char[is.na(predictions_numeric)] <- NA

if (length(res_char) != nrow(test_dt)) {
  stop("Predictions (", length(res_char), ") != test cells (", nrow(test_dt), ")")
}

outfile <- file.path(args$output_dir, paste0(args$name, "_predicted_labels.txt"))
write.table(
  file      = outfile,
  x         = res_char,
  col.names = FALSE,
  row.names = FALSE,
  quote     = FALSE,
  na        = '""'
)

message("\nPredictions saved to: ", outfile)
message("Unique labels: ", length(unique(predictions_numeric)))
message("\nSuccess!")
