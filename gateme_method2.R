#!/usr/bin/env Rscript

# GateMeClass wrapper for omnibenchmark - Method 2
# - Separate training and annotation steps
# - Drops technical markers before processing
# - Applies arcsinh transformation
# - Logs marker table between steps for transparency

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
# Argument parser
# ---------------------------------------------------------------------
parser <- ArgumentParser(description = "GateMeClass wrapper - Method 2 (separate train/annotate)")

parser$add_argument('--train.data.matrix', dest = "train_matrix", type = "character")
parser$add_argument('--labels_train', dest = "train_labels", type = "character")
parser$add_argument('--test.data.matrix', dest = "test_matrix", type = "character")
parser$add_argument('--labels_test', dest = "test_labels", type = "character")
parser$add_argument('--seed', type = "integer", default = 42)
parser$add_argument("--output_dir", "-o", dest = "output_dir", type = "character", default = getwd())
parser$add_argument("--name", "-n", dest = "name", type = "character")
parser$add_argument("--GMM_parameterization", type = "character", default = "V")
parser$add_argument("--sampling", type = "double", default = 0.3)
parser$add_argument("--k", type = "integer", default = 20)
parser$add_argument("--cofactor", type = "double", default = 5)

args <- parser$parse_args()
set.seed(args$seed)

message("=== GateMeClass Configuration (Method 2) ===")
message("Train matrix: ", args$train_matrix)
message("Train labels: ", args$train_labels)
message("Test matrix:  ", args$test_matrix)
message("GMM: ", args$GMM_parameterization, ", sampling: ", args$sampling, ", k: ", args$k)
message("Cofactor: ", args$cofactor)
message("============================================\n")

# ---------------------------------------------------------------------
# Load training data
# ---------------------------------------------------------------------
message("Loading training data...")
train_dt <- fread(args$train_matrix)
train_labels_numeric <- fread(args$train_labels, header = FALSE)[[1]]

message("  Original: ", nrow(train_dt), " cells × ", ncol(train_dt), " columns")

if (nrow(train_dt) != length(train_labels_numeric)) {
  stop("Training matrix rows don't match label count!")
}

# Standardize column names
names(train_dt) <- gsub("[^A-Za-z0-9_]", "_", names(train_dt))
names(train_dt) <- gsub("_+", "_", names(train_dt))

if ("col" %in% names(train_dt)) {
  train_dt <- train_dt[, !"col", with = FALSE]
}

# Drop empty labels
train_labels_char <- as.character(train_labels_numeric)
valid_idx <- !is.na(train_labels_char) & train_labels_char != "" & train_labels_char != '""'

if (sum(!valid_idx) > 0) {
  message("  Dropping ", sum(!valid_idx), " cells with empty/NA labels")
  train_dt <- train_dt[valid_idx, ]
  train_labels_char <- train_labels_char[valid_idx]
}

# ---------------------------------------------------------------------
# Load test data
# ---------------------------------------------------------------------
message("\nLoading test data...")
test_dt <- fread(args$test_matrix)
message("  Original: ", nrow(test_dt), " cells × ", ncol(test_dt), " columns")

names(test_dt) <- gsub("[^A-Za-z0-9_]", "_", names(test_dt))
names(test_dt) <- gsub("_+", "_", names(test_dt))

if ("col" %in% names(test_dt)) {
  test_dt <- test_dt[, !"col", with = FALSE]
}

# ---------------------------------------------------------------------
# Drop technical markers BEFORE transformation
# ---------------------------------------------------------------------
technical_markers <- c("Time", "Cell_length", "DNA1", "DNA2", "Viability", "event_number")

drop_train <- intersect(names(train_dt), technical_markers)
drop_test <- intersect(names(test_dt), technical_markers)

if (length(drop_train) > 0) {
  message("\n*** Dropping technical markers from TRAIN: ", paste(drop_train, collapse = ", "))
  train_dt <- train_dt[, !names(train_dt) %in% drop_train, with = FALSE]
}

if (length(drop_test) > 0) {
  message("*** Dropping technical markers from TEST: ", paste(drop_test, collapse = ", "))
  test_dt <- test_dt[, !names(test_dt) %in% drop_test, with = FALSE]
}

# Sanity check
if (!all(names(train_dt) == names(test_dt))) {
  stop("Marker names differ between train and test after preprocessing!")
}

message("\nFinal dimensions:")
message("  Train: ", nrow(train_dt), " cells × ", ncol(train_dt), " BIOLOGICAL markers")
message("  Test:  ", nrow(test_dt), " cells × ", ncol(test_dt), " BIOLOGICAL markers")

# ---------------------------------------------------------------------
# Map labels to cell type names
# ---------------------------------------------------------------------
unique_ids <- sort(unique(train_labels_char))
id_to_name <- setNames(paste0("CellType_", unique_ids), unique_ids)
name_to_id <- setNames(unique_ids, paste0("CellType_", unique_ids))

train_labels_celltype <- id_to_name[train_labels_char]
names(train_labels_celltype) <- NULL

message("\nCell types (", length(unique(train_labels_celltype)), "): ",
        paste(sort(unique(train_labels_celltype)), collapse = ", "))

# ---------------------------------------------------------------------
# Transform with arcsinh
# ---------------------------------------------------------------------
message("\nApplying arcsinh transformation (cofactor = ", args$cofactor, ")...")

train_matrix_raw <- as.matrix(train_dt)
test_matrix_raw <- as.matrix(test_dt)

train_matrix_transformed <- asinh(train_matrix_raw / args$cofactor)
test_matrix_transformed <- asinh(test_matrix_raw / args$cofactor)

message("  Transformed range: [",
        round(min(train_matrix_transformed), 2), ", ",
        round(max(train_matrix_transformed), 2), "]")

# Transpose for GateMeClass (markers × cells)
train_matrix <- t(train_matrix_transformed)
test_matrix <- t(test_matrix_transformed)

message("  GateMeClass format:")
message("    Train: ", nrow(train_matrix), " markers × ", ncol(train_matrix), " cells")
message("    Test:  ", nrow(test_matrix), " markers × ", ncol(test_matrix), " cells")

# =====================================================================
# STEP 1: TRAINING - Extract marker table
# =====================================================================
message("\n=== STEP 1: GateMeClass_train ===")
message("Extracting marker table from training data...\n")

marker_table <- GateMeClass_train(
  reference = train_matrix,
  labels = train_labels_celltype,
  GMM_parameterization = args$GMM_parameterization,
  verbose = TRUE,
  seed = args$seed
)

message("\n=== Training Complete ===")

# CRITICAL: Verify marker table format
message("\nMarker table structure:")
message("  Dimensions: ", nrow(marker_table), " rows × ", ncol(marker_table), " columns")
message("  Column names: ", paste(names(marker_table), collapse = ", "))

# GateMeClass_train returns NARROW format: "Cell" and "Gate" columns
# This requires narrow_marker_table = TRUE in annotate
is_narrow <- all(names(marker_table) %in% c("Cell", "Gate")) && ncol(marker_table) == 2

if (!is_narrow) {
  warning("UNEXPECTED: marker_table is NOT narrow format (Cell, Gate)!")
  warning("  This may cause annotation to fail.")
  warning("  Columns found: ", paste(names(marker_table), collapse = ", "))
}

message("\nMarker table preview:")
print(head(marker_table, 10))
message("\nMarker table format: ", if(is_narrow) "NARROW (Cell, Gate)" else "WIDE (multiple marker columns)")

# =====================================================================
# STEP 2: ANNOTATION - Apply marker table to test data
# =====================================================================
message("\n=== STEP 2: GateMeClass_annotate ===")
message("Annotating test data using extracted marker table...\n")

# GateMeClass_train returns narrow format (Cell, Gate), so narrow_marker_table = TRUE
# If training returned wide format, this would need to be FALSE
res <- GateMeClass_annotate(
  exp_matrix = test_matrix,
  marker_table = marker_table,
  GMM_parameterization = args$GMM_parameterization,
  reject_option = FALSE,
  sampling = args$sampling,
  k = args$k,
  verbose = TRUE,
  narrow_marker_table = TRUE,  # GateMeClass_train returns narrow format
  seed = args$seed
)

message("\n=== Annotation Complete ===")

# ---------------------------------------------------------------------
# Save results
# ---------------------------------------------------------------------
predictions_celltype <- res$labels
predictions_numeric <- as.numeric(name_to_id[predictions_celltype])

# Format: "k.0" for numeric, "" for NA
res_char <- paste0(predictions_numeric, ".0")
res_char[is.na(predictions_numeric)] <- NA

if (length(res_char) != nrow(test_dt)) {
  stop("Predictions (", length(res_char), ") != test cells (", nrow(test_dt), ")")
}

outfile <- file.path(args$output_dir, paste0(args$name, "_predicted_labels.txt"))

write.table(
  file = outfile,
  x = res_char,
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE,
  na = '""'
)

message("\nPredictions saved to: ", outfile)
message("Unique labels: ", length(unique(predictions_numeric)))
message("Label distribution:")
print(table(predictions_celltype))

message("\n=== SUCCESS ===")
