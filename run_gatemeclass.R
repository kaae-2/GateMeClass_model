#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
})

# ---------------------------------------------------------------------
# Install / load GateMeClass
# ---------------------------------------------------------------------
if (!requireNamespace("GateMeClass", quietly = TRUE)) {
  message("GateMeClass not found, installing from GitHub...")
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org")
  }
  tryCatch(
    remotes::install_github("simo1c/GateMeClass", upgrade = "never"),
    error = function(e) stop("GateMeClass install failed: ", conditionMessage(e))
  )
}
suppressPackageStartupMessages(library(GateMeClass))

# Safety: if unqualified train() is used internally, point it to caret::train
if (requireNamespace("caret", quietly = TRUE)) {
  assign("train", get("train", envir = asNamespace("caret")), envir = .GlobalEnv)
}

# ---------------------------------------------------------------------
# Args
# ---------------------------------------------------------------------
parser <- ArgumentParser(description="GateMeClass wrapper (Method 3, robust marker naming)")

parser$add_argument("--train.data.matrix", type="character", required=TRUE)
parser$add_argument("--labels_train",      type="character", required=TRUE)
parser$add_argument("--test.data.matrix",  type="character", required=TRUE)
parser$add_argument("--labels_test",       type="character", required=FALSE, default=NULL)

parser$add_argument("--seed", type="integer", default=42)
parser$add_argument("--output_dir", "-o", type="character", default=getwd())
parser$add_argument("--name", "-n", type="character", required=TRUE)

parser$add_argument("--GMM_parameterization", type="character", default="V")
parser$add_argument("--sampling", type="double", default=0.3)
parser$add_argument("--k", type="integer", default=0)   # default safe
parser$add_argument("--cofactor", type="double", default=5)

args <- parser$parse_args()
set.seed(args$seed)

message("=== GateMeClass Configuration ===")
message("Train matrix: ", args$`train.data.matrix`)
message("Train labels: ", args$labels_train)
message("Test matrix:  ", args$`test.data.matrix`)
message("GMM: ", args$GMM_parameterization, ", sampling: ", args$sampling, ", k: ", args$k)
message("Cofactor: ", args$cofactor)
message("Output dir: ", args$output_dir)
message("Name: ", args$name)
message("=================================\n")

dir.create(args$output_dir, recursive=TRUE, showWarnings=FALSE)

# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------
sanitize_names <- function(x) {
  x <- gsub("[^A-Za-z0-9_]", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

keep_numeric_only <- function(dt) {
  is_num <- vapply(dt, is.numeric, logical(1))
  cols <- names(dt)[is_num]
  dropped <- names(dt)[!is_num]
  if (length(dropped) > 0) {
    message("  Dropping non-numeric columns (metadata): ", paste(dropped, collapse=", "))
  }
  dt[, cols, with = FALSE]
}

# ---------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------
message("Loading training data...")
train_dt <- fread(args$`train.data.matrix`)
train_y  <- fread(args$labels_train, header=FALSE)[[1]]

message("Loading test data...")
test_dt <- fread(args$`test.data.matrix`)

if (nrow(train_dt) != length(train_y)) {
  stop(sprintf("Training rows (%d) != train labels (%d)", nrow(train_dt), length(train_y)))
}

# Sanitize names and remove accidental 'col'
names(train_dt) <- sanitize_names(names(train_dt))
names(test_dt)  <- sanitize_names(names(test_dt))

if ("col" %in% names(train_dt)) { train_dt <- train_dt[, !"col", with=FALSE]; message("  Removed 'col' from train") }
if ("col" %in% names(test_dt))  { test_dt  <- test_dt[,  !"col", with=FALSE]; message("  Removed 'col' from test") }

# Keep numeric only
train_dt <- keep_numeric_only(train_dt)
test_dt  <- keep_numeric_only(test_dt)

# Align markers by name (and reorder test columns to match train)
missing <- setdiff(union(names(train_dt), names(test_dt)), intersect(names(train_dt), names(test_dt)))
if (length(missing) > 0) stop(paste("Marker mismatch between train/test:", paste(missing, collapse=", ")))
test_dt <- test_dt[, names(train_dt), with=FALSE]

message("Train: ", nrow(train_dt), " cells x ", ncol(train_dt), " markers")
message("Test:  ", nrow(test_dt),  " cells x ", ncol(test_dt),  " markers")

# Clean train labels
train_y_char <- as.character(train_y)
valid <- !is.na(train_y_char) & train_y_char != "" & train_y_char != '""'
if (sum(!valid) > 0) {
  message("Dropping ", sum(!valid), " unlabeled cells from TRAIN")
  train_dt <- train_dt[valid, ]
  train_y_char <- train_y_char[valid]
}
train_labels <- train_y_char
names(train_labels) <- NULL
message("Training labels (unique): ", paste(sort(unique(train_labels)), collapse=", "))

# ---------------------------------------------------------------------
# Marker renaming (robust): rename BOTH by position
# ---------------------------------------------------------------------
if (ncol(train_dt) != ncol(test_dt)) {
  stop(sprintf("Train/Test have different number of markers: %d vs %d",
               ncol(train_dt), ncol(test_dt)))
}

message("Renaming markers to simple names (M1..Mn) by POSITION (robust).")

train_old <- names(train_dt)
test_old  <- names(test_dt)
simple_markers <- paste0("M", seq_len(ncol(train_dt)))

# Save mapping (train/test originals may differ)
mapping_file <- file.path(args$output_dir, paste0(args$name, "_marker_mapping.tsv"))
fwrite(
  data.table(Simplified = simple_markers,
             TrainOriginal = train_old,
             TestOriginal  = test_old),
  mapping_file,
  sep = "\t"
)
message("Saved marker mapping: ", mapping_file)

setnames(train_dt, train_old, simple_markers)
setnames(test_dt,  test_old,  simple_markers)

# ---------------------------------------------------------------------
# Transform + transpose (GateMeClass expects markers x cells)
# ---------------------------------------------------------------------
message("Applying arcsinh transformation (cofactor=", args$cofactor, ")...")

train_m <- as.matrix(train_dt)
test_m  <- as.matrix(test_dt)

train_m <- asinh(train_m / args$cofactor)
test_m  <- asinh(test_m  / args$cofactor)

train_m <- t(train_m)
test_m  <- t(test_m)

rownames(train_m) <- simple_markers
rownames(test_m)  <- simple_markers

# ---------------------------------------------------------------------
# Method 3: train + annotate in one step
# ---------------------------------------------------------------------
message("\n=== Running GateMeClass_annotate (Method 3) ===")

k_to_use <- if (is.null(args$k) || args$k <= 0) NULL else args$k

res <- GateMeClass_annotate(
  exp_matrix = test_m,
  marker_table = NULL,
  train_parameters = list(
    reference = train_m,
    labels = train_labels
  ),
  GMM_parameterization = args$GMM_parameterization,
  reject_option = FALSE,
  sampling = args$sampling,
  k = k_to_use,
  verbose = TRUE,
  seed = args$seed
)

pred <- res$labels
message("Predicted labels: ", length(pred), " cells; unique=", length(unique(pred)))

# ---------------------------------------------------------------------
# Save predictions: numeric + ".0", NA -> ""
# ---------------------------------------------------------------------
outfile <- file.path(args$output_dir, paste0(args$name, "_predicted_labels.txt"))

pred_num <- suppressWarnings(as.numeric(pred))
out <- paste0(pred_num, ".0")
out[is.na(pred_num)] <- NA

write.table(
  x = out,
  file = outfile,
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE,
  na = '""'
)

message("Saved predictions: ", outfile)
message("\nDone.")


