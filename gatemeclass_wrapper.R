#!/usr/bin/env Rscript

if (!requireNamespace("argparse", quietly = TRUE)) {
  stop("Missing R package 'argparse'. Update the gatemeclass conda env.")
}
if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Missing R package 'data.table'. Update the gatemeclass conda env.")
}

script_path <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_path)
script_dir <- if (length(script_path) > 0) {
  dirname(normalizePath(script_path))
} else {
  getwd()
}
local_lib <- file.path(script_dir, ".r_libs")
dir.create(local_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(local_lib, .libPaths()))

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(dplyr)
  library(stringi)
  library(stringr)
  library(caret)
  library(mclust)
  library(batchelor)
  library(plyr)
})

# ---------------------------------------------------------------------
# Load GateMeClass from local source
# ---------------------------------------------------------------------
comboGrid <- function(x, y, repetition = FALSE) {
  grid <- expand.grid(x, y, stringsAsFactors = FALSE)
  if (!repetition) {
    grid <- grid[grid[[1]] != grid[[2]], , drop = FALSE]
  }
  as.matrix(grid)
}

local_src <- file.path(script_dir, "gatemeclass_model", "R", "GateMeClass.R")
if (!file.exists(local_src)) {
  stop("Missing GateMeClass source at ", local_src)
}
source(local_src)

if (requireNamespace("caret", quietly = TRUE)) {
  assign("train", get("train", envir = asNamespace("caret")), envir = .GlobalEnv)
}

# ---------------------------------------------------------------------
# Args - matching omnibenchmark pipeline
# ---------------------------------------------------------------------
parser <- ArgumentParser(description="GateMeClass wrapper for omnibenchmark pipeline")

# Pipeline inputs (ALL are tar.gz archives in omnibenchmark!)
parser$add_argument("--data.train_matrix", type="character", required=TRUE,
                    help="tar.gz archive containing training matrix CSV")
parser$add_argument("--data.train_labels", type="character", required=TRUE,
                    help="tar.gz archive containing training labels CSV")
parser$add_argument("--data.test_matrix",  type="character", required=TRUE,
                    help="tar.gz archive containing test matrix CSVs")

# Optional label key (not used by default pipeline, but useful)
parser$add_argument("--data.label_key", type="character", required=FALSE, default=NULL,
                    help="Optional: JSON.gz with id_to_label mapping")

parser$add_argument("--seed", type="integer", default=42)
parser$add_argument("--output_dir", "-o", type="character", default=getwd())
parser$add_argument("--name", "-n", type="character", required=TRUE)

# GateMeClass parameters (tutorial-recommended defaults)
parser$add_argument("--GMM_parameterization", type="character", default="V",
                    help="GMM variance: 'V' (Variable) or 'E' (Equal)")
parser$add_argument("--sampling", type="double", default=1.0,
                    help="Fraction of cells to use (0.0-1.0)")
parser$add_argument("--k", type="integer", default=20,
                    help="k parameter for KNN refinement")
parser$add_argument("--cofactor", type="double", default=5,
                    help="Cofactor for arcsinh transformation")

args <- parser$parse_args()
set.seed(args$seed)

message("GateMeClass: starting")

log_ts <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(sprintf("GateMeClass: %s %s", timestamp, msg))
}

# Convert output_dir to absolute path IMMEDIATELY
output_dir_abs <- normalizePath(args$output_dir, mustWork = FALSE)
dir.create(output_dir_abs, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

list_csv_members <- function(tar_path) {
  members <- untar(tar_path, list = TRUE)
  csv_members <- members[grepl("\\.csv(\\.gz)?$", members)]
  sort(csv_members)
}

extract_member <- function(tar_path, member, extract_dir) {
  untar(tar_path, exdir = extract_dir, files = member)
  file.path(extract_dir, member)
}

clean_member_name <- function(file_name) {
  clean <- basename(file_name)
  clean <- gsub("\\.csv\\.gz$", "", clean)
  clean <- gsub("\\.csv$", "", clean)
  clean <- gsub("\\.matrix$", "", clean)
  clean <- gsub("\\.labels$", "", clean)
  clean
}

get_sample_number <- function(file_name, fallback) {
  base <- basename(file_name)
  base <- gsub("\\.csv(\\.gz)?$", "", base)
  m <- regexpr("[0-9]+(?!.*[0-9])", base, perl = TRUE)
  if (m[1] == -1) {
    return(as.character(fallback))
  }
  substr(base, m[1], m[1] + attr(m, "match.length") - 1)
}

#' Load label key from JSON.gz
load_label_key <- function(path) {
  if (is.null(path) || !file.exists(path)) {
    return(NULL)
  }
  
  if (grepl("\\.gz$", path)) {
    con <- gzfile(path, "rt")
    json_text <- paste(readLines(con, warn = FALSE), collapse = "")
    close(con)
  } else {
    json_text <- paste(readLines(path, warn = FALSE), collapse = "")
  }
  
  parsed <- jsonlite::fromJSON(json_text)
  
  # Handle {"id_to_label": {"1": "CD4", "2": "CD8", ...}}
  if ("id_to_label" %in% names(parsed)) {
    id_to_label <- parsed$id_to_label
    mapping <- list()
    for (id in names(id_to_label)) {
      mapping[[id]] <- id_to_label[[id]]
    }
    return(mapping)
  }
  
  return(NULL)
}

# ---------------------------------------------------------------------
# Load data from tar.gz archives
# ---------------------------------------------------------------------
log_ts("Loading data")

tmp_root <- file.path(tempdir(), paste0("gatemeclass_", Sys.getpid()))
unlink(tmp_root, recursive = TRUE)
dir.create(tmp_root, showWarnings = FALSE, recursive = TRUE)

# Train matrix (tar.gz with single CSV inside)
train_members <- list_csv_members(args$`data.train_matrix`)
if (length(train_members) == 0) {
  stop("No CSV files found in archive: ", args$`data.train_matrix`)
}
train_extract_dir <- file.path(tmp_root, "train_matrix")
dir.create(train_extract_dir, showWarnings = FALSE, recursive = TRUE)
train_start <- Sys.time()
train_matrix_path <- extract_member(args$`data.train_matrix`, train_members[[1]], train_extract_dir)
train_dt <- fread(train_matrix_path, header = FALSE)
log_ts(sprintf("Loaded train matrix in %.2fs", as.numeric(difftime(Sys.time(), train_start, units = "secs"))))

# Train labels (tar.gz with single CSV inside)
label_members <- list_csv_members(args$`data.train_labels`)
if (length(label_members) == 0) {
  stop("No CSV files found in archive: ", args$`data.train_labels`)
}
train_labels_extract_dir <- file.path(tmp_root, "train_labels")
dir.create(train_labels_extract_dir, showWarnings = FALSE, recursive = TRUE)
labels_start <- Sys.time()
train_labels_path <- extract_member(args$`data.train_labels`, label_members[[1]], train_labels_extract_dir)
train_labels_dt <- fread(train_labels_path, header = FALSE)
train_y <- train_labels_dt[[1]]
log_ts(sprintf("Loaded train labels in %.2fs", as.numeric(difftime(Sys.time(), labels_start, units = "secs"))))

# Test matrices (tar.gz with multiple CSVs)
test_members <- list_csv_members(args$`data.test_matrix`)
if (length(test_members) == 0) {
  stop("No CSV files found in archive: ", args$`data.test_matrix`)
}
test_sample_names <- vapply(test_members, clean_member_name, character(1))

# Label key (optional)
label_key <- load_label_key(args$`data.label_key`)

# ---------------------------------------------------------------------
# Validate and prepare training data
# ---------------------------------------------------------------------
if (nrow(train_dt) != length(train_y)) {
  stop(sprintf("Training rows (%d) != train labels (%d)", nrow(train_dt), length(train_y)))
}

log_ts("Preparing training data")

# Remove unlabeled cells (label = 0, -1, or NA) from training
unlabeled_mask <- is.na(train_y) | train_y == 0 | train_y == -1
n_unlabeled <- sum(unlabeled_mask)

if (n_unlabeled > 0) {
  train_dt <- train_dt[!unlabeled_mask, ]
  train_y <- train_y[!unlabeled_mask]
}

# Convert integer labels to character for GateMeClass
if (!is.null(label_key)) {
  train_labels <- sapply(train_y, function(y) {
    key <- as.character(y)
    if (key %in% names(label_key)) label_key[[key]] else paste0("Type_", y)
  })
} else {
  train_labels <- paste0("Type_", train_y)
}


# ---------------------------------------------------------------------
# Marker naming
# ---------------------------------------------------------------------
n_markers <- ncol(train_dt)
simple_markers <- paste0("M", seq_len(n_markers))
setnames(train_dt, names(train_dt), simple_markers)

# ---------------------------------------------------------------------
# Transform data (arcsinh)
# ---------------------------------------------------------------------
log_ts("Transforming data")

train_m <- as.matrix(train_dt)
train_m <- asinh(train_m / args$cofactor)

# ---------------------------------------------------------------------
# Transpose for GateMeClass (expects markers x cells)
# ---------------------------------------------------------------------
log_ts("Running predictions")

train_m <- t(train_m)
rownames(train_m) <- simple_markers

invisible(train_m)

# ---------------------------------------------------------------------
# Run GateMeClass on each test sample
# ---------------------------------------------------------------------
k_to_use <- if (is.null(args$k) || args$k <= 0) 20 else args$k

test_extract_dir <- file.path(tmp_root, "test_samples")
dir.create(test_extract_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------
# Convert predictions back to integers and save as tar.gz
# ---------------------------------------------------------------------
message("GateMeClass: writing archive")

# Create reverse mapping: label -> integer
if (!is.null(label_key)) {
  label_to_id <- setNames(as.integer(names(label_key)), unlist(label_key))
} else {
  # Extract integer from "Type_X" format
  label_to_id <- NULL
}

# Create temp directory for tar contents
tmp_pred_dir <- file.path(tempdir(), paste0("predictions_", Sys.getpid()))
unlink(tmp_pred_dir, recursive = TRUE)
dir.create(tmp_pred_dir, showWarnings = FALSE, recursive = TRUE)

# Save each test sample's predictions IN THE SAME ORDER as input
for (idx in seq_along(test_members)) {
  test_member <- test_members[[idx]]
  test_name <- test_sample_names[[idx]]
  sample_start <- Sys.time()

  test_path <- extract_member(args$`data.test_matrix`, test_member, test_extract_dir)
  test_dt <- fread(test_path, header = FALSE)
  unlink(test_path)

  if (ncol(test_dt) != n_markers) {
    stop(sprintf("Test sample '%s' has %d markers, expected %d",
                 test_name, ncol(test_dt), n_markers))
  }
  setnames(test_dt, names(test_dt), simple_markers)

  test_m <- as.matrix(test_dt)
  test_m <- asinh(test_m / args$cofactor)
  test_m <- t(test_m)
  rownames(test_m) <- simple_markers

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
    verbose = FALSE,
    seed = args$seed
  )

  pred_labels <- res$labels
  
  # Convert back to integers
  if (!is.null(label_to_id)) {
    pred_int <- label_to_id[pred_labels]
    pred_int[is.na(pred_int)] <- NA
  } else {
    # Extract from "Type_X" format
    pred_int <- as.integer(gsub("^Type_", "", pred_labels))
  }
  
  # Format as "X.0" to match DG-CyTOF output format (flow_metrics expects this)
  # Empty string for NA/unknown
  out_labels <- ifelse(is.na(pred_int), "", sprintf("%.1f", as.numeric(pred_int)))
  
  sample_number <- get_sample_number(test_name, idx)
  tmp_file <- file.path(tmp_pred_dir, sprintf("%s-prediction-%s.csv", args$name, sample_number))

  writeLines(out_labels, tmp_file)

  log_ts(sprintf("Processed %s in %.2fs", test_name, as.numeric(difftime(Sys.time(), sample_start, units = "secs"))))
  
}

# Create the tar.gz archive
pred_archive <- file.path(output_dir_abs, paste0(args$name, "_predicted_labels.tar.gz"))

if (file.exists(pred_archive)) {
  file.remove(pred_archive)
}

old_wd <- getwd()
setwd(tmp_pred_dir)

# List files in sorted order (must match test labels order!)
files_to_tar <- sort(list.files(".", pattern = "\\.csv$"))

tar(pred_archive, files = files_to_tar, compression = "gzip")

setwd(old_wd)
unlink(tmp_pred_dir, recursive = TRUE)
unlink(tmp_root, recursive = TRUE)

if (!file.exists(pred_archive)) {
  stop("Failed to create predictions archive!")
}

message("GateMeClass: done")
