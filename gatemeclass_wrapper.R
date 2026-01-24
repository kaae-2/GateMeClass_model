#!/usr/bin/env Rscript

if (!requireNamespace("argparse", quietly = TRUE)) {
  install.packages("argparse", repos = "https://cloud.r-project.org")
}

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
})

# ---------------------------------------------------------------------
# Install / load GateMeClass
# ---------------------------------------------------------------------
if (!requireNamespace("GateMeClass", quietly = TRUE)) {
  conda_env <- Sys.getenv("CONDA_DEFAULT_ENV")
  local_src <- file.path(script_dir, "gatemeclass_model")
  if (identical(conda_env, "gatemeclass") && dir.exists(local_src)) {
    message("GateMeClass: installing from local source at ", local_src)
    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes", repos = "https://cloud.r-project.org")
    }
    remotes::install_local(local_src, upgrade = "never")
  } else if (identical(conda_env, "gatemeclass")) {
    stop(
      "GateMeClass is missing. Clone it to ",
      local_src,
      " and rerun, or install in the env with: ",
      "R -e 'remotes::install_github(\"simo1c/GateMeClass\")'"
    )
  }
  if (!requireNamespace("GateMeClass", quietly = TRUE)) {
    message("GateMeClass: installing from GitHub...")
    message("GateMeClass not found, installing from GitHub...")
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    if (!requireNamespace("batchelor", quietly = TRUE)) {
      BiocManager::install("batchelor", ask = FALSE, update = FALSE)
    }
    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes", repos = "https://cloud.r-project.org")
    }
    tryCatch(
      remotes::install_github("simo1c/GateMeClass", upgrade = "never"),
      error = function(e) stop("GateMeClass install failed: ", conditionMessage(e))
    )
  }
}
suppressPackageStartupMessages(library(GateMeClass))

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

# Convert output_dir to absolute path IMMEDIATELY
output_dir_abs <- normalizePath(args$output_dir, mustWork = FALSE)
dir.create(output_dir_abs, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

#' Extract CSVs from tar.gz archive
#' Returns named list of data.tables, sorted by filename
#' If single_file=TRUE, returns just the first data.table (for train)
extract_from_tar <- function(tar_path, single_file = FALSE) {
  # silent extraction
  
  tmp_dir <- tempdir()
  extract_dir <- file.path(tmp_dir, paste0("extract_", basename(tar_path), "_", Sys.getpid()))
  unlink(extract_dir, recursive = TRUE)
  dir.create(extract_dir, showWarnings = FALSE, recursive = TRUE)
  
  untar(tar_path, exdir = extract_dir)
  
  # Find CSV files (may be .csv or .csv.gz)
  csv_files <- list.files(extract_dir, pattern = "\\.csv(\\.gz)?$", 
                          full.names = TRUE, recursive = TRUE)
  
  if (length(csv_files) == 0) {
    stop("No CSV files found in archive: ", tar_path)
  }
  
  # Sort for consistent ordering
  csv_files <- sort(csv_files)
  # no per-file logging
  
  # Read all files (no header - pipeline outputs have no header)
  result <- lapply(csv_files, function(f) {
    fread(f, header = FALSE)
  })
  
  # Clean names
  clean_names <- basename(csv_files)
  clean_names <- gsub("\\.csv\\.gz$", "", clean_names)
  clean_names <- gsub("\\.csv$", "", clean_names)
  clean_names <- gsub("\\.matrix$", "", clean_names)
  clean_names <- gsub("\\.labels$", "", clean_names)
  names(result) <- clean_names
  
  unlink(extract_dir, recursive = TRUE)
  
  if (single_file) {
    # Return just the data.table for single-file archives (train)
    return(result[[1]])
  }
  
  return(result)
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
message("GateMeClass: loading data")

# Train matrix (tar.gz with single CSV inside)
train_dt <- extract_from_tar(args$`data.train_matrix`, single_file = TRUE)

# Train labels (tar.gz with single CSV inside)
train_labels_dt <- extract_from_tar(args$`data.train_labels`, single_file = TRUE)
train_y <- train_labels_dt[[1]]

# Test matrices (tar.gz with multiple CSVs)
test_list <- extract_from_tar(args$`data.test_matrix`, single_file = FALSE)
# Store original order for output
test_sample_names <- names(test_list)

# Label key (optional)
label_key <- load_label_key(args$`data.label_key`)

# ---------------------------------------------------------------------
# Validate and prepare training data
# ---------------------------------------------------------------------
if (nrow(train_dt) != length(train_y)) {
  stop(sprintf("Training rows (%d) != train labels (%d)", nrow(train_dt), length(train_y)))
}

message("GateMeClass: preparing training data")

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

for (test_name in test_sample_names) {
  if (ncol(test_list[[test_name]]) != n_markers) {
    stop(sprintf("Test sample '%s' has %d markers, expected %d",
                 test_name, ncol(test_list[[test_name]]), n_markers))
  }
  setnames(test_list[[test_name]], names(test_list[[test_name]]), simple_markers)
}

# ---------------------------------------------------------------------
# Transform data (arcsinh)
# ---------------------------------------------------------------------
message("GateMeClass: transforming data")

train_m <- as.matrix(train_dt)
train_m <- asinh(train_m / args$cofactor)

test_matrices <- lapply(test_list, function(dt) {
  m <- as.matrix(dt)
  asinh(m / args$cofactor)
})

# ---------------------------------------------------------------------
# Transpose for GateMeClass (expects markers x cells)
# ---------------------------------------------------------------------
message("GateMeClass: running predictions")

train_m <- t(train_m)
rownames(train_m) <- simple_markers

test_matrices <- lapply(test_matrices, function(m) {
  m <- t(m)
  rownames(m) <- simple_markers
  m
})

invisible(train_m)

# ---------------------------------------------------------------------
# Run GateMeClass on each test sample
# ---------------------------------------------------------------------
k_to_use <- if (is.null(args$k) || args$k <= 0) 20 else args$k

all_predictions <- list()

for (test_name in test_sample_names) {
  test_m <- test_matrices[[test_name]]
  
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
  
  pred_labels <- res$labels
  all_predictions[[test_name]] <- pred_labels
}

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
for (idx in seq_along(test_sample_names)) {
  test_name <- test_sample_names[[idx]]
  pred_labels <- all_predictions[[test_name]]
  
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

if (!file.exists(pred_archive)) {
  stop("Failed to create predictions archive!")
}

message("GateMeClass: done")
