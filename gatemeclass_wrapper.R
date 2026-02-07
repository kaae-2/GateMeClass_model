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
Sys.setenv(R_LIBS_USER = local_lib)

thread_envs <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")
if (all(Sys.getenv(thread_envs) == "")) {
  cores <- parallel::detectCores(logical = TRUE)
  if (is.na(cores) || cores < 1) {
    cores <- 1
  }
  cores <- min(cores, 3)
  Sys.setenv(
    OMP_NUM_THREADS = cores,
    OPENBLAS_NUM_THREADS = cores,
    MKL_NUM_THREADS = cores,
    VECLIB_MAXIMUM_THREADS = cores
  )
}

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

gate_source <- file.path(script_dir, "gatemeclass_model", "R", "GateMeClass.R")
if (file.exists(gate_source)) {
  source(gate_source)
} else {
  stop("Missing GateMeClass source at ", gate_source, ". Vendor the module into the repo.")
}

orig_parse_marker_table <- parse_marker_table
parse_marker_table <- function(marker_table, narrow_marker_table, extended_marker_table) {
  if (!is.null(marker_table) && "Cell" %in% names(marker_table)) {
    marker_table$Cell <- make.unique(as.character(marker_table$Cell))
    if (any(duplicated(marker_table$Cell))) {
      marker_table$Cell <- sprintf("%s__%d", as.character(marker_table$Cell), seq_len(nrow(marker_table)))
    }
  }
  orig_parse_marker_table(marker_table, narrow_marker_table, extended_marker_table)
}

assign("parse_marker_table", parse_marker_table, envir = environment(GateMeClass_train))
assign("parse_marker_table", parse_marker_table, envir = environment(GateMeClass_annotate))

if (requireNamespace("caret", quietly = TRUE)) {
  assign("train", get("train", envir = asNamespace("caret")), envir = .GlobalEnv)
}

parser <- ArgumentParser(description = "GateMeClass wrapper for omnibenchmark pipeline")
parser$add_argument("--data.train_matrix", type = "character", required = TRUE,
                    help = "tar.gz archive containing training matrix CSV")
parser$add_argument("--data.train_labels", type = "character", required = TRUE,
                    help = "tar.gz archive containing training labels CSV")
parser$add_argument("--data.test_matrix", type = "character", required = TRUE,
                    help = "tar.gz archive containing test matrix CSVs")
parser$add_argument("--data.label_key", type = "character", required = FALSE, default = NULL,
                    help = "Optional: JSON.gz with id_to_label mapping")
parser$add_argument("--seed", type = "integer", default = 42)
parser$add_argument("--output_dir", "-o", type = "character", default = getwd())
parser$add_argument("--name", "-n", type = "character", required = TRUE)
parser$add_argument("--GMM_parameterization", type = "character", default = "V",
                    help = "GMM variance: 'V' (Variable) or 'E' (Equal)")
parser$add_argument("--sampling", type = "double", default = 0.1,
                    help = "Fraction of cells to use (0.0-1.0)")
parser$add_argument("--k", type = "integer", default = 20,
                    help = "k parameter for KNN refinement")
parser$add_argument("--cofactor", type = "double", default = 5,
                    help = "Cofactor for arcsinh transformation")
parser$add_argument("--sampling_imp_vars", type = "double", default = -1,
                    help = "Fraction of training cells for variable-importance step (<=0 defaults to --sampling / 10)")

args <- parser$parse_args()
set.seed(args$seed)

message("GateMeClass: starting")

output_dir_abs <- normalizePath(args$output_dir, mustWork = FALSE)
dir.create(output_dir_abs, recursive = TRUE, showWarnings = FALSE)

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

sanitize_matrix_dt <- function(dt, sample_name) {
  na_count <- sum(is.na(dt))
  if (na_count > 0) {
    message(sprintf("GateMeClass: replacing %d missing values with 0 in '%s'", na_count, sample_name))
    dt[is.na(dt)] <- 0
  }

  m <- as.matrix(dt)
  non_finite <- !is.finite(m)
  non_finite_count <- sum(non_finite)
  if (non_finite_count > 0) {
    message(sprintf("GateMeClass: replacing %d non-finite values with 0 in '%s'", non_finite_count, sample_name))
    m[non_finite] <- 0
    dt <- as.data.table(m)
  }

  dt
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

message("GateMeClass: loading data")

tmp_root <- file.path(tempdir(), paste0("gatemeclass_", Sys.getpid()))
unlink(tmp_root, recursive = TRUE)
dir.create(tmp_root, showWarnings = FALSE, recursive = TRUE)

train_members <- list_csv_members(args$`data.train_matrix`)
if (length(train_members) == 0) {
  stop("No CSV files found in archive: ", args$`data.train_matrix`)
}
train_extract_dir <- file.path(tmp_root, "train_matrix")
dir.create(train_extract_dir, showWarnings = FALSE, recursive = TRUE)
train_matrix_path <- extract_member(args$`data.train_matrix`, train_members[[1]], train_extract_dir)
train_dt <- fread(train_matrix_path, header = FALSE)
train_dt <- sanitize_matrix_dt(train_dt, "train")

label_members <- list_csv_members(args$`data.train_labels`)
if (length(label_members) == 0) {
  stop("No CSV files found in archive: ", args$`data.train_labels`)
}
train_labels_extract_dir <- file.path(tmp_root, "train_labels")
dir.create(train_labels_extract_dir, showWarnings = FALSE, recursive = TRUE)
train_labels_path <- extract_member(args$`data.train_labels`, label_members[[1]], train_labels_extract_dir)
train_labels_dt <- fread(train_labels_path, header = FALSE)
train_y <- train_labels_dt[[1]]

test_members <- list_csv_members(args$`data.test_matrix`)
if (length(test_members) == 0) {
  stop("No CSV files found in archive: ", args$`data.test_matrix`)
}
test_sample_names <- vapply(test_members, clean_member_name, character(1))

label_key <- load_label_key(args$`data.label_key`)

if (nrow(train_dt) != length(train_y)) {
  stop(sprintf("Training rows (%d) != train labels (%d)", nrow(train_dt), length(train_y)))
}

message("GateMeClass: preparing training data")

unlabeled_mask <- is.na(train_y) | train_y == 0 | train_y == -1
if (sum(unlabeled_mask) > 0) {
  train_dt <- train_dt[!unlabeled_mask, ]
  train_y <- train_y[!unlabeled_mask]
}

if (!is.null(label_key)) {
  train_labels <- sapply(train_y, function(y) {
    key <- as.character(y)
    if (key %in% names(label_key)) label_key[[key]] else paste0("Type_", y)
  })
} else {
  train_labels <- paste0("Type_", train_y)
}

n_markers <- ncol(train_dt)
simple_markers <- paste0("M", seq_len(n_markers))
setnames(train_dt, names(train_dt), simple_markers)

message("GateMeClass: transforming data")
train_m <- as.matrix(train_dt)
train_m <- asinh(train_m / args$cofactor)
train_m <- t(train_m)
rownames(train_m) <- simple_markers

message("GateMeClass: running predictions")
k_to_use <- if (is.null(args$k) || args$k <= 0) 20 else args$k
sampling_imp_vars_to_use <- args$sampling_imp_vars
if (is.na(sampling_imp_vars_to_use) || sampling_imp_vars_to_use <= 0) {
  sampling_imp_vars_to_use <- args$sampling / 10
}
sampling_imp_vars_to_use <- max(min(sampling_imp_vars_to_use, 1), 1e-06)

message("GateMeClass: training marker table")
marker_table <- GateMeClass_train(
  reference = train_m,
  labels = train_labels,
  GMM_parameterization = args$GMM_parameterization,
  sampling_imp_vars = sampling_imp_vars_to_use,
  seed = args$seed,
  verbose = FALSE
)

test_extract_dir <- file.path(tmp_root, "test_samples")
dir.create(test_extract_dir, showWarnings = FALSE, recursive = TRUE)

process_sample <- function(idx) {
  test_member <- test_members[[idx]]
  test_name <- test_sample_names[[idx]]

  tryCatch({
    test_path <- extract_member(args$`data.test_matrix`, test_member, test_extract_dir)
    test_dt <- fread(test_path, header = FALSE)
    unlink(test_path)

    if (ncol(test_dt) != n_markers) {
      stop(sprintf("Test sample '%s' has %d markers, expected %d",
                   test_name, ncol(test_dt), n_markers))
    }
    test_dt <- sanitize_matrix_dt(test_dt, test_name)
    setnames(test_dt, names(test_dt), simple_markers)

    test_m <- as.matrix(test_dt)
    test_m <- asinh(test_m / args$cofactor)
    test_m <- t(test_m)
    rownames(test_m) <- simple_markers

    res <- GateMeClass_annotate(
      exp_matrix = test_m,
      marker_table = marker_table,
      GMM_parameterization = args$GMM_parameterization,
      reject_option = TRUE,
      sampling = args$sampling,
      k = k_to_use,
      verbose = FALSE,
      seed = args$seed
    )

    list(ok = TRUE, name = test_name, labels = res$labels)
  }, error = function(e) {
    list(ok = FALSE, name = test_name, error = conditionMessage(e))
  })
}

cores <- parallel::detectCores(logical = TRUE)
if (is.na(cores) || cores < 1) {
  cores <- 1
}
cores <- min(cores, 3)
if (length(test_members) > 1 && cores > 1) {
  results <- parallel::mclapply(seq_along(test_members), process_sample, mc.cores = cores, mc.preschedule = FALSE)
} else {
  results <- lapply(seq_along(test_members), process_sample)
}

all_predictions <- list()
for (res in results) {
  if (is.null(res$ok) || !isTRUE(res$ok)) {
    stop(sprintf("GateMeClass failed for sample '%s': %s", res$name, res$error))
  }
  all_predictions[[res$name]] <- res$labels
}

message("GateMeClass: writing archive")

if (!is.null(label_key)) {
  label_to_id <- setNames(as.integer(names(label_key)), unlist(label_key))
} else {
  label_to_id <- NULL
}

tmp_pred_dir <- file.path(tempdir(), paste0("predictions_", Sys.getpid()))
unlink(tmp_pred_dir, recursive = TRUE)
dir.create(tmp_pred_dir, showWarnings = FALSE, recursive = TRUE)

for (idx in seq_along(test_sample_names)) {
  test_name <- test_sample_names[[idx]]
  pred_labels <- all_predictions[[test_name]]

  if (!is.null(label_to_id)) {
    pred_int <- label_to_id[pred_labels]
    pred_int[is.na(pred_int)] <- NA
  } else {
    pred_int <- as.integer(gsub("^Type_", "", pred_labels))
  }

  out_labels <- ifelse(is.na(pred_int), "", sprintf("%.1f", as.numeric(pred_int)))

  sample_number <- get_sample_number(test_name, idx)
  tmp_file <- file.path(tmp_pred_dir, sprintf("%s-prediction-%s.csv", args$name, sample_number))
  writeLines(out_labels, tmp_file)
}

pred_archive <- file.path(output_dir_abs, paste0(args$name, "_predicted_labels.tar.gz"))

if (file.exists(pred_archive)) {
  file.remove(pred_archive)
}

old_wd <- getwd()
setwd(tmp_pred_dir)
files_to_tar <- sort(list.files(".", pattern = "\\.csv$"))
tar(pred_archive, files = files_to_tar, compression = "gzip")
setwd(old_wd)

unlink(tmp_pred_dir, recursive = TRUE)
unlink(tmp_root, recursive = TRUE)

if (!file.exists(pred_archive)) {
  stop("Failed to create predictions archive!")
}

message("GateMeClass: done")
