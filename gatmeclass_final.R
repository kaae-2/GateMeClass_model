#!/usr/bin/env Rscript
#!/usr/bin/env Rscript

# GateMeClass wrapper for omnibenchmark
# Handles unlabeled cells in training data

library(argparse)

# Auto-install GateMeClass if not available
if (!require("GateMeClass", quietly = TRUE)) {
  message("GateMeClass not found, installing from GitHub...")
  if (!require("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org")
  }
  remotes::install_github("simo1c/GateMeClass")
}


library(GateMeClass)
library(data.table)


read_matrix <- function(path) {
  message("Reading matrix: ", path)
  dt <- fread(path, header = TRUE, data.table = FALSE)
  as.matrix(dt)
}

read_labels <- function(path) {
  message("Reading labels: ", path)
  dt <- fread(path, header = FALSE, data.table = FALSE)
  if (ncol(dt) != 1L) {
    stop("Expected label file with 1 column in ", path,
         ", got ", ncol(dt), " columns.")
  }
  as.vector(dt[[1L]])
}

main <- function() {
  parser <- ArgumentParser(
    description = "GateMeClass method for Omnibenchmark CyTOF benchmark"
  )

  ## Standard omnibenchmark arguments
  parser$add_argument(
    "--output_dir", "-o",
    dest = "output_dir",
    type = "character",
    required = TRUE,
    help = "Output directory where prediction file will be written."
  )
  parser$add_argument(
    "--name", "-n",
    dest = "name",
    type = "character",
    required = TRUE,
    help = "Dataset name (used in output filename)."
  )

  ## Inputs from the 'preprocessing' stage
  # IDs must match the YAML input IDs; dest names can be anything
  parser$add_argument(
    "--test.data.matrix",
    dest = "test_data_matrix",
    type = "character",
    required = TRUE,
    help = "Path to test data matrix (.gz)."
  )
  parser$add_argument(
    "--labels_test",
    dest = "labels_test",
    type = "character",
    required = TRUE,
    help = "Path to test labels (.gz). Not used for training."
  )
  parser$add_argument(
    "--train.data.matrix",
    dest = "train_data_matrix",
    type = "character",
    required = TRUE,
    help = "Path to train data matrix (.gz)."
  )
  parser$add_argument(
    "--labels_train",
    dest = "labels_train",
    type = "character",
    required = TRUE,
    help = "Path to train labels (.gz)."
  )

  ## GateMeClass hyperparameters (wired from YAML `parameters`)
  parser$add_argument(
    "--GMM_parameterization",
    dest   = "gmm_param",
    type   = "character",
    default = "V",
    help   = "GMM variance parameterization: 'V' (variable) or 'E' (equal)."
  )
  parser$add_argument(
    "--sampling",
    dest   = "sampling",
    type   = "double",
    default = 0.1,
    help   = "Fraction of cells used internally by GateMeClass."
  )
  parser$add_argument(
    "--k",
    dest   = "k",
    type   = "integer",
    default = 20L,
    help   = "k for k-NN/MNN label refinement."
  )
  parser$add_argument(
    "--seed",
    dest   = "seed",
    type   = "integer",
    default = 1L,
    help   = "Random seed."
  )

  args <- parser$parse_args()
  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)
  set.seed(args$seed)

  ## 1) Load training data
  train_x <- read_matrix(args$train_data_matrix)
  train_y <- read_labels(args$labels_train)

  if (nrow(train_x) != length(train_y)) {
    stop(
      "Mismatch between training matrix rows (", nrow(train_x),
      ") and training labels length (", length(train_y), ")."
    )
  }

  message("Fitting GateMeClass marker table on training data ...")
  marker_table <- GateMeClass_train(
    reference             = train_x,
    labels                = train_y,
    GMM_parameterization  = args$gmm_param,
    verbose               = TRUE,
    seed                  = args$seed
  )

  ## 2) Load test data
  test_x <- read_matrix(args$test_data_matrix)

  message("Annotating test data with GateMeClass ...")
  res <- GateMeClass_annotate(
    exp_matrix           = test_x,
    marker_table         = marker_table,
    reject_option        = FALSE,
    GMM_parameterization = args$gmm_param,
    k                    = args$k,
    sampling             = args$sampling,
    narrow_marker_table  = FALSE,
    verbose              = TRUE,
    seed                 = args$seed
  )

  pred_labels <- as.character(res$labels)

  if (length(pred_labels) != nrow(test_x)) {
    stop(
      "GateMeClass returned ", length(pred_labels),
      " labels for ", nrow(test_x), " test cells (",
      nrow(test_x), ")."
    )
  }

  ## 3) Write predictions in the format expected by the benchmark
  out_file <- file.path(
    args$output_dir,
    paste0(args$name, "_predicted_labels.txt")
  )

  message("Writing predictions to: ", out_file)
  write.table(
    pred_labels,
    file      = out_file,
    quote     = FALSE,
    sep       = "\n",
    row.names = FALSE,
    col.names = FALSE
  )

  message("Done.")
}

if (identical(environment(), globalenv())) {
  tryCatch(
    main(),
    error = function(e) {
      message("Fatal error in GateMeClass module: ", conditionMessage(e))
      quit(status = 1)
    }
  )
}
