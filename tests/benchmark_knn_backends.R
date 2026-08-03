#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
test_dir <- dirname(normalizePath(script_path))

suppressPackageStartupMessages({
  library(BiocNeighbors)
  library(caret)
})

source(file.path(test_dir, "..", "gatemeclass_model", "R", "GateMeClass.R"))

set.seed(42)
n_training <- 10000L
n_query <- 20000L
n_markers <- 10L
training_x <- matrix(rnorm(n_training * n_markers), ncol = n_markers)
query <- matrix(rnorm(n_query * n_markers), ncol = n_markers)
colnames(training_x) <- colnames(query) <- paste0("M", seq_len(n_markers))
training_labels <- factor(sample(c("A", "B", "C", "D"), n_training, replace = TRUE))

results <- list()
reference <- NULL
for (backend in gatemeclass_knn_backends) {
  set.seed(101)
  fit_elapsed <- system.time({
    fit <- gatemeclass_fit_knn(training_x, training_labels, 20, backend)
  })[["elapsed"]]
  predict_elapsed <- system.time({
    prediction <- as.character(gatemeclass_predict_knn(fit, query, 5000L))
  })[["elapsed"]]
  if (is.null(reference)) {
    reference <- prediction
  }
  results[[backend]] <- data.frame(
    backend = backend,
    fit_seconds = fit_elapsed,
    predict_seconds = predict_elapsed,
    disagreements_with_caret = sum(prediction != reference)
  )
}

print(do.call(rbind, results), row.names = FALSE)
