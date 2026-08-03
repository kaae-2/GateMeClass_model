#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
test_dir <- dirname(normalizePath(script_path))

suppressPackageStartupMessages({
  library(BiocNeighbors)
  library(caret)
  library(data.table)
  library(dplyr)
  library(mclust)
  library(plyr)
  library(stringi)
  library(stringr)
})

source(file.path(test_dir, "..", "gatemeclass_model", "R", "GateMeClass.R"))

failures <- character(0)

expect_true <- function(value, description) {
  if (!isTRUE(value)) {
    failures <<- c(failures, description)
  }
}

predict_backend <- function(backend, training_x, training_labels, query, k, chunk_size = 17L) {
  set.seed(101)
  fit <- gatemeclass_fit_knn(training_x, training_labels, k, backend)
  as.character(gatemeclass_predict_knn(fit, query, chunk_size))
}

predict_legacy_caret <- function(training_x, training_labels, query, k) {
  training_set <- data.frame(labels = as.character(training_labels), training_x)
  set.seed(101)
  fit <- caret::train(
    labels ~ .,
    data = training_set,
    method = "knn",
    trControl = caret::trainControl(method = "none"),
    tuneGrid = data.frame(k = k)
  )
  as.character(predict(fit, newdata = query))
}

set.seed(7)
class_names <- c("Alpha / mature", "Beta cells", "Gamma")
centers <- rbind(rep(-3, 10), rep(0, 10), rep(3, 10))
training_labels <- factor(rep(class_names, each = 200))
training_x <- do.call(rbind, lapply(seq_len(3), function(index) {
  matrix(rnorm(2000, centers[index, ], 0.25), ncol = 10, byrow = TRUE)
}))
colnames(training_x) <- paste0("M", seq_len(ncol(training_x)))
query <- rbind(
  training_x[c(1:25, 201:225, 401:425), ],
  matrix(rnorm(250, 0, 0.25), ncol = 10)
)
colnames(query) <- colnames(training_x)

predictions <- lapply(gatemeclass_knn_backends, function(backend) {
  predict_backend(backend, training_x, training_labels, query, k = 20)
})
names(predictions) <- gatemeclass_knn_backends

expect_true(
  all(vapply(predictions, identical, logical(1), predictions[["caret"]])),
  "backends differ on the no-tie k=20 fixture"
)
expect_true(
  identical(sort(unique(predictions[["kmknn"]])), sort(class_names)),
  "non-syntactic class labels were not preserved"
)

set.seed(101)
kmknn_fit <- gatemeclass_fit_knn(training_x, training_labels, 20, "kmknn")
set.seed(202)
kmknn_one_chunk <- as.character(gatemeclass_predict_knn(kmknn_fit, query, nrow(query)))
set.seed(202)
kmknn_many_chunks <- as.character(gatemeclass_predict_knn(kmknn_fit, query, 7L))
expect_true(
  identical(kmknn_one_chunk, kmknn_many_chunks),
  "KMKNN query chunking changed no-tie predictions"
)

permutation <- rev(seq_len(nrow(query)))
permuted <- predict_backend("kmknn", training_x, training_labels, query[permutation, ], 20)
expect_true(
  identical(permuted[order(permutation)], predictions[["kmknn"]]),
  "KMKNN changed query-row ordering"
)

small_x <- training_x[1:5, , drop = FALSE]
small_labels <- factor(rep("Alpha / mature", 5))
for (backend in gatemeclass_knn_backends) {
  set.seed(303)
  fit <- gatemeclass_fit_knn(small_x, small_labels, 20, backend)
  prediction <- as.character(gatemeclass_predict_knn(fit, query[1:3, ], 2L))
  expect_true(identical(fit$k, 5L), paste(backend, "did not cap k to training rows"))
  expect_true(
    identical(prediction, rep("Alpha / mature", 3)),
    paste(backend, "failed after capping k")
  )
}

tie_x <- rbind(c(-1, 0), c(1, 0), c(0, -1), c(0, 1))
colnames(tie_x) <- c("M1", "M2")
tie_labels <- factor(c("A", "A", "B", "B"))
tie_query <- matrix(c(0, 0), nrow = 1, dimnames = list(NULL, colnames(tie_x)))
for (backend in gatemeclass_knn_backends) {
  first <- predict_backend(backend, tie_x, tie_labels, tie_query, k = 4)
  second <- predict_backend(backend, tie_x, tie_labels, tie_query, k = 4)
  expect_true(identical(first, second), paste(backend, "is not reproducible for a vote tie"))
  expect_true(first %in% c("A", "B"), paste(backend, "returned an invalid tied class"))
}
expect_true(
  identical(
    predict_backend("caret", tie_x, tie_labels, tie_query, k = 4),
    predict_legacy_caret(tie_x, tie_labels, tie_query, k = 4)
  ),
  "optimized caret changed legacy vote-tie behavior"
)

exp_matrix <- rbind(
  M1 = c(10.0, 9.8, 0.0, 0.2, 10.2, 9.9, 0.1, 0.3),
  M2 = c(0.0, 0.2, 10.0, 9.8, 10.2, 9.9, 0.1, 0.3)
)
marker_table <- data.frame(
  Cell = c("A", "B"),
  M1 = c("+", "*"),
  M2 = c("*", "+"),
  check.names = FALSE
)
annotate_backend <- function(backend, return_cell_signatures) {
  GateMeClass_annotate(
    exp_matrix = exp_matrix,
    marker_table = marker_table,
    GMM_parameterization = "E",
    reject_option = FALSE,
    k = 1,
    sampling = 1,
    narrow_marker_table = FALSE,
    verbose = FALSE,
    seed = 17,
    diagnostics = TRUE,
    knn_backend = backend,
    knn_query_chunk_size = 2L,
    return_cell_signatures = return_cell_signatures
  )
}
annotated <- lapply(gatemeclass_knn_backends, annotate_backend, return_cell_signatures = FALSE)
annotated_with_signatures <- lapply(
  gatemeclass_knn_backends,
  annotate_backend,
  return_cell_signatures = TRUE
)
names(annotated) <- gatemeclass_knn_backends
names(annotated_with_signatures) <- gatemeclass_knn_backends
for (backend in gatemeclass_knn_backends) {
  expect_true(
    identical(annotated[[backend]]$labels, annotated_with_signatures[[backend]]$labels),
    paste(backend, "labels-only mode changed annotation labels")
  )
}
expect_true(
  all(vapply(annotated, function(result) is.null(result$cell_signatures), logical(1))),
  "annotation labels-only mode retained cell signatures"
)
expect_true(
  identical(
    unname(vapply(annotated, function(result) result$diagnostics$knn_backend, character(1))),
    gatemeclass_knn_backends
  ),
  "annotation diagnostics do not identify the selected backend"
)

if (length(failures)) {
  stop(paste(failures, collapse = "\n"), call. = FALSE)
}

message("GateMeClass KNN backend tests passed")
