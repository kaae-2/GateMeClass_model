#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
test_dir <- dirname(normalizePath(script_path))

suppressPackageStartupMessages({
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

annotate <- function(diagnostics) {
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
    diagnostics = diagnostics,
    knn_backend = "caret"
  )
}

without_messages <- capture.output(without_diagnostics <- annotate(FALSE), type = "message")
diagnostic_messages <- capture.output(with_diagnostics <- annotate(TRUE), type = "message")

expect_true(
  identical(without_diagnostics$labels, with_diagnostics$labels),
  "diagnostics changed predictions"
)
expect_true(is.null(without_diagnostics$diagnostics), "disabled diagnostics were returned")
expect_true(is.list(with_diagnostics$diagnostics), "enabled diagnostics were not returned")
expected_messages <- c(
  "GateMeClass diagnostics - marker-expression GMM start",
  "GateMeClass diagnostics - marker-expression GMM end",
  "GateMeClass diagnostics - cell_classification start",
  "GateMeClass diagnostics - cell_classification end",
  "GateMeClass diagnostics - KNN fit start",
  "GateMeClass diagnostics - KNN fit end",
  "GateMeClass diagnostics - KNN predict start",
  "GateMeClass diagnostics - KNN predict end"
)
expect_true(!length(without_messages), "disabled diagnostics emitted phase messages")
expect_true(identical(diagnostic_messages, expected_messages), "diagnostic phase messages are invalid")

diagnostics <- with_diagnostics$diagnostics
expected_fields <- c(
  "total_cells",
  "requested_sampling",
  "initially_sampled_cells",
  "marker_expression_gmm_elapsed_seconds",
  "cell_classification_elapsed_seconds",
  "directly_classified_cells",
  "directly_unclassified_cells",
  "knn_entered",
  "knn_backend",
  "knn_training_rows",
  "knn_prediction_rows",
  "knn_fit_elapsed_seconds",
  "knn_predict_elapsed_seconds",
  "final_unclassified_cells",
  "total_annotate_elapsed_seconds"
)
expect_true(identical(names(diagnostics), expected_fields), "diagnostic fields differ from the schema")
expect_true(identical(diagnostics$total_cells, 8L), "total cell count is invalid")
expect_true(identical(diagnostics$requested_sampling, 1), "requested sampling is invalid")
expect_true(identical(diagnostics$initially_sampled_cells, 8L), "sampled cell count is invalid")
expect_true(
  diagnostics$directly_classified_cells + diagnostics$directly_unclassified_cells == diagnostics$total_cells,
  "direct classification counters do not cover all cells"
)
expect_true(isTRUE(diagnostics$knn_entered), "sampling=1 unclassified cells did not enter KNN")
expect_true(identical(diagnostics$knn_backend, "caret"), "KNN backend was not recorded")
expect_true(
  diagnostics$knn_training_rows == diagnostics$directly_classified_cells,
  "KNN training row count is invalid"
)
expect_true(
  diagnostics$knn_prediction_rows == diagnostics$directly_unclassified_cells,
  "KNN prediction row count is invalid"
)

timing_fields <- c(
  "marker_expression_gmm_elapsed_seconds",
  "cell_classification_elapsed_seconds",
  "knn_fit_elapsed_seconds",
  "knn_predict_elapsed_seconds",
  "total_annotate_elapsed_seconds"
)
expect_true(
  all(vapply(diagnostics[timing_fields], function(value) {
    is.numeric(value) && length(value) == 1 && is.finite(value) && value >= 0
  }, logical(1))),
  "timing fields must be finite non-negative numeric scalars"
)
expect_true(
  is.numeric(diagnostics$final_unclassified_cells) &&
    length(diagnostics$final_unclassified_cells) == 1 &&
    diagnostics$final_unclassified_cells >= 0 &&
    diagnostics$final_unclassified_cells <= diagnostics$total_cells,
  "final unclassified count is invalid"
)

without_signatures <- GateMeClass_annotate(
  exp_matrix = exp_matrix,
  marker_table = marker_table,
  GMM_parameterization = "E",
  reject_option = FALSE,
  k = 1,
  sampling = 1,
  narrow_marker_table = FALSE,
  verbose = FALSE,
  seed = 17,
  diagnostics = FALSE,
  knn_backend = "caret",
  return_cell_signatures = FALSE
)
expect_true(
  identical(without_diagnostics$labels, without_signatures$labels),
  "labels-only mode changed predictions"
)
expect_true(is.null(without_signatures$cell_signatures), "labels-only mode retained signatures")

jsonl_record <- c(list(sample_name = "sample-1", test_member = "sample-1.csv"), diagnostics)
jsonl_line <- jsonlite::toJSON(jsonl_record, auto_unbox = TRUE, null = "null", na = "null")
jsonl_decoded <- jsonlite::fromJSON(jsonl_line)
expect_true(length(jsonl_line) == 1 && jsonlite::validate(jsonl_line), "diagnostics JSONL record is invalid")
expect_true(
  identical(jsonl_decoded$sample_name, "sample-1") &&
    identical(jsonl_decoded$test_member, "sample-1.csv"),
  "diagnostics JSONL record lost sample identity"
)

wrapper_source <- paste(readLines(file.path(test_dir, "..", "gatemeclass_wrapper.R")), collapse = "\n")
wrapper_patterns <- c(
  "_gatemeclass_diagnostics.partial.jsonl",
  "unlink(diagnostics_partial_path)",
  "writeLines(jsonlite::toJSON(record",
  "flush(con)",
  "GateMeClass diagnostics - sample start:",
  "GateMeClass diagnostics - sample end:",
  "append_partial_diagnostics(results[[idx]]$diagnostics)",
  "_gatemeclass_diagnostics.json",
  "--knn-backend",
  "--knn-query-chunk-size",
  "--workers",
  "return_cell_signatures = FALSE"
)
expect_true(
  all(vapply(wrapper_patterns, grepl, logical(1), x = wrapper_source, fixed = TRUE)),
  "wrapper diagnostics progress source checks failed"
)

if (length(failures)) {
  stop(paste(failures, collapse = "\n"), call. = FALSE)
}

message("GateMeClass diagnostics tests passed")
