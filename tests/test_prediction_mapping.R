#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
test_dir <- dirname(normalizePath(script_path))
source(file.path(test_dir, "..", "prediction_mapping.R"))

failures <- character(0)

expect_equal <- function(actual, expected, description) {
  if (!identical(unname(actual), expected)) {
    failures <<- c(failures, sprintf(
      "%s: expected %s, got %s",
      description,
      paste(expected, collapse = ","),
      paste(actual, collapse = ",")
    ))
  }
}

expect_error <- function(expression, patterns, description) {
  error <- tryCatch({
    force(expression)
    NULL
  }, error = identity)
  matches <- !is.null(error) && all(vapply(
    patterns,
    grepl,
    logical(1),
    x = conditionMessage(error),
    fixed = TRUE
  ))
  if (!matches) {
    message <- if (is.null(error)) "no error" else conditionMessage(error)
    failures <<- c(failures, sprintf("%s: got %s", description, message))
  }
}

label_to_id <- c(Monocyte = 4L, B_cell = 9L)
expect_equal(
  map_prediction_labels(c("Monocyte", "Unclassified"), label_to_id, "sample-5.csv"),
  c(4L, 0L),
  "metadata names and rejection label map to benchmark ids"
)
expect_error(
  map_prediction_labels(c("Monocyte", "Mystery"), label_to_id, "sample-5.csv"),
  c("sample-5.csv", "Mystery"),
  "unexpected metadata label fails contextually"
)

expect_equal(
  map_prediction_labels(c("Type_7", "Unclassified"), NULL, "sample-8.csv"),
  c(7L, 0L),
  "fallback labels and rejection label map to benchmark ids"
)
expect_error(
  map_prediction_labels(c("Type_7", "Other"), NULL, "sample-8.csv"),
  c("sample-8.csv", "Other"),
  "non-parseable fallback label fails contextually"
)

if (length(failures)) {
  stop(paste(failures, collapse = "\n"), call. = FALSE)
}

message("prediction mapping tests passed")
