map_prediction_labels <- function(pred_labels, label_to_id, sample_name) {
  pred_labels <- as.character(pred_labels)
  pred_int <- rep(NA_integer_, length(pred_labels))
  rejected <- !is.na(pred_labels) & pred_labels == "Unclassified"
  pred_int[rejected] <- 0L
  classify <- !rejected

  if (!is.null(label_to_id)) {
    pred_int[classify] <- unname(label_to_id[pred_labels[classify]])
  } else {
    fallback_labels <- pred_labels[classify]
    parseable <- !is.na(fallback_labels) & grepl("^Type_-?[0-9]+$", fallback_labels)
    parsed <- rep(NA_integer_, length(fallback_labels))
    parsed[parseable] <- suppressWarnings(as.integer(sub("^Type_", "", fallback_labels[parseable])))
    pred_int[classify] <- parsed
  }

  unmapped <- unique(pred_labels[is.na(pred_int)])
  if (length(unmapped)) {
    unmapped[is.na(unmapped)] <- "<NA>"
    stop(sprintf(
      "Cannot map GateMeClass label(s) for sample '%s': %s",
      sample_name,
      paste(unmapped, collapse = ", ")
    ), call. = FALSE)
  }

  pred_int
}
