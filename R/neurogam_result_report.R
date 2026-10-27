#' Prepare a common report for neurogam analysis results
#'
#' @description
#' Internal reporting helper used by the 1D, 2D, and spatiotemporal print and
#' summary methods. Operates only on the standardized result tables and P2.2
#' metadata; never refits a model or recomputes posterior quantities.
#'
#' @param x A recognized neurogam result.
#' @return A list containing `metadata`, `clusters`, `n_clusters`, `model_class`,
#'   and `n_draws`.
#' @keywords internal
neurogam_report_data <- function (x) {

    # Allow reporting results without an attached fitted model.
    # This modifies only the local copy of the result.
    if (!"model" %in% names(x)) {
        x["model"] <- list(NULL)
    }

    view <- as_neurogam_result_view(x)
    metadata <- view$metadata

    model_class <- if (is.null(view$model)) {
        NULL
    } else {
        class(view$model)[1L]
    }

    # The metadata stores the requested number of posterior draws, when known.
    # Only query brms when a real brmsfit is supplied: test fixtures and older
    # results may have missing or non-brms models.
    n_draws <- metadata$n_post_samples
    if (inherits(view$model, "brmsfit") &&
        requireNamespace("brms", quietly = TRUE)) {
        actual_draws <- tryCatch(
            brms::ndraws(view$model),
            error = function(e) NULL
        )
        if (!is.null(actual_draws)) {
            n_draws <- actual_draws
        }
    }

    list(
        metadata = metadata,
        clusters = view$cluster_summary,
        n_clusters = nrow(view$cluster_summary),
        model_class = model_class,
        n_draws = n_draws
    )
}

#' Format one metadata value for console output
#' @keywords internal
neurogam_report_value <- function(x) {
    if (is.null(x) || length(x) == 0L || all(is.na(x))) {
        return("not recorded")
    }

    paste(as.character(x), collapse = ", ")
}

#' Print a labelled metadata item
#' @keywords internal
neurogam_report_field <- function(label, value, show_missing = FALSE) {
    if (is.null(value) && !show_missing) {
        return(invisible(NULL))
    }

    cat("  ", label, ": ", neurogam_report_value(value), "\n", sep = "")
    invisible(NULL)
}

#' Provide a readable cluster table without mutating the source
#' @keywords internal
neurogam_report_cluster_table <- function(clusters, digits) {
    if (!is.data.frame(clusters) || nrow(clusters) == 0L) {
        return(clusters)
    }

    out <- as.data.frame(clusters)
    numeric_cols <- names(out)[vapply(out, is.numeric, logical(1))]
    rounded_cols <- setdiff(numeric_cols, c("id", "n_points", "n_sensors"))

    for (col in rounded_cols) {
        out[[col]] <- round(out[[col]], digits = digits)
    }

    preferred <- c("participant", "sign", "id", "onset", "offset", "duration")
    preferred <- c(
        preferred,
        grep("(_min|_max|^n_points$|^n_sensors$)", names(out), value = TRUE)
    )
    preferred <- unique(preferred[preferred %in% names(out)])
    out[c(preferred, setdiff(names(out), preferred))]
}

#' Validate console report formatting options
#' @keywords internal
neurogam_validate_report_options <- function(digits, max_clusters) {
    if (!is.numeric(digits) || length(digits) != 1L || is.na(digits) ||
        !is.finite(digits) || digits < 0 || digits != as.integer(digits)) {
        stop("`digits` must be a nonnegative integer.", call. = FALSE)
    }

    if (!is.numeric(max_clusters) || length(max_clusters) != 1L ||
        is.na(max_clusters) || max_clusters < 1 ||
        (!is.infinite(max_clusters) &&
         (!is.finite(max_clusters) || max_clusters != as.integer(max_clusters)))) {
        stop("`max_clusters` must be a positive integer or Inf.", call. = FALSE)
    }

    invisible(TRUE)
}

#' Print a common neurogam analysis report
#' @keywords internal
neurogam_print_report <- function(
    x,
    digits = 3L,
    max_clusters = 12L,
    detailed = FALSE
) {
    neurogam_validate_report_options(digits, max_clusters)
    report <- neurogam_report_data(x)
    metadata <- report$metadata

    display_type <- switch(
        metadata$analysis_type,
        temporal = "Temporal (1D)",
        temporal_generalization = "Temporal generalization (2D)",
        spatiotemporal = "Spatiotemporal EEG (sensor x time)"
    )

    cat("\n==== Time-resolved GAM results ================================\n")
    cat("  Analysis: ", display_type, "\n", sep = "")
    neurogam_report_field("Model coordinates", metadata$model_coordinates)
    neurogam_report_field("Cluster coordinates", metadata$cluster_coordinates)
    neurogam_report_field("Contrast", metadata$contrast)
    neurogam_report_field("Posterior odds mode", metadata$odds_mode, show_missing = TRUE)
    neurogam_report_field("Odds threshold", metadata$threshold)
    cat("  Clusters found: ", report$n_clusters, "\n", sep = "")

    if (detailed) {
        cat("\nModel and inference:\n")
        neurogam_report_field("Model type", metadata$multilevel)
        neurogam_report_field("Backend", report$model_class)
        neurogam_report_field("Posterior draws", report$n_draws)
        neurogam_report_field("Predictor", metadata$predictor_id)
        neurogam_report_field("Contrast levels", metadata$contrast_levels)
        neurogam_report_field("Threshold direction", metadata$threshold_type)
        neurogam_report_field("Null value", metadata$null_value)
        neurogam_report_field("ROPE half-width", metadata$rope)
        neurogam_report_field("Baseline window", metadata$baseline_window)
        neurogam_report_field("Credible interval", metadata$credible_interval)
        neurogam_report_field("Sensor adjacency", metadata$adjacency_source)
        neurogam_report_field("Prediction locations", metadata$n_prediction_locations)
        neurogam_report_field("Cluster member points", metadata$n_cluster_members)
    }

    clusters <- report$clusters
    if (report$n_clusters == 0L) {
        cat("\nNo clusters exceeded the threshold.\n")
        cat("=================================================================\n")
        return(invisible(x))
    }

    if (detailed) {
        cat("\nCluster statistics:\n")
        if ("sign" %in% names(clusters)) {
            sign_counts <- table(as.character(clusters$sign))
            cat("  By sign: ", paste(
                paste(names(sign_counts), as.integer(sign_counts), sep = " = "),
                collapse = "; "
            ), "\n", sep = "")
        }

        if ("n_points" %in% names(clusters)) {
            sizes <- clusters$n_points
            sizes <- sizes[is.finite(sizes)]
            if (length(sizes)) {
                cat("  Total cluster points: ", sum(sizes), "\n", sep = "")
                cat("  Median cluster size: ",
                    round(stats::median(sizes), digits), "\n", sep = "")
            }
        }

        duration <- if ("duration" %in% names(clusters)) {
            clusters$duration
        } else if (all(c("onset", "offset") %in% names(clusters))) {
            clusters$offset - clusters$onset
        } else {
            numeric(0)
        }
        duration <- duration[is.finite(duration)]
        if (length(duration)) {
            cat("  Mean duration: ", round(mean(duration), digits), "\n", sep = "")
            cat("  Median duration: ",
                round(stats::median(duration), digits), "\n", sep = "")
        }

        # For a 2D temporal generalization matrix, summarize both axes
        # using the coordinate names stored by P2.2 rather than hard-coding
        # time1/time2 or reconstructing clusters from legacy cell rows.
        if (identical(metadata$analysis_type, "temporal_generalization")) {
            for (coordinate in metadata$cluster_coordinates) {
                min_col <- paste0(coordinate, "_min")
                max_col <- paste0(coordinate, "_max")
                if (all(c(min_col, max_col) %in% names(clusters))) {
                    spans <- clusters[[max_col]] - clusters[[min_col]]
                    spans <- spans[is.finite(spans)]
                    if (length(spans)) {
                        cat("  Median span (", coordinate, "): ",
                            round(stats::median(spans), digits), "\n", sep = "")
                    }
                }
            }
        }

        if ("n_sensors" %in% names(clusters)) {
            values <- clusters$n_sensors
            values <- values[is.finite(values)]
            if (length(values)) {
                cat("  Maximum sensors in one cluster: ", max(values), "\n", sep = "")
            }
        }
    }

    count_to_show <- if (is.infinite(max_clusters)) {
        report$n_clusters
    } else {
        min(report$n_clusters, as.integer(max_clusters))
    }

    cat("\nCluster table:\n\n")
    table <- neurogam_report_cluster_table(
        utils::head(clusters, count_to_show),
        digits = digits
    )
    print(table, row.names = FALSE)

    if (count_to_show < report$n_clusters) {
        cat("  ... ", report$n_clusters - count_to_show,
            " additional cluster(s) not shown.\n", sep = "")
    }
    cat("=================================================================\n")
    invisible(x)
}
