#' Print method for 1D \code{clusters_results} objects
#'
#' This method provides a concise console representation of the output from
#' \code{\link{testing_through_time}}, including the number of detected
#' clusters and a compact table summarising each cluster's onset, offset,
#' and duration. Values are rounded for readability.
#'
#' @param x An object of class \code{"clusters_results"} as returned by
#'   \code{\link{testing_through_time}}.
#' @param digits Integer; number of decimal places used when printing numeric
#'   values (default: \code{3}).
#' @param ... Additional arguments (currently ignored).
#'
#' @details
#' The printed cluster table includes:
#' \itemize{
#'   \item \code{cluster_id}: numeric identifier of the cluster;
#'   \item \code{cluster_onset}: estimated temporal onset of the cluster;
#'   \item \code{cluster_offset}: estimated temporal offset of the cluster;
#'   \item \code{duration}: duration of the cluster, computed as
#'     \code{cluster_offset - cluster_onset}.
#' }
#'
#' If no clusters exceed the posterior odds threshold, an informative message
#' is displayed and no table is printed.
#'
#' @return The input object \code{x}, returned invisibly.
#'
#' @seealso \code{\link{summary.clusters_results_1d}},
#'   \code{\link{testing_through_time}}
#'
#' @export
print.clusters_results_1d <- function (x, digits = 3, ...) {

    cat("\n==== Time-resolved GAM results ================================\n\n")

    clusters <- x$clusters
    n_clust <- if (is.null(clusters) ) 0L else nrow(clusters)

    # number of clusters
    cat("Clusters found: ", n_clust, "\n\n", sep = "")

    # if no clusters, stop early
    if (n_clust == 0) {

        cat("\nNo clusters exceed the threshold.\n\n")
        cat("=================================================================\n")

        return (invisible(x) )

    }

    # optional participant column
    has_participant <- "participant" %in% colnames(clusters)

    # helper function to round safely
    myround <- function (v) round(v, digits)

    # prepare 1D cluster table
    if ("participant" %in% colnames(x$clusters) ) {

        clust_tbl <- x$clusters |>
            dplyr::mutate(
                onset = myround(.data$onset),
                offset = myround(.data$offset),
                duration = myround(.data$offset - .data$onset)
                ) |>
            dplyr::select(
                .data$participant, .data$sign, .data$id,
                .data$onset, .data$offset, .data$duration
                )

    } else {

        clust_tbl <- x$clusters |>
            dplyr::mutate(
                onset = myround(.data$onset),
                offset = myround(.data$offset),
                duration = myround(.data$offset - .data$onset)
                ) |>
            dplyr::select(
                .data$sign, .data$id, .data$onset,
                .data$offset, .data$duration
                )

    }

    # print nicely
    print(data.frame(clust_tbl), row.names = FALSE)
    cat("\n=================================================================\n")

    return (invisible(x) )

}

#' @export
print.clusters_results_2d <- function (x, digits = 3, ...) {

    cat("\n==== Time-resolved GAM results ================================\n\n")

    clusters <- x$clusters
    n_clust <- if (is.null(clusters) ) 0L else max(clusters$id)

    # number of clusters
    cat("Clusters found: ", n_clust, "\n\n", sep = "")

    # helper function to round safely
    myround <- function (v) round(v, digits)

    # if no clusters, stop early
    if (n_clust == 0) {

        cat("\nNo clusters exceed the threshold.\n\n")
        cat("=================================================================\n")

        return (invisible(x) )

    }

    # pointwise 2D format with arbitrary time column names
    # identify the two time columns as the two numeric columns not in bookkeeping.
    cn <- colnames(clusters)
    bookkeeping <- c("id", "sign", "value", "n_points", "participant")
    candidate_cols <- setdiff(cn, bookkeeping)

    # keep only numeric candidates
    is_num <- vapply(clusters[, candidate_cols, drop = FALSE], is.numeric, logical(1) )
    time_cols <- candidate_cols[is_num]
    is_2d <- length(time_cols) == 2L

    # prepare 2D cluster table
    has_n_points <- "n_points" %in% names(clusters)

    # summarize each cluster id by its min/max on both axes
    clust_tbl <- clusters |>
        dplyr::group_by(.data$id, .data$sign, dplyr::across(dplyr::any_of("participant") ) ) |>
        dplyr::summarise(
            onset_time1 = min(.data$time1, na.rm = TRUE),
            offset_time1 = max(.data$time1, na.rm = TRUE),
            onset_time2 = min(.data$time2, na.rm = TRUE),
            offset_time2 = max(.data$time2, na.rm = TRUE),
            # n_points = dplyr::if_any(dplyr::all_of("n_points"), ~ TRUE) |> {
            #     # if n_points exists, it is constant per cluster in your examples;
            #     # keep the first non-NA; otherwise compute as number of rows per id.
            #     if ("n_points" %in% cn) dplyr::first(.data$n_points) else dplyr::n()
            # },
            n_points = if (has_n_points) dplyr::first(stats::na.omit(.data$n_points) ) else dplyr::n(),
            value_min = min(.data$value, na.rm = TRUE),
            value_max = max(.data$value, na.rm = TRUE),
            .groups = "drop"
            ) |>
        dplyr::mutate(
            span_time1 = .data$offset_time1 - .data$onset_time1,
            span_time2 = .data$offset_time2 - .data$onset_time2
            ) |>
        dplyr::mutate(
            onset_time1 = myround(.data$onset_time1),
            offset_time1 = myround(.data$offset_time1),
            span_time1 = myround(.data$span_time1),
            onset_time2 = myround(.data$onset_time2),
            offset_time2 = myround(.data$offset_time2),
            span_time2 = myround(.data$span_time2),
            value_min = myround(.data$value_min),
            value_max = myround(.data$value_max)
            )

    # optional participant column
    has_participant <- "participant" %in% colnames(clusters)

    # order columns
    if (has_participant) {

        clust_tbl <- clust_tbl |>
            dplyr::select(
                .data$participant, .data$sign, .data$id,
                dplyr::starts_with("onset_"), dplyr::starts_with("offset_"), dplyr::starts_with("span_"),
                .data$n_points,
                dplyr::any_of(c("value_min", "value_max") )
                )

    } else {

        clust_tbl <- clust_tbl |>
            dplyr::select(
                .data$sign, .data$id,
                dplyr::starts_with("onset_"), dplyr::starts_with("offset_"), dplyr::starts_with("span_"),
                .data$n_points,
                dplyr::any_of(c("value_min", "value_max") )
                )

    }

    print(data.frame(clust_tbl), row.names = FALSE)
    cat("\n=================================================================\n")

    return (invisible(x) )

}

#' Summary method for 1D \code{clusters_results} objects
#'
#' Produces a detailed textual summary of a time-resolved Bayesian GAMM
#' analysis, including model metadata, the number of clusters detected, and
#' descriptive statistics of cluster durations. A rounded cluster table is
#' also printed.
#'
#' @param object An object of class \code{"clusters_results"} created by
#'   \code{\link{testing_through_time}}.
#' @param digits Integer; number of decimal places used when printing numeric
#'   values (default: \code{3}).
#' @param ... Additional arguments (currently ignored).
#'
#' @details
#' The summary prints:
#' \itemize{
#'   \item the model type used (\code{"full"}, \code{"summary"},
#'     or \code{"group"});
#'   \item the class of the underlying \pkg{brms} model and the number of
#'     posterior draws;
#'   \item the number of clusters detected by the posterior odds threshold;
#'   \item descriptive statistics for cluster durations (minimum, maximum,
#'     mean, median, and total duration);
#'   \item a neatly formatted table listing each cluster's onset, offset,
#'     and duration.
#' }
#'
#' If no clusters were detected, the function prints a message and returns
#' invisibly.
#'
#' @return The object \code{object}, returned invisibly.
#'
#' @seealso \code{\link{print.clusters_results_1d}},
#'   \code{\link{testing_through_time}}
#'
#' @export
summary.clusters_results_1d <- function (object, digits = 3, ...) {

    cat("\n==== Time-resolved GAM results ================================\n\n")

    # model type
    if (!is.null(object$multilevel) ) {

        cat("Model type: ", object$multilevel, "\n", sep = "")

    }

    # class of backend model
    if (!is.null(object$model) ) {

        cat("Backend model: ", class(object$model)[1], "\n", sep = "")
        cat("Posterior draws: ", brms::ndraws(object$model), "\n", sep = "")

    }

    # number of clusters
    n_clust <- nrow(object$clusters)

    # no clusters, simple summary
    if (n_clust == 0) {

        cat("\nNo clusters exceeded the threshold.\n\n")

        return (invisible(object) )

    }

    # compute durations
    cl <- object$clusters |> dplyr::mutate(duration = .data$offset - .data$onset)

    # basic cluster stats
    cat("\nCluster statistics:\n")
    cat("  Mean cluster duration: ",
        round(mean(cl$duration), digits), "\n", sep = "")
    cat("  Median cluster duration: ",
        round(stats::median(cl$duration), digits), "\n", sep = "")
    cat("  Min cluster duration: ",
        round(min(cl$duration), digits), "\n", sep = "")
    cat("  Max cluster duration: ",
        round(max(cl$duration), digits), "\n", sep = "")

    # prepare cluster table
    if ("participant" %in% colnames(object$clusters) ) {

        cl_print <- cl |>
            dplyr::mutate(
                onset = round(.data$onset,  digits),
                offset = round(.data$offset, digits),
                duration = round(.data$duration, digits)
                ) |>
            dplyr::select(
                .data$participant, .data$sign, .data$id,
                .data$onset, .data$offset, .data$duration
                )

    } else {

        cl_print <- cl |>
            dplyr::mutate(
                onset = round(.data$onset,  digits),
                offset = round(.data$offset, digits),
                duration = round(.data$duration, digits)
                ) |>
            dplyr::select(
                .data$sign, .data$id,
                .data$onset, .data$offset, .data$duration
                )

    }

    cat("\nCluster table:\n\n")
    print(data.frame(cl_print), row.names = FALSE)
    cat("\n=================================================================\n")

    return (invisible(object) )

}

#' Summary method for 2D \code{clusters_results} objects
#'
#' Produces a detailed textual summary of a time-resolved Bayesian GAMM
#' analysis, including model metadata, the number of clusters detected, and
#' descriptive statistics of cluster durations. A rounded cluster table is
#' also printed.
#'
#' @param object An object of class \code{"clusters_results"} created by
#'   \code{\link{testing_through_time}}.
#' @param digits Integer; number of decimal places used when printing numeric
#'   values (default: \code{3}).
#' @param ... Additional arguments (currently ignored).
#'
#' @return The object \code{object}, returned invisibly.
#'
#' @export
summary.clusters_results_2d <- function (object, digits = 3, ...) {

    cat("\n==== Time-resolved GAM results ================================\n\n")

    # model type
    if (!is.null(object$multilevel) ) {

        cat("Model type: ", object$multilevel, "\n", sep = "")

    }

    # class of backend model
    if (!is.null(object$model) ) {

        cat("Backend model: ", class(object$model)[1], "\n", sep = "")
        cat("Posterior draws: ", brms::ndraws(object$model), "\n", sep = "")
    }

    clusters <- object$clusters

    # no clusters object
    if (is.null(clusters) || nrow(clusters) == 0) {

        cat("\nNo clusters exceeded the threshold.\n\n")

        return (invisible(object) )

    }

    # identify the two time columns as numeric columns not in bookkeeping
    cn <- colnames(clusters)
    bookkeeping <- c("id", "sign", "value", "n_points", "participant")
    candidate_cols <- setdiff(cn, bookkeeping)

    is_num <- vapply(
        clusters[, candidate_cols, drop = FALSE],
        is.numeric,
        logical(1)
        )

    time_cols <- candidate_cols[is_num]

    if (length(time_cols) != 2L) {

        stop (
            "Could not identify exactly two numeric time columns for 2D clusters. ",
            "Found: ", paste(time_cols, collapse = ", "), "."
            )

    }

    t1 <- time_cols[1]
    t2 <- time_cols[2]

    has_participant <- "participant" %in% cn
    has_n_points <- "n_points" %in% cn
    has_value <- "value" %in% cn

    # count clusters robustly (id may repeat across participants)
    if (has_participant) {

        n_clust <- dplyr::n_distinct(
            paste(clusters$participant, clusters$id, clusters$sign, sep = "|")
            )

    } else {

        n_clust <- dplyr::n_distinct(
            paste(clusters$id, clusters$sign, sep = "|")
            )

    }

    cat("Clusters found: ", n_clust, "\n", sep = "")

    # summarise pointwise 2D clusters into one row per cluster (and participant if present)
    cl_sum <- clusters |>
        dplyr::group_by(
          .data$id,
          .data$sign,
          dplyr::across(dplyr::any_of("participant") )
          ) |>
        dplyr::summarise(
          onset_t1 = min(.data[[t1]], na.rm = TRUE),
          offset_t1 = max(.data[[t1]], na.rm = TRUE),
          onset_t2 = min(.data[[t2]], na.rm = TRUE),
          offset_t2 = max(.data[[t2]], na.rm = TRUE),
          n_points = if (has_n_points) {
              dplyr::first(stats::na.omit(.data$n_points) )
              } else {
                  dplyr::n()
              },
          value_min = if (has_value) min(.data$value, na.rm = TRUE) else NA_real_,
          value_max = if (has_value) max(.data$value, na.rm = TRUE) else NA_real_,
          .groups = "drop"
          ) |>
        dplyr::mutate(
            span_t1 = .data$offset_t1 - .data$onset_t1,
            span_t2 = .data$offset_t2 - .data$onset_t2
            )

    # if n_points exists but was all-NA for a cluster, fall back to row count
    if (has_n_points) {

        cl_sum <- cl_sum |>
            dplyr::mutate(
                n_points = dplyr::if_else(
                    is.na(.data$n_points),
                    as.numeric(NA), # placeholder; will be replaced next line
                    .data$n_points
                    )
                )

        # replace remaining NAs with counts by id (+ participant if present)
        if (anyNA(cl_sum$n_points) ) {

            counts <- clusters |>
                dplyr::group_by(
                    .data$id, .data$sign,
                    dplyr::across(dplyr::any_of("participant") )
                    ) |>
                dplyr::summarise(n_points_fallback = dplyr::n(), .groups = "drop")

          cl_sum <- cl_sum |>
              dplyr::left_join(counts, by = intersect(names(cl_sum), names(counts) ) ) |>
              dplyr::mutate(
                  n_points = dplyr::if_else(
                      is.na(.data$n_points),
                      as.numeric(.data$n_points_fallback),
                      as.numeric(.data$n_points)
                      )
                  ) |>
              dplyr::select(-.data$n_points_fallback)

        }

    }

    # cluster statistics (2D analogue of 1D durations)
    cat("\nCluster statistics:\n")
    cat("  Mean span (", t1, "): ", round(mean(cl_sum$span_t1), digits), "\n", sep = "")
    cat("  Median span (", t1, "): ", round(stats::median(cl_sum$span_t1), digits), "\n", sep = "")
    cat("  Min span (", t1, "): ", round(min(cl_sum$span_t1), digits), "\n", sep = "")
    cat("  Max span (", t1, "): ", round(max(cl_sum$span_t1), digits), "\n", sep = "")

    cat("  Mean span (", t2, "): ", round(mean(cl_sum$span_t2), digits), "\n", sep = "")
    cat("  Median span (", t2, "): ", round(stats::median(cl_sum$span_t2), digits), "\n", sep = "")
    cat("  Min span (", t2, "): ", round(min(cl_sum$span_t2), digits), "\n", sep = "")
    cat("  Max span (", t2, "): ", round(max(cl_sum$span_t2), digits), "\n", sep = "")

    cat("  Mean cluster size (n_points): ", round(mean(cl_sum$n_points), digits), "\n", sep = "")
    cat("  Median cluster size (n_points): ", round(stats::median(cl_sum$n_points), digits), "\n", sep = "")
    cat("  Min cluster size (n_points): ", round(min(cl_sum$n_points), digits), "\n", sep = "")
    cat("  Max cluster size (n_points): ", round(max(cl_sum$n_points), digits), "\n", sep = "")
    cat("  Total points across clusters: ", sum(cl_sum$n_points), "\n", sep = "")

    if (has_value) {

        cat("  Value range across clusters: [",
            round(min(cl_sum$value_min, na.rm = TRUE), digits),
            ", ",
            round(max(cl_sum$value_max, na.rm = TRUE), digits),
            "]\n",
            sep = ""
            )

    }

    # prepare rounded cluster table (with original time names in column headers)
    cl_print <- cl_sum |>
        dplyr::mutate(
            dplyr::across(dplyr::where(is.numeric), ~ round(.x, digits) )
            )

    # order columns (parallel to print.clusters_results_2d)
    if (has_participant) {

        cl_print <- cl_print |>
            dplyr::select(
                .data$participant, .data$sign, .data$id,
                dplyr::starts_with("onset_"),
                dplyr::starts_with("offset_"),
                dplyr::starts_with("span_"),
                .data$n_points,
                dplyr::any_of(c("value_min", "value_max") )
                )

    } else {

        cl_print <- cl_print |>
            dplyr::select(
                .data$sign, .data$id,
                dplyr::starts_with("onset_"),
                dplyr::starts_with("offset_"),
                dplyr::starts_with("span_"),
                .data$n_points,
                dplyr::any_of(c("value_min", "value_max") )
                )

    }

    cat("\nCluster table:\n\n")
    print(data.frame(cl_print), row.names = FALSE)
    cat("\n=================================================================\n")

    return (invisible(object) )

}
