#' Identify a neurogam analysis result
#'
#' @description
#' Internal type detection shared by the read-only result accessors. The
#' existing concrete S3 classes remain authoritative for dispatch.
#'
#' @param x A neurogam result object.
#'
#' @return One of `"temporal"`, `"temporal_generalization"`, or
#'   `"spatiotemporal"`.
#' @keywords internal
neurogam_result_type <- function(x) {
    if (inherits(x, "clusters_results_st")) {
        return("spatiotemporal")
    }

    if (inherits(x, "clusters_results_2d")) {
        return("temporal_generalization")
    }

    if (inherits(x, "clusters_results_1d")) {
        return("temporal")
    }

    stop(
        "`x` must inherit from `clusters_results_1d`, ",
        "`clusters_results_2d`, or `clusters_results_st`.",
        call. = FALSE
    )
}

#' Validate a neurogam analysis result
#'
#' @description
#' Validate the common result contract without requiring a fitted Stan model.
#' This checks structure, not the statistical correctness of the estimates.
#' Existing concrete S3 classes and legacy fields are not changed.
#'
#' @param x A result returned by `testing_through_time()` or
#'   `analyse_bgam_spatiotemporal()`.
#' @param check_membership If `TRUE`, also verify that cluster membership
#'   counts agree with `cluster_summary$n_points`. Cluster identity is defined
#'   by sign and ID, and additionally by participant when applicable.
#'
#' @return Invisibly returns `TRUE` or raises an informative error.
#' @export
validate_neurogam_result <- function(x, check_membership = FALSE) {
    neurogam_result_type(x)

    required <- c(
        "predictions", "cluster_summary", "cluster_membership", "model"
    )
    missing <- setdiff(required, names(x))
    if (length(missing)) {
        stop(
            "Result lacks required fields: ",
            paste(missing, collapse = ", "),
            call. = FALSE
        )
    }

    tables <- c("predictions", "cluster_summary", "cluster_membership")
    invalid <- tables[!vapply(x[tables], is.data.frame, logical(1))]
    if (length(invalid)) {
        stop(
            "These result fields must be data frames: ",
            paste(invalid, collapse = ", "),
            call. = FALSE
        )
    }

    if (!"prob_ratio" %in% names(x$predictions)) {
        stop("`predictions` must contain `prob_ratio`.", call. = FALSE)
    }

    for (nm in c("cluster_summary", "cluster_membership")) {
        if (nrow(x[[nm]]) > 0L &&
            !all(c("id", "sign") %in% names(x[[nm]]))) {
            stop(
                "Nonempty `", nm, "` must have `id` and `sign` columns.",
                call. = FALSE
            )
        }
    }

    if (!is.logical(check_membership) || length(check_membership) != 1L ||
        is.na(check_membership)) {
        stop("`check_membership` must be TRUE or FALSE.", call. = FALSE)
    }

    if (check_membership && nrow(x$cluster_summary) > 0L) {
        summary <- x$cluster_summary
        membership <- x$cluster_membership

        if (!"n_points" %in% names(summary)) {
            stop("`cluster_summary` must contain `n_points`.", call. = FALSE)
        }
        if (nrow(membership) == 0L) {
            stop("Nonempty cluster summaries lack memberships.", call. = FALSE)
        }

        group_cols <- c(
            intersect("participant", intersect(names(summary), names(membership))),
            "id", "sign"
        )
        size <- stats::aggregate(
            rep.int(1L, nrow(membership)),
            by = membership[group_cols],
            FUN = length
        )
        names(size)[ncol(size)] <- ".count"
        joined <- merge(summary, size, by = group_cols, all.x = TRUE)

        if (nrow(joined) != nrow(summary) ||
            anyNA(joined$.count) ||
            any(joined$.count != joined$n_points)) {
            stop(
                "Cluster membership counts disagree with `cluster_summary`.",
                call. = FALSE
            )
        }
    }

    invisible(TRUE)
}

#' Access standardized neurogam result metadata
#'
#' @description
#' Derives consistent metadata for temporal (1D), cross-temporal (2D), and
#' sensor-by-time results without mutating the original result. Metadata absent
#' from older results are kept `NULL` rather than guessed. In particular, the
#' posterior-odds convention is not inferred from the model or result class.
#'
#' @inheritParams validate_neurogam_result
#'
#' @return A named list with analysis type, coordinate names, and any recorded
#'   inference or adjacency settings.
#' @export
neurogam_result_metadata <- function(x) {
    validate_neurogam_result(x)

    analysis_type <- neurogam_result_type(x)
    columns <- if (is.list(x$input_columns)) x$input_columns else list()

    time_id <- if (identical(analysis_type, "spatiotemporal")) {
        if (!is.null(columns$time_id)) columns$time_id else "time"
    } else if (identical(analysis_type, "temporal_generalization")) {
        c("time1", "time2")
    } else {
        "time"
    }

    sensor_id <- if (identical(analysis_type, "spatiotemporal")) {
        if (!is.null(columns$sensor_id)) {
            columns$sensor_id
        } else if ("sensor" %in% names(x$predictions)) {
            "sensor"
        } else if ("channel" %in% names(x$predictions)) {
            "channel"
        } else {
            NULL
        }
    } else {
        NULL
    }

    space_id <- if (identical(analysis_type, "spatiotemporal")) {
        if (!is.null(columns$space_id)) columns$space_id else {
            intersect(c("xproj", "yproj"), names(x$predictions))
        }
    } else {
        NULL
    }

    model_coordinates <- if (identical(analysis_type, "spatiotemporal")) {
        c(space_id, time_id)
    } else {
        time_id
    }

    cluster_coordinates <- if (identical(analysis_type, "spatiotemporal")) {
        c(sensor_id, time_id)
    } else {
        time_id
    }

    metadata <- list(
        analysis_type = analysis_type,
        model_coordinates = model_coordinates,
        cluster_coordinates = cluster_coordinates,
        model_ndim = if (identical(analysis_type, "spatiotemporal")) 3L else length(time_id),
        sensor_id = sensor_id,
        time_id = time_id,
        space_id = space_id,
        n_prediction_locations = nrow(x$predictions),
        n_clusters = nrow(x$cluster_summary),
        n_cluster_members = nrow(x$cluster_membership),
        multilevel = x$multilevel,
        threshold = x$threshold,
        threshold_type = x$threshold_type,
        null_value = x$null_value,
        rope = x$rope,
        odds_mode = x$odds_mode,
        contrast = x$contrast,
        contrast_levels = x$contrast_levels,
        adjacency_source = x$adjacency_source
    )

    # P2.2: use recorded settings when present, preserving the P2.1 fallback
    # for historical objects that have no attached metadata.
    recorded <- attr(x, "neurogam_metadata", exact = TRUE)
    if (!is.null(recorded)) {
        if (!is.list(recorded)) {
            stop("Stored neurogam metadata must be a list.", call. = FALSE)
        }
        for (name in names(recorded)) {
            metadata[name] <- list(recorded[[name]])
        }
    }

    # Counts describe the object as it exists now; do not cache them as
    # construction-time metadata that may become stale after subsetting.
    metadata$n_prediction_locations <- nrow(x$predictions)
    metadata$n_clusters <- nrow(x$cluster_summary)
    metadata$n_cluster_members <- nrow(x$cluster_membership)

    metadata
}

#' Create a standardized read-only view of neurogam results
#'
#' @description
#' Creates a common interface to predictions and detected clusters without
#' changing the input object's class, list fields, or legacy `clusters` table.
#' Use this as an intermediate compatibility layer before introducing unified
#' S3 `print()`, `summary()`, and `plot()` methods.
#'
#' @inheritParams validate_neurogam_result
#'
#' @return A list of class `neurogam_result_view` containing `predictions`,
#'   `cluster_summary`, `cluster_membership`, `model`, `summary_data`, and
#'   `metadata`. It holds the existing objects, not recomputed values.
#' @export
as_neurogam_result_view <- function(x, check_membership = FALSE) {
    validate_neurogam_result(x, check_membership = check_membership)

    result <- list(
        predictions = x$predictions,
        cluster_summary = x$cluster_summary,
        cluster_membership = x$cluster_membership,
        model = x$model,
        summary_data = x$summary_data,
        metadata = neurogam_result_metadata(x)
    )

    class(result) <- c("neurogam_result_view", "list")
    result
}
