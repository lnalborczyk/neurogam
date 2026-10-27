#' Construct metadata for a neurogam result
#'
#' @description
#' Records the analysis specification that is known when a result is created.
#' This function does not infer settings from a fitted model or posterior
#' summaries. A `NULL` field represents a setting not recorded for that result.
#'
#' @param analysis_type One of `"temporal"`, `"temporal_generalization"`, or
#'   `"spatiotemporal"`.
#' @param time_id Names of time coordinates in the resulting predictions.
#' @param space_id Names of projected spatial coordinates, for spatiotemporal
#'   results only.
#' @param sensor_id Name of the sensor identifier column, for spatiotemporal
#'   results only.
#' @param participant_id,outcome_id,predictor_id Original input column names.
#' @param multilevel Model aggregation mode.
#' @param threshold,threshold_type Threshold used for cluster detection.
#' @param null_value Null reference for the effect or contrast.
#' @param rope Fixed ROPE half-width, when one was specified.
#' @param baseline_window Baseline interval defining a data-dependent ROPE,
#'   when applicable. The resulting numerical ROPE width is not stored here.
#' @param credible_interval Credible interval width.
#' @param odds_mode Posterior-odds strategy actually used by the analysis.
#' @param contrast Description of the estimated effect or difference.
#' @param contrast_levels Explicit condition labels, when available.
#' @param adjacency_source One of `"auto_knn"`, `"provided"`, or `NULL`.
#' @param participant_clusters Whether participant-level clusters were computed.
#' @param n_post_samples Number of posterior draws requested or used.
#'
#' @return A named metadata list.
#' @keywords internal
make_neurogam_metadata <- function(
    analysis_type,
    time_id,
    space_id = NULL,
    sensor_id = NULL,
    participant_id = NULL,
    outcome_id = NULL,
    predictor_id = NULL,
    multilevel = NULL,
    threshold = NULL,
    threshold_type = NULL,
    null_value = NULL,
    rope = NULL,
    baseline_window = NULL,
    credible_interval = NULL,
    odds_mode = NULL,
    contrast = NULL,
    contrast_levels = NULL,
    adjacency_source = NULL,
    participant_clusters = NULL,
    n_post_samples = NULL
) {
    analysis_type <- match.arg(
        analysis_type,
        c("temporal", "temporal_generalization", "spatiotemporal")
    )

    expected_time_dimensions <- if (identical(analysis_type, "temporal_generalization")) {
        2L
    } else {
        1L
    }

    if (!is.character(time_id) || length(time_id) != expected_time_dimensions ||
        anyNA(time_id) || any(!nzchar(time_id)) || anyDuplicated(time_id)) {
        stop("`time_id` is invalid for the selected analysis type.", call. = FALSE)
    }

    if (identical(analysis_type, "spatiotemporal")) {
        if (!is.character(space_id) || length(space_id) != 2L ||
            anyNA(space_id) || any(!nzchar(space_id))) {
            stop("Spatiotemporal metadata requires two spatial coordinates.", call. = FALSE)
        }
        if (!is.character(sensor_id) || length(sensor_id) != 1L ||
            is.na(sensor_id) || !nzchar(sensor_id)) {
            stop("Spatiotemporal metadata requires a sensor identifier.", call. = FALSE)
        }
        if (anyDuplicated(c(space_id, sensor_id, time_id))) {
            stop("Sensor, spatial, and temporal columns must be distinct.", call. = FALSE)
        }
    } else if (!is.null(space_id) || !is.null(sensor_id)) {
        stop("Spatial coordinates are only valid for spatiotemporal metadata.",
             call. = FALSE)
    }

    if (!is.null(odds_mode)) {
        odds_mode <- match.arg(
            odds_mode,
            #c("legacy_one_sided", "legacy_directional", "spatiotemporal")
            c(
                "directional_rope",
                "legacy_one_sided",
                "legacy_directional",
                "spatiotemporal"
            )
        )
    }

    if (!is.null(threshold) &&
        (!is.numeric(threshold) || length(threshold) != 1L ||
         is.na(threshold) || !is.finite(threshold) || threshold <= 1)) {
        stop("`threshold` must be a finite number greater than 1.", call. = FALSE)
    }

    if (!is.null(threshold_type)) {
        threshold_type <- match.arg(threshold_type, c("both", "above", "below"))
    }

    if (!is.null(rope) &&
        (!is.numeric(rope) || length(rope) != 1L ||
         is.na(rope) || !is.finite(rope) || rope < 0)) {
        stop("`rope` must be a nonnegative finite number.", call. = FALSE)
    }

    list(
        metadata_version = 1L,
        analysis_type = analysis_type,
        model_coordinates = c(space_id, time_id),
        cluster_coordinates = if (identical(analysis_type, "spatiotemporal")) {
            c(sensor_id, time_id)
        } else {
            time_id
        },
        model_ndim = as.integer(length(c(space_id, time_id))),
        sensor_id = sensor_id,
        time_id = time_id,
        space_id = space_id,
        participant_id = participant_id,
        outcome_id = outcome_id,
        predictor_id = predictor_id,
        multilevel = multilevel,
        threshold = threshold,
        threshold_type = threshold_type,
        null_value = null_value,
        rope = rope,
        baseline_window = baseline_window,
        credible_interval = credible_interval,
        odds_mode = odds_mode,
        contrast = contrast,
        contrast_levels = contrast_levels,
        adjacency_source = adjacency_source,
        participant_clusters = participant_clusters,
        n_post_samples = n_post_samples
    )
}

#' Attach analysis metadata without altering result fields or S3 classes
#'
#' @description
#' Stores versioned metadata as the `neurogam_metadata` attribute. This keeps
#' legacy list field names and classes unchanged and survives `saveRDS()`.
#' Existing metadata are replaced by fields supplied in the new metadata list.
#'
#' @param x A recognized neurogam result.
#' @param metadata Named metadata returned by `make_neurogam_metadata()`.
#'
#' @return The same result with a metadata attribute.
#' @keywords internal
set_neurogam_result_metadata <- function(x, metadata) {
    analysis_type <- neurogam_result_type(x)

    if (!is.list(metadata) || is.null(names(metadata)) ||
        anyNA(names(metadata)) || any(!nzchar(names(metadata))) ||
        anyDuplicated(names(metadata))) {
        stop("`metadata` must be a uniquely named list.", call. = FALSE)
    }

    recorded_type <- metadata$analysis_type
    if (!is.null(recorded_type) && !identical(recorded_type, analysis_type)) {
        stop("Metadata analysis type conflicts with the result class.",
             call. = FALSE)
    }

    existing <- attr(x, "neurogam_metadata", exact = TRUE)
    if (is.null(existing)) {
        existing <- list()
    }
    if (!is.list(existing)) {
        stop("Existing result metadata are malformed.", call. = FALSE)
    }

    for (name in names(metadata)) {
        # Use single-bracket assignment to retain explicit NULL values.
        existing[name] <- list(metadata[[name]])
    }
    existing$analysis_type <- analysis_type
    attr(x, "neurogam_metadata") <- existing
    x
}
