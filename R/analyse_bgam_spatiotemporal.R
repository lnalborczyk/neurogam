#' Summarise posterior draws at arbitrary locations
#'
#' Summarises a matrix of paired posterior contrasts or one-sample
#' expectations (draws in rows, observation locations in columns).
#' Posterior odds are directional and can use a symmetric ROPE.
#'
#' @param draws Numeric matrix, posterior draws by locations.
#' @param null_value Null value against which draws are assessed.
#' @param rope Non-negative half-width of a ROPE around `null_value`.
#' @param credible_interval Width of equal-tailed credible intervals.
#'
#' @return Data frame with posterior estimates, probabilities, and odds.
#'
#' @author Ladislas Nalborczyk \email{ladislas.nalborczyk@@cnrs.fr}.
#'
#' @keywords internal
summarise_spatiotemporal_draws <- function(
    draws,
    null_value = 0,
    rope = 0,
    credible_interval = 0.95,
    odds_mode = "legacy_one_sided"
) {
    summarise_posterior_matrix(
        draws = draws,
        null_value = null_value,
        rope = rope,
        credible_interval = credible_interval,
        odds_mode = odds_mode
    )
}

#' Identify posterior-odds clusters in sensor-by-time BGAMM predictions
#'
#' Uses an existing univariate `brmsfit` object and explicit prediction grids.
#' When sensor adjacency is omitted, an approximate adjacency graph is
#' inferred from projected sensor coordinates. An explicitly supplied physical
#' sensor topology is preferable for confirmatory analyses.
#'
#' @param model Fitted univariate `brmsfit` model.
#' @param newdata Data frame containing exactly one row per sensor and time.
#'   Additional model predictor columns, including grouping factors, must
#'   be present as required by `brms::posterior_epred()`.
#' @param sensor_adjacency Named adjacency matrix, two-column sensor edge list,
#'   or `NULL` (default) to infer adjacency from projected sensor positions.
#' @param newdata_b Optional second condition grid with exactly the same rows
#'   and sensor/time coordinates as `newdata`. Inference then concerns A - B.
#' @param sensor_id Name of the electrode identifier column.
#' @param time_id Name of the time coordinate column.
#' @param spatial_id Names of two projected coordinate columns (for plotting).
#' @param null_value Null value for the expected outcome or contrast.
#' @param rope Nonnegative half-width of a fixed ROPE around `null_value`.
#' @param threshold Posterior odds threshold, strictly greater than one.
#' @param threshold_type Whether to detect above/below/both directions.
#' @param credible_interval Width of equal-tailed posterior intervals.
#' @param draw_ids Optional posterior draw indices. When omitted, use all draws.
#'   Explicitly supplying indices ensures pairing across conditions and chunks.
#' @param chunk_size Number of sensor-time locations processed at once.
#' @param re_formula `NA` for population-level predictions (default), `NULL`
#'   to include fitted group-level deviations.
#' @param max_gap Optional `list(time = ...)` passed to `find_clusters_nd()`.
#' @param adjacency_k Number of nearest neighbors considered when adjacency
#'   is generated automatically.
#' @param adjacency_distance_factor Maximum auto-generated edge length
#'   relative to the median nearest-neighbor distance. Set to `NULL` to
#'   disable distance pruning.
#'
#' @return List with `predictions`, `cluster_summary`, `cluster_membership`,
#'   `model`, `sensor_adjacency`, `adjacency_source`, and metadata.
#' @seealso [make_sensor_adjacency()], [plot_sensor_adjacency()]
#'
#' @author Ladislas Nalborczyk \email{ladislas.nalborczyk@@cnrs.fr}.
#'
#' @export
analyse_bgam_spatiotemporal <- function(
    model,
    newdata,
    sensor_adjacency = NULL,
    newdata_b = NULL,
    sensor_id = "sensor",
    time_id = "time",
    spatial_id = c("xproj", "yproj"),
    null_value = 0,
    rope = 0,
    threshold = 10,
    threshold_type = c("both", "above", "below"),
    credible_interval = 0.95,
    draw_ids = NULL,
    chunk_size = 256L,
    re_formula = NA,
    max_gap = NULL,
    adjacency_k = 4L,
    adjacency_distance_factor = 1.8,
    odds_mode = "directional_rope"
    ) {

    threshold_type <- match.arg(threshold_type)

    if (!inherits(model, "brmsfit") ) {

        stop ("`model` must be a fitted `brmsfit` object.")

    }

    if (!is.data.frame(newdata) || nrow(newdata) == 0L) {

        stop ("`newdata` must be a nonempty data frame.")

    }

    if (length(sensor_id) != 1L || length(time_id) != 1L ||
        length(spatial_id) != 2L ||
        anyDuplicated(c(sensor_id, time_id, spatial_id) ) ) {

        stop ("Provide four distinct sensor, time, and spatial column names.")

    }

    position_columns <- c(sensor_id, time_id, spatial_id)

    if (!all(position_columns %in% names(newdata) ) ) {

        stop ("`newdata` must contain: ", paste(position_columns, collapse = ", "))

    }

    if (anyNA(newdata[position_columns]) ) {

        stop ("Sensor and coordinate columns must not have missing values.")

    }

    if (!is.numeric(newdata[[time_id]]) ||
        any(!is.finite(newdata[[time_id]]) ) ||
        !all(vapply(newdata[spatial_id], is.numeric, logical(1) ) ) ||
        any(!is.finite(as.matrix(newdata[spatial_id]) ) ) ) {

        stop ("Time and projected spatial coordinates must be numeric.")

    }

    if (anyDuplicated(newdata[c(sensor_id, time_id)]) ) {

        stop ("Each sensor/time combination must occur exactly once.")

    }

    layout <- unique(newdata[c(sensor_id, spatial_id)])

    if (anyDuplicated(layout[sensor_id]) ) {

        stop ("Each sensor must have only one pair of spatial coordinates.")

    }

    # Derive adjacency once before any expensive posterior prediction calls
    adjacency_source <- if (is.null(sensor_adjacency) ) {

        "auto_knn"

    } else {

        "provided"

    }

    if (is.null(sensor_adjacency) ) {

        sensor_adjacency <- make_sensor_adjacency(
            data = layout,
            sensor_id = sensor_id,
            spatial_id = spatial_id,
            k = adjacency_k,
            distance_factor = adjacency_distance_factor
            )

    }

    if (!is.null(newdata_b) ) {

        if (!is.data.frame(newdata_b) || nrow(newdata_b) != nrow(newdata) ||
            !all(position_columns %in% names(newdata_b) ) ) {

            stop ("`newdata_b` must match the rows and spatial columns in `newdata`.")

        }

        matched <- vapply(position_columns, function (column) {
            isTRUE(all.equal(
                newdata[[column]],
                newdata_b[[column]],
                check.attributes = FALSE
                ) )
        }, logical(1) )
        if (!all(matched) ) {

            stop ("`newdata` and `newdata_b` must have identical, ordered locations.")

        }

    }

    if (!is.numeric(chunk_size) || length(chunk_size) != 1L ||
        !is.finite(chunk_size) || chunk_size < 1 ||
        chunk_size != as.integer(chunk_size) ) {

        stop ("`chunk_size` must be a positive integer.")

    }

    if (!is.null(draw_ids) &&
        (!is.numeric(draw_ids) || !length(draw_ids) ||
         anyNA(draw_ids) || any(draw_ids < 1) ||
         any(draw_ids != as.integer(draw_ids) ) || anyDuplicated(draw_ids) ) ) {

        stop ("`draw_ids` must be distinct positive integer indices.")

    }

    summaries <- chunked_epred(
        model = model,
        newdata = newdata,
        newdata_b = newdata_b,
        draw_ids = draw_ids,
        chunk_size = chunk_size,
        re_formula = re_formula,
        summarise_fun = function(draws) {
            summarise_spatiotemporal_draws(
                draws,
                null_value = null_value,
                rope = rope,
                credible_interval = credible_interval
            )
        }
    )

    predictions <- cbind(
        newdata[position_columns],
        do.call(rbind, summaries)
        )

    rownames(predictions) <- NULL

    # Clustering uses either the supplied or automatically inferred topology.
    cluster_data <- predictions[c(sensor_id, time_id, "prob_ratio")]
    names(cluster_data)[3L] <- "value"

    if (identical(odds_mode, "directional_rope")) {

        cluster_data$value <- directional_rope_cluster_score(
            odds_above = predictions$odds_above,
            odds_below = predictions$odds_below,
            threshold = threshold
        )

    }

    components <- find_clusters_nd(
        data = cluster_data,
        coordinates = time_id,
        value = "value",
        sensor_id = sensor_id,
        sensor_adjacency = sensor_adjacency,
        threshold = threshold,
        threshold_type = threshold_type,
        max_gap = max_gap
        )

    result <- list(
        predictions = predictions,
        cluster_summary = components$clusters,
        cluster_membership = components$membership,
        model = model,
        null_value = null_value,
        rope = rope,
        threshold = threshold,
        threshold_type = threshold_type,
        contrast = if (is.null(newdata_b) ) "one_sample" else "A_minus_B",
        sensor_adjacency = sensor_adjacency,
        adjacency_source = adjacency_source,
        sensor_positions = layout
        )

    class(result) <- "clusters_results_st"

    result <- set_neurogam_result_metadata(
        result,
        make_neurogam_metadata(
            analysis_type = "spatiotemporal",
            time_id = time_id,
            space_id = spatial_id,
            sensor_id = sensor_id,
            threshold = threshold,
            threshold_type = threshold_type,
            null_value = null_value,
            rope = rope,
            credible_interval = credible_interval,
            odds_mode = "legacy_one_sided",
            contrast = if (is.null(newdata_b)) "one_sample" else "A_minus_B",
            adjacency_source = adjacency_source,
            n_post_samples = if (is.null(draw_ids)) NULL else length(draw_ids)
        )
    )

    return(result)

}
