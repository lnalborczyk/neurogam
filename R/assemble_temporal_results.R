#' Detect temporal clusters from posterior summaries
#'
#' @description
#' Internal compatibility boundary between posterior inference and clustering.
#' Preserves the historical 1D/2D cluster output alongside the standardized
#' cluster summary and membership tables.
#'
#' @param predictions Data frame containing coordinates and `prob_ratio`.
#' @param is_2d_time Whether predictions use `time1` and `time2` instead of `time`.
#' @param participant_clusters Whether predictions are grouped by participant.
#' @param threshold Posterior-odds threshold.
#' @param threshold_type One of `"both"`, `"above"`, or `"below"`.
#'
#' @return A list containing `clusters`, `cluster_summary`, and
#'   `cluster_membership`.
#' @keywords internal
assemble_temporal_clusters <- function (
    predictions,
    is_2d_time,
    participant_clusters,
    threshold,
    threshold_type,
    odds_mode = "directional_rope"
    ) {

    coordinates <- if (is_2d_time) c("time1", "time2") else "time"
    cluster_group <- if (participant_clusters) "participant" else NULL
    cluster_columns <- c(coordinates, cluster_group)

    if (!all(c(cluster_columns, "prob_ratio") %in% names(predictions) ) ) {

        stop (
            "Posterior summaries lack the required clustering columns.",
            call. = FALSE
            )

    }

    # cluster_data <- as.data.frame(
    #     predictions[, c(cluster_columns, "prob_ratio"), drop = FALSE]
    #     )
    #
    # names(cluster_data)[names(cluster_data) == "prob_ratio"] <- "value"
    #
    # if (identical(odds_mode, "directional_rope")) {
    #
    #     cluster_data$value <- directional_rope_cluster_score(
    #         odds_above = predictions$odds_above,
    #         odds_below = predictions$odds_below,
    #         threshold = threshold
    #     )
    #
    # }

    cluster_data <- as.data.frame(
        predictions[, cluster_columns, drop = FALSE]
    )

    if (identical(odds_mode, "directional_rope")) {

        cluster_data$value <- directional_rope_cluster_score(
            odds_above = predictions$odds_above,
            odds_below = predictions$odds_below,
            threshold = threshold
        )

    } else {

        cluster_data$value <- predictions$prob_ratio

    }

    cluster_data <- cluster_data[
        stats::complete.cases(cluster_data),
        ,
        drop = FALSE
        ]

    if (!nrow(cluster_data) ) {

        stop (
            "No complete posterior-odds observations available for clustering.",
            call. = FALSE
            )

    }

    cluster_nd <- find_clusters_nd(
        data = cluster_data,
        coordinates = coordinates,
        value = "value",
        group = cluster_group,
        threshold = threshold,
        threshold_type = threshold_type
        )

    return (
        list(
            clusters = cluster_nd_to_legacy(
                result = cluster_nd,
                time_id = coordinates,
                group = cluster_group
                ),
            cluster_summary = cluster_nd$clusters,
            cluster_membership = cluster_nd$membership
            )
        )

}

#' Construct legacy-compatible temporal results
#'
#' @description
#' Internal constructor. Maintains the existing field order and S3 class names.
#'
#' @param clustering Output from `assemble_temporal_clusters()`.
#' @param predictions Posterior summary data frame.
#' @param model Fitted brms model or test fixture.
#' @param summary_data Data used for fitting.
#' @param multilevel Multilevel mode.
#' @param is_2d_time Whether this is a 2D temporal analysis.
#'
#' @return A `clusters_results_1d` or `clusters_results_2d` object.
#' @keywords internal
assemble_temporal_results <- function (
    clustering,
    predictions,
    model,
    summary_data,
    multilevel,
    is_2d_time,
    metadata = NULL
) {

    results <- list(
        clusters = clustering$clusters,
        cluster_summary = clustering$cluster_summary,
        cluster_membership = clustering$cluster_membership,
        predictions = predictions,
        model = model,
        summary_data = summary_data,
        multilevel = multilevel
        )

    class(results) <- if (is_2d_time) {

        "clusters_results_2d"

    } else {

        "clusters_results_1d"

    }

    if (!is.null(metadata)) {
        results <- set_neurogam_result_metadata(results, metadata)
    }

    return (results)

}
