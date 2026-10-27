#' Prepare a neurogam result for model-independent plotting
#'
#' @description
#' Internal adapter. A saved result may lack its model field, which is not
#' needed for plots built from saved prediction/cluster tables. P2.1 validation
#' remains strict; any repair is applied to this function's local copy only.
#' @keywords internal
neurogam_plot_view <- function(x) {
    if (!"model" %in% names(x)) {
        x["model"] <- list(NULL)
    }
    as_neurogam_result_view(x)
}

#' Check a typed neurogam plotting request
#' @keywords internal
neurogam_plot_type <- function(type, choices) {
    if (!is.character(type) || length(type) != 1L ||
        is.na(type) || !nzchar(type)) {
        stop("`type` must be a single nonempty string.", call. = FALSE)
    }
    if (!type %in% choices) {
        stop(
            "Unsupported `type`: ", type, ". Available types: ",
            paste(choices, collapse = ", "), ".",
            call. = FALSE
        )
    }
    type
}

#' Reject plotting arguments that cannot be forwarded to a typed plot
#' @keywords internal
neurogam_check_plot_dots <- function(...) {
    dots <- list(...)
    if (length(dots)) {
        names_dots <- names(dots)
        if (is.null(names_dots)) {
            names_dots <- rep("", length(dots))
        }
        names_dots[!nzchar(names_dots)] <- "<unnamed>"
        stop(
            "Unused plotting argument(s) for this `type`: ",
            paste(names_dots, collapse = ", "), ".",
            call. = FALSE
        )
    }
    invisible(NULL)
}

#' Choose the ordinate stored in a neurogam prediction table
#' @keywords internal
neurogam_prediction_value <- function(predictions) {
    options <- c("post_prob", "estimate")
    available <- intersect(options, names(predictions))
    if (!length(available)) {
        stop(
            "`predictions` must include `post_prob` or `estimate` ",
            "for a predictions plot.",
            call. = FALSE
        )
    }
    available[1L]
}

#' Explicit plotting options for temporal results
#'
#' @description
#' Called only when the user supplies `type` to `plot()` for a 1D result.
#' The historical no-`type` plot remains in R/plot.R, unchanged.
#' @keywords internal
neurogam_plot_typed_1d <- function(
    x,
    type,
    null_value = 0,
    clusters_colour = "black",
    theme = theme_neurogam(),
    ...
) {
    type <- neurogam_plot_type(type, c("predictions", "clusters"))
    neurogam_check_plot_dots(...)
    # if (!ggplot2::is.theme(theme)) {
    #     stop("`theme` must be a ggplot2 theme.", call. = FALSE)
    # }
    if (!ggplot2::is_theme(theme)) {
        stop("`theme` must be a ggplot2 theme.", call. = FALSE)
    }
    view <- neurogam_plot_view(x)
    time_id <- view$metadata$time_id
    if (!is.character(time_id) || length(time_id) != 1L) {
        stop("Expected one temporal coordinate in result metadata.", call. = FALSE)
    }

    if (identical(type, "predictions")) {
        predictions <- view$predictions
        value_id <- neurogam_prediction_value(predictions)
        if (!time_id %in% names(predictions)) {
            stop("Missing temporal coordinate in `predictions`.", call. = FALSE)
        }
        p <- ggplot2::ggplot(
            predictions,
            ggplot2::aes(x = .data[[time_id]], y = .data[[value_id]])
        ) + ggplot2::geom_hline(
            yintercept = null_value, linetype = 2
        )
        if (all(c("lower", "upper") %in% names(predictions))) {
            p <- p + ggplot2::geom_ribbon(
                ggplot2::aes(
                    ymin = .data$lower, ymax = .data$upper, y = NULL
                ),
                fill = clusters_colour, alpha = 0.18
            )
        }
        p <- p + ggplot2::geom_line(colour = clusters_colour, linewidth = 0.8)
        if ("participant" %in% names(predictions)) {
            p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$participant))
        }
        return(p + theme + ggplot2::labs(
            x = time_id, y = "Posterior estimate"
        ))
    }

    clusters <- view$cluster_summary
    if (!all(c("onset", "offset") %in% names(clusters))) {
        if (nrow(clusters) > 0L) {
            stop(
                "1D `cluster_summary` must have `onset` and `offset`.",
                call. = FALSE
            )
        }
        clusters$onset <- numeric(0)
        clusters$offset <- numeric(0)
    }
    if (!"sign" %in% names(clusters)) {
        clusters$sign <- character(nrow(clusters))
    }
    p <- ggplot2::ggplot(clusters) +
        ggplot2::geom_segment(
            ggplot2::aes(
                x = .data$onset, xend = .data$offset,
                y = .data$sign, yend = .data$sign,
                colour = .data$sign
            ),
            linewidth = 3, lineend = "butt", show.legend = FALSE
        ) +
        ggplot2::labs(x = time_id, y = "Cluster direction") + theme
    if ("participant" %in% names(clusters)) {
        p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$participant))
    }
    p
}

#' Explicit plotting options for temporal-generalisation results
#' @keywords internal
neurogam_plot_typed_2d <- function(
    x,
    type,
    theme = theme_neurogam(),
    palette = "vik",
    midpoint = 0.5,
    axes_labels = c("Training time (s)", "Testing time (s)"),
    ...
) {
    type <- neurogam_plot_type(type, c("predictions", "clusters"))
    neurogam_check_plot_dots(...)
    # if (!ggplot2::is.theme(theme)) {
    #     stop("`theme` must be a ggplot2 theme.", call. = FALSE)
    # }
    if (!ggplot2::is_theme(theme)) {
        stop("`theme` must be a ggplot2 theme.", call. = FALSE)
    }
    if (!palette %in% scico::scico_palette_names()) {
        stop("`palette` is not a scico palette.", call. = FALSE)
    }
    if (!is.numeric(midpoint) || length(midpoint) != 1L ||
        !is.finite(midpoint)) {
        stop("`midpoint` must be finite and numeric.", call. = FALSE)
    }
    if (!is.character(axes_labels) || length(axes_labels) != 2L ||
        anyNA(axes_labels)) {
        stop("`axes_labels` must contain two labels.", call. = FALSE)
    }

    view <- neurogam_plot_view(x)
    coordinates <- view$metadata$time_id
    if (length(coordinates) != 2L) {
        stop("Expected two temporal coordinates in result metadata.", call. = FALSE)
    }
    field <- if (identical(type, "predictions")) {
        view$predictions
    } else {
        view$cluster_membership
    }
    if (!all(coordinates %in% names(field))) {
        stop(
            "Result table lacks the recorded temporal coordinates.",
            call. = FALSE
        )
    }

    p <- ggplot2::ggplot(
        field, ggplot2::aes(
            x = .data[[coordinates[1L]]],
            y = .data[[coordinates[2L]]]
        )
    )
    if (identical(type, "predictions")) {
        value_id <- neurogam_prediction_value(field)
        colours <- scico::scico(3L, palette = palette)
        p <- p + ggplot2::geom_tile(
            ggplot2::aes(fill = .data[[value_id]])
        ) + ggplot2::scale_fill_gradient2(
            low = colours[1L], mid = colours[2L], high = colours[3L],
            midpoint = midpoint, name = "Posterior estimate"
        )
    } else {
        if (!"sign" %in% names(field)) {
            field$sign <- character(nrow(field))
            p <- ggplot2::ggplot(
                field, ggplot2::aes(
                    x = .data[[coordinates[1L]]],
                    y = .data[[coordinates[2L]]]
                )
            )
        }
        p <- p + ggplot2::geom_tile(
            ggplot2::aes(fill = .data$sign)
        ) + ggplot2::scale_fill_manual(
            values = c(positive = "#BC5241", negative = "#376DA4"),
            na.value = "grey80", name = "Cluster direction"
        )
    }
    if ("participant" %in% names(field)) {
        p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$participant))
    }
    p + ggplot2::coord_equal() + theme + ggplot2::labs(
        x = axes_labels[1L], y = axes_labels[2L]
    )
}

#' Normalize a requested spatiotemporal plot
#' @keywords internal
neurogam_st_plot_request <- function(type, dots = list()) {
    type <- neurogam_plot_type(
        type, c("topography", "trajectories", "predictions", "clusters")
    )
    if (!is.list(dots)) {
        stop("`dots` must be a list.", call. = FALSE)
    }
    if (length(dots)) {
        if (is.null(names(dots)) || anyNA(names(dots)) ||
            any(!nzchar(names(dots))) || anyDuplicated(names(dots))) {
            stop("All additional spatiotemporal plotting arguments must have unique names.",
                 call. = FALSE)
        }
        reserved <- intersect(names(dots), c("results", "observed_data", "observed_data_b"))
        if (length(reserved)) {
            stop("Do not repeat `results`, `observed_data`, or `observed_data_b` in `...`.",
                 call. = FALSE)
        }
    }
    if (identical(type, "predictions")) {
        if (!is.null(dots$cluster_style) &&
            !identical(dots$cluster_style, "none")) {
            stop(
                "`type = 'predictions'` requires `cluster_style = 'none'`.",
                call. = FALSE
            )
        }
        dots$cluster_style <- "none"
        type <- "trajectories"
    } else if (identical(type, "clusters")) {
        if (is.null(dots$cluster_style)) {
            dots$cluster_style <- "segments"
        }
        type <- "trajectories"
    }
    list(type = type, dots = dots)
}

#' Plot spatiotemporal neurogam results with a unified interface
#'
#' @description
#' S3 method for 3D sensor-by-time results. Delegates to the existing
#' [plot_spatiotemporal_topography()] and [plot_sensor_trajectories()]
#' implementations; no predictions or clusters are recomputed.
#'
#' @param x A `clusters_results_st` object.
#' @param type One of `"topography"` (default), `"trajectories"`,
#'   `"predictions"`, or `"clusters"`. The latter two use trajectories with
#'   cluster highlighting switched off or on, respectively.
#' @param observed_data Original long-form EEG observations required by the
#'   existing spatiotemporal plotting helpers.
#' @param observed_data_b Optional observations for the second condition.
#' @param ... Named parameters passed to the selected specialized plotter.
#'
#' @return A `ggplot` object.
#' @examples
#' \dontrun{
#' plot(results_st, type = "topography", observed_data = eeg_data2)
#' plot(results_st, type = "trajectories", observed_data = eeg_data2,
#'      layout = "facets", cluster_style = "segments")
#' }
#' @export
plot.clusters_results_st <- function(
    x,
    type = "topography",
    observed_data = NULL,
    observed_data_b = NULL,
    ...
) {
    request <- neurogam_st_plot_request(type, list(...))
    if (is.null(observed_data) || !is.data.frame(observed_data)) {
        stop(
            "`observed_data` must be a data frame for spatiotemporal ",
            "plots. Supply the original long-form EEG observations.",
            call. = FALSE
        )
    }
    if (!is.null(observed_data_b) && !is.data.frame(observed_data_b)) {
        stop("`observed_data_b` must be a data frame.", call. = FALSE)
    }
    plotting_function <- switch(
        request$type,
        topography = plot_spatiotemporal_topography,
        trajectories = plot_sensor_trajectories
    )
    do.call(
        plotting_function,
        c(
            list(
                results = x,
                observed_data = observed_data,
                observed_data_b = observed_data_b
            ),
            request$dots
        )
    )
}
