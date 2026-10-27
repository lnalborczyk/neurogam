#' Compare observed and predicted EEG time courses across sensors
#'
#' Plots observed sensor-wise means (black) against posterior median model
#' predictions (blue), optionally with credible ribbons. Trajectories can be
#' placed in approximate scalp positions, sharing one amplitude scale, or
#' displayed as ordinary facets. Detected clusters can appear as translucent
#' rectangles or as thicker, opaque portions of the *predicted* trajectories.
#' Cluster segments respect sign, cluster ID and gaps in the time grid.
#'
#' @param results Output from [analyse_bgam_spatiotemporal()].
#' @param observed_data Original long-form observations, or condition A.
#' @param observed_data_b Observed condition B for A - B contrasts.
#' @param outcome_id Observed outcome column.
#' @param sensor_id Sensor identifier column.
#' @param time_id Numeric time column.
#' @param participant_id Optional participant ID for an equal-participant-
#'   weighted observed average. With `NULL`, observations are pooled.
#' @param spatial_id Two projected sensor coordinates.
#' @param layout `"scalp"` for mini-trajectories at electrode positions or
#'   `"facets"` for one independent panel per sensor.
#' @param cluster_style `"rectangles"` for transparent backgrounds,
#'   `"segments"` for highlighted prediction lines, or `"none"`.
#' @param cluster_alpha Opacity of cluster rectangle backgrounds.
#' @param highlighted_linewidth Width of cluster-highlighted segments.
#' @param max_gap Maximum time gap for joining cluster members. Defaults to
#'   1.5 times the median sampling interval.
#' @param show_ci Plot posterior credible interval ribbons.
#' @param ncol Number of facet columns when `layout = "facets"`.
#' @param free_y Use independent amplitude limits in faceted plots.
#' @param sensors Optional vector of sensors to display.
#' @param sensor_order Optional vector specifying the order of facet panels.
#' @param mini_width Width of each scalp mini-trace in projected units;
#'   defaults to 80 percent of median nearest-sensor spacing.
#' @param mini_height Height of each scalp mini-trace in projected units;
#'   defaults to 75 percent of median nearest-sensor spacing.
#' @param amplitude_limits Optional common amplitude range for scalp traces.
#' @param show_labels Display sensor labels in scalp layout.
#' @param head_outline Add a schematic head circle and nose in scalp layout.
#'
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' plot_sensor_trajectories(
#'     results_st, eeg_data2, outcome_id = "voltage",
#'     sensor_id = "channel", layout = "scalp",
#'     cluster_style = "segments"
#' )
#' }
#' @importFrom rlang .data
#'
#' @author Ladislas Nalborczyk \email{ladislas.nalborczyk@@cnrs.fr}.
#'
#' @export
plot_sensor_trajectories <- function(
    results,
    observed_data,
    observed_data_b = NULL,
    outcome_id = "voltage",
    sensor_id = "channel",
    time_id = "time",
    spatial_id = c("xproj", "yproj"),
    participant_id = NULL,
    layout = c("scalp", "facets"),
    cluster_style = c("rectangles", "segments", "none"),
    cluster_alpha = 0.18,
    highlighted_linewidth = 1.1,
    max_gap = NULL,
    show_ci = TRUE,
    ncol = 8L,
    free_y = FALSE,
    sensors = NULL,
    sensor_order = NULL,
    mini_width = NULL,
    mini_height = NULL,
    amplitude_limits = NULL,
    show_labels = TRUE,
    head_outline = TRUE
    ) {

    layout <- match.arg(layout)
    cluster_style <- match.arg(cluster_style)
    if (!is.numeric(cluster_alpha) || length(cluster_alpha) != 1L ||
        !is.finite(cluster_alpha) || cluster_alpha < 0 || cluster_alpha > 1) {
        stop("`cluster_alpha` must be between zero and one.", call. = FALSE)
    }
    if (!is.numeric(highlighted_linewidth) ||
        length(highlighted_linewidth) != 1L ||
        !is.finite(highlighted_linewidth) || highlighted_linewidth <= 0) {
        stop("`highlighted_linewidth` must be positive.", call. = FALSE)
    }

    dat <- prepare_spatiotemporal_plot_data(
        results, observed_data, observed_data_b,
        outcome_id, sensor_id, time_id, spatial_id, participant_id
    )
    full_head <- st_head_geometry(dat)
    if (!is.null(sensors)) {
        missing_sensors <- setdiff(as.character(sensors), unique(dat$sensor))
        if (length(missing_sensors)) {
            stop("Unknown sensor(s): ", paste(missing_sensors, collapse = ", "),
                 call. = FALSE)
        }
        dat <- dat[dat$sensor %in% as.character(sensors), , drop = FALSE]
    }
    if (!nrow(dat)) {
        stop("No sensors are available for plotting.", call. = FALSE)
    }
    settings <- st_time_settings(dat, max_gap)
    intervals <- make_cluster_intervals_per_sensor(dat, settings$max_gap)
    segments <- make_cluster_segments_per_sensor(dat, settings$max_gap)
    members <- dat[!is.na(dat$cluster_id), , drop = FALSE]
    palette <- c(positive = "#BC5241", negative = "#376DA4")

    if (layout == "facets") {
        order_sensors <- if (is.null(sensor_order)) {
            sort(unique(dat$sensor))
        } else {
            c(intersect(as.character(sensor_order), unique(dat$sensor)),
              setdiff(sort(unique(dat$sensor)), as.character(sensor_order)))
        }
        dat$sensor <- factor(dat$sensor, levels = unique(order_sensors))
        intervals$sensor <- factor(intervals$sensor, levels = unique(order_sensors))
        segments$sensor <- factor(segments$sensor, levels = unique(order_sensors))
        members$sensor <- factor(members$sensor, levels = unique(order_sensors))

        p <- ggplot2::ggplot()
        if (cluster_style == "rectangles" && nrow(intervals)) {
            p <- p + ggplot2::geom_rect(
                data = intervals,
                ggplot2::aes(
                    xmin = .data$xmin, xmax = .data$xmax,
                    ymin = -Inf, ymax = Inf, fill = .data$sign
                ),
                alpha = cluster_alpha, inherit.aes = FALSE
            )
        }
        if (show_ci) {
            p <- p + ggplot2::geom_ribbon(
                data = dat,
                ggplot2::aes(
                    x = .data$time, ymin = .data$lower,
                    ymax = .data$upper, group = .data$sensor
                ),
                fill = "#3973A6", alpha = 0.13
            )
        }
        p <- p +
            ggplot2::geom_line(
                data = dat,
                ggplot2::aes(
                    x = .data$time, y = .data$observed,
                    group = .data$sensor
                ),
                color = "#202020", linewidth = 0.45, na.rm = TRUE
            ) +
            ggplot2::geom_line(
                data = dat,
                ggplot2::aes(
                    x = .data$time, y = .data$predicted,
                    group = .data$sensor
                ),
                color = "#3275B2", linewidth = 0.58,
                linetype = "dashed"
            )
        if (cluster_style == "segments" && nrow(members)) {
            if (nrow(segments)) {
                p <- p + ggplot2::geom_segment(
                    data = segments,
                    ggplot2::aes(
                        x = .data$time, y = .data$predicted,
                        xend = .data$timeend, yend = .data$predictedend,
                        color = .data$cluster_sign
                    ),
                    linewidth = highlighted_linewidth, alpha = 1
                )
            }
            p <- p + ggplot2::geom_point(
                data = members,
                ggplot2::aes(
                    x = .data$time, y = .data$predicted,
                    color = .data$cluster_sign
                ),
                size = 0.7, alpha = 0.85
            )
        }
        if (cluster_style == "rectangles" && nrow(intervals)) {
            p <- p + ggplot2::scale_fill_manual(
                values = palette, name = "Cluster"
            )
        }
        if (cluster_style == "segments" && nrow(members)) {
            p <- p + ggplot2::scale_color_manual(
                values = palette, name = "Cluster"
            )
        }
        return(
            p +
                ggplot2::facet_wrap(
                    ggplot2::vars(.data$sensor), ncol = ncol,
                    scales = if (free_y) "free_y" else "fixed"
                ) +
                ggplot2::labs(
                    x = "Time (s)", y = outcome_id,
                    subtitle = "Observed mean: black solid; predicted: blue dashed"
                ) +
                ggplot2::theme_minimal(base_size = 10) +
                #theme_neurogam() +
                ggplot2::theme(
                    #panel.border = ggplot2::element_blank(),
                    panel.grid.minor = ggplot2::element_blank(),
                    strip.text = ggplot2::element_text(size = 8)
                )
        )
    }

    # Scalp layout: transform individual temporal and amplitude axes into
    # small, uniformly scaled traces positioned at the sensor coordinates.
    positions <- unique(dat[c("sensor", "x", "y")])
    n_sensors <- nrow(positions)
    if (n_sensors > 1L) {
        distances <- as.matrix(stats::dist(positions[c("x", "y")]))
        diag(distances) <- Inf
        spacing <- stats::median(apply(distances, 1L, min))
    } else {
        spacing <- 0.5
    }
    if (is.null(mini_width)) {
        mini_width <- 0.8 * spacing
    }
    if (is.null(mini_height)) {
        mini_height <- 0.75 * spacing
    }
    if (any(!is.finite(c(mini_width, mini_height))) ||
        any(c(mini_width, mini_height) <= 0)) {
        stop("Mini-trace width and height must be positive finite values.",
             call. = FALSE)
    }

    range_time <- range(dat$time)
    span_time <- diff(range_time)
    if (span_time == 0) {
        stop("Scalp mini-trajectories require at least two time points.",
             call. = FALSE)
    }
    range_amp <- if (is.null(amplitude_limits)) {
        range(c(dat$observed, dat$predicted,
                if (show_ci) c(dat$lower, dat$upper)),
              na.rm = TRUE)
    } else {
        if (!is.numeric(amplitude_limits) || length(amplitude_limits) != 2L ||
            any(!is.finite(amplitude_limits)) ||
            amplitude_limits[2L] <= amplitude_limits[1L]) {
            stop("`amplitude_limits` must be two increasing numbers.",
                 call. = FALSE)
        }
        amplitude_limits
    }
    if (diff(range_amp) == 0) {
        range_amp <- range_amp + c(-0.5, 0.5)
    }
    amp_mid <- mean(range_amp)
    amp_span <- diff(range_amp)
    st_x <- function(time, x) {
        x + ((time - range_time[1L]) / span_time - 0.5) * mini_width
    }
    st_y <- function(value, y) {
        y + (value - amp_mid) / amp_span * mini_height
    }

    dat$plot_x <- st_x(dat$time, dat$x)
    dat$plot_observed <- st_y(dat$observed, dat$y)
    dat$plot_predicted <- st_y(dat$predicted, dat$y)
    dat$plot_lower <- st_y(dat$lower, dat$y)
    dat$plot_upper <- st_y(dat$upper, dat$y)

    p <- ggplot2::ggplot()
    if (head_outline) {
        geom <- full_head
        p <- p +
            ggplot2::geom_path(
                data = geom$outline,
                ggplot2::aes(x = .data$x, y = .data$y),
                color = "grey55", linewidth = 0.55
            ) +
            ggplot2::geom_path(
                data = geom$nose,
                ggplot2::aes(x = .data$x, y = .data$y),
                color = "grey55", linewidth = 0.55
            )
    }

    if (cluster_style == "rectangles" && nrow(intervals)) {
        anchors <- merge(intervals, positions, by = "sensor", sort = FALSE)
        anchors$plot_xmin <- st_x(anchors$xmin, anchors$x)
        anchors$plot_xmax <- st_x(anchors$xmax, anchors$x)
        anchors$plot_ymin <- anchors$y - mini_height / 2
        anchors$plot_ymax <- anchors$y + mini_height / 2
        p <- p + ggplot2::geom_rect(
            data = anchors,
            ggplot2::aes(
                xmin = .data$plot_xmin, xmax = .data$plot_xmax,
                ymin = .data$plot_ymin, ymax = .data$plot_ymax,
                fill = .data$sign
            ),
            inherit.aes = FALSE, alpha = cluster_alpha
        )
    }

    # Baseline at zero activity when zero belongs to the displayed range.
    if (range_amp[1L] <= 0 && range_amp[2L] >= 0) {
        baselines <- positions
        baselines$xmin <- baselines$x - mini_width / 2
        baselines$xmax <- baselines$x + mini_width / 2
        baselines$baseline <- st_y(0, baselines$y)
        p <- p + ggplot2::geom_segment(
            data = baselines,
            ggplot2::aes(
                x = .data$xmin, xend = .data$xmax,
                y = .data$baseline, yend = .data$baseline
            ),
            color = "grey80", linewidth = 0.25
        )
    }
    if (show_ci) {
        p <- p + ggplot2::geom_ribbon(
            data = dat,
            ggplot2::aes(
                x = .data$plot_x, ymin = .data$plot_lower,
                ymax = .data$plot_upper, group = .data$sensor
            ),
            fill = "#3973A6", alpha = 0.12
        )
    }
    p <- p +
        ggplot2::geom_line(
            data = dat,
            ggplot2::aes(
                x = .data$plot_x, y = .data$plot_observed,
                group = .data$sensor
            ),
            color = "#202020", linewidth = 0.35, na.rm = TRUE
        ) +
        ggplot2::geom_line(
            data = dat,
            ggplot2::aes(
                x = .data$plot_x, y = .data$plot_predicted,
                group = .data$sensor
            ),
            color = "#3275B2", linewidth = 0.46,
            linetype = "dashed"
        )

    if (cluster_style == "segments" && nrow(members)) {
        if (nrow(segments)) {
            segments$plot_x <- st_x(segments$time, segments$x)
            segments$plot_y <- st_y(segments$predicted, segments$y)
            segments$plot_xend <- st_x(segments$timeend, segments$x)
            segments$plot_yend <- st_y(segments$predictedend, segments$y)
            p <- p + ggplot2::geom_segment(
                data = segments,
                ggplot2::aes(
                    x = .data$plot_x, y = .data$plot_y,
                    xend = .data$plot_xend, yend = .data$plot_yend,
                    color = .data$cluster_sign
                ),
                linewidth = highlighted_linewidth, alpha = 1
            )
        }
        members$plot_x <- st_x(members$time, members$x)
        members$plot_y <- st_y(members$predicted, members$y)
        p <- p + ggplot2::geom_point(
            data = members,
            ggplot2::aes(
                x = .data$plot_x, y = .data$plot_y,
                color = .data$cluster_sign
            ),
            size = 0.45, alpha = 0.85
        )
    }
    if (show_labels) {
        positions$label_y <- positions$y + mini_height * 0.6
        p <- p + ggplot2::geom_text(
            data = positions,
            ggplot2::aes(
                x = .data$x, y = .data$label_y,
                label = .data$sensor
            ),
            size = 2.3, check_overlap = TRUE
        )
    }

    if (cluster_style == "rectangles" && nrow(intervals)) {
        p <- p + ggplot2::scale_fill_manual(
            values = palette, name = "Cluster"
        )
    }
    if (cluster_style == "segments" && nrow(members)) {
        p <- p + ggplot2::scale_color_manual(
            values = palette, name = "Cluster"
        )
    }
    p +
        ggplot2::coord_equal(clip = "off") +
        ggplot2::labs(
            title = "Observed and predicted EEG time courses",
            subtitle = "Black: observed mean; blue dashed: posterior median expectation",
            caption = paste0(
                "All mini-traces share the same time and amplitude scales. ",
                "Activity limits: [",
                paste(signif(range_amp, 3), collapse = ", "), "]"
            ),
            x = NULL, y = NULL
        ) +
        ggplot2::theme_void(base_size = 11) +
        #theme_neurogam() +
        ggplot2::theme(
            #panel.border = ggplot2::element_blank(),
            plot.margin = grid::unit(rep(12, 4), "pt"),
            legend.position = "right"
        )
}
