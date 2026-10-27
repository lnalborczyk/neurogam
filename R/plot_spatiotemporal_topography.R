#' Compare observed, predicted and residual EEG scalp maps through time
#'
#' Plots interpolated maps of observed sensor-level means, population-level
#' posterior median expected responses, and residuals (`observed - predicted`).
#' Inferences and clusters remain at the *original sensor locations*, not at
#' interpolated pixels. Interpolation is for visualisation only and uses
#' inverse-distance weighting (IDW) without additional package dependencies.
#'
#' @param results Object returned by [analyse_bgam_spatiotemporal()].
#' @param observed_data Original long-form observations, or condition A data.
#' @param observed_data_b Observed condition B if `results` is an A - B contrast.
#' @param outcome_id Observed outcome column.
#' @param sensor_id Sensor identifier shared by observed data and predictions.
#' @param time_id Numeric time column.
#' @param participant_id Optional participant ID for an equal-participant-
#'   weighted observed average. With `NULL`, observations are pooled.
#' @param spatial_id Two projected spatial coordinate columns.
#' @param times Optional requested times; closest existing samples are used.
#' @param n_times Number of evenly distributed panels if `times = NULL`.
#' @param what One or more of `"observed"`, `"predicted"`, `"residual"`.
#' @param grid_n Number of interpolated pixels per grid axis.
#' @param idw_power Positive IDW inverse-distance exponent.
#' @param limits Symmetric range half-width of a shared color scale, or `NULL`
#'   for automatic scaling across all selected panels.
#' @param contour Whether to draw interpolated iso-value contours.
#' @param show_sensors Show original electrode locations.
#' @param show_labels Show sensor names (can be busy with dense montages).
#' @param head_outline Add schematic scalp circle and nose.
#' @param cluster_overlay `"predicted"` highlights inference clusters only on
#'   predicted maps; `"all"` on all maps; `"none"` disables highlights.
#' @param cluster_alpha Transparency of cluster-marking rings.
#' @param cluster_size Size of cluster-marking rings.
#'
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' plot_spatiotemporal_topography(
#'     results_st, eeg_data2, outcome_id = "voltage",
#'     sensor_id = "channel", n_times = 8
#' )
#' }
#'
#' @importFrom rlang .data
#'
#' @author Ladislas Nalborczyk \email{ladislas.nalborczyk@@cnrs.fr}.
#'
#' @export
plot_spatiotemporal_topography <- function(
    results,
    observed_data,
    observed_data_b = NULL,
    outcome_id = "voltage",
    sensor_id = "channel",
    time_id = "time",
    spatial_id = c("xproj", "yproj"),
    participant_id = NULL,
    times = NULL,
    n_times = 8L,
    what = c("observed", "predicted", "residual"),
    grid_n = 85L,
    idw_power = 2,
    limits = NULL,
    contour = TRUE,
    show_sensors = TRUE,
    show_labels = FALSE,
    head_outline = TRUE,
    cluster_overlay = c("predicted", "all", "none"),
    cluster_alpha = 0.9,
    cluster_size = 2.8
    ) {
    cluster_overlay <- match.arg(cluster_overlay)
    what <- match.arg(
        what,
        choices = c("observed", "predicted", "residual"),
        several.ok = TRUE
    )
    if (!is.numeric(n_times) || length(n_times) != 1L ||
        !is.finite(n_times) || n_times < 1 || n_times != floor(n_times)) {
        stop("`n_times` must be a positive integer.", call. = FALSE)
    }
    if (!is.numeric(grid_n) || length(grid_n) != 1L ||
        !is.finite(grid_n) || grid_n < 10 || grid_n != floor(grid_n)) {
        stop("`grid_n` must be an integer of at least 10.", call. = FALSE)
    }
    if (!is.numeric(idw_power) || length(idw_power) != 1L ||
        !is.finite(idw_power) || idw_power <= 0) {
        stop("`idw_power` must be positive.", call. = FALSE)
    }
    if (!is.numeric(cluster_alpha) || length(cluster_alpha) != 1L ||
        !is.finite(cluster_alpha) || cluster_alpha < 0 || cluster_alpha > 1) {
        stop("`cluster_alpha` must be between zero and one.", call. = FALSE)
    }

    dat <- prepare_spatiotemporal_plot_data(
        results, observed_data, observed_data_b,
        outcome_id, sensor_id, time_id, spatial_id, participant_id
    )
    unique_times <- sort(unique(dat$time))
    if (is.null(times)) {
        indices <- unique(round(seq(
            1, length(unique_times),
            length.out = min(n_times, length(unique_times))
        )))
        selected <- unique_times[indices]
    } else {
        if (!is.numeric(times) || !length(times) || any(!is.finite(times))) {
            stop("`times` must contain finite numeric values.", call. = FALSE)
        }
        if (any(times < min(unique_times) | times > max(unique_times))) {
            stop("Requested `times` must be inside the prediction time range.",
                 call. = FALSE)
        }
        selected <- sort(unique(vapply(times, function(z) {
            unique_times[which.min(abs(unique_times - z))]
        }, numeric(1))))
    }

    dat <- dat[dat$time %in% selected, , drop = FALSE]
    labels <- sprintf("%.4f s", selected)
    label_for_time <- function(z) {
        factor(sprintf("%.4f s", z), levels = labels)
    }

    columns <- list()
    for (view in what) {
        d <- dat[c("sensor", "time", "x", "y", view)]
        names(d)[5L] <- "value"
        d$view <- view
        d$time_label <- label_for_time(d$time)
        columns[[view]] <- d
    }
    long <- do.call(rbind, columns)
    rownames(long) <- NULL
    long$view <- factor(long$view, levels = what)

    if (is.null(limits)) {
        lim <- max(abs(long$value), na.rm = TRUE)
        if (!is.finite(lim) || lim == 0) {
            lim <- 1
        }
    } else {
        if (!is.numeric(limits) || length(limits) != 1L ||
            !is.finite(limits) || limits <= 0) {
            stop("`limits` must be a positive finite half-range.",
                 call. = FALSE)
        }
        lim <- limits
    }

    positions <- unique(dat[c("sensor", "x", "y")])
    if (nrow(positions) < 3L ||
        qr(cbind(1, positions$x, positions$y))$rank < 3L) {
        stop("Topographies require at least three non-collinear sensors.",
             call. = FALSE)
    }
    geom <- st_head_geometry(dat)
    side <- seq(-geom$radius, geom$radius, length.out = grid_n)
    target <- expand.grid(
        x = geom$center[1L] + side,
        y = geom$center[2L] + side
    )
    inside <- (target$x - geom$center[1L])^2 +
        (target$y - geom$center[2L])^2 <= geom$radius^2
    groups <- split(
        seq_len(nrow(long)),
        interaction(long$view, long$time_label, drop = TRUE)
    )
    surfaces <- lapply(groups, function(index) {
        d <- long[index, , drop = FALSE]
        out <- target
        out$value <- interpolate_st_topography(
            d, out$x, out$y, power = idw_power
        )
        out$value[!inside] <- NA_real_
        out$view <- d$view[1L]
        out$time_label <- d$time_label[1L]
        out
    })
    surface <- do.call(rbind, surfaces)
    rownames(surface) <- NULL
    surface$view <- factor(surface$view, levels = what)
    surface$time_label <- factor(surface$time_label, levels = labels)

    p <- ggplot2::ggplot(
        surface,
        ggplot2::aes(x = .data$x, y = .data$y)
    ) +
        ggplot2::geom_raster(
            ggplot2::aes(fill = .data$value), na.rm = TRUE
        ) +
        ggplot2::scale_fill_gradient2(
            low = "#346AA1", mid = "#FAFAFA", high = "#C7584D",
            midpoint = 0, limits = c(-lim, lim),
            oob = function(x, range) pmin(pmax(x, range[1L]), range[2L]),
            na.value = NA_character_, name = "Activity"
        )

    if (contour && any(is.finite(surface$value))) {
        p <- p + ggplot2::geom_contour(
            ggplot2::aes(z = .data$value),
            color = "grey30", linewidth = 0.2,
            bins = 7, na.rm = TRUE
        )
    }

    panel_keys <- unique(long[c("view", "time_label")])
    if (head_outline) {
        head <- merge(geom$outline, panel_keys, by = NULL)
        nose <- merge(geom$nose, panel_keys, by = NULL)
        p <- p +
            ggplot2::geom_path(
                data = head,
                ggplot2::aes(x = .data$x, y = .data$y),
                inherit.aes = FALSE, linewidth = 0.45
            ) +
            ggplot2::geom_path(
                data = nose,
                ggplot2::aes(x = .data$x, y = .data$y),
                inherit.aes = FALSE, linewidth = 0.45
            )
    }

    sensors <- unique(long[c("sensor", "time_label", "x", "y")])
    sensors <- merge(sensors, data.frame(view = what), by = NULL)
    sensors$view <- factor(sensors$view, levels = what)
    if (show_sensors) {
        p <- p + ggplot2::geom_point(
            data = sensors,
            ggplot2::aes(x = .data$x, y = .data$y),
            inherit.aes = FALSE, color = "#202020", size = 0.6
        )
    }
    if (show_labels) {
        p <- p + ggplot2::geom_text(
            data = sensors,
            ggplot2::aes(x = .data$x, y = .data$y, label = .data$sensor),
            inherit.aes = FALSE, size = 2.1, vjust = -0.7
        )
    }
    if (cluster_overlay != "none") {
        marks <- dat[!is.na(dat$cluster_id), , drop = FALSE]
        if (nrow(marks)) {
            views <- if (cluster_overlay == "all") what else {
                intersect("predicted", what)
            }
            if (length(views)) {
                marks <- merge(marks, data.frame(view = views), by = NULL)
                marks$view <- factor(marks$view, levels = what)
                marks$time_label <- label_for_time(marks$time)
                p <- p +
                    ggplot2::geom_point(
                        data = marks,
                        ggplot2::aes(
                            x = .data$x, y = .data$y,
                            color = .data$cluster_sign
                        ),
                        inherit.aes = FALSE, shape = 1,
                        alpha = cluster_alpha, size = cluster_size,
                        stroke = 1.05
                    ) +
                    ggplot2::scale_color_manual(
                        values = c(positive = "#B53E32", negative = "#225AA0"),
                        name = "Posterior cluster"
                    )
            }
        }
    }

    p +
        ggplot2::facet_grid(
            rows = ggplot2::vars(.data$view),
            cols = ggplot2::vars(.data$time_label),
            drop = FALSE
        ) +
        ggplot2::coord_equal(expand = FALSE) +
        ggplot2::labs(x = NULL, y = NULL) +
        ggplot2::theme_void(base_size = 11) +
        ggplot2::theme(
            strip.text = ggplot2::element_text(size = 9),
            strip.background = ggplot2::element_rect(
                fill = "grey95", color = NA
            ),
            panel.spacing = grid::unit(0.35, "lines"),
            legend.position = "right"
        )
}
