#' Plot spatio-temporal EEG data as topographies or 2D surfaces
#'
#' Create facetted EEG maps across multiple time points from either
#' i) a \pkg{brms} spatio-temporal model (predicted values) or
#' ii) raw/summarised EEG data in long format (observed values).
#'
#' The function can render either a scalp topography (\code{type = "topo"})
#' with a head outline and optional sensor markers/labels, or a simple 2D
#' surface plot (\code{type = "surface"}) without a head outline.
#'
#' @param x Either a data frame containing EEG data (raw or summarised) or a
#'   fitted \code{brmsfit} object. If a data frame is provided, it must contain
#'   columns specified by \code{time_col}, \code{x_col}, \code{y_col}, and
#'   \code{value_col}. If a \code{brmsfit} is provided, predictions are computed
#'   from \code{x$data} and \code{x} must include a compatible spatio-temporal
#'   smooth (e.g., \code{gp(time, xproj, yproj, ...)} or \code{t2(time, xproj, yproj, ...)}).
#' @param type Character; plot type. \code{"topo"} draws a scalp outline and masks
#'   values outside the head. \code{"surface"} draws only the interpolated/predicted
#'   field.
#'
#' @param times Numeric vector of time points to plot. If \code{NULL}, all unique
#'   times found in the input are used (or a subset if \code{n_times} is provided).
#' @param n_times Optional integer; if not \code{NULL}, selects \code{n_times}
#'   approximately equally spaced time points from \code{times}.
#'
#' @param sensors Optional data frame of sensor positions used to draw the head
#'   outline, sensor points, and sensor labels. Must contain columns \code{xproj}
#'   and \code{yproj}. If \code{NULL}, sensor coordinates are inferred from the
#'   unique \code{x_col}/\code{y_col} pairs in \code{x} (labels will then be unavailable).
#'
#' @param value_col Character; name of the value column in raw data (e.g., \code{"voltage"}).
#' @param time_col,x_col,y_col Character; column names in \code{x} giving time and 2D sensor
#'   coordinates (projected).
#'
#' @param grid_res Integer; grid resolution used to evaluate predictions and/or
#'   interpolate data to a regular grid (approximately \code{grid_res^2} pixels per facet).
#' @param head_expand Numeric; expansion factor used to compute the head radius
#'   from sensor coordinates (values > 1 add padding).
#'
#' @param re_formula Passed to \code{\link[brms:fitted.brmsfit]{brms::fitted}} when \code{x}
#'   is a \code{brmsfit}. Use \code{NULL} for default behaviour, or \code{NA} to exclude
#'   group-level terms.
#' @param probs Numeric vector of length 2 giving the quantiles returned by
#'   \code{brms::fitted(..., probs = probs)}.
#' @param ndraws Optional integer; number of posterior draws used by \code{brms::fitted}.
#'
#' @param show_sensors Logical; if \code{TRUE} and \code{type = "topo"}, plot sensor points.
#' @param sensor_size Numeric; point size for sensors.
#' @param sensor_labels Logical; if \code{TRUE}, add sensor labels (requires \code{sensors}
#'   to contain \code{sensor_label_col}).
#' @param sensor_label_col Character; column name in \code{sensors} providing sensor names.
#' @param sensor_label_size Numeric; text size for sensor labels.
#' @param sensor_label_repel Logical; if \code{TRUE}, use \pkg{ggrepel} for labels.
#'
#' @param contours Logical; if \code{TRUE}, add contour lines.
#' @param contour_bins Integer; number of contour bins.
#' @param palette Character; palette name passed to \code{\link[ggplot2]{scale_fill_distiller}}.
#'
#' @param facet_nrow,facet_ncol Integers; layout for \code{\link[ggplot2]{facet_wrap}}.
#' @param facet_scales Character; facet scaling, passed to \code{facet_wrap(scales = ...)}.
#' @param facet_label_prefix Character; prefix used when generating facet labels (e.g., \code{"Time: "}).
#' @param facet_unit Character; unit suffix used when generating facet labels (e.g., \code{"s"}).
#'
#' @param fill_limits Either a character string or numeric vector.
#'   If character, one of \code{"global_quantile"}, \code{"global"}, or \code{"none"}.
#'   If numeric, a vector of length 2 specifying explicit colour scale limits
#'   (e.g., \code{c(-5, 5)}).
#' @param limit_quantiles Numeric vector of length 2; quantiles used when
#'   \code{fill_limits = "global_quantile"}.
#'
#' @param theme A \code{\link[ggplot2:theme]{theme}} object
#'   modifying the appearance of the plots.
#'
#' @return A \code{ggplot} object.
#'
#' @details
#' For \code{type = "topo"}, values are masked outside the head radius computed from
#' sensor coordinates and a head outline (circle + nose + ears) is drawn.
#'
#' For raw data, values are interpolated to a regular grid per time point using
#' \code{st_interp_to_grid()}, which relies on \pkg{akima} and can optionally fill
#' missing grid cells using nearest-neighbour (\pkg{FNN}).
#'
#' For \code{brmsfit} inputs, predictions are obtained via \code{brms::fitted()} on a
#' regular grid. When \code{type = "topo"}, missing values inside the head can be
#' filled using \code{st_fill_head()} (nearest-neighbour).
#'
#' This function relies on internal helpers such as \code{st_take_n_times()},
#' \code{st_make_grid()}, \code{st_head_radius()}, \code{st_head_outline()},
#' \code{st_build_time_labels()}, \code{st_order_time_factor()},
#' \code{st_compute_limits()}, \code{st_predict_brms()}, \code{st_fill_head()},
#' and \code{st_interp_to_grid()}.
#'
#' @author Ladislas Nalborczyk \email{ladislas.nalborczyk@@cnrs.fr}
#'
#' @examples
#' \dontrun{
#' # Summarised EEG data (long format)
#' plot_eeg(
#'   eeg_data_summary,
#'   type = "topo",
#'   sensors = sensors,
#'   times = c(0, 0.1, 0.2, 0.3),
#'   grid_res = 80,
#'   contours = FALSE,
#'   facet_nrow = 2
#'   )
#'
#' # brms model predictions
#' plot_eeg(
#'   spatio_temporal_gam,
#'   type = "topo",
#'   sensors = sensors,
#'   times = c(0, 0.1, 0.2, 0.3),
#'   ndraws = 200,
#'   # modifying default ggplot2 theme
#'   theme = theme_bw(base_size = 12, base_family = "Open Sans")
#'   )
#' }
#'
#' @export
plot_eeg <- function (
        # either raw EEG data or a brms spatio-temporal model
        x,
        type = c("topo", "surface"),
        times = NULL,
        # optional, if NULL, choose N equally spaced from available times
        n_times = NULL,
        # required for topo head radius and sensor points
        sensors = NULL,
        # raw data value column
        value_col = "voltage",
        time_col = "time",
        x_col = "xproj",
        y_col = "yproj",

        # head grid
        grid_res = 100,
        head_expand = 1.1,

        # prediction controls (brms only)
        re_formula = NULL,
        probs = c(0.025, 0.975),
        ndraws = NULL,

        # display options
        show_sensors = TRUE,
        sensor_size = 1.5,
        sensor_labels = FALSE,
        sensor_label_col = "channel",
        sensor_label_size = 3,
        sensor_label_repel = TRUE,

        contours = TRUE,
        contour_bins = 10,
        palette = "RdBu",

        facet_nrow = NULL,
        facet_ncol = NULL,
        facet_scales = "fixed",
        facet_label_prefix = "Time: ",
        facet_unit = "s",

        # colour scaling
        fill_limits = "global_quantile",
        limit_quantiles = c(0.01, 0.99),

        # ggplot2 theme
        theme = ggplot2::theme_void()
        ) {

    type <- match.arg(type)

    if (!is.numeric(fill_limits) ) {

        fill_limits <- match.arg(fill_limits, choices = c("global_quantile", "global", "none") )

    } else {

        if (length(fill_limits) != 2 || any(!is.finite(fill_limits) ) ) {

            stop ("`fill_limits` must be a numeric vector of length 2 with finite values.", call. = FALSE)

        }

    }

    is_brms <- inherits(x, "brmsfit")
    is_df <- is.data.frame(x)

    if (!is_brms && !is_df) stop ("`x` must be a data.frame (raw) or a brmsfit.", call. = FALSE)

    if (!ggplot2::is.theme(theme) ) {

        stop ("Argument 'theme' should be a 'theme' object.")

    }

    # get data for ranges / available times
    base_dat <- if (is_brms) x$data else x

    # rename columns internally (so downstream is simpler)
    dat_std <- base_dat |>
        dplyr::transmute(
            time = .data[[time_col]],
            xproj = .data[[x_col]],
            yproj = .data[[y_col]],
            value = if (!is_brms) .data[[value_col]] else NA_real_
            )

    # choose times
    if (is.null(times) ) {

        # default = unique times in data/model data
        times <- sort(unique(dat_std$time) )

    } else {

        times <- as.numeric(times)

    }

    if (!is.null(n_times) ) {

        times <- st_take_n_times(sort(unique(times) ), n_times)

    }

    # if sensors not provided, try to build from data (unique coords)
    if (is.null(sensors) ) {

        # NB: labels won't be available in this case
        sensors <- dat_std |> dplyr::distinct(.data$xproj, .data$yproj)

    } else {

        # standardize required cols
        sensors <- sensors

        if (!all(c("xproj","yproj") %in% names(sensors) ) ) {

            stop ("`sensors` must contain columns xproj and yproj.", call. = FALSE)

        }

    }

    # x/y ranges
    x_range <- range(dat_std$xproj, na.rm = TRUE)
    y_range <- range(dat_std$yproj, na.rm = TRUE)

    # build prediction grid (for brms) or data slice (for raw)
    if (is_brms) {

        if (type == "topo") {

            r_head <- st_head_radius(sensors, head_expand = head_expand)
            grid_xy <- st_make_grid(c(-r_head, r_head), c(-r_head, r_head), grid_res = grid_res)

        } else {

            grid_xy <- st_make_grid(x_range, y_range, grid_res = grid_res)

        }

        newdata <- do.call(
            rbind,
            lapply(times, function (tt) {
                cbind(grid_xy, time = tt, time_facet = tt)
            })
        )

        pred_df <- st_predict_brms(
            fit = x,
            newdata = newdata,
            re_formula = re_formula,
            probs = probs,
            ndraws = ndraws
            )

        plot_df <- pred_df |>
            dplyr::transmute(
                xproj = .data$xproj,
                yproj = .data$yproj,
                time = .data$time_facet,
                Estimate = .data$Estimate,
                Q2.5 = .data$Q2.5,
                Q97.5 = .data$Q97.5
                )

        # fills the head
        if (type == "topo") {

            r_head <- st_head_radius(sensors, head_expand = head_expand)
            plot_df <- st_fill_head(plot_df, r_head = r_head, fill_na = "nearest")

        }

    } else {

        # raw: filter to requested times (tolerant matching via nearest)
        # if times not exactly present, pick nearest in data
        avail <- sort(unique(dat_std$time) )
        nearest <- function (t0) avail[which.min(abs(avail - t0) )]
        times_use <- vapply(times, nearest, numeric(1) )

        # assume already one value per (time, xproj, yproj)
        plot_df <- dat_std |>
            dplyr::filter(.data$time %in% times_use) |>
            dplyr::rename(Estimate = .data$value)

        # regular grid for raster plotting
        # interpolate each time slice to the grid
        # grid_xy <- st_make_grid(x_range, y_range, grid_res = grid_res)
        # r_head <- if (type == "topo") st_head_radius(sensors, head_expand = head_expand) else NULL

        if (type == "topo") {

            r_head <- st_head_radius(sensors, head_expand = head_expand)
            grid_xy <- st_make_grid(c(-r_head, r_head), c(-r_head, r_head), grid_res = grid_res)

        } else {

            grid_xy <- st_make_grid(x_range, y_range, grid_res = grid_res)

        }

        plot_df <- plot_df |>
            dplyr::group_by(.data$time) |>
            dplyr::group_modify(~st_interp_to_grid(
                df_t = .x,
                grid_xy = grid_xy,
                value = "Estimate",
                r_head = r_head,
                # fills the head
                fill_na = "nearest",
                extrap = FALSE
                ) ) |>
            dplyr::ungroup()

    }

    # facet labels
    time_labels <- st_build_time_labels(
        times = sort(unique(plot_df$time) ),
        unit = facet_unit,
        prefix = facet_label_prefix
        )

    lab_map <- data.frame(time = sort(unique(plot_df$time) ), time_lab = time_labels)

    plot_df <- plot_df |>
        dplyr::left_join(lab_map, by = "time") |>
        dplyr::mutate(time_lab = st_order_time_factor(.data$time_lab) )

    # topo masking (only for topo)
    if (type == "topo") {

        r_head <- st_head_radius(sensors, head_expand = head_expand)
        plot_df <- plot_df |> dplyr::filter(sqrt(.data$xproj^2 + .data$yproj^2) <= r_head)
        head_df <- st_head_outline(r_head)

    } else {

        head_df <- NULL
        r_head <- NULL

    }

    # limits for colour scale
    if (is.numeric(fill_limits) ) {

        lim <- fill_limits

    } else {

        lim <- st_compute_limits(
            plot_df,
            z = "Estimate",
            limits = fill_limits,
            limit_quantiles = limit_quantiles
            )

    }

    # build base plot
    p <- plot_df |>
        ggplot2::ggplot(ggplot2::aes(x = .data$xproj, y = .data$yproj) ) +
        ggplot2::geom_raster(ggplot2::aes(fill = .data$Estimate), interpolate = TRUE) +
        ggplot2::coord_equal() +
        theme +
        ggplot2::scale_fill_distiller(
            palette = palette,
            direction = -1,
            limits = lim,
            oob = scales::squish
            ) +
        ggplot2::labs(fill = if (is_brms) "Predicted" else "Observed")

    # optional contours
    if (contours) {

        p <- p + ggplot2::geom_contour(
            ggplot2::aes(z = .data$Estimate),
            bins = contour_bins,
            linewidth = 0.5,
            alpha = 0.8,
            color = "black"
            )

    }

    # head outline + sensors for topo
    if (type == "topo") {

        p <- p +
            ggplot2::geom_path(
                data = head_df,
                ggplot2::aes(x = .data$x, y = .data$y, group = .data$part),
                inherit.aes = FALSE,
                linewidth = 0.5
                )

        if (show_sensors) {

            sensors_plot <- merge(
                sensors,
                data.frame(time_lab = unique(plot_df$time_lab), stringsAsFactors = FALSE),
                by = NULL
                )

            p <- p + ggplot2::geom_point(
                data = sensors_plot,
                ggplot2::aes(x = .data$xproj, y = .data$yproj),
                inherit.aes = FALSE,
                size = sensor_size,
                na.rm = TRUE
                )

            if (sensor_labels) {

                if (!sensor_label_col %in% names(sensors) ) {

                    stop (sprintf("`sensor_label_col` (%s) not found in `sensors`.", sensor_label_col), call. = FALSE)

                }

                sensors_plot$label <- as.character(sensors_plot[[sensor_label_col]])

                if (sensor_label_repel) {

                    p <- p + ggrepel::geom_text_repel(
                        data = sensors_plot,
                        ggplot2::aes(x = .data$xproj, y = .data$yproj, label = .data$label),
                        inherit.aes = FALSE,
                        size = sensor_label_size,
                        max.overlaps = Inf,
                        min.segment.length = 0,
                        na.rm = TRUE
                        )

                } else {

                    p <- p + ggplot2::geom_text(
                        data = sensors_plot,
                        ggplot2::aes(x = .data$xproj, y = .data$yproj, label = .data$label),
                        inherit.aes = FALSE,
                        size = sensor_label_size,
                        vjust = -0.6,
                        na.rm = TRUE
                        )

                }

            }

        }

    }

    # facets
    p <- p + ggplot2::facet_wrap(
        ~time_lab,
        nrow = facet_nrow,
        ncol = facet_ncol,
        scales = facet_scales
        )

    return (p)

}
