#' Plot EEG sensor positions on a 2D scalp
#'
#' Visualise the 2D layout of EEG sensors using projected coordinates (e.g.,
#' \code{xproj} and \code{yproj}). Optionally draws a head outline, sensor
#' markers, and/or sensor labels. A subset of sensors can be highlighted (by
#' name or by row index), with non-highlighted sensors optionally dimmed.
#'
#' @param sensors A data frame containing sensor coordinates. Must include the
#'   columns specified by \code{x_col} and \code{y_col}. If \code{show_labels = TRUE}
#'   and labels are desired, it must also include \code{label_col}.
#' @param x_col,y_col Character scalars giving the column names in \code{sensors}
#'   corresponding to the x and y projected sensor coordinates.
#' @param label_col Character scalar giving the column name in \code{sensors} used
#'   as sensor labels (e.g., \code{"channel"}).
#'
#' @param show_head Logical; if \code{TRUE}, draw a head outline (circle + nose + ears).
#' @param head_expand Numeric scalar; multiplier applied to the maximum sensor radius
#'   when computing the head radius. Values > 1 add padding around sensors.
#'
#' @param show_points Logical; if \code{TRUE}, draw sensor markers for non-highlighted sensors.
#' @param point_size Numeric; size of non-highlighted sensor markers.
#'
#' @param show_labels Logical; if \code{TRUE}, add sensor labels.
#' @param label_size Numeric; text size for sensor labels.
#' @param label_repel Logical; if \code{TRUE}, use \pkg{ggrepel} to reduce label overlap.
#'   If \code{FALSE}, uses \code{\link[ggplot2]{geom_text}}.
#' @param label_only_highlight Logical; if \code{TRUE} and \code{highlight} is not \code{NULL},
#'   show labels only for highlighted sensors.
#'
#' @param highlight Optional vector specifying sensors to highlight. If numeric,
#'   treated as row indices of \code{sensors}. Otherwise treated as sensor names
#'   that are matched against \code{sensors[[label_col]]}.
#' @param highlight_col Character; colour used for highlighted sensors.
#' @param highlight_size Numeric; marker size for highlighted sensors.
#' @param highlight_shape Numeric; point shape for highlighted sensors (as in \code{ggplot2}).
#' @param highlight_alpha Numeric in [0, 1]; alpha transparency for highlighted sensors.
#' @param dim_others Logical; if \code{TRUE}, apply \code{other_alpha} to non-highlighted
#'   sensor markers when \code{show_points = TRUE}.
#' @param other_alpha Numeric in [0, 1]; alpha transparency for non-highlighted sensors when
#'   \code{dim_others = TRUE}.
#'
#' @param xlim,ylim Optional numeric vectors of length 2 giving x and y limits.
#'   If \code{NULL}, limits are computed automatically from the head radius (if
#'   \code{show_head = TRUE}) or from the range of sensor coordinates otherwise.
#' @param theme A \code{\link[ggplot2:theme]{theme}} object
#'   modifying the appearance of the plots.
#'
#' @return A \code{ggplot} object.
#'
#' @details
#' The head outline is computed from the maximum sensor radius (multiplied by
#' \code{head_expand}) using \code{st_head_radius()} and \code{st_head_outline()}.
#' These helper functions must be available in the calling environment.
#'
#' Highlighting works in two modes:
#' \itemize{
#'   \item If \code{highlight} is numeric, indices are used to mark rows of \code{sensors}.
#'   \item If \code{highlight} is character, values are matched against \code{sensors[[label_col]]}.
#' }
#'
#' @author Ladislas Nalborczyk \email{ladislas.nalborczyk@@cnrs.fr}
#'
#' @examples
#' \dontrun{
#' # Basic sensor layout with head outline and labels
#' plot_sensors(sensors, show_points = TRUE, show_labels = TRUE)
#'
#' # Highlight a few channels by name
#' plot_sensors(
#'   sensors,
#'   show_points = TRUE,
#'   show_labels = TRUE,
#'   highlight = c("CZ", "C1", "C2"),
#'   label_only_highlight = TRUE
#'   )
#'
#' # Highlight by row index
#' plot_sensors(
#'   sensors,
#'   show_points = TRUE,
#'   highlight = c(10, 25, 42),
#'   highlight_col = "dodgerblue"
#'   )
#' }
#'
#' @export
plot_sensors <- function (
        sensors,
        x_col = "xproj",
        y_col = "yproj",
        label_col = "channel",

        show_head = TRUE,
        head_expand = 1.1,

        show_points = FALSE,
        point_size = 2,

        show_labels = TRUE,
        label_size = 3,
        label_repel = TRUE,
        label_only_highlight = FALSE,

        # vector of sensor names or row indices
        highlight = NULL,
        highlight_col = "orangered",
        highlight_size = 3,
        highlight_shape = 16,
        highlight_alpha = 1,
        dim_others = TRUE,
        other_alpha = 0.3,

        xlim = NULL,
        ylim = NULL,
        theme = ggplot2::theme_void()
        ) {

    if (!all(c(x_col, y_col) %in% names(sensors) ) ) {

        stop (
            sprintf(
                "`sensors` must contain columns '%s' and '%s'.",
                x_col, y_col
                ),
            call. = FALSE
            )

    }

    if (!ggplot2::is.theme(theme) ) {

        stop ("Argument 'theme' should be a 'theme' object.")

    }

    # standardise names for plotting
    sens <- sensors
    sens$xproj <- sens[[x_col]]
    sens$yproj <- sens[[y_col]]

    # identify highlighted sensors
    sens$.highlight <- FALSE

    if (!is.null(highlight) ) {

        if (is.numeric(highlight) ) {

            sens$.highlight[highlight] <- TRUE

        } else {

            if (!label_col %in% names(sensors) ) {

                stop ("Highlighting by name requires `label_col`.", call. = FALSE)

            }

            sens$.highlight <- sens[[label_col]] %in% highlight

        }

    }

    # head outline
    head_df <- NULL
    r_head <- NULL

    if (show_head) {

        r_head <- st_head_radius(
            sensors = data.frame(xproj = sens$xproj, yproj = sens$yproj),
            head_expand = head_expand
            )

        head_df <- st_head_outline(r_head)

    }

    p <- ggplot2::ggplot() +
        ggplot2::coord_equal() +
        theme

    if (show_head) {

        p <- p + ggplot2::geom_path(
            data = head_df,
            ggplot2::aes(x = .data$x, y = .data$y, group = .data$part),
            inherit.aes = FALSE,
            linewidth = 0.5
            )

    }

    if (show_points) {

        p <- p + ggplot2::geom_point(
            data = sens[!sens$.highlight, ],
            ggplot2::aes(x = .data$xproj, y = .data$yproj),
            inherit.aes = FALSE,
            size = point_size,
            alpha = if (dim_others) other_alpha else 1
            )

    }

    # highlighted sensors
    if (any(sens$.highlight) ) {

        p <- p + ggplot2::geom_point(
            data = sens[sens$.highlight, ],
            ggplot2::aes(x = .data$xproj, y = .data$yproj),
            inherit.aes = FALSE,
            size = highlight_size,
            color = highlight_col,
            shape = highlight_shape,
            alpha = highlight_alpha
            )

    }

    # labels
    if (show_labels) {

        lab_df <- sens

        if (label_only_highlight && any(sens$.highlight) ) {

            lab_df <- sens[sens$.highlight, ]

        }

        lab_df$label <- as.character(lab_df[[label_col]])

        if (label_repel) {

            p <- p + ggrepel::geom_text_repel(
                data = lab_df,
                ggplot2::aes(x = .data$xproj, y = .data$yproj, label = .data$label),
                inherit.aes = FALSE,
                size = label_size,
                max.overlaps = Inf,
                min.segment.length = 0
                )

        } else {

            p <- p + ggplot2::geom_text(
                data = lab_df,
                ggplot2::aes(x = .data$xproj, y = .data$yproj, label = .data$label),
                inherit.aes = FALSE,
                size = label_size,
                vjust = -0.6
                )

        }

    }

    # sensible limits
    if (is.null(xlim) || is.null(ylim) ) {

        if (show_head && !is.null(r_head) ) {

            if (is.null(xlim)) xlim <- c(-r_head, r_head)
            if (is.null(ylim)) ylim <- c(-r_head, 1.15 * r_head)

        } else {

            if (is.null(xlim) ) xlim <- range(sens$xproj, na.rm = TRUE)
            if (is.null(ylim) ) ylim <- range(sens$yproj, na.rm = TRUE)

        }

    }

    final_plot <- p + ggplot2::xlim(xlim) + ggplot2::ylim(ylim)

    return (final_plot)

}
