#' Plot EEG sensor adjacency with ggplot2
#'
#' Displays electrode locations and the edges of a supplied adjacency matrix.
#' Sensor coordinates are matched to the matrix by electrode name, so the
#' order of rows in `data` need not match the matrix order. The returned plot
#' can be customized with regular ggplot2 layers and themes.
#'
#' @param data Sensor coordinate table or a sensor-by-time data frame. Repeated
#'   sensor rows are allowed when their coordinates are consistent.
#' @param adjacency Symmetric, named logical or numeric adjacency matrix.
#'   Nonzero off-diagonal entries indicate edges. The row and column names
#'   must contain exactly the sensor identifiers in `data`.
#' @param sensor_id Name of the sensor identifier column.
#' @param spatial_id Character vector of two projected coordinate columns.
#' @param show_labels Whether to display sensor names.
#' @param edge_colour Colour of the adjacency edges.
#' @param edge_linewidth Width of adjacency edges in mm.
#' @param point_colour Colour of sensor markers.
#' @param point_size Size of sensor markers.
#' @param label_size Size of sensor labels in ggplot2 text units.
#' @param check_overlap If `TRUE`, suppress overlapping sensor labels using
#'   ggplot2's `geom_text()` overlap detection.
#' @param show_axes Whether to show coordinate axes and their labels.
#'
#' @return A ggplot2 plot showing sensor positions and adjacency edges.
#' @examples
#' sensors <- data.frame(
#'     sensor = c("Fz", "Cz", "Pz"),
#'     xproj = c(0, 0, 0),
#'     yproj = c(1, 0, -1)
#' )
#' adjacency <- make_sensor_adjacency(sensors, k = 1, distance_factor = 1.5)
#' plot_sensor_adjacency(sensors, adjacency)
#' @importFrom rlang .data
#' @export
plot_sensor_adjacency <- function (
    data,
    adjacency,
    sensor_id = "sensor",
    spatial_id = c("xproj", "yproj"),
    show_labels = TRUE,
    edge_colour = "grey65",
    edge_linewidth = 0.6,
    point_colour = "#264E70",
    point_size = 2.4,
    label_size = 3,
    check_overlap = TRUE,
    show_axes = TRUE
    ) {

    positions <- prepare_sensor_positions(
        data,
        sensor_id = sensor_id,
        spatial_id = spatial_id
        )

    sensors <- as.character(positions[[sensor_id]])

    if (!is.matrix(adjacency) ||
        !(is.logical(adjacency) || is.numeric(adjacency) ) ||
        is.null(rownames(adjacency) ) ||
        is.null(colnames(adjacency) ) ||
        anyDuplicated(rownames(adjacency) ) ||
        anyDuplicated(colnames(adjacency) ) ||
        !setequal(rownames(adjacency), sensors) ||
        !setequal(colnames(adjacency), sensors) ) {

        stop (
            "`adjacency` must be a named logical or numeric matrix ",
            "covering all sensors.",
            call. = FALSE
            )

    }

    adjacency <- adjacency[sensors, sensors, drop = FALSE]

    if (anyNA(adjacency) ||
        (is.numeric(adjacency) && any(!is.finite(adjacency) ) ) ||
        !isTRUE(all.equal(unname(adjacency), unname(t(adjacency) ) ) ) ) {

        stop (
            "`adjacency` must be symmetric with finite, nonmissing entries.",
            call. = FALSE
            )

    }

    position_data <- data.frame(
        sensor = sensors,
        x = positions[[spatial_id[1L]]],
        y = positions[[spatial_id[2L]]],
        stringsAsFactors = FALSE
        )

    # Only the upper triangle is used, so every undirected edge appears once.
    edge_indices <- which(
        upper.tri(adjacency) & adjacency != 0,
        arr.ind = TRUE
        )

    edge_data <- data.frame(
        x = position_data$x[edge_indices[, 1L]],
        y = position_data$y[edge_indices[, 1L]],
        xend = position_data$x[edge_indices[, 2L]],
        yend = position_data$y[edge_indices[, 2L]]
        )

    plot <- ggplot2::ggplot() +
        ggplot2::labs(
            title = "EEG sensor adjacency",
            x = spatial_id[1L],
            y = spatial_id[2L]
            ) +
        ggplot2::coord_equal() +
        ggplot2::theme(panel.grid = ggplot2::element_blank() )

    if (nrow(edge_data) > 0L) {

        plot <- plot +
            ggplot2::geom_segment(
                data = edge_data,
                mapping = ggplot2::aes(
                    x = .data$x,
                    y = .data$y,
                    xend = .data$xend,
                    yend = .data$yend
                    ),
                colour = edge_colour,
                linewidth = edge_linewidth
                )

    }

    plot <- plot +
        ggplot2::geom_point(
            data = position_data,
            mapping = ggplot2::aes(x = .data$x, y = .data$y),
            colour = point_colour,
            size = point_size
            )

    if (isTRUE(show_labels) ) {

        plot <- plot +
            ggplot2::geom_text(
                data = position_data,
                mapping = ggplot2::aes(
                    x = .data$x,
                    y = .data$y,
                    label = .data$sensor
                    ),
                vjust = -0.8,
                size = label_size,
                check_overlap = check_overlap
                )

    }

    if (!isTRUE(show_axes) ) {

        plot <- plot +
            ggplot2::theme(
                axis.title = ggplot2::element_blank(),
                axis.text = ggplot2::element_blank(),
                axis.ticks = ggplot2::element_blank()
                )

    }

    plot

}
