#' Prepare a sensor coordinate table
#'
#' @param data Data frame with one or more rows per sensor.
#' @param sensor_id Name of the sensor identifier column.
#' @param spatial_id Names of two numeric projected spatial columns.
#'
#' @return One row per sensor, sorted by sensor identifier.
#' @keywords internal
prepare_sensor_positions <- function (
    data,
    sensor_id = "sensor",
    spatial_id = c("xproj", "yproj")
    ) {
    if (!is.data.frame(data)) {
        stop("`data` must be a data frame.", call. = FALSE)
    }

    if (!is.character(sensor_id) || length(sensor_id) != 1L ||
        is.na(sensor_id) || !nzchar(sensor_id) ||
        !is.character(spatial_id) || length(spatial_id) != 2L ||
        anyNA(spatial_id) || any(!nzchar(spatial_id)) ||
        anyDuplicated(c(sensor_id, spatial_id))) {
        stop("Supply one sensor column and two distinct spatial columns.",
             call. = FALSE)
    }

    columns <- c(sensor_id, spatial_id)
    if (!all(columns %in% names(data))) {
        stop(
            "`data` must contain: ", paste(columns, collapse = ", "),
            call. = FALSE
        )
    }
    if (!nrow(data)) {
        stop("`data` must contain at least one sensor.", call. = FALSE)
    }

    if (!is.character(data[[sensor_id]]) &&
        !is.factor(data[[sensor_id]]) &&
        !is.numeric(data[[sensor_id]])) {
        stop("The sensor identifier must be character, factor, or numeric.",
             call. = FALSE)
    }
    ids <- as.character(data[[sensor_id]])
    if (anyNA(ids) || any(!nzchar(trimws(ids)))) {
        stop("Sensor identifiers must be nonmissing and nonempty.",
             call. = FALSE)
    }

    if (!all(vapply(data[spatial_id], is.numeric, logical(1))) ||
        any(!is.finite(as.matrix(data[spatial_id])))) {
        stop("Spatial coordinates must be finite numeric values.",
             call. = FALSE)
    }

    positions <- data.frame(
        sensor_name = ids,
        x = data[[spatial_id[1L]]],
        y = data[[spatial_id[2L]]],
        stringsAsFactors = FALSE
    )
    positions <- unique(positions)
    if (anyDuplicated(positions$sensor_name)) {
        stop("Each sensor must have one consistent coordinate pair.",
             call. = FALSE)
    }
    if (anyDuplicated(positions[c("x", "y")])) {
        stop("Different sensors cannot share identical projected positions.",
             call. = FALSE)
    }

    positions <- positions[order(positions$sensor_name), , drop = FALSE]
    names(positions) <- c(sensor_id, spatial_id)
    rownames(positions) <- NULL
    positions
}

#' Construct approximate EEG sensor adjacency from projected positions
#'
#' Builds a symmetric, named logical adjacency matrix using spatial
#' k-nearest neighbors. To avoid connecting distant sensors, candidate edges
#' longer than `distance_factor` times the median nearest-neighbor distance
#' are omitted. The result is an *approximation* based on 2D positions;
#' inspect the graph before scientific use.
#'
#' @param data Data frame containing sensor names and projected locations.
#'   Repeated rows per sensor (e.g., a sensor-by-time prediction grid) are OK
#'   provided each sensor has a single consistent position.
#' @param sensor_id Sensor identifier column (default `"sensor"`).
#' @param spatial_id Two numeric spatial coordinate columns.
#' @param k Number of nearest neighbors considered per electrode. The
#'   resulting undirected graph may have more than `k` edges per sensor.
#' @param distance_factor Maximum candidate edge length divided by the
#'   median nearest-neighbor distance. Use `NULL` to disable filtering.
#' @param warn_disconnected Issue a warning if the inferred graph has
#'   more than one connected component.
#'
#' @return Symmetric logical matrix with sensor names on both axes and
#'   `FALSE` on the diagonal.
#' @export
make_sensor_adjacency <- function(
    data,
    sensor_id = "sensor",
    spatial_id = c("xproj", "yproj"),
    k = 4L,
    distance_factor = 1.8,
    warn_disconnected = TRUE
) {
    positions <- prepare_sensor_positions(
        data,
        sensor_id = sensor_id,
        spatial_id = spatial_id
    )

    if (!is.numeric(k) || length(k) != 1L || is.na(k) ||
        !is.finite(k) || k < 1 || k != floor(k)) {
        stop("`k` must be a positive integer.", call. = FALSE)
    }
    if (!is.null(distance_factor) &&
        (!is.numeric(distance_factor) || length(distance_factor) != 1L ||
         is.na(distance_factor) || !is.finite(distance_factor) ||
         distance_factor <= 0)) {
        stop("`distance_factor` must be NULL or a positive finite number.",
             call. = FALSE)
    }
    if (!is.logical(warn_disconnected) ||
        length(warn_disconnected) != 1L || is.na(warn_disconnected)) {
        stop("`warn_disconnected` must be TRUE or FALSE.", call. = FALSE)
    }

    ids <- as.character(positions[[sensor_id]])
    n <- length(ids)
    adjacency <- matrix(
        FALSE,
        nrow = n,
        ncol = n,
        dimnames = list(ids, ids)
    )
    if (n == 1L) {
        return(adjacency)
    }

    xy <- as.matrix(positions[spatial_id])
    distances <- as.matrix(stats::dist(xy))
    diag(distances) <- Inf

    nearest_distance <- apply(distances, 1L, min)
    max_distance <- if (is.null(distance_factor)) {
        Inf
    } else {
        distance_factor * stats::median(nearest_distance)
    }

    for (i in seq_len(n)) {
        neighbors <- order(distances[i, ], seq_len(n))
        neighbors <- neighbors[seq_len(min(as.integer(k), n - 1L))]
        neighbors <- neighbors[distances[i, neighbors] <= max_distance]
        if (length(neighbors)) {
            adjacency[i, neighbors] <- TRUE
        }
    }

    # kNN need not be symmetric; make all selected pairs undirected.
    adjacency <- adjacency | t(adjacency)
    diag(adjacency) <- FALSE

    if (warn_disconnected && requireNamespace("igraph", quietly = TRUE)) {
        n_components <- igraph::components(
            igraph::graph_from_adjacency_matrix(adjacency, mode = "undirected")
        )$no
        if (n_components > 1L) {
            warning(
                "Automatic sensor adjacency has ", n_components,
                " disconnected components. Inspect the montage and consider ",
                "increasing `distance_factor`, adjusting `k`, or providing ",
                "an explicit adjacency matrix.",
                call. = FALSE
            )
        }
    }

    return (adjacency)

}
