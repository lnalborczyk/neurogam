#' Find connected threshold-exceeding clusters on arbitrary coordinates
#'
#' @param data Data frame containing unique coordinate combinations and a statistic.
#' @param coordinates Character vector of regular-grid coordinate names. For
#'   sensor-time data, normally `"time"` (sensor topology is separate).
#' @param value Column containing positive posterior odds (default `"value"`).
#' @param threshold Posterior-odds threshold, strictly greater than one.
#' @param threshold_type `"both"`, `"above"`, or `"below"`.
#' @param group Optional grouping column; connected components never cross groups.
#' @param sensor_id Optional column with sensor identifiers.
#' @param sensor_adjacency Named logical/numeric matrix with sensor row/column
#'   names, or two-column data frame with sensor pairs.
#' @param max_gap Optional named list of maximum allowed differences for numeric
#'   grid coordinates (e.g. `list(time = 0.004)`). By default adjacent observed
#'   coordinate levels connect, even when physically far apart.
#' @return List with `clusters` (one row per component), `membership` (one
#'   row per threshold-exceeding sample), and `edges` (full-domain graph edges).
#'
#' @author Ladislas Nalborczyk \email{ladislas.nalborczyk@@cnrs.fr}.
#'
#' @export
find_clusters_nd <- function (
        data, coordinates = "time", value = "value",
        threshold = 10, threshold_type = c("both", "above", "below"),
        group = NULL, sensor_id = NULL, sensor_adjacency = NULL, max_gap = NULL
        ) {

    threshold_type <- match.arg(threshold_type)
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Install 'igraph' to use find_clusters_nd().", call. = FALSE)
  }
  if (!is.data.frame(data)) stop("`data` must be a data frame.")
  if (length(threshold) != 1L || !is.numeric(threshold) ||
      !is.finite(threshold) || threshold <= 1) {
    stop("`threshold` must be a finite posterior-odds threshold > 1.")
  }
  cols <- c(coordinates, sensor_id, group, value)
  if (!length(coordinates) && is.null(sensor_id)) stop("Supply coordinates or sensor_id.")
  if (anyNA(cols) || anyDuplicated(cols) || !all(cols %in% names(data))) {
    stop("Coordinate, sensor, group and value columns must exist and be distinct.")
  }
  if (!is.numeric(data[[value]]) || anyNA(data[[value]]) ||
      any(data[[value]] < 0)) {
    stop("`value` must contain nonnegative posterior odds without missing values.")
  }
  if (anyNA(data[c(coordinates, sensor_id, group)])) {
    stop("Missing coordinate, sensor or group identifiers are not supported.")
  }
  if (!is.null(sensor_id) && is.null(sensor_adjacency)) {
    stop("Provide `sensor_adjacency` when `sensor_id` is supplied.")
  }
  if (!is.null(max_gap)) {
    if (!is.list(max_gap) || is.null(names(max_gap)) ||
        !all(names(max_gap) %in% coordinates) ||
        any(!vapply(max_gap, function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0, logical(1)))) {
      stop("`max_gap` must be a named list of positive numeric tolerances.")
    }
  }
  all_outputs <- list()
  groups <- if (is.null(group)) list(data) else split(data, data[[group]], drop = TRUE)
  for (g in seq_along(groups)) {
    d <- as.data.frame(groups[[g]])
    rownames(d) <- NULL
    domain_cols <- c(coordinates, sensor_id)
    if (anyDuplicated(d[domain_cols])) stop("Duplicate coordinate combinations within a group.")
    n <- nrow(d)
    if (!n) next
    # Map each sample to integer indices on each axis, independent of row order.
    axes <- lapply(coordinates, function(nm) sort(unique(d[[nm]])))
    names(axes) <- coordinates
    indices <- lapply(coordinates, function(nm) match(d[[nm]], axes[[nm]]))
    names(indices) <- coordinates
    if (!is.null(sensor_id)) {
      sensors <- sort(unique(as.character(d[[sensor_id]])))
      indices[[sensor_id]] <- match(as.character(d[[sensor_id]]), sensors)
    }
    idx <- as.data.frame(indices, check.names = FALSE)
    keys <- do.call(paste, c(unname(idx), sep = "|"))
    if (anyDuplicated(keys)) stop("Internal duplicate grid keys.")
    make_keys <- function(x) do.call(paste, c(unname(x), sep = "|"))
    edge_from <- integer(0)
    edge_to <- integer(0)
    # One-step neighbors along each regular coordinate axis (von Neumann topology).
    for (nm in coordinates) {
      source <- which(idx[[nm]] < length(axes[[nm]]))
      if (!length(source)) next
      if (!is.null(max_gap[[nm]])) {
        levels <- axes[[nm]]
        if (!is.numeric(levels)) stop("max_gap requires a numeric axis: ", nm)
        keep <- (levels[idx[[nm]][source] + 1L] - levels[idx[[nm]][source]]) <= max_gap[[nm]]
        source <- source[keep]
      }
      if (!length(source)) next
      target_idx <- idx[source, , drop = FALSE]
      target_idx[[nm]] <- target_idx[[nm]] + 1L
      target <- match(make_keys(target_idx), keys)
      valid <- !is.na(target)
      edge_from <- c(edge_from, source[valid])
      edge_to <- c(edge_to, target[valid])
    }
    # Edges between sensor neighbors at precisely identical other coordinates.
    if (!is.null(sensor_id)) {
      if (is.matrix(sensor_adjacency)) {
        a <- sensor_adjacency
        if (is.null(rownames(a)) || is.null(colnames(a)) ||
            anyDuplicated(rownames(a)) || anyDuplicated(colnames(a)) ||
            !setequal(rownames(a), colnames(a)) ||
            !all(sensors %in% rownames(a))) {
          stop("Sensor adjacency matrix must have matching unique row/column names covering all sensors.")
        }
        a <- a[sensors, sensors, drop = FALSE]
        if (anyNA(a) || !isTRUE(all.equal(unname(a), unname(t(a))))) {
          stop("Sensor adjacency matrix must be symmetric and nonmissing.")
        }
        diag(a) <- 0
        pairs <- which(upper.tri(a) & a != 0, arr.ind = TRUE)
        pair_names <- if (nrow(pairs)) cbind(sensors[pairs[, 1]], sensors[pairs[, 2]]) else matrix(character(), ncol = 2)
      } else if (is.data.frame(sensor_adjacency) && ncol(sensor_adjacency) == 2L) {
        pair_names <- as.matrix(sensor_adjacency)
        if (anyNA(pair_names) || !all(pair_names %in% sensors)) {
          stop("Sensor edge list references missing or unknown sensors.")
        }
        pair_names <- matrix(as.character(pair_names), ncol = 2)
      } else stop("`sensor_adjacency` must be a named matrix or two-column edge list.")
      for (i in seq_len(nrow(pair_names))) {
        source <- which(as.character(d[[sensor_id]]) == pair_names[i, 1])
        if (!length(source)) next
        target_idx <- idx[source, , drop = FALSE]
        target_idx[[sensor_id]] <- match(pair_names[i, 2], sensors)
        target <- match(make_keys(target_idx), keys)
        valid <- !is.na(target)
        edge_from <- c(edge_from, source[valid])
        edge_to <- c(edge_to, target[valid])
      }
    }
    edges <- data.frame(from = edge_from, to = edge_to)
    graph <- igraph::make_empty_graph(n = n, directed = FALSE)
    if (nrow(edges)) graph <- igraph::add_edges(graph, as.vector(t(as.matrix(edges))))
    # Calculate components separately for each sign, never mixing them.
    signs <- switch(threshold_type, both = c("positive", "negative"),
                    above = "positive", below = "negative")
    cluster_rows <- list()
    member_rows <- list()
    next_id <- 0L
    for (sign in signs) {
      hit <- if (sign == "positive") d[[value]] >= threshold else d[[value]] <= 1 / threshold
      active <- which(hit)
      if (!length(active)) next
      subgraph <- igraph::induced_subgraph(graph, vids = active)
      comp <- igraph::components(subgraph)$membership
      for (c_id in sort(unique(comp))) {
        members <- active[which(comp == c_id)]
        next_id <- next_id + 1L
        m <- d[members, , drop = FALSE]
        m$id <- next_id
        m$sign <- sign
        member_rows[[length(member_rows) + 1L]] <- m
        summary <- data.frame(id = next_id, sign = sign,
                              n_points = length(members), stringsAsFactors = FALSE)
        if (length(coordinates)) {
          for (nm in coordinates) {
            summary[[paste0(nm, "_min")]] <- min(m[[nm]])
            summary[[paste0(nm, "_max")]] <- max(m[[nm]])
          }
        }
        if (!is.null(sensor_id)) summary$n_sensors <- length(unique(m[[sensor_id]]))
        if ("time" %in% coordinates && is.numeric(m$time)) {
          summary$onset <- min(m$time)
          summary$offset <- max(m$time)
          summary$duration <- summary$offset - summary$onset
        }
        cluster_rows[[length(cluster_rows) + 1L]] <- summary
      }
    }
    empty_members <- d[FALSE, , drop = FALSE]
    empty_members$id <- integer()
    empty_members$sign <- character()
    empty_clusters <- data.frame(id = integer(), sign = character(), n_points = integer())
    # Add domain-specific columns to zero-row result for schema consistency.
    for (nm in coordinates) {
      empty_clusters[[paste0(nm, "_min")]] <- d[[nm]][FALSE]
      empty_clusters[[paste0(nm, "_max")]] <- d[[nm]][FALSE]
    }
    if (!is.null(sensor_id)) empty_clusters$n_sensors <- integer()
    if ("time" %in% coordinates && is.numeric(d$time)) {
      empty_clusters$onset <- numeric()
      empty_clusters$offset <- numeric()
      empty_clusters$duration <- numeric()
    }
    cl <- if (length(cluster_rows)) do.call(rbind, cluster_rows) else empty_clusters
    mb <- if (length(member_rows)) do.call(rbind, member_rows) else empty_members
    if (!is.null(group)) cl[[group]] <- rep(d[[group]][1], nrow(cl))
    all_outputs[[length(all_outputs) + 1L]] <- list(clusters = cl, membership = mb, edges = edges)
  }
  if (!length(all_outputs)) stop("Input must contain at least one observation.")
  combine <- function(nm) {
    z <- do.call(rbind, lapply(all_outputs, `[[`, nm))
    rownames(z) <- NULL
    z
  }
  list(clusters = combine("clusters"), membership = combine("membership"),
       edges = if (length(all_outputs) == 1L) all_outputs[[1]]$edges else NULL)
}
