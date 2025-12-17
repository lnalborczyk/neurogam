#' Find contiguous clusters in a time series
#'
#' Identify contiguous clusters of time points in a time series where a
#' variable exceeds a positive threshold or falls below a negative threshold.
#' Clusters are defined as consecutive time points satisfying the thresholding
#' condition.
#'
#' By default, the function detects \strong{both positive and negative clusters}
#' in a single call and returns a column indicating the cluster sign.
#'
#' If a grouping variable is provided (e.g., \code{"participant"}), clusters
#' are detected independently within each group level.
#'
#' @param data A data frame containing at least two columns named
#'   \code{"time"} and \code{"value"}. If \code{group} is not \code{NULL},
#'   \code{data} must also contain that grouping column.
#' @param threshold Numeric scalar specifying the (positive) threshold used to
#'   define clusters. Positive clusters are defined where
#'   \code{value >= threshold}, and negative clusters where
#'   \code{value <= -threshold}. Must be non-negative when
#'   \code{both_signs = TRUE}.
#' @param group Optional grouping column name (character scalar) used to find
#'   clusters independently within each group level (e.g.,
#'   \code{"participant"}). Set to \code{NULL} (default) to ignore grouping.
#' @param above_threshold Logical scalar used only when
#'   \code{both_signs = FALSE}. If \code{TRUE} (default), clusters are formed
#'   where \code{value >= threshold}; if \code{FALSE}, clusters are formed
#'   where \code{value <= threshold}.
#' @param both_signs Logical scalar indicating whether to detect both positive
#'   and negative clusters in a single call (default: \code{TRUE}). When
#'   \code{TRUE}, \code{above_threshold} is ignored.
#'
#' @return A data frame with one row per detected cluster and columns:
#'   \itemize{
#'     \item \code{cluster_id}: integer cluster index (starting at 1). When
#'       \code{group} is provided, \code{cluster_id} restarts at 1 within each
#'       group level;
#'     \item \code{cluster_onset}: time of the first point in the cluster;
#'     \item \code{cluster_offset}: time of the last point in the cluster;
#'     \item \code{n_points}: number of time points in the cluster;
#'     \item \code{sign}: character indicating cluster sign
#'       (\code{"positive"} or \code{"negative"}).
#'   }
#'
#'   If \code{group} is not \code{NULL}, the returned data frame also contains
#'   the grouping column (named as in \code{group}).
#'
#'   If no clusters are found, an empty data frame with the same column structure
#'   is returned.
#'
#' @details
#' The function assumes that the \code{time} variable is numeric and that
#' consecutive rows correspond to consecutive time points (within each group
#' if grouping is used). Internally, the data are:
#' \enumerate{
#'   \item filtered to remove rows with missing values;
#'   \item arranged by \code{time} (and by \code{group} then \code{time}, if used);
#'   \item thresholded to identify positive and/or negative excursions;
#'   \item segmented into runs of consecutive threshold-exceeding values,
#'         which define clusters.
#' }
#'
#' @importFrom rlang .data
#' @importFrom dplyr all_of
#'
#' @author Ladislas Nalborczyk \email{ladislas.nalborczyk@@cnrs.fr}
#'
#' @examples
#' set.seed(666)
#' df <- data.frame(
#'   time = seq(0, 1, length.out = 100),
#'   value = c(
#'     rnorm(30, 0, 1),
#'     rnorm(20,  4, 1), # positive cluster
#'     rnorm(20, -4, 1), # negative cluster
#'     rnorm(30, 0, 1)
#'     )
#'   )
#'
#' # Detect both positive and negative clusters
#' find_clusters(data = df, threshold = 3)
#'
#' # One-sided detection (positive only)
#' find_clusters(data = df, threshold = 3, both_signs = FALSE, above_threshold = TRUE)
#'
#' # Grouped example (e.g., per participant)
#' df_g <- rbind(
#'   transform(df, participant = "P01"),
#'   transform(df, participant = "P02")
#'   )
#'
#' find_clusters(
#'   data = df_g,
#'   threshold = 3,
#'   group = "participant"
#'   )
#'
#' @export
find_clusters <- function (
        data,
        threshold = 10,
        group = NULL,
        above_threshold = TRUE,
        both_signs = TRUE
        ) {

    stopifnot("`data` must be a data frame." = is.data.frame(data) )

    if (!is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold) ) {

        stop ("`threshold` must be a finite numeric scalar.", call. = FALSE)

    }

    if (both_signs && threshold < 0) {

        stop ("`threshold` must be >= 0 when `both_signs = TRUE`.", call. = FALSE)

    }

    if (!is.null(group) ) {

        if (!is.character(group) || length(group) != 1L || group == "") {

            stop ("`group` must be NULL or a single non-empty character string.", call. = FALSE)

        }

        if (!group %in% names(data) ) {

            stop ("Grouping column '", group, "' not found in `data`.", call. = FALSE)

        }

    }

    required_columns <- if (is.null(group) ) c("time", "value") else c("time", "value", group)

    missing_cols <- setdiff(required_columns, names(data) )

    if (length(missing_cols) > 0L) {

        stop(
            "Data is missing required column(s): ",
            paste(missing_cols, collapse = ", "),
            call. = FALSE
            )

    }

    data <- data |>
        dplyr::select(dplyr::all_of(required_columns) ) |>
        dplyr::filter(
            !is.na(.data$time),
            !is.na(.data$value),
            if (is.null(group) ) TRUE else !is.na(.data[[group]])
            )

    data <- if (is.null(group) ) {

        dplyr::arrange(data, .data$time)

    } else {

        dplyr::arrange(data, .data[[group]], .data$time)

    }

    # compute clusters given a predicate evaluated within the current data/group
    compute_clusters <- function (dat, direction = c("positive", "negative"), sign_label) {

        direction <- match.arg(direction)

        if (is.null(group) ) {

            out <- dat |>
                dplyr::mutate(
                    hit = if (direction == "positive") .data$value >=  threshold else .data$value <= -threshold,
                    change = dplyr::lag(.data$hit, default = FALSE) != .data$hit,
                    cluster_id = cumsum(.data$change & .data$hit)
                    ) |>
                dplyr::filter(.data$hit) |>
                dplyr::group_by(.data$cluster_id) |>
                dplyr::summarise(
                    cluster_onset = dplyr::first(.data$time),
                    cluster_offset = dplyr::last(.data$time),
                    n_points = dplyr::n(),
                    .groups = "drop"
                    ) |>
                dplyr::mutate(sign = sign_label) |>
                data.frame()

        } else {

            out <- dat |>
                dplyr::group_by(.data[[group]]) |>
                dplyr::mutate(
                    hit = if (direction == "positive") .data$value >=  threshold else .data$value <= -threshold,
                    change = dplyr::lag(.data$hit, default = FALSE) != .data$hit,
                    cluster_id = cumsum(.data$change & .data$hit)
                    ) |>
                dplyr::filter(.data$hit) |>
                dplyr::group_by(.data[[group]], .data$cluster_id) |>
                dplyr::summarise(
                    cluster_onset = dplyr::first(.data$time),
                    cluster_offset = dplyr::last(.data$time),
                    n_points = dplyr::n(),
                    .groups = "drop"
                    ) |>
                dplyr::mutate(sign = sign_label) |>
                data.frame()

        }

        return (out)

    }

    if (both_signs) {

        clusters <- dplyr::bind_rows(
            compute_clusters(data, direction = "positive", sign_label = "positive"),
            compute_clusters(data, direction = "negative", sign_label = "negative")
            )

    } else {

        if (!above_threshold) {

            data <- data |> dplyr::mutate(value = -.data$value)
            threshold <- -threshold

        }

        clusters <- if (is.null(group) ) {
            data |>
                dplyr::mutate(
                    hit = .data$value >= threshold,
                    change = dplyr::lag(.data$hit, default = FALSE) != .data$hit,
                    cluster_id = cumsum(.data$change & .data$hit)
                    ) |>
                dplyr::filter(.data$hit) |>
                dplyr::group_by(.data$cluster_id) |>
                dplyr::summarise(
                    cluster_onset = dplyr::first(.data$time),
                    cluster_offset = dplyr::last(.data$time),
                    n_points = dplyr::n(),
                    .groups = "drop"
                    ) |>
                dplyr::mutate(sign = if (above_threshold) "positive" else "negative") |>
                data.frame()

        } else {

            data |>
                dplyr::group_by(.data[[group]]) |>
                dplyr::mutate(
                    hit = .data$value >= threshold,
                    change = dplyr::lag(.data$hit, default = FALSE) != .data$hit,
                    cluster_id = cumsum(.data$change & .data$hit)
                    ) |>
                dplyr::filter(.data$hit) |>
                dplyr::group_by(.data[[group]], .data$cluster_id) |>
                dplyr::summarise(
                    cluster_onset = dplyr::first(.data$time),
                    cluster_offset = dplyr::last(.data$time),
                    n_points = dplyr::n(),
                    .groups = "drop"
                    ) |>
                dplyr::mutate(sign = if (above_threshold) "positive" else "negative") |>
                data.frame()

        }

    }

    return (clusters)

}
