#' Find contiguous clusters in a time series
#'
#' Identify contiguous clusters of time points in a time series where a
#' variable exceeds a positive threshold, falls below a negative threshold,
#' or both. Clusters are defined as consecutive time points satisfying the
#' thresholding condition.
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
#'   \code{value <= 1/threshold}. Must be non-negative when
#'   \code{threshold_type = "both"}.
#' @param group Optional grouping column name (character scalar) used to find
#'   clusters independently within each group level (e.g.,
#'   \code{"participant"}). Set to \code{NULL} (default) to ignore grouping.
#' @param threshold_type Character scalar controlling which clusters are
#'   detected. Must be one of \code{"above"}, \code{"below"}, or \code{"both"}
#'   (default). When \code{"above"}, clusters are formed where
#'   \code{value >= threshold}. When \code{"below"}, clusters are formed where
#'   \code{value <= 1/threshold}. When \code{"both"}, both types are detected
#'   and the returned data include a \code{sign} column.
#'
#' @return A data frame with one row per detected cluster and columns:
#'   \itemize{
#'     \item \code{id}: integer cluster index (starting at 1). When
#'       \code{group} is provided, \code{id} restarts at 1 within each
#'       group level;
#'     \item \code{onset}: time of the first point in the cluster;
#'     \item \code{offset}: time of the last point in the cluster;
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
#' find_clusters(data = df, threshold = 3, threshold_type = "both")
#'
#' # One-sided detection (positive only)
#' find_clusters(data = df, threshold = 3, threshold_type = "above")
#'
#' # One-sided detection (negative only)
#' find_clusters(data = df, threshold = 3, threshold_type = "below")
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
#'   group = "participant",
#'   threshold_type = "both"
#'   )
#'
#' @export
find_clusters <- function (
        data,
        threshold = 10,
        group = NULL,
        threshold_type = c("both", "above", "below")
        ) {

    stopifnot("`data` must be a data frame." = is.data.frame(data) )
    threshold_type <- match.arg(threshold_type)

    if (!is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold) ) {

        stop ("`threshold` must be a finite numeric scalar.", call. = FALSE)

    }

    if (threshold_type == "both" && threshold < 0) {

        stop ("`threshold` must be >= 0 when `threshold_type = \"both\"`.", call. = FALSE)

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

        stop (
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

    compute_clusters <- function (dat, direction = c("positive", "negative"), sign_label) {

        direction <- match.arg(direction)

        if (is.null(group) ) {

            out <- dat |>
                dplyr::mutate(
                    hit = if (direction == "positive") {
                        .data$value >= threshold
                    } else {
                        .data$value <= 1 / threshold
                    },
                    change = dplyr::lag(.data$hit, default = FALSE) != .data$hit,
                    id = cumsum(.data$change & .data$hit)
                    ) |>
                dplyr::filter(.data$hit) |>
                dplyr::group_by(.data$id) |>
                dplyr::summarise(
                    onset = dplyr::first(.data$time),
                    offset = dplyr::last(.data$time),
                    n_points = dplyr::n(),
                    .groups = "drop"
                    ) |>
                dplyr::mutate(sign = sign_label) |>
                data.frame()

        } else {

            out <- dat |>
                dplyr::group_by(.data[[group]]) |>
                dplyr::mutate(
                    hit = if (direction == "positive") {
                        .data$value >= threshold
                    } else {
                        .data$value <= 1 / threshold
                    },
                    change = dplyr::lag(.data$hit, default = FALSE) != .data$hit,
                    id = cumsum(.data$change & .data$hit)
                    ) |>
                dplyr::filter(.data$hit) |>
                dplyr::group_by(.data[[group]], .data$id) |>
                dplyr::summarise(
                    onset = dplyr::first(.data$time),
                    offset = dplyr::last(.data$time),
                    n_points = dplyr::n(),
                    .groups = "drop"
                    ) |>
                dplyr::mutate(sign = sign_label) |>
                data.frame()

        }

        return (out)

    }

    clusters <- switch (
        threshold_type,
        both = {
            out <- dplyr::bind_rows(
                compute_clusters(dat = data, direction = "positive", sign_label = "positive"),
                compute_clusters(dat = data, direction = "negative", sign_label = "negative")
                )
            if (!is.null(group) ) out <- dplyr::arrange(out, .data[[group]])
            out
        },
        above = {
            out <- compute_clusters(dat = data, direction = "positive", sign_label = "positive")
            if (!is.null(group) ) out <- dplyr::arrange(out, .data[[group]])
            out
        },
        below = {
            out <- compute_clusters(dat = data, direction = "negative", sign_label = "negative")
            if (!is.null(group) ) out <- dplyr::arrange(out, .data[[group]])
            out
        }

    )

    return (clusters)

}
