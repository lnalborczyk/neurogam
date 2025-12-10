#' Find contiguous clusters in a time series
#'
#' Identify contiguous clusters of time points in a time series where a
#' variable exceeds (or falls below) a given threshold. Clusters are defined
#' as consecutive time points where the thresholding condition holds.
#'
#' @param data A data frame containing at least two columns named
#'   \code{"time"} and \code{"value"}. Additional columns are ignored.
#' @param threshold Numeric scalar; threshold used to define clusters.
#'   When \code{above_threshold = TRUE}, clusters are defined where
#'   \code{value >= threshold}. When \code{above_threshold = FALSE},
#'   clusters are defined where \code{value <= threshold}.
#' @param above_threshold Logical scalar; if \code{TRUE} (default), clusters
#'   are formed where \code{value >= threshold}. If \code{FALSE}, clusters
#'   are formed where \code{value <= threshold}.
#'
#' @return A data frame with one row per detected cluster and columns:
#'   \itemize{
#'     \item \code{cluster_id}: integer cluster index (starting at 1);
#'     \item \code{cluster_onset}: time of the first point in the cluster;
#'     \item \code{cluster_offset}: time of the last point in the cluster;
#'     \item \code{n_points}: number of time points in the cluster.
#'   }
#'   If no clusters are found, an empty data frame with these columns is
#'   returned.
#'
#' @details
#' The function assumes that the \code{time} variable is numeric and that
#' consecutive rows correspond to consecutive time points. Internally, the
#' data are first arranged by \code{time} and any rows with \code{NA} in
#' \code{time} or \code{value} are removed.
#'
#' @importFrom rlang .data
#'
#' @author Ladislas Nalborczyk \email{ladislas.nalborczyk@@cnrs.fr}
#'
#' @examples
#' \dontrun{
#' set.seed(666)
#' df <- data.frame(
#'   time  = seq(0, 1, length.out = 100),
#'   value = c(rnorm(40, 0, 1), rnorm(20, 5, 1), rnorm(40, 0, 1))
#'   )
#'
#' # find clusters where value >= 3
#' cl <- find_clusters(data = df, threshold = 3, above_threshold = TRUE)
#' cl
#' }
#'
#' @export

find_clusters <- function (data, threshold = 10, above_threshold = TRUE) {

    # basic checks
    stopifnot("`data` must be a data frame." = is.data.frame(data) )

    if (!is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold) ) {

        stop ("`threshold` must be a finite numeric scalar.", call. = FALSE)

    }

    if (!is.logical(above_threshold) || length(above_threshold) != 1L || is.na(above_threshold) ) {

        stop ("`above_threshold` must be a single non-NA logical value.", call. = FALSE)

    }

    # required columns
    required_columns <- c("time", "value")
    assertthat::assert_that(
        identical(sort(colnames(data) ), sort(required_columns) ),
        msg = paste(
            "Data must have exactly two columns named 'time' and 'value'. Found:",
            paste(colnames(data), collapse = ", ")
            )
        )

    # if (isFALSE(above_threshold) ) {
    #
    #     threshold <- -threshold
    #     data <- data |> dplyr::mutate(value = if (above_threshold) .data$value else -.data$value)
    #
    # }

    # keep only what we need and ensure proper ordering / no NAs
    data <- data |>
        dplyr::select(.data$time, .data$value) |>
        dplyr::arrange(.data$time) |>
        dplyr::filter(!is.na(.data$time), !is.na(.data$value) )

    # handle "below threshold" case by flipping sign
    # if we want value <= threshold, we can equivalently look for (-value) >= (-threshold).
    if (!above_threshold) {

        data <- data |> dplyr::mutate(value = -.data$value)
        threshold <- -threshold

    }

    # clusters <- data |>
    #     dplyr::mutate(
    #         above = .data$value >= threshold,
    #         change = dplyr::lag(.data$above, default = FALSE) != .data$above,
    #         cluster_id = cumsum(.data$change & .data$above)
    #         ) |>
    #     dplyr::filter(.data$above) |>
    #     dplyr::group_by(.data$cluster_id) |>
    #     dplyr::summarise(
    #         cluster_onset = dplyr::first(.data$time),
    #         cluster_offset = dplyr::last(.data$time),
    #         .groups = "drop"
    #         ) |>
    #     data.frame()

    # identify contiguous runs above threshold
    clusters <- data |>
        dplyr::mutate(
            above = .data$value >= threshold,
            change = dplyr::lag(.data$above, default = FALSE) != .data$above,
            cluster_id = cumsum(.data$change & .data$above)
            ) |>
        dplyr::filter(.data$above) |>
        dplyr::group_by(.data$cluster_id) |>
        dplyr::summarise(
            cluster_onset = dplyr::first(.data$time),
            cluster_offset = dplyr::last(.data$time),
            n_points = dplyr::n(),
            .groups = "drop"
            ) |>
        # in case no clusters found, ensure consistent structure
        data.frame()

    return (clusters)

}
