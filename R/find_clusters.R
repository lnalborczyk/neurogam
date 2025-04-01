#' Finding clusters (groups of timesteps) in time-series
#'
#' Finding clusters (groups of timesteps) in time-series that exceeds some threshold.
#'
#' @param data Dataframe, containing time-series of posterior probability ratios (in long format).
#' @param threshold Numeric, indicates the threshold for finding clusters.
#'
#' @return The identified clusters (i.e., onset and offset).
#'
#' @importFrom rlang .data
#'
#' @author Ladislas Nalborczyk \email{ladislas.nalborczyk@@gmail.com}.
#'
#' @export

find_clusters <- function (data, threshold) {

    # checking required column names
    required_columns <- c("time", "prob_ratio")
    assertthat::assert_that(
        all(required_columns %in% colnames(data) ),
        msg = paste(
            "Missing columns:",
            paste(setdiff(required_columns, colnames(data) ), collapse = ", ")
            )
        )

    clusters <- data |>
        tidyr::pivot_longer(cols = -.data$time) |>
        dplyr::group_by(.data$name) |>
        dplyr::mutate(above_threshold = .data$value >= threshold) |>
        dplyr::mutate(cluster_change = c(TRUE, diff(.data$above_threshold) != 0) ) |>
        dplyr::filter(.data$above_threshold) |>
        dplyr::mutate(cluster_id = cumsum(.data$cluster_change) ) |>
        dplyr::group_by(.data$name, .data$cluster_id) |>
        dplyr::summarise(cluster_onset = dplyr::first(.data$time), cluster_offset = dplyr::last(.data$time) ) |>
        dplyr::ungroup() |>
        # dplyr::select(.data$cluster_onset, .data$cluster_offset) |>
        dplyr::select(-.data$name) |>
        data.frame()

    return (clusters)

}
