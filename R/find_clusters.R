#' Finding clusters (groups of timesteps) in time-series
#'
#' Finding clusters (groups of timesteps) in time-series that exceeds some threshold.
#'
#' @param data Dataframe, containing exactly two columns called "time" and "value" (in long format).
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

    required_columns <- c("time", "value")
    assertthat::assert_that(
        identical(sort(colnames(data) ), sort(required_columns) ),
        msg = paste(
            "Data must have exactly two columns named 'time' and 'value'. Found:",
            paste(colnames(data), collapse = ", ")
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
        dplyr::select(-.data$name) |>
        data.frame()

    return (clusters)

}
