#' Finding clusters (groups of timesteps) in time-series
#'
#' Finding clusters (groups of timesteps) in time-series that exceeds some threshold.
#'
#' @param data Dataframe, with reverse correlation data (in long format).
#' @param threshold Numeric/Character/Factor, column in data specifying the participant ID.
#'
#' @return The original data augmented with the kernel.
#'
#' @importFrom rlang .data
#'
#' @author Ladislas Nalborczyk \email{ladislas.nalborczyk@@gmail.com}.
#'
#' @export

find_clusters <- function (data, threshold = 0.05) {

    data |>
        dplyr::pivot_longer(cols = -.data$time) |>
        dplyr::group_by(.data$name) |>
        dplyr::mutate(above_threshold = .data$value <= threshold) |>
        dplyr::mutate(cluster_change = c(TRUE, diff(.data$above_threshold) != 0) ) |>
        dplyr::filter(.data$above_threshold) |>
        dplyr::mutate(cluster_id = cumsum(.data$cluster_change) ) |>
        dplyr::group_by(.data$name, .data$cluster_id) |>
        dplyr::summarise(cluster_onset = first(.data$time), cluster_offset = last(.data$time) ) |>
        dplyr::select(.data$cluster_onset, .data$cluster_offset)

}
