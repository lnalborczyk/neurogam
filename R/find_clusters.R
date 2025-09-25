#' Finding clusters (groups of timesteps) in time-series
#'
#' Finding clusters (groups of timesteps) in time-series that exceeds some threshold.
#'
#' @param data Dataframe, containing exactly two columns called "time" and "value" (in long format).
#' @param threshold Numeric, indicates the threshold for finding clusters.
#' @param above_threshold Logical, indicates whether we should threshold above (or below) value.
#'
#' @return The identified clusters (i.e., onset and offset).
#'
#' @importFrom rlang .data
#'
#' @author Ladislas Nalborczyk \email{ladislas.nalborczyk@@cnrs.fr}.
#'
#' @export

find_clusters <- function (data, threshold = 10, above_threshold = TRUE) {

    stopifnot("data must be a dataframe..." = is.data.frame(data) )
    stopifnot("threshold must be a numeric..." = is.numeric(threshold) )
    stopifnot("above_threshold must be a logical..." = is.logical(above_threshold) )

    required_columns <- c("time", "value")
    assertthat::assert_that(
        identical(sort(colnames(data) ), sort(required_columns) ),
        msg = paste(
            "Data must have exactly two columns named 'time' and 'value'. Found:",
            paste(colnames(data), collapse = ", ")
            )
        )

    if (isFALSE(above_threshold) ) {

        threshold <- -threshold
        data <- data |> dplyr::mutate(value = if (above_threshold) .data$value else -.data$value)

    }

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
            .groups = "drop"
            ) |>
        data.frame()

    return (clusters)

}
