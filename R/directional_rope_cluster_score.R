#' Construct a temporary score for directional ROPE clustering
#'
#' @keywords internal
directional_rope_cluster_score <- function (
    odds_above,
    odds_below,
    threshold
    ) {

    if (
        !is.numeric(threshold) ||
        length(threshold) != 1L ||
        !is.finite(threshold) ||
        threshold <= 1
        ) {

        stop ("`threshold` must be greater than 1.", call. = FALSE)

    }

    if (
        !is.numeric(odds_above) ||
        !is.numeric(odds_below) ||
        length(odds_above) != length(odds_below) ||
        anyNA(odds_above) ||
        anyNA(odds_below) ||
        any(odds_above < 0) ||
        any(odds_below < 0)
        ) {

        stop ("Invalid directional posterior odds.", call. = FALSE)

    }

    positive <- odds_above >= threshold
    negative <- odds_below >= threshold

    if (any(positive & negative) ) {

        stop (
            "Positive and negative evidence masks overlap.",
            call. = FALSE
            )

    }

    # This score is used ONLY for cluster detection.
    # It does not replace the continuous directional odds.
    score <- rep(1, length(odds_above) )
    score[positive] <- odds_above[positive]
    score[negative] <- 1 / odds_below[negative]

    return (score)

}
