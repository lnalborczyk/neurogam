#' Select posterior draw indices for prediction
#'
#' @param n_available Number of draws in a fitted model.
#' @param n_requested NULL for all draws, otherwise a positive integer.
#' @return NULL when all draws are requested, otherwise integer draw indices.
#' @keywords internal
select_prediction_draws <- function (n_available, n_requested = NULL) {

    if (length(n_available) != 1L || !is.numeric(n_available) ||

        is.na(n_available) || !is.finite(n_available) ||
        n_available < 1 || n_available != as.integer(n_available)) {

        stop ("`n_available` must be a positive integer.", call. = FALSE)

    }

    if (is.null(n_requested) ) {

        return (NULL)

    }

    if (length(n_requested) != 1L || !is.numeric(n_requested) ||
        is.na(n_requested) || !is.finite(n_requested) ||
        n_requested < 1 || n_requested != as.integer(n_requested) ||
        n_requested > n_available) {

        stop (
            "`n_requested` must be a positive integer no larger than `n_available`.",
            call. = FALSE
            )

    }

    if (n_requested == 1L) {

        return (1L)

    }

    return (as.integer(round(seq.int(1L, n_available, length.out = n_requested) ) ) )

}

#' Split observation positions into consecutive prediction chunks
#'
#' @param n Number of prediction locations.
#' @param chunk_size Maximum number of locations in each chunk.
#' @return List of integer position vectors.
#' @keywords internal
prediction_chunks <- function (n, chunk_size = 256L) {

    valid <- function (x) {

        is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x) && x >= 1 && x == as.integer(x)

    }

    if (!valid(n) || !valid(chunk_size)) {

        stop ("`n` and `chunk_size` must be positive integers.", call. = FALSE)

    }

    positions <- seq_len(n)
    split(positions, ceiling(positions / chunk_size) )

}

#' Retrieve chunked posterior expectation draws with paired contrasts
#'
#' @param model A fitted model, passed to `predict_fun`.
#' @param newdata Prediction grid for condition A.
#' @param newdata_b Optional matched prediction grid for condition B.
#' @param draw_ids Optional indices, identical for every chunk and condition.
#' @param chunk_size Number of prediction locations per chunk.
#' @param re_formula Group-level effects specification.
#' @param predict_fun Function accepting model, newdata, re_formula, draw_ids.
#'   Defaults to `brms::posterior_epred`.
#' @param summarise_fun Optional function applied to each draws-by-location
#'   chunk. If NULL, return a list of draw matrices.
#' @return List of chunk outputs in prediction-grid order.
#' @keywords internal
chunked_epred <- function (
    model,
    newdata,
    newdata_b = NULL,
    draw_ids = NULL,
    chunk_size = 256L,
    re_formula = NA,
    predict_fun = brms::posterior_epred,
    summarise_fun = NULL
    ) {

    if (!is.data.frame(newdata) || nrow(newdata) == 0L) {

        stop ("`newdata` must be a nonempty data frame.", call. = FALSE)

    }

    if (!is.null(newdata_b) &&
        (!is.data.frame(newdata_b) || nrow(newdata_b) != nrow(newdata) ) ) {

        stop (
            "`newdata_b` must have the same number of rows as `newdata`.",
            call. = FALSE
            )

    }

    if (!is.null(draw_ids) &&
        (!is.numeric(draw_ids) || !length(draw_ids) || anyNA(draw_ids) ||
         any(!is.finite(draw_ids)) || any(draw_ids < 1) ||
         any(draw_ids != as.integer(draw_ids)) || anyDuplicated(draw_ids))) {
        stop("`draw_ids` must contain distinct positive integer indices.",
             call. = FALSE)
    }

    chunks <- prediction_chunks(nrow(newdata), chunk_size)
    out <- vector("list", length(chunks))

    for (i in seq_along(chunks) ) {

        rows <- chunks[[i]]
        draws <- predict_fun(
            model,
            newdata = newdata[rows, , drop = FALSE],
            re_formula = re_formula,
            draw_ids = draw_ids
            )

        if (!is.matrix(draws) || !is.numeric(draws) ||
            ncol(draws) != length(rows) || nrow(draws) < 1L) {

            stop (
                "`predict_fun` must return a numeric draws-by-locations matrix.",
                call. = FALSE
                )

        }

        if (!is.null(newdata_b) ) {

            draws_b <- predict_fun(
                model,
                newdata = newdata_b[rows, , drop = FALSE],
                re_formula = re_formula,
                draw_ids = draw_ids
                )

            if (!is.matrix(draws_b) || !identical(dim(draws), dim(draws_b) ) ) {

                stop (
                    "Condition predictions have incompatible dimensions.",
                    call. = FALSE
                    )

            }

            draws <- draws - draws_b

        }

        out[[i]] <- if (is.null(summarise_fun) ) {

            draws

        } else {

            summarise_fun(draws)

        }

    }

    return (out)

}
