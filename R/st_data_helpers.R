#' Create participant-level EEG summaries for spatial-temporal fitting
#'
#' @param data Trial-level or previously summarized data frame.
#' @param participant_id, outcome_id, outcome_sem, time_id, sensor_id,
#'   space_id, predictor_id Source column names; `predictor_id` may be `NA`.
#' @param use_se Whether the model will use observation-level SEM.
#'
#' @return Data frame with standardized model columns.
#' @importFrom rlang .data
#' @keywords internal
make_spatiotemporal_summary_data <- function(
    data,
    participant_id,
    outcome_id,
    outcome_sem = NULL,
    time_id = "time",
    sensor_id = "channel",
    space_id = c("xproj", "yproj"),
    predictor_id = NA_character_,
    use_se = TRUE
) {
    has_predictor <- length(predictor_id) == 1L && !is.na(predictor_id)
    required <- c(participant_id, outcome_id, time_id, sensor_id, space_id)
    if (has_predictor) {
        required <- c(required, predictor_id)
    }
    if (!is.null(outcome_sem)) {
        required <- c(required, outcome_sem)
    }
    if (!all(required %in% names(data))) {
        stop("Missing data columns: ",
             paste(setdiff(required, names(data)), collapse = ", "),
             call. = FALSE)
    }
    if (length(space_id) != 2L || anyDuplicated(required)) {
        stop("Participant, response, sensor and coordinate columns must be distinct.",
             call. = FALSE)
    }

    x <- data.frame(
        participant = as.character(data[[participant_id]]),
        sensor = as.character(data[[sensor_id]]),
        xproj = as.numeric(data[[space_id[1L]]]),
        yproj = as.numeric(data[[space_id[2L]]]),
        time = as.numeric(data[[time_id]]),
        response = as.numeric(data[[outcome_id]]),
        stringsAsFactors = FALSE
    )
    if (has_predictor) {
        original <- data[[predictor_id]]
        if (is.numeric(original)) {
            stop("Continuous spatiotemporal predictors are not yet supported.",
                 call. = FALSE)
        }
        x$predictor <- factor(original)
        if (nlevels(x$predictor) != 2L) {
            stop("Spatiotemporal condition comparisons require exactly two levels.",
                 call. = FALSE)
        }
    }
    if (!is.null(outcome_sem)) {
        x$outcome_sem <- as.numeric(data[[outcome_sem]])
    }
    if (anyNA(x) || any(!is.finite(x$time)) ||
        any(!is.finite(x$xproj)) || any(!is.finite(x$yproj)) ||
        any(!is.finite(x$response)) ||
        (!is.null(outcome_sem) && any(!is.finite(x$outcome_sem)))) {
        stop("Spatiotemporal input contains missing or non-finite values.",
             call. = FALSE)
    }
    positions <- unique(x[c("sensor", "xproj", "yproj")])
    if (anyDuplicated(positions["sensor"])) {
        stop("Each sensor must have exactly one pair of projected coordinates.",
             call. = FALSE)
    }

    keys <- c("participant", "sensor", "xproj", "yproj", "time")
    if (has_predictor) {
        keys <- c(keys, "predictor")
    }
    if (is.null(outcome_sem)) {
        summarized <- x |>
            dplyr::group_by(dplyr::across(dplyr::all_of(keys))) |>
            dplyr::summarise(
                outcome_mean = mean(.data$response),
                outcome_sem = if (dplyr::n() >= 2L) {
                    stats::sd(.data$response) / sqrt(dplyr::n())
                } else {
                    NA_real_
                },
                .groups = "drop"
            )
        if (isTRUE(use_se) && anyNA(summarized$outcome_sem)) {
            stop("`use_se = TRUE` requires at least two trials per participant/sensor/time/condition, or a supplied `outcome_sem`. Alternatively use `use_se = FALSE`.",
                 call. = FALSE)
        }
    } else {
        if (anyDuplicated(x[keys])) {
            stop("With `outcome_sem` provided, input must already contain one row per participant/sensor/time/condition.",
                 call. = FALSE)
        }
        x$outcome_mean <- x$response
        summarized <- x[c(keys, "outcome_mean", "outcome_sem")]
    }
    if (!isTRUE(use_se)) {
        summarized$outcome_sem[is.na(summarized$outcome_sem)] <- 0
    }
    summarized$participant <- factor(summarized$participant)
    summarized$sensor <- as.character(summarized$sensor)
    if (has_predictor) {
        summarized$predictor <- factor(
            summarized$predictor,
            levels = levels(x$predictor)
        )
    }
    as.data.frame(summarized)
}

#' Build paired EEG sensor-time prediction grids
#'
#' @param summary_data Output of `make_spatiotemporal_summary_data()`.
#' @param contrast_levels Optional two-element character vector defining A - B.
#' @return List with `a`, `b` grids and a readable contrast label.
#' @keywords internal
make_spatiotemporal_prediction_grid <- function(summary_data, contrast_levels = NULL) {
    positions <- unique(summary_data[c("sensor", "xproj", "yproj")])
    times <- data.frame(time = sort(unique(summary_data$time)))
    grid <- merge(positions, times, by = NULL, sort = FALSE)
    grid <- grid[order(grid$time, grid$sensor), , drop = FALSE]
    rownames(grid) <- NULL
    grid$participant <- factor(
        as.character(summary_data$participant[1L]),
        levels = levels(summary_data$participant)
    )
    if (!"predictor" %in% names(summary_data)) {
        if (!is.null(contrast_levels)) {
            stop("`contrast_levels` requires a categorical predictor.", call. = FALSE)
        }
        return(list(a = grid, b = NULL, contrast = "one_sample"))
    }
    lev <- levels(summary_data$predictor)
    if (is.null(contrast_levels)) {
        contrast_levels <- rev(lev)  # coefficient direction: level 2 minus level 1
    }
    if (!is.character(contrast_levels) || length(contrast_levels) != 2L ||
        identical(contrast_levels[1L], contrast_levels[2L]) ||
        !setequal(contrast_levels, lev)) {
        stop("`contrast_levels` must name the two condition levels in A - B order.",
             call. = FALSE)
    }
    a <- grid
    b <- grid
    a$predictor <- factor(contrast_levels[1L], levels = lev)
    b$predictor <- factor(contrast_levels[2L], levels = lev)
    list(
        a = a,
        b = b,
        contrast = paste(contrast_levels, collapse = " - ")
    )
}
