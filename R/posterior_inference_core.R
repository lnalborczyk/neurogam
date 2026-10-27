#' Summarise posterior draws at arbitrary prediction locations
#'
#' Each column is one prediction location, and each row is a paired posterior
#' draw. This low-level implementation does not know about EEG coordinates,
#' sensors, conditions, or time. Those are handled by adapters.
#'
#' The `legacy_one_sided` convention is the default across neurogam.
#' Other conventions remain available explicitly for reproducibility.
#' Changing the default is an intentional inference-policy change.
#'
#' @param draws Finite numeric matrix, posterior draws by locations.
#' @param null_value Point-null value. For `legacy_directional`, the null is
#'   fixed at zero and `rope` is the half-width around it.
#' @param rope Non-negative half-width of the symmetric ROPE.
#' @param credible_interval Width of the equal-tailed credible interval.
#' @param odds_mode One of `legacy_one_sided`, `legacy_directional`, or
#'   `spatiotemporal`.
#' @param odds_cap Finite-sample cap for legacy odds. Defaults to the number
#'   of posterior draws supplied; ignored for `spatiotemporal`.
#'
#' @return Data frame with one row per prediction location.
#'
#' @keywords internal
summarise_posterior_matrix <- function (
    draws,
    null_value = 0,
    rope = 0,
    credible_interval = 0.95,
    odds_mode = c("legacy_one_sided", "legacy_directional", "spatiotemporal"),
    odds_cap = NULL
    ) {

    odds_mode <- match.arg(odds_mode)

    if (!is.matrix(draws) || !is.numeric(draws) ||
        nrow(draws) < 1L || ncol(draws) < 1L || any(!is.finite(draws))) {
        stop("`draws` must be a finite numeric matrix with positive dimensions.",
             call. = FALSE)
    }

    scalar_number <- function(x) {
        is.numeric(x) && length(x) == 1L &&
            !is.na(x) && is.finite(x)
    }

    if (!scalar_number(null_value)) {
        stop("`null_value` must be a finite numeric scalar.", call. = FALSE)
    }
    if (!scalar_number(rope) || rope < 0) {
        stop("`rope` must be a nonnegative finite numeric scalar.", call. = FALSE)
    }
    if (!scalar_number(credible_interval) ||
        credible_interval <= 0 || credible_interval >= 1) {
        stop("`credible_interval` must be strictly between zero and one.",
             call. = FALSE)
    }
    if (is.null(odds_cap)) {
        odds_cap <- nrow(draws)
    }
    if (odds_mode != "spatiotemporal" &&
        (!scalar_number(odds_cap) || odds_cap < 1)) {
        stop("`odds_cap` must be a positive numeric scalar (>= 1) in legacy modes.",
             call. = FALSE)
    }

    if (odds_mode == "legacy_directional") {
        if (null_value != 0) {
            stop("`legacy_directional` compares against zero; use `rope`.",
                 call. = FALSE)
        }
        prob_above <- colMeans(draws > rope)
        prob_below <- colMeans(draws < -rope)
        prob_inside <- colMeans(abs(draws) <= rope)
    } else {
        prob_above <- colMeans(draws > (null_value + rope))
        prob_below <- colMeans(draws < (null_value - rope))
        prob_inside <- 1 - prob_above - prob_below
    }

    odds_above <- prob_above / (1 - prob_above)
    odds_below <- prob_below / (1 - prob_below)

    if (odds_mode == "legacy_one_sided") {
        prob_ratio <- pmax(pmin(odds_above, odds_cap), 1 / odds_cap)
    } else if (odds_mode == "legacy_directional") {
        odds_above <- pmin(odds_above, odds_cap)
        odds_below <- pmin(odds_below, odds_cap)
        prob_ratio <- rep(1, ncol(draws))
        positive <- prob_above > 0.5
        negative <- !positive & prob_below > 0.5
        prob_ratio[positive] <- odds_above[positive]
        prob_ratio[negative] <- 1 / pmax(
            odds_below[negative],
            1 / odds_cap
        )
    } else {
        # No finite-sample odds cap in the existing spatiotemporal API.
        prob_ratio <- rep(1, ncol(draws))
        positive <- prob_above > 0.5
        negative <- !positive & prob_below > 0.5
        prob_ratio[positive] <- odds_above[positive]
        prob_ratio[negative] <- 1 / odds_below[negative]
    }

    alpha <- (1 - credible_interval) / 2
    interval <- vapply(
        seq_len(ncol(draws)),
        function(j) {
            stats::quantile(
                draws[, j],
                probs = c(alpha, 1 - alpha),
                names = FALSE
            )
        },
        numeric(2L)
    )

    return (
        data.frame(
            estimate = apply(draws, 2L, stats::median),
            lower = interval[1L, ],
            upper = interval[2L, ],
            prob_above = prob_above,
            prob_below = prob_below,
            prob_inside = prob_inside,
            odds_above = odds_above,
            odds_below = odds_below,
            prob_ratio = prob_ratio,
            row.names = NULL
            )
        )

}

#' Group long-format posterior draws and use the shared inference core
#'
#' @param data Long data frame with coordinate columns and one draw value.
#' @param groups Coordinate/group columns identifying a prediction location.
#' @param value Column containing posterior draw values.
#' @param ... Arguments passed to `summarise_posterior_matrix()`.
#' @return A data frame in coordinate order with one row per location.
#'
#' @keywords internal
summarise_posterior_long <- function (data, groups, value = ".epred", ...) {
    if (!is.data.frame(data) || !length(groups) ||
        anyDuplicated(groups) ||
        !all(c(groups, value) %in% names(data))) {
        stop("`data` must contain the grouping columns and value column.",
             call. = FALSE)
    }
    if (!nrow(data) || anyNA(data[groups])) {
        stop("`data` must contain nonmissing prediction locations.",
             call. = FALSE)
    }

    grouped <- dplyr::group_by(
        data,
        dplyr::across(dplyr::all_of(groups))
    )
    group_keys <- dplyr::group_keys(grouped)
    pieces <- dplyr::group_split(grouped, .keep = TRUE)

    summaries <- lapply(pieces, function(piece) {
        summarise_posterior_matrix(
            matrix(piece[[value]], ncol = 1L),
            ...
        )
    })

    result <- dplyr::bind_cols(
        group_keys,
        dplyr::bind_rows(summaries)
        )

    return (as.data.frame(result) )

}

#' Create paired A/B posterior contrasts in long format
#'
#' @param post_draws Long predictions with `.draw`, `.epred`, and `predictor`.
#' @param coordinates Time or spatial coordinates; participant optional.
#' @param predictor_type `categorical` or `continuous`.
#' @return A data frame with `epred_diff = level2 - level1` by draw/location.
#'   Categorical factor levels determine order; character values are sorted;
#'   numeric continuous levels are sorted numerically after rounding.
#'
#' @keywords internal
pair_posterior_contrasts <- function(
    post_draws,
    coordinates,
    predictor_type = c("categorical", "continuous")
    ) {
    predictor_type <- match.arg(predictor_type)
    required <- c(coordinates, "predictor", ".epred", ".draw")

    if (!is.data.frame(post_draws) || !length(coordinates) ||
        anyDuplicated(coordinates) ||
        !all(required %in% names(post_draws))) {
        stop("Missing columns or invalid coordinates in posterior contrast data.",
             call. = FALSE)
    }

    # Always create this local object before working with `draws`.
    draws <- as.data.frame(post_draws[, required, drop = FALSE])

    if (anyNA(draws[c(coordinates, "predictor", ".draw", ".epred")])) {
        stop("Posterior contrast data must not contain missing values.",
             call. = FALSE)
    }

    if (predictor_type == "continuous") {
        if (!is.numeric(draws$predictor)) {
            stop("Continuous predictor must be numeric.", call. = FALSE)
        }

        # Sort *numerically* before string conversion: 2 < 10.
        values <- round(draws$predictor, digits = 3L)
        conditions <- as.character(sort(unique(values)))
        draws$predictor <- as.character(values)
    } else if (is.factor(draws$predictor)) {
        # An explicit factor-level order takes precedence over input row order.
        conditions <- levels(droplevels(draws$predictor))
        draws$predictor <- as.character(draws$predictor)
    } else {
        # Character predictors have a stable, alphabetical direction.
        draws$predictor <- as.character(draws$predictor)
        conditions <- sort(unique(draws$predictor))
    }

    if (length(conditions) != 2L || anyNA(conditions)) {
        stop("Exactly two predictor levels are required.", call. = FALSE)
    }

    key_columns <- c(coordinates, ".draw", "predictor")
    if (anyDuplicated(draws[key_columns])) {
        stop("Duplicate posterior draws within a condition/location.",
             call. = FALSE)
    }

    paired <- tidyr::pivot_wider(
        draws,
        id_cols = dplyr::all_of(c(coordinates, ".draw")),
        names_from = "predictor",
        values_from = ".epred"
    )

    if (!all(conditions %in% names(paired)) ||
        anyNA(paired[conditions])) {
        stop("Posterior draws must be paired across the two conditions.",
             call. = FALSE)
    }

    paired$epred_diff <- paired[[conditions[2L]]] - paired[[conditions[1L]]]

    return (as.data.frame(paired) )

}
