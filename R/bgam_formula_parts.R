#' Internal helpers for constructing temporal BGAMM formulas
#'
#' These functions preserve the established 1D/2D formula conventions.
#' In particular, the order of terms and the use of factor smooths remain
#' unchanged. They are not part of the public API.
#'
#' @keywords internal
#' @name bgam_formula_parts
NULL

bgam_formula_response_lhs <- function(
    family_name,
    multilevel,
    use_se,
    include_ar_term
) {
    if (identical(family_name, "binomial")) {
        return("success | trials(trials)")
    }

    if (identical(multilevel, "group")) {
        return("outcome_mean")
    }

    if (isTRUE(use_se) || isTRUE(include_ar_term)) {
        return("outcome_mean | se(outcome_sem, sigma = TRUE)")
    }

    "outcome_mean"
}

bgam_formula_temporal_smooth_terms <- function(
    is_2d_time,
    predictor_type,
    within_between,
    kvalue,
    bs,
    t2_full,
    include_by_smooth,
    include_by_varying_smooth
) {
    wb <- within_between[1L]
    half_k <- floor(kvalue / 2)

    if (!is_2d_time) {
        if (identical(predictor_type, "none")) {
            smooth <- glue::glue("s(time, bs = '{bs}', k = {kvalue})")
            varying <- glue::glue(
                "s(time, participant, bs = 'fs', m = 1, k = {half_k})"
            )
        } else if (isTRUE(include_by_smooth)) {
            smooth <- glue::glue(
                "s(time, bs = '{bs}', k = {kvalue}, by = predictor)"
            )
            if (identical(wb, "within-subject") &&
                isTRUE(include_by_varying_smooth)) {
                varying <- glue::glue(
                    "s(time, participant, bs = 'fs', m = 1, k = {half_k}, by = predictor)"
                )
            } else {
                varying <- glue::glue(
                    "s(time, participant, bs = 'fs', m = 1, k = {half_k})"
                )
            }
        } else {
            smooth <- glue::glue("s(time, bs = '{bs}', k = {kvalue})")
            varying <- glue::glue(
                "s(time, participant, bs = 'fs', m = 1, k = {kvalue})"
            )
        }
    } else {
        if (identical(predictor_type, "none")) {
            smooth <- glue::glue(
                "t2(time1, time2, bs = '{bs}', k = c({kvalue}, {kvalue}), full = {t2_full})"
            )
            varying <- glue::glue(
                "t2(time1, time2, participant, bs = 'fs', k = {half_k})"
            )
        } else if (isTRUE(include_by_smooth)) {
            smooth <- glue::glue(
                "t2(time1, time2, by = predictor, bs = '{bs}', k = c({kvalue}, {kvalue}), full = {t2_full})"
            )
            if (identical(wb, "within-subject")) {
                varying <- glue::glue(
                    "t2(time1, time2, participant, by = predictor, bs = 'fs', k = {half_k})"
                )
            } else {
                varying <- glue::glue(
                    "t2(time1, time2, participant, bs = 'fs', k = {half_k})"
                )
            }
        } else {
            smooth <- glue::glue(
                "t2(time1, time2, bs = '{bs}', k = c({kvalue}, {kvalue}), full = {t2_full})"
            )
            varying <- glue::glue(
                "t2(time1, time2, bs = '{bs}', k = c({half_k}, {half_k}), full = {t2_full})"
            )
        }
    }

    list(
        smooth = as.character(smooth),
        varying = as.character(varying)
    )
}

bgam_formula_temporal_rhs <- function(
    multilevel,
    predictor_type,
    within_between,
    varying_smooth,
    smooth,
    varying
) {
    if (identical(multilevel, "group")) {
        terms <- if (identical(predictor_type, "none")) {
            c("1", smooth)
        } else {
            c("1", "predictor", smooth)
        }
        return(paste(terms, collapse = " + "))
    }

    if (identical(predictor_type, "none")) {
        terms <- c("1", "(1 | participant)", smooth)
        if (isTRUE(varying_smooth)) {
            terms <- c(terms, varying)
        }
        return(paste(terms, collapse = " + "))
    }

    random_term <- if (identical(within_between[1L], "within-subject")) {
        "(1 + predictor | participant)"
    } else {
        "(1 | participant)"
    }

    terms <- c("1", "predictor", random_term, smooth)
    if (isTRUE(varying_smooth)) {
        # Preserve the historical convention: factor smooths replace the
        # explicit varying intercept/slope for models with predictors.
        terms <- c(terms[-3L], varying)
    }

    paste(terms, collapse = " + ")
}
