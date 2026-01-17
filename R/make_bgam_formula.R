#' Build a \pkg{brms} formula for time-resolved BGAMMs
#'
#' Construct the appropriate \code{\link[brms]{brmsformula}} used internally by
#' \code{\link{testing_through_time}}, depending on the response distribution
#' (Gaussian vs. Binomial), whether a categorical predictor is provided, the
#' multilevel mode (\code{"summary"} vs. \code{"group"}), and the chosen random/varying
#' effects structure.
#'
#' The function assumes that the data passed to \code{brm()} has already been
#' reshaped to use the internal column names expected by \code{testing_through_time}:
#' \code{time}, \code{participant}, \code{predictor} (optional), and one of
#' \code{outcome_mean} / \code{outcome_sd} (Gaussian) or \code{success} / \code{trials}
#' (Binomial).
#'
#' @param family A family object describing the response distribution. Must be
#'   either \code{gaussian()} or \code{binomial()} (as used in \code{testing_through_time}).
#' @param multilevel Character; which model family to build. One of
#'   \code{"summary"} (participant-level varying effects) or \code{"group"}
#'   (population-level GAM without varying effects).
#' @param predictor_type Character; type of predictor. One of
#'   \code{"none"}, \code{"categorical"}, or \code{"continuous"}.
#' @param within_between Character; only used when \code{predictor_type = "categorical"}
#'   and \code{multilevel = "summary"}. Must be one of \code{"within-subject"}
#'   or \code{"between-subject"} to determine whether to include a varying slope
#'   \code{(1 + predictor | participant)} or only \code{(1 | participant)}.
#' @param time_id Character vector specifying the temporal variable(s) used in
#'   the smooth term. By default \code{"time"} (1D). If a character vector of
#'   length 2 is provided (e.g., \code{c("train_time", "test_time")}), the
#'   function constructs a 2D tensor-product smooth using
#'   \code{t2(train_time, test_time, ...)} (rather than \code{s(time, ...)}).
#'
#'   Currently, 2D \code{time_id} is supported for \code{multilevel = "group"}
#'   models (population-level smooth surface). Participant-level varying smooths
#'   for 2D surfaces are not implemented.
#' @param t2_full Logical; If TRUE, then there is a separate penalty for each
#' combination of null space column and range space, see \code{\link[mgcv]{t2}}.
#' Only use when fitting 2D temporal models (i.e., when \code{time_id} contains
#' two temporal variables).
#' @param kvalue Numeric; basis dimension \code{k} used in \code{s(time, ..., k = kvalue)}.
#' @param bs Character; spline basis used in \code{s(time, bs = bs, ...)} (e.g., \code{"tp"}).
#' @param include_ar_term Logical; if \code{TRUE}, adds an AR(1) autocorrelation
#'   structure within participant via
#'.  \code{autocor = brms::ar(time = "time", gr = "participant", p = 1, cov = TRUE)}.
#' @param varying_smooth Logical; if \code{TRUE} (default), add a factor-smooth interaction
#'   \code{s(participant, time, bs = "fs", ...)} in \code{multilevel = "summary"} models.
#'   If \code{FALSE}, only include the intercept/slope varying effects.
#' @param include_by_smooth Logical; if \code{TRUE} (default) and
#'   \code{predictor_type = "categorical"}, the time smooth is specified with
#'   \code{by = predictor}.
#' @param use_se Logical; whether to include known or internally computed
#'   measurement error via \code{y | se(outcome_sd)} in the model formula.
#'
#' @return A \code{\link[brms]{brmsformula}} object.
#'
#' @examples
#' \dontrun{
#' # Gaussian, summary, no predictor
#' make_bgam_formula(
#'   family = gaussian(),
#'   multilevel = "summary",
#'   predictor_type = "none",
#'   kvalue = 20, bs = "tp"
#'   )
#'
#' # Binomial, summary, within-subject predictor
#' make_bgam_formula(
#'   family = binomial(),
#'   multilevel = "summary",
#'   predictor_type = "categorical",
#'   within_between = "within-subject",
#'   kvalue = 20, bs = "tp"
#'   )
#' }
#'
#' @author Ladislas Nalborczyk \email{ladislas.nalborczyk@@cnrs.fr}.
#'
#' @export
make_bgam_formula <- function (
        family,
        multilevel = c("summary", "group"),
        predictor_type = c("none", "categorical", "continuous"),
        within_between = c(NA, "within-subject", "between-subject"),
        time_id = "time",
        t2_full = FALSE,
        kvalue = 20,
        bs = "tp",
        include_ar_term = FALSE,
        varying_smooth = TRUE,
        include_by_smooth = TRUE,
        use_se = TRUE
        ) {

    multilevel <- match.arg(multilevel)
    predictor_type <- match.arg(predictor_type)
    include_ar_term <- isTRUE(include_ar_term)

    # validate time_id (1D or 2D)
    if (!is.character(time_id) || !(length(time_id) %in% c(1L, 2L) ) ) {

        stop (
            "`time_id` must be a character vector of length 1 or 2 ",
            "(e.g., \"time\" or c(\"train_time\", \"test_time\")).",
            call. = FALSE
            )

    }

    is_2d_time <- length(time_id) == 2L

    # restrict families the same way testing_through_time() does
    fam_name <- tryCatch ({

        if (is.list(family) && !is.null(family$family) ) {

            as.character(family$family)

        } else {

            stop ("not a family object")

        }

    }, error = function (e) {

        stop ("`family` must be a valid family object (see ?brm).", call. = FALSE)

    })

    allowed_families <- c("gaussian", "binomial")

    if (!fam_name %in% allowed_families) {

        stop (
            "Unsupported `family`: '", fam_name, "'. ",
            "Currently supported families are: gaussian() and binomial().",
            call. = FALSE
            )

    }

    # optional AR(1) autocorrelation term
    autocor_term <- NULL

    if (isTRUE(include_ar_term) ) {

        if (multilevel == "summary") {

            autocor_term <- ~ar(time = time, gr = ar_series, p = 1, cov = TRUE)

        } else {

            autocor_term <- ~ar(time = time, gr = ar_series, p = 1, cov = FALSE)

        }


    }

    # LHS
    formula_lhs <- if (identical(fam_name, "binomial") ) {

        # expects columns: success, trials
        "success | trials(trials)"

    } else {

        if (multilevel == "group") {

            # expects column outcome_mean
            "outcome_mean"

        } else {

            if (isTRUE(include_ar_term) ) {

                # expects columns: outcome_mean, outcome_sd
                "outcome_mean | se(outcome_sd, sigma = TRUE)"

            } else {

                if (use_se) {

                    # expects columns: outcome_mean, outcome_sd
                    "outcome_mean | se(outcome_sd)"

                } else {

                    "outcome_mean"

                }

            }

        }

    }

    # predictor_type == "categorical"
    wb <- within_between[1]

    # smooth terms (1D vs 2D)
    if (!is_2d_time) {

        # smooth terms
        if (predictor_type == "none") {

            smooth_term <- glue::glue("s(time, bs = '{bs}', k = {kvalue})")
            varying_smooth_term <- glue::glue("s(participant, time, bs = 'fs', m = 1, k = {kvalue})")

        } else {

            if (include_by_smooth) {

                smooth_term <- glue::glue("s(time, bs = '{bs}', k = {kvalue}, by = predictor)")

                if (wb == "within-subject") {

                    # smooth_term <- glue::glue("s(time, bs = '{bs}', k = {kvalue}, by = predictor)")
                    varying_smooth_term <- glue::glue("s(participant, time, bs = 'fs', m = 1, k = {kvalue}, by = predictor)")

                } else {

                    # smooth_term <- glue::glue("s(time, bs = '{bs}', k = {kvalue})")
                    varying_smooth_term <- glue::glue("s(participant, time, bs = 'fs', m = 1, k = {kvalue})")

                }

            } else {

                smooth_term <- glue::glue("s(time, bs = '{bs}', k = {kvalue})")
                varying_smooth_term <- glue::glue("s(participant, time, bs = 'fs', m = 1, k = {kvalue})")

            }

        }

    } else {

        # 2D temporal
        # Use t2(time1, time2, ...) and mirror kvalue across both dimensions

        if (predictor_type == "none") {

            smooth_term <- glue::glue("t2(time1, time2, bs = '{bs}', k = c({kvalue}, {kvalue}), full = {t2_full})")
            varying_smooth_term <- glue::glue("t2(time1, time2, participant, bs = 'fs', k = {kvalue})")

        } else {

            if (include_by_smooth) {

                smooth_term <- glue::glue("t2(time1, time2, by = predictor, bs = '{bs}', k = c({kvalue}, {kvalue}), full = {t2_full})")

                if (wb == "within-subject") {

                    varying_smooth_term <- glue::glue("t2(time1, time2, participant, by = predictor, bs = 'fs', k = {kvalue})")

                } else {

                    varying_smooth_term <- glue::glue("t2(time1, time2, participant, bs = 'fs', k = {kvalue})")

                }

            } else {

                smooth_term <- glue::glue("t2(time1, time2, bs = '{bs}', k = c({kvalue}, {kvalue}), full = {t2_full})")
                varying_smooth_term <- glue::glue("t2(time1, time2, bs = '{bs}', k = c({kvalue}, {kvalue}), full = {t2_full})")

            }

        }

    }

    # group-level model (no varying effect)
    if (multilevel == "group") {

        rhs <- if (predictor_type == "none") {

            glue::glue("1 + {smooth_term}")

        } else {

            glue::glue("1 + predictor + {smooth_term}")

        }

        formula_str <- glue::glue("{formula_lhs} ~ {rhs}")

        if (include_ar_term) {

            return (brms::bf(formula = as.character(formula_str), autocor = autocor_term) )

        } else {

            return (brms::bf(formula = as.character(formula_str) ) )

        }

    }

    # summary-level model (varying effects)
    if (predictor_type == "none") {

        rhs_parts <- c(
            "1",
            "(1 | participant)",
            as.character(smooth_term)
            )

        if (varying_smooth) {

            rhs_parts <- c(rhs_parts, as.character(varying_smooth_term) )

        }

        formula_str <- paste0(formula_lhs, " ~ ", paste(rhs_parts, collapse = " + ") )

        if (isTRUE(include_ar_term) ) {

            return (brms::bf(formula_str, autocor = autocor_term) )

        } else {

            return (brms::bf(formula_str) )

        }

    }

    rand_eff <- if (identical(wb, "within-subject") ) {

        "(1 + predictor | participant)"

    } else {

        "(1 | participant)"

    }

    rhs_parts <- c(
        "1",
        "predictor",
        rand_eff,
        as.character(smooth_term)
        )

    if (isTRUE(varying_smooth) ) {

        rhs_parts <- c(rhs_parts, as.character(varying_smooth_term) )

    }

    formula_str <- paste0(formula_lhs, " ~ ", paste(rhs_parts, collapse = " + ") )

    if (isTRUE(include_ar_term) ) {

        return (brms::bf(formula_str, autocor = autocor_term) )

    } else {

        return (brms::bf(formula_str) )

    }

}
