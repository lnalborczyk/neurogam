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
#' \code{outcome_mean} / \code{outcome_sem} (Gaussian) or \code{success} / \code{trials}
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
#' @param include_by_varying_smooth Logical; if \code{TRUE} (non default) and
#'   \code{predictor_type = "categorical"}, the time varying smooth is specified with
#'   \code{by = predictor}.
#' @param use_se Logical; whether to include known or internally computed
#'   measurement error via \code{y | se(outcome_sem)} in the model formula.
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
#' @param space_id Optional vector of two spatial coordinate names. If given,
#'   build a 3D spatial-temporal `t2()` smooth instead of the 1D/2D smooth.
#' @param tensor_d Tensor marginal dimensions (currently `c(2, 1)` for
#'   two spatial coordinates jointly and one temporal coordinate).
#' @export
make_bgam_formula <- function(
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
        include_by_varying_smooth = FALSE,
        use_se = TRUE,
        space_id = NULL,
        tensor_d = c(2L, 1L)
        ) {

    # STAGE 3 DISPATCH: leave all existing 1D and 2D paths unchanged.
    if (!is.null(space_id) ) {

        return (
            make_bgam_formula_spatiotemporal(
                family = family,
                multilevel = match.arg(multilevel),
                predictor_type = match.arg(predictor_type),
                within_between = within_between[1L],
                time_id = time_id,
                space_id = space_id,
                tensor_d = tensor_d,
                kvalue = kvalue,
                bs = bs,
                t2_full = t2_full,
                varying_smooth = varying_smooth,
                include_ar_term = include_ar_term,
                include_by_smooth = include_by_smooth,
                use_se = use_se
                )
            )

    }

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

    # Compose the temporal model from independently testable components.
    # All formulas and term ordering retain the pre-refactor conventions.
    formula_lhs <- bgam_formula_response_lhs(
        family_name = fam_name,
        multilevel = multilevel,
        use_se = use_se,
        include_ar_term = include_ar_term
    )

    smooth_terms <- bgam_formula_temporal_smooth_terms(
        is_2d_time = is_2d_time,
        predictor_type = predictor_type,
        within_between = within_between,
        kvalue = kvalue,
        bs = bs,
        t2_full = t2_full,
        include_by_smooth = include_by_smooth,
        include_by_varying_smooth = include_by_varying_smooth
    )

    formula_rhs <- bgam_formula_temporal_rhs(
        multilevel = multilevel,
        predictor_type = predictor_type,
        within_between = within_between,
        varying_smooth = varying_smooth,
        smooth = smooth_terms$smooth,
        varying = smooth_terms$varying
    )

    formula_str <- paste0(formula_lhs, " ~ ", formula_rhs)
    if (isTRUE(include_ar_term)) {
        return(brms::bf(formula_str, autocor = autocor_term))
    }

    brms::bf(formula_str)

}
