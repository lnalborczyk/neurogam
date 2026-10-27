#' Construct a spatiotemporal Bayesian GAMM formula
#'
#' Internal formula builder for one temporal and two spatial coordinates.
#' The spatial coordinates form one marginal of a tensor-product smooth,
#' and time forms the other marginal.
#'
#' @param family A Gaussian family object.
#' @param multilevel Must be `"summary"`.
#' @param predictor_type `"none"` or `"categorical"`.
#' @param within_between `"within-subject"` or `"between-subject"` when
#'   a predictor is supplied.
#' @param time_id One temporal variable name.
#' @param space_id Two spatial variable names.
#' @param kvalue Basis dimensions for the two marginal bases.
#' @param bs Basis families for the two marginal bases.
#' @param tensor_d Dimensions per marginal, currently `c(2, 1)`.
#' @param t2_full Whether to use the full `mgcv::t2()` decomposition.
#' @param varying_smooth Participant-specific 3D smooths are not supported.
#' @param include_ar_term AR residual structure is not supported here.
#' @param include_by_smooth Include a condition-specific smooth.
#' @param use_se Use `se(outcome_sem, sigma = TRUE)` on the left-hand side.
#'
#' @return A `brmsformula` object.
#' @keywords internal
make_bgam_formula_spatiotemporal <- function(
    family,
    multilevel = "summary",
    predictor_type = "none",
    within_between = NA_character_,
    time_id = "time",
    space_id = c("xproj", "yproj"),
    kvalue = c(15L, 10L),
    bs = c("ps", "ps"),
    tensor_d = c(2L, 1L),
    t2_full = FALSE,
    varying_smooth = FALSE,
    include_ar_term = FALSE,
    include_by_smooth = TRUE,
    use_se = TRUE
) {
    if (!identical(multilevel, "summary")) {
        stop("Spatiotemporal formulas currently require `multilevel = 'summary'`.",
             call. = FALSE)
    }
    if (!identical(family$family, "gaussian")) {
        stop("Spatiotemporal formulas currently support Gaussian models only.",
             call. = FALSE)
    }
    if (!predictor_type %in% c("none", "categorical")) {
        stop("Spatiotemporal formulas support no predictor or a binary categorical predictor.",
             call. = FALSE)
    }
    if (isTRUE(varying_smooth)) {
        stop("Participant-specific 3D smooths are not yet supported; set `varying_smooth = FALSE`.",
             call. = FALSE)
    }
    if (isTRUE(include_ar_term)) {
        stop("AR residual structure is not yet supported in the 3D interface.",
             call. = FALSE)
    }
    coordinates <- c(space_id, time_id)
    if (!is.character(space_id) || length(space_id) != 2L ||
        !is.character(time_id) || length(time_id) != 1L ||
        anyNA(coordinates) || any(!nzchar(coordinates)) ||
        anyDuplicated(coordinates) ||
        any(!grepl("^[A-Za-z.][A-Za-z0-9._]*$", coordinates))) {
        stop("Supply three distinct syntactically valid coordinate names.",
             call. = FALSE)
    }
    if (!is.numeric(tensor_d) || length(tensor_d) != 2L ||
        anyNA(tensor_d) || any(tensor_d != c(2, 1))) {
        stop("For two spatial coordinates and one time variable, use `tensor_d = c(2, 1)`.",
             call. = FALSE)
    }
    if (length(kvalue) == 1L) {
        kvalue <- rep(kvalue, 2L)
    }
    if (length(bs) == 1L) {
        bs <- rep(bs, 2L)
    }
    if (!is.numeric(kvalue) || length(kvalue) != 2L ||
        anyNA(kvalue) || any(!is.finite(kvalue)) ||
        any(kvalue < 2L) || any(kvalue != as.integer(kvalue))) {
        stop("`kvalue` must contain two positive integer basis dimensions (>= 2).",
             call. = FALSE)
    }
    if (!is.character(bs) || length(bs) != 2L ||
        anyNA(bs) || any(!nzchar(bs))) {
        stop("`bs` must contain one basis name for each marginal.", call. = FALSE)
    }
    if (!is.logical(t2_full) || length(t2_full) != 1L || is.na(t2_full)) {
        stop("`t2_full` must be TRUE or FALSE.", call. = FALSE)
    }
    if (identical(predictor_type, "categorical") &&
        !within_between %in% c("within-subject", "between-subject")) {
        stop("Specify `within_between` for categorical predictors.", call. = FALSE)
    }

    smooth_args <- c(
        coordinates,
        sprintf("bs = c(%s)", paste(sprintf('"%s"', bs), collapse = ", ")),
        sprintf("k = c(%s)", paste(as.integer(kvalue), collapse = ", ")),
        "d = c(2, 1)",
        sprintf("full = %s", if (t2_full) "TRUE" else "FALSE")
    )
    if (identical(predictor_type, "categorical") && isTRUE(include_by_smooth)) {
        smooth_args <- c(smooth_args, "by = predictor")
    }
    smooth <- sprintf("t2(%s)", paste(smooth_args, collapse = ", "))

    rhs <- c("1", smooth)
    if (identical(predictor_type, "categorical")) {
        rhs <- c("1", "predictor", smooth)
        random_term <- if (identical(within_between, "within-subject")) {
            "(1 + predictor | participant)"
        } else {
            "(1 | participant)"
        }
    } else {
        random_term <- "(1 | participant)"
    }
    rhs <- c(rhs, random_term)
    lhs <- if (isTRUE(use_se)) {
        "outcome_mean | se(outcome_sem, sigma = TRUE)"
    } else {
        "outcome_mean"
    }
    brms::bf(stats::as.formula(paste(lhs, "~", paste(rhs, collapse = " + "))))
}
