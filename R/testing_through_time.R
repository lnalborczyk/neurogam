#' Time-resolved testing based on BGAMMs
#'
#' Fits time-resolved Bayesian generalised additive (multilevel) models (BGAMMs)
#' using \pkg{brms}, and computes posterior odds for an effect at each time
#' point. The effect can be either i) a deviation of the outcome from a
#' reference value (e.g., zero or a chance level), or ii) a difference between two
#' groups/conditions.
#'
#' @param data A data frame in long format containing time-resolved data.
#' @param participant_id Character; name of the column in \code{data}
#' specifying participant IDs.
#' @param outcome_id Character; name of the column in \code{data} containing
#' the outcome values (e.g., M/EEG amplitude, decoding accuracy).
#' @param time_id Character; name of the column in \code{data}
#' containing time information (e.g., in seconds or samples).
#' @param predictor_id Character; name of the column in \code{data}
#'   containing either:
#'   \itemize{
#'     \item A \emph{binary} categorical predictor (e.g., group or condition),
#'       in which case the function tests, at each time point, whether the
#'       difference between the two levels exceeds
#'       \code{chance_level + sesoi};
#'     \item A \emph{continuous} numeric predictor, in which case the function
#'       tests, at each time point, whether the \emph{slope} of the outcome
#'       with respect to the predictor differs from \code{chance_level + sesoi}
#'       (typically with \code{chance_level = 0}).
#'     \item If \code{predictor_id = NA}, the function tests whether the outcome differs
#'       from \code{chance_level + sesoi} over time (useful for decoding accuracies,
#'       for instance).
#'   }
#' @param family A \pkg{brms} family object describing the response
#'   distribution to be used in the model (defaults to \code{gaussian()}).
#' @param kvalue Numeric; basis dimension \code{k} passed to the smooth term
#'   \code{s(time, ..., k = kvalue)}.
#' @param bs Character; Character scalar; type of spline basis to be used by \pkg{brms}
#'   (passed to \code{s()}, e.g., \code{"tp"} for thin-plate splines).
#' @param multilevel Character; which model to fit. One of
#'   \itemize{
#'     \item \code{"summary"}: GAMM fitted to participant-level summary
#'       statistics (mean outcome and its standard deviation);
#'     \item \code{"group"}: Group-level GAM fitted to participant-averaged
#'       data (no random/varying effects).
#'   }
#' @param by_ppt Logical; should we return clusters at the participant-level.
#' @param warmup Numeric; number of warm-up iterations per chain.
#' @param iter Numeric; total number of iterations per chain (including warmup).
#' @param chains Numeric; number of MCMCs.
#' @param cores Numeric; number of parallel cores to use.
#' @param backend Character; package to use as the backend for fitting the
#'   Stan model. One of \code{"cmdstanr"} (default) or \code{"rstan"}.
#' @param threshold Numeric; threshold on the posterior odds
#'   (\code{prob_ratio}) used to define contiguous temporal clusters. Values
#'   greater than 1 favour the hypothesis that the effect exceeds
#'   \code{chance_level + sesoi}.
#' @param n_post_samples Numeric; number of posterior draws used to compute
#'   posterior probabilities. If \code{NULL} (default), all available draws
#'   from the fitted model are used.
#' @param chance_level Numeric; reference value for the outcome (e.g., 0.5 for
#'   decoding accuracy). Only used when testing against a constant (i.e.,
#'   when there is no \code{predictor_id} or when the effect is a difference
#'   from chance).
#' @param sesoi Numeric; smallest effect size of interest (SESOI). The
#'   posterior probability is computed for the effect being strictly larger
#'   than \code{chance_level + sesoi}.
#' @param credible_interval Numeric; width of the credible (quantile) interval.
#'
#' @return An object of class \code{"clusters_results"}, which is a list with
#'   elements:
#'   \itemize{
#'     \item \code{clusters}: a data frame with one row per detected cluster
#'       (e.g., \code{cluster_onset}, \code{cluster_offset}, \code{duration});
#'     \item \code{predictions}: a data frame with time-resolved posterior
#'       summaries (posterior median, credible interval, posterior
#'       probabilities, and odds \code{prob_ratio});
#'     \item \code{data}: data used to fit the \pkg{brms} model
#'       (possibly summarised);
#'     \item \code{model}: the fitted \pkg{brms} model object;
#'     \item \code{multilevel}: the value of the \code{multilevel} argument.
#'   }
#'
#'   The object has an associated \code{plot()} method for visualising the
#'   smoothed time course and detected clusters, as well as \code{print()} and
#'   \code{summary()} methods.
#'
#' @details
#' Internally, the function:
#' \enumerate{
#'   \item builds a formula with a smooth term over time (optionally by group);
#'   \item fits a \pkg{brms} model according to \code{multilevel};
#'   \item uses \pkg{tidybayes} to extract posterior predictions over time;
#'   \item computes, at each time point, the posterior probability that the
#'     effect (or condition difference) exceeds \code{chance_level + sesoi};
#'   \item converts this into posterior odds (\code{prob_ratio}) and applies
#'     a clustering procedure (\code{find_clusters()}) over time.
#' }
#'
#' @importFrom rlang .data
#' @importFrom stats gaussian
#'
#' @examples
#' \dontrun{
#' # import some simulated EEG data
#' data(eeg_data)
#' head(eeg_data)
#'
#' # fit the BGAMM to identify clusters
#' results <- testing_through_time(data = eeg_data)
#'
#' # display the identified clusters
#' print(results$clusters)
#'
#' # plot the GAM-smoothed signal and identified clusters
#' plot(results)
#' }
#'
#' @author Ladislas Nalborczyk \email{ladislas.nalborczyk@@cnrs.fr}.
#'
#' @seealso \code{\link[brms]{brm}}
#'
#' @export
testing_through_time <- function (
        data,
        participant_id = "participant", outcome_id = "eeg",
        time_id = "time", predictor_id = "condition",
        family = gaussian(), kvalue = 20, bs = "tp",
        multilevel = c("summary", "group"), by_ppt = FALSE,
        warmup = 1000, iter = 2000, chains = 4, cores = 4,
        backend = "cmdstanr",
        threshold = 10, n_post_samples = NULL,
        chance_level = 0, sesoi = 0, credible_interval = 0.95
        ) {

    # some tests for variable types
    stopifnot("data must be a dataframe..." = is.data.frame(data) )
    stopifnot("kvalue must be a numeric..." = is.numeric(kvalue) )
    stopifnot("warmup must be a numeric..." = is.numeric(warmup) )
    stopifnot("iter must be a numeric..." = is.numeric(iter) )
    stopifnot("chains must be a numeric..." = is.numeric(chains) )
    stopifnot("cores must be a numeric..." = is.numeric(cores) )
    stopifnot("threshold must be a numeric..." = is.numeric(threshold) )
    stopifnot("chance_level must be a numeric..." = is.numeric(chance_level) )
    stopifnot("sesoi must be a numeric..." = is.numeric(sesoi) )
    stopifnot("bs must be a character..." = is.character(bs) )

    # multilevel should be one of above
    multilevel <- match.arg(multilevel)

    # some more tests
    if (!backend %in% c("cmdstanr", "rstan") ) {

        stop ("`backend` must be either 'cmdstanr' or 'rstan'.", call. = FALSE)

    }

    if (warmup >= iter) {

        stop ("`iter` must be strictly larger than `warmup`.", call. = FALSE)

    }

    if (!is.null(n_post_samples) ) {

        if (!is.numeric(n_post_samples) || length(n_post_samples) != 1L || n_post_samples <= 0) {

            stop ("`n_post_samples` must be NULL or a positive numeric scalar.", call. = FALSE)

        }

    }

    # checking required column names
    required_columns <- c(participant_id, outcome_id, time_id)

    if (!is.na(predictor_id) ) {

        required_columns <- c(required_columns, predictor_id)

    }

    assertthat::assert_that(
        all(required_columns %in% colnames(data) ),
        msg = paste(
            "Missing columns:",
            paste(setdiff(required_columns, colnames(data) ), collapse = ", ")
            )
        )

    # checking predictor type
    predictor_type <- "none"

    if (!is.na(predictor_id) ) {

        pred_vec <- data[[predictor_id]]

        if (is.numeric(pred_vec) ) {

            predictor_type <- "continuous"

            stop ("Continuous predictors are not yet supported, please use brms.", call. = FALSE)

        } else {

            predictor_type <- "categorical"

            # optional: enforce 2 levels for categorical case
            if (length(unique(pred_vec) ) != 2) {

                stop ("For categorical `predictor_id`, there must be exactly 2 levels.", call. = FALSE)

            }

        }

    }

    if (multilevel == "summary") {

        # construct the smooth term dynamically
        if (is.na(predictor_id) ) {

            # reshape and summarising the data
            summary_data <- data |>
                # reshape the original variables
                dplyr::mutate(participant = .data[[participant_id]]) |>
                dplyr::mutate(outcome = .data[[outcome_id]]) |>
                dplyr::mutate(time = .data[[time_id]]) |>
                # summarise per participant
                dplyr::summarise(
                    outcome_mean = mean(.data$outcome),
                    outcome_sd = stats::sd(.data$outcome),
                    .by = c(.data$participant, .data$time)
                    )

            # define the smooth terms
            smooth_term <- glue::glue("s(time, bs = '{bs}', k = {kvalue})")
            varying_smooth <- glue::glue("s(participant, time, bs = 'fs', m = 1, k = {kvalue})")

        } else {

            # reshape and summarising the data
            summary_data <- data |>
                # reshape the original variables
                dplyr::mutate(predictor = as.factor(.data[[predictor_id]]) ) |>
                dplyr::mutate(participant = .data[[participant_id]]) |>
                dplyr::mutate(outcome = .data[[outcome_id]]) |>
                dplyr::mutate(time = .data[[time_id]]) |>
                # summarise per participant
                dplyr::summarise(
                    outcome_mean = mean(.data$outcome),
                    outcome_sd = stats::sd(.data$outcome),
                    .by = c(.data$participant, .data$predictor, .data$time)
                    )

            # test whether predictor is varying within or between participant
            within_between <- check_within_between(
                data = summary_data,
                participant = "participant",
                predictor = "predictor"
                )

            # define the smooth terms
            smooth_term <- glue::glue("s(time, bs = '{bs}', k = {kvalue}, by = predictor)")
            varying_smooth <- glue::glue("s(participant, time, bs = 'fs', m = 1, k = {kvalue})")

        }

        # full formula
        if (is.na(predictor_id) ) {

            # varying_intercept <- glue::glue("s(participant, bs = 're')")
            # varying_slope <- glue::glue("s(participant, time, bs = 're')")
            # varying_smooth <- glue::glue("s(participant, time, bs = 'fs', m = 1, k = {kvalue})")
            # formula_str <- glue::glue("outcome_mean | se(outcome_sd) ~ 1 + {smooth_term} + {varying_intercept} + {varying_slope}")
            # formula_str <- glue::glue("outcome_mean | se(outcome_sd) ~ 1 + {smooth_term} + {varying_intercept} + {varying_slope} + {varying_smooth}")
            formula_str <- glue::glue("outcome_mean | se(outcome_sd) ~ 1 + (1 | participant) + {smooth_term} + {varying_smooth}")

        } else {

            if (within_between$classification == "between-subject") {

                formula_str <- glue::glue(
                    "outcome_mean | se(outcome_sd) ~ 1 + predictor + (1 | participant) + {smooth_term} + {varying_smooth}"
                    )

            } else if (within_between$classification == "within-subject") {

                formula_str <- glue::glue(
                    "outcome_mean | se(outcome_sd) ~ 1 + predictor + (1 + predictor | participant) + {smooth_term} + {varying_smooth}"
                    )

            }

        }

        # convert to formula
        formula_obj <- brms::bf(formula_str)

        # fit the model
        brms_gam <- brms::brm(
            formula = formula_obj,
            data = summary_data,
            family = family,
            warmup = warmup,
            iter = iter,
            chains = chains,
            cores = cores,
            backend = backend,
            # speedup
            # running on 2x4 cores
            # threads = brms::threading(threads = 2),
            # stan_model_args = list(
            #     stanc_options = list("O1"),
            #     cpp_options = list(stan_threads = TRUE)
            #     )
            stan_model_args = list(stanc_options = list("O1") )
            )

        if (is.null(n_post_samples) ) {

            n_post_samples <- brms::ndraws(brms_gam)

        }

        # computing the posterior odds over time
        if (is.na(predictor_id) ) {

            if (by_ppt) {

                # newdata grid over time
                newdata_grid <- tidyr::crossing(
                    time = sort(unique(brms_gam$data$time) ),
                    participant = sort(unique(brms_gam$data$participant) )
                    ) |>
                    # dummy outcome_sd to satisfy | se(outcome_sd)
                    dplyr::mutate(outcome_sd = 1)

                # retrieving posterior predictions (draws)
                post_draws <- tidybayes::add_epred_draws(
                    object = brms_gam,
                    newdata = newdata_grid,
                    ndraws = n_post_samples,
                    re_formula = NULL
                    ) |>
                    data.frame()

            } else {

                # newdata grid over time
                newdata_grid <- tidyr::crossing(
                    time = sort(unique(brms_gam$data$time) )
                    ) |>
                    # dummy outcome_sd to satisfy | se(outcome_sd)
                    dplyr::mutate(outcome_sd = 1) |>
                    # arbitrary participant to satisfy validate_data()
                    # dplyr::mutate(participant = brms_gam$data$participant[1])
                    # setting participant = NA to satisfy validate_data()
                    dplyr::mutate(participant = NA)

                # retrieving posterior predictions (draws)
                post_draws <- tidybayes::add_epred_draws(
                    object = brms_gam,
                    newdata = newdata_grid,
                    ndraws = n_post_samples,
                    re_formula = NA
                    ) |>
                    data.frame()

            }

        } else if (predictor_type == "categorical") {

            if (by_ppt) {

                if (within_between$classification == "within-subject") {

                    # TO-DO

                } else if (within_between$classification == "between-subject") {

                    cat(
                        "We can not estimate clusters at the participant level when the predictor varies across participants...\nSwitching to by_ppt = FALSE\n"
                        )

                    by_ppt <- FALSE

                }

            } else {

                # newdata grid over time and predictor
                newdata_grid <- tidyr::crossing(
                    time = sort(unique(brms_gam$data$time) ),
                    predictor = levels(brms_gam$data$predictor)
                    ) |>
                    # dummy outcome_sd to satisfy | se(outcome_sd)
                    dplyr::mutate(outcome_sd = 1) |>
                    # arbitrary participant to satisfy validate_data()
                    # dplyr::mutate(participant = brms_gam$data$participant[1])
                    dplyr::mutate(participant = NA)

                # retrieving posterior predictions (draws)
                post_draws <- tidybayes::add_epred_draws(
                    object = brms_gam,
                    newdata = newdata_grid,
                    ndraws = n_post_samples,
                    re_formula = NA
                    ) |>
                    data.frame()

            }

        } else if (predictor_type == "continuous") {

            # TO-DO...

        }

    } else if (multilevel == "group") {

        # construct the smooth term dynamically
        if (is.na(predictor_id) ) {

            # reshaping and summarising the data
            summary_data <- data |>
                # reshaping the original variables
                dplyr::mutate(participant = .data[[participant_id]]) |>
                dplyr::mutate(outcome = .data[[outcome_id]]) |>
                dplyr::mutate(time = .data[[time_id]]) |>
                # summarising per participant
                dplyr::summarise(
                    outcome_mean = mean(.data$outcome),
                    outcome_sd = stats::sd(.data$outcome),
                    .by = c(.data$participant, .data$time)
                    )

            # define the smooth term
            smooth_term <- glue::glue("s(time, bs = '{bs}', k = {kvalue})")

        } else if (predictor_type == "categorical") {

            # reshape and summarising the data
            summary_data <- data |>
                # reshape the original variables
                dplyr::mutate(predictor = as.factor(.data[[predictor_id]]) ) |>
                dplyr::mutate(participant = .data[[participant_id]]) |>
                dplyr::mutate(outcome = .data[[outcome_id]]) |>
                dplyr::mutate(time = .data[[time_id]]) |>
                # summarise per participant
                dplyr::summarise(
                    outcome_mean = mean(.data$outcome),
                    outcome_sd = stats::sd(.data$outcome),
                    .by = c(.data$participant, .data$predictor, .data$time)
                    )

            # define the smooth term
            smooth_term <- glue::glue("s(time, bs = '{bs}', k = {kvalue}, by = predictor)")

        } else if (predictor_type == "continuous") {

            # TO-DO...

        }

        # full formula
        if (is.na(predictor_id) ) {

            formula_str <- glue::glue("outcome_mean ~ 1 + {smooth_term}")

        } else {

            formula_str <- glue::glue("outcome_mean ~ 1 + predictor + {smooth_term}")

        }

        # convert to formula
        formula_obj <- brms::bf(formula_str)

        # fit the model
        brms_gam <- brms::brm(
            formula = formula_obj,
            data = summary_data,
            family = family,
            warmup = warmup,
            iter = iter,
            chains = chains,
            cores = cores,
            backend = backend,
            # speedup
            # running on 2x4 cores
            # threads = brms::threading(threads = 2),
            # stan_model_args = list(
            #     stanc_options = list("O1"),
            #     cpp_options = list(stan_threads = TRUE)
            #     )
            stan_model_args = list(stanc_options = list("O1") )
            )

        if (is.null(n_post_samples) ) {

            n_post_samples <- brms::ndraws(brms_gam)

        }

        # computing the posterior odds over time
        if (is.na(predictor_id) ) {

            # newdata grid over time
            newdata_grid <- tidyr::crossing(time = sort(unique(brms_gam$data$time) ) )

            # retrieving posterior predictions (draws)
            post_draws <- tidybayes::add_epred_draws(
                object = brms_gam,
                newdata = newdata_grid,
                ndraws = n_post_samples,
                re_formula = NA
                ) |>
                data.frame()

        } else if (predictor_type == "categorical") {

            # newdata grid over time and predictor
            newdata_grid <- tidyr::crossing(
                time = sort(unique(brms_gam$data$time) ),
                predictor = levels(brms_gam$data$predictor)
                )

            # retrieving posterior predictions (draws)
            post_draws <- tidybayes::add_epred_draws(
                object = brms_gam,
                newdata = newdata_grid,
                ndraws = n_post_samples,
                re_formula = NA
                ) |>
                data.frame()

        } else if (predictor_type == "continuous") {

            # TO-DO...

        }

    }

    # compute the posterior odds
    if (is.na(predictor_id) ) {

        prob_y_above <- .compute_one_sample_prob(
            post_draws = post_draws,
            by_ppt = by_ppt,
            threshold = chance_level + sesoi,
            n_post_samples = n_post_samples,
            credible_interval = credible_interval
            )

    } else if (predictor_type == "categorical") {

        prob_y_above <- .compute_two_sample_prob(
            post_draws = post_draws,
            by_ppt = by_ppt,
            threshold = chance_level + sesoi,
            n_post_samples = n_post_samples,
            credible_interval = credible_interval
            )

    } else if (predictor_type == "continuous") {

        # TO-DO...

    }

    if (by_ppt) {

        # find the clusters
        clusters <- find_clusters(
            data = prob_y_above |> dplyr::select(.data$time, .data$participant, value = .data$prob_ratio),
            group = "participant",
            threshold = threshold
            )

    } else {

        # find the clusters
        clusters <- find_clusters(
            data = prob_y_above |> dplyr::select(.data$time, value = .data$prob_ratio),
            group = NULL,
            threshold = threshold
            )

    }

    # combine the results in a list
    clusters_results <- list(
        clusters = clusters,
        predictions = prob_y_above,
        model = brms_gam,
        multilevel = multilevel
        )

    # assign a new class to the list
    class(clusters_results) <- "clusters_results"

    # return the clusters and posterior probabilities
    return (clusters_results)

}

.compute_one_sample_prob <- function (post_draws, threshold, by_ppt, n_post_samples, credible_interval) {

    # validate credible interval
    if (!is.numeric(credible_interval) ||

        length(credible_interval) != 1L ||
        credible_interval <= 0 || credible_interval >= 1) {

        stop ("`credible_interval` must be a number strictly between 0 and 1.")

    }

    alpha <- (1 - credible_interval) / 2
    lower_q <- alpha
    upper_q <- 1 - alpha

    if (by_ppt) {

        prob_y_above <- post_draws |>
            dplyr::group_by(.data$time, .data$participant) |>
            dplyr::summarise(prob_above = mean(.data$.epred > threshold) ) |>
            dplyr::mutate(prob_ratio = .data$prob_above / (1 - .data$prob_above) ) |>
            # dplyr::ungroup() |>
            dplyr::mutate(
                prob_ratio = pmin(.data$prob_ratio, n_post_samples),
                prob_ratio = pmax(.data$prob_ratio, 1 / n_post_samples)
                # prob_ratio = min(.data$prob_ratio, n_post_samples),
                # prob_ratio = max(.data$prob_ratio, 1 / n_post_samples)
                ) |>
            dplyr::ungroup() |>
            data.frame()

        post_prob_slope <- post_draws |>
            dplyr::group_by(.data$time, .data$participant) |>
            dplyr::summarise(
                post_prob = stats::quantile(.data$.epred, probs = 0.5, na.rm = TRUE),
                lower = stats::quantile(.data$.epred, probs = lower_q, na.rm = TRUE),
                upper = stats::quantile(.data$.epred, probs = upper_q, na.rm = TRUE)
                ) |>
            dplyr::ungroup()

        results <- dplyr::left_join(prob_y_above, post_prob_slope, by = c("time", "participant") )

    } else {

        prob_y_above <- post_draws |>
            dplyr::group_by(.data$time) |>
            dplyr::summarise(
                prob_above = mean(.data$.epred > threshold)
                ) |>
            dplyr::mutate(prob_ratio = .data$prob_above / (1 - .data$prob_above) ) |>
            dplyr::ungroup() |>
            dplyr::mutate(
                prob_ratio = pmin(.data$prob_ratio, n_post_samples),
                prob_ratio = pmax(.data$prob_ratio, 1 / n_post_samples)
                ) |>
            data.frame()

        post_prob_slope <- post_draws |>
            dplyr::group_by(.data$time) |>
            dplyr::summarise(
                post_prob = stats::quantile(.data$.epred, probs = 0.5, na.rm = TRUE),
                lower = stats::quantile(.data$.epred, probs = lower_q, na.rm = TRUE),
                upper = stats::quantile(.data$.epred, probs = upper_q, na.rm = TRUE)
                ) |>
            dplyr::ungroup()

        results <- dplyr::left_join(prob_y_above, post_prob_slope, by = "time")

    }

    return (results)

}

.compute_two_sample_prob <- function (post_draws, threshold, by_ppt, n_post_samples, credible_interval) {

    # validate credible interval
    if (!is.numeric(credible_interval) ||

        length(credible_interval) != 1L ||
        credible_interval <= 0 || credible_interval >= 1) {

        stop ("`credible_interval` must be a number strictly between 0 and 1.")

    }

    alpha <- (1 - credible_interval) / 2
    lower_q <- alpha
    upper_q <- 1 - alpha

    cond1 <- unique(post_draws$predictor)[1]
    cond2 <- unique(post_draws$predictor)[2]

    if (by_ppt) {

        # To-be checked...
        post_diff <- post_draws |>
            dplyr::select(.data$time, .data$predictor, .data$participant, .data$.epred, .data$.draw) |>
            tidyr::pivot_wider(names_from = .data$predictor, values_from = .data$.epred) |>
            dplyr::mutate(epred_diff = .data[[cond2]] - .data[[cond1]])

        prob_y_above <- post_diff |>
            dplyr::group_by(.data$time, .data$participant) |>
            dplyr::summarise(
                prob_above = mean(.data$epred_diff > threshold)
                ) |>
            dplyr::mutate(prob_ratio = .data$prob_above / (1 - .data$prob_above) ) |>
            dplyr::ungroup() |>
            dplyr::mutate(
                prob_ratio = min(.data$prob_ratio, n_post_samples),
                prob_ratio = max(.data$prob_ratio, 1 / n_post_samples)
                ) |>
            data.frame()

        post_prob_slope <- post_diff |>
            dplyr::group_by(.data$time, .data$participant) |>
            dplyr::summarise(
                post_prob = stats::quantile(.data$epred_diff, probs = 0.5),
                lower = stats::quantile(.data$epred_diff, probs = lower_q),
                upper = stats::quantile(.data$epred_diff, probs = upper_q)
                ) |>
            dplyr::ungroup()

        results <- dplyr::left_join(prob_y_above, post_prob_slope, by = c("time", "participant") )

    } else {

        post_diff <- post_draws |>
            dplyr::select(.data$time, .data$predictor, .data$.epred, .data$.draw) |>
            tidyr::pivot_wider(names_from = .data$predictor, values_from = .data$.epred) |>
            dplyr::mutate(epred_diff = .data[[cond2]] - .data[[cond1]])

        prob_y_above <- post_diff |>
            dplyr::group_by(.data$time) |>
            dplyr::summarise(
                prob_above = mean(.data$epred_diff > threshold)
                ) |>
            dplyr::mutate(prob_ratio = .data$prob_above / (1 - .data$prob_above) ) |>
            dplyr::ungroup() |>
            dplyr::mutate(
                prob_ratio = pmin(.data$prob_ratio, n_post_samples),
                prob_ratio = pmax(.data$prob_ratio, 1 / n_post_samples)
                ) |>
            data.frame()

        post_prob_slope <- post_diff |>
            dplyr::group_by(.data$time) |>
            dplyr::summarise(
                post_prob = stats::quantile(.data$epred_diff, probs = 0.5),
                lower = stats::quantile(.data$epred_diff, probs = lower_q),
                upper = stats::quantile(.data$epred_diff, probs = upper_q)
                ) |>
            dplyr::ungroup()

        results <- dplyr::left_join(prob_y_above, post_prob_slope, by = "time")

    }

    return (results)

}

#' @export
plot.clusters_results <- function (x, clusters_y = -Inf, clusters_colour = "black", lineend = "butt", ...) {

    # retrieve the empirical data
    emp_data <- x$model$data

    # group-level or participant-level clusters?
    # group_level <- ifelse(test = ncol(x$clusters) == 5, yes = FALSE, no = TRUE)
    group_level <- ifelse(
        test = "participant" %in% colnames(x$clusters),
        yes = FALSE, no = TRUE
        )

    # reconstruct a raw data time course when possible
    if (!is.null(x$multilevel) && "outcome_mean" %in% names(emp_data) ) {

        if ("predictor" %in% names(emp_data) ) {

            cond1 <- levels(emp_data$predictor)[1]
            cond2 <- levels(emp_data$predictor)[2]

            reshaped_data <- emp_data |>
                dplyr::summarise(
                    outcome_mean = mean(.data$outcome_mean),
                    .by = c(.data$time, .data$predictor)
                    ) |>
                tidyr::pivot_wider(
                    names_from  = .data$predictor,
                    values_from = .data$outcome_mean
                    ) |>
                dplyr::mutate(outcome_mean = .data[[cond2]] - .data[[cond1]])

        } else {

            if (group_level) {

                reshaped_data <- emp_data |>
                    dplyr::summarise(
                        outcome_mean = mean(.data$outcome_mean),
                        .by = .data$time
                        )

            } else {

                reshaped_data <- emp_data |>
                    dplyr::summarise(
                        outcome_mean = mean(.data$outcome_mean),
                        .by = c(.data$participant, .data$time)
                        )

            }

        }

    } else {

        reshaped_data <- NULL

    }

    p <- ggplot2::ggplot(
        data = if (is.null(reshaped_data) ) x$predictions else reshaped_data,
        ggplot2::aes(
            x = .data$time,
            y = if (!is.null(reshaped_data) ) .data$outcome_mean else .data$post_prob)
        ) +
        ggplot2::geom_hline(yintercept = 0.0, linetype = 2) +
        ggplot2::geom_ribbon(
            data = x$predictions,
            ggplot2::aes(x = .data$time, y = NULL, ymin = .data$lower, ymax = .data$upper),
            fill = clusters_colour, alpha = 0.2
            ) +
        ggplot2::geom_line(
            data = x$predictions,
            ggplot2::aes(x = .data$time, y = .data$post_prob),
            colour = clusters_colour,
            linewidth = 1
            )

    if (!is.null(reshaped_data) ) {

        p <- p + ggplot2::geom_line(linewidth = 0.5)

    }

    p <- p +
        ggplot2::geom_segment(
            data = x$clusters,
            ggplot2::aes(
                x = .data$cluster_onset,
                xend = .data$cluster_offset,
                y = clusters_y,
                yend = clusters_y
                ),
            colour = clusters_colour,
            inherit.aes = FALSE,
            lineend = lineend,
            linewidth = 5
            ) +
        ggplot2::theme_bw() +
        ggplot2::labs(x = "Time", y = "Observed and predicted effect")

    # if clusters are available at the participant level
    if ("participant" %in% colnames(x$clusters) ) {

        p + ggplot2::facet_wrap(~participant, scales = "free")

    } else {

        p

    }

}

#' Print method for \code{clusters_results} objects
#'
#' This method provides a concise console representation of the output from
#' \code{\link{testing_through_time}}, including the number of detected
#' clusters and a compact table summarising each cluster's onset, offset,
#' and duration. Values are rounded for readability.
#'
#' @param x An object of class \code{"clusters_results"} as returned by
#'   \code{\link{testing_through_time}}.
#' @param digits Integer; number of decimal places used when printing numeric
#'   values (default: \code{3}).
#' @param ... Additional arguments (currently ignored).
#'
#' @details
#' The printed cluster table includes:
#' \itemize{
#'   \item \code{cluster_id}: numeric identifier of the cluster;
#'   \item \code{cluster_onset}: estimated temporal onset of the cluster;
#'   \item \code{cluster_offset}: estimated temporal offset of the cluster;
#'   \item \code{duration}: duration of the cluster, computed as
#'     \code{cluster_offset - cluster_onset}.
#' }
#'
#' If no clusters exceed the posterior odds threshold, an informative message
#' is displayed and no table is printed.
#'
#' @return The input object \code{x}, returned invisibly.
#'
#' @seealso \code{\link{summary.clusters_results}},
#'   \code{\link{testing_through_time}}
#'
#' @export
print.clusters_results <- function (x, digits = 3, ...) {

    cat("\n==== Time-resolved GAMM results ===============================\n\n")

    # number of clusters
    n_clust <- nrow(x$clusters)
    # cat("Clusters found: ", n_clust, "\n\n", sep = "")
    cat("Clusters found: ", "\n\n", sep = "")

    # if no clusters, stop early
    if (n_clust == 0) {

        cat("\nNo clusters exceed the threshold.\n\n")

        return (invisible(x) )

    }

    # prepare cluster table
    if ("participant" %in% colnames(x$clusters) ) {

        clust_tbl <- x$clusters |>
            dplyr::mutate(
                cluster_onset = round(.data$cluster_onset,  digits),
                cluster_offset = round(.data$cluster_offset, digits),
                duration = round(.data$cluster_offset - .data$cluster_onset, digits)
                ) |>
            dplyr::select(.data$participant, .data$cluster_id, .data$cluster_onset, .data$cluster_offset, .data$duration)

    } else {

        clust_tbl <- x$clusters |>
            dplyr::mutate(
                cluster_onset = round(.data$cluster_onset,  digits),
                cluster_offset = round(.data$cluster_offset, digits),
                duration = round(.data$cluster_offset - .data$cluster_onset, digits)
                ) |>
            dplyr::select(.data$cluster_id, .data$cluster_onset, .data$cluster_offset, .data$duration)

    }

    # print nicely
    print(data.frame(clust_tbl), row.names = FALSE)
    cat("\n=================================================================\n")
    invisible(x)

}

#' Summary method for \code{clusters_results} objects
#'
#' Produces a detailed textual summary of a time-resolved Bayesian GAMM
#' analysis, including model metadata, the number of clusters detected, and
#' descriptive statistics of cluster durations. A rounded cluster table is
#' also printed.
#'
#' @param object An object of class \code{"clusters_results"} created by
#'   \code{\link{testing_through_time}}.
#' @param digits Integer; number of decimal places used when printing numeric
#'   values (default: \code{3}).
#' @param ... Additional arguments (currently ignored).
#'
#' @details
#' The summary prints:
#' \itemize{
#'   \item the model type used (\code{"full"}, \code{"summary"},
#'     or \code{"group"});
#'   \item the class of the underlying \pkg{brms} model and the number of
#'     posterior draws;
#'   \item the number of clusters detected by the posterior odds threshold;
#'   \item descriptive statistics for cluster durations (minimum, maximum,
#'     mean, median, and total duration);
#'   \item a neatly formatted table listing each cluster's onset, offset,
#'     and duration.
#' }
#'
#' If no clusters were detected, the function prints a message and returns
#' invisibly.
#'
#' @return The object \code{object}, returned invisibly.
#'
#' @seealso \code{\link{print.clusters_results}},
#'   \code{\link{testing_through_time}}
#'
#' @export
summary.clusters_results <- function (object, digits = 3, ...) {

    cat("\n==== Time-resolved GAMM results ===============================\n\n")

    # model type
    if (!is.null(object$multilevel) ) {

        cat("Model type: ", object$multilevel, "\n", sep = "")

    }

    # class of backend model
    if (!is.null(object$model) ) {

        cat("Backend model: ", class(object$model)[1], "\n", sep = "")
        cat("Posterior draws: ", brms::ndraws(object$model), "\n", sep = "")

    }

    # number of clusters
    n_clust <- nrow(object$clusters)

    # no clusters, simple summary
    if (n_clust == 0) {

        cat("\nNo clusters exceeded the threshold.\n\n")

        return (invisible(object) )

    }

    # compute durations
    cl <- object$clusters |>
        dplyr::mutate(duration = .data$cluster_offset - .data$cluster_onset)

    # basic cluster stats
    cat("\nCluster statistics:\n")
    cat("  Mean cluster duration: ",
        round(mean(cl$duration), digits), "\n", sep = "")
    cat("  Median cluster duration: ",
        round(stats::median(cl$duration), digits), "\n", sep = "")
    cat("  Min cluster duration: ",
        round(min(cl$duration), digits), "\n", sep = "")
    cat("  Max cluster duration: ",
        round(max(cl$duration), digits), "\n", sep = "")

    # prepare cluster table
    if ("participant" %in% colnames(object$clusters) ) {

        cl_print <- cl |>
            dplyr::mutate(
                cluster_onset = round(.data$cluster_onset,  digits),
                cluster_offset = round(.data$cluster_offset, digits),
                duration = round(.data$duration, digits)
                ) |>
            dplyr::select(
                .data$participant, .data$cluster_id, .data$cluster_onset,
                .data$cluster_offset, .data$duration
                )

    } else {

        cl_print <- cl |>
            dplyr::mutate(
                cluster_onset = round(.data$cluster_onset,  digits),
                cluster_offset = round(.data$cluster_offset, digits),
                duration = round(.data$duration, digits)
                ) |>
            dplyr::select(
                .data$cluster_id, .data$cluster_onset,
                .data$cluster_offset, .data$duration
                )

    }

    cat("\nCluster table:\n\n")
    print(data.frame(cl_print), row.names = FALSE)
    cat("\n=================================================================\n")
    invisible(object)

}

#' Posterior predictive checks
#'
#' Generate posterior predictive checks (PPCs) from a fitted Bayesian
#' time-resolved GAMM stored in a \code{clusters_results} object.
#' PPCs can be produced either at the group level or separately for each participant.
#'
#' At the group level, predictions are obtained by simulating from the posterior
#' using \code{\link[brms]{posterior_predict}} with \code{re_formula = NA},
#' after collapsing the original data across participants (by time).
#' At the participant level, PPCs are generated using
#' \code{\link[brms]{pp_check}} with grouped ribbons.
#'
#' @param object A \code{clusters_results} object containing a fitted
#'   \code{\link[brms]{brmsfit}} model in \code{object$model}.
#' @param ppc_type Character string specifying the type of PPC to generate.
#'   Either \code{"group"} (default) for group-level PPCs (ignoring participant
#'   identity) or \code{"participant"} for participant-wise PPCs.
#' @param ndraws Integer specifying the number of posterior draws to use for
#'   the PPC. Defaults to 500.
#' @param group_var Optional character; name of the grouping variable to use for
#' grouped PPCs at the group level. If NULL (default), the function uses
#' "predictor" when present in model$data and binary (two levels).
#' @param cores Numeric; number of parallel cores to use (only used when \code{ppc_type = "participant"}).
#' @param ... Currently unused. Included for future extensions.
#'
#' @details
#' \itemize{
#'   \item \strong{Group-level PPCs} are computed by averaging numeric variables
#'   across participants at each time point, and simulating posterior predictive
#'   draws with random effects excluded (\code{re_formula = NA}).
#'   This provides a marginal, population-level posterior predictive check.
#'
#'   \item \strong{Participant-level PPCs} are computed using grouped ribbon
#'   plots, showing posterior predictive distributions separately for each
#'   participant.
#' }
#'
#' The returned object is a \code{ggplot2} object produced by
#' \code{\link[bayesplot]{ppc_ribbon}} or \code{\link[brms]{pp_check}},
#' depending on the selected \code{ppc_type}.
#'
#' @return
#' A \code{ggplot} object visualising the posterior predictive check.
#' The plot is printed to the active graphics device and also returned invisibly.
#'
#' @seealso
#' \code{\link[brms]{pp_check}},
#' \code{\link[brms]{posterior_predict}},
#' \code{\link[bayesplot]{ppc_ribbon}}
#'
#' @examples
#' \dontrun{
#' # Group-level PPC
#' ppc(object = res, ppc_type = "group")
#'
#' # Participant-level PPC
#' ppc(object = res, ppc_type = "participant")
#' }
#'
#' @export
ppc <- function (
        object,
        ppc_type = c("group", "participant"),
        ndraws = 500,
        group_var = NULL,
        cores = 4,
        ...
        ) {

    fit <- object$model
    ppc_type <- match.arg(ppc_type)

    if (!inherits(object, "clusters_results") ) {

        stop ("`object` must be of class 'clusters_results'.", call. = FALSE)

    }

    if (is.null(fit) || !inherits(fit, "brmsfit") ) {

        stop ("`object$model` must be a valid 'brmsfit' object.", call. = FALSE)

    }

    # determine grouping variable (if any and if NULL)
    data_fit <- fit$data

    if (is.null(group_var) ) {

        if ("predictor" %in% names(data_fit) ) {

            g <- data_fit[["predictor"]]

            # coerce to factor for level checking
            if (!is.factor(g) ) {

                g <- factor(g)

            }

            if (nlevels(g) == 2) {

                group_var <- "predictor"

            } else {

                group_var <- NULL

            }

        } else {

            group_var <- NULL

        }

    } else {

        if (!group_var %in% names(data_fit) ) {

            stop (
                "Specified `group_var` '", group_var,
                "' not found in model$data.",
                call. = FALSE
                )

        }

    }

    # PPC per group
    if (ppc_type == "group") {

        if (is.null(group_var) ) {

            # grid for group-level prediction
            newdata <- data_fit |>
                dplyr::summarise(dplyr::across(dplyr::where(is.numeric), mean, na.rm = TRUE), .by = .data$time) |>
                dplyr::mutate(participant = NA)

            # simulate from posterior at the group level
            yrep <- brms::posterior_predict(
                object = fit,
                newdata = newdata,
                re_formula = NA,
                ndraws = ndraws
                )

            # observed y on that same grid
            y_obs <- newdata$outcome_mean
            x_time <- newdata$time

            # ribbon PPC
            ppc_plot <- bayesplot::ppc_ribbon(
                x = x_time,
                y = y_obs,
                yrep = yrep,
                prob = 0.5,
                prob_outer = 0.5,
                alpha = 0.5
                ) +
                ggplot2::theme_bw() +
                ggplot2::labs(x = "Time (s)")

        } else {

            # grid for group-level prediction
            newdata <- data_fit |>
                dplyr::summarise(
                    dplyr::across(dplyr::where(is.numeric), mean, na.rm = TRUE),
                    .by = c(.data$predictor, .data$time)
                    ) |>
                dplyr::mutate(participant = NA) |>
                dplyr::arrange(.data$predictor, .data$time)

            # simulate from posterior at the group level
            yrep <- brms::posterior_predict(
                object = fit,
                newdata = newdata,
                re_formula = NA,
                ndraws = ndraws
                )

            # observed y on that same grid
            y_obs <- data_fit |>
                dplyr::group_by(.data$time, .data$predictor) |>
                dplyr::summarise(y = mean(.data$outcome_mean, na.rm = TRUE), .groups = "drop") |>
                dplyr::arrange(.data$predictor, .data$time) |>
                dplyr::pull(.data$y)

            # x-axis timesteps
            x_time <- newdata$time

            # ribbon PPC
            ppc_plot <- bayesplot::ppc_ribbon_grouped(
                x = x_time,
                y = y_obs,
                group = newdata$predictor,
                yrep = yrep,
                prob = 0.5,
                prob_outer = 0.5,
                alpha = 0.5
                ) +
                ggplot2::theme_bw() +
                ggplot2::labs(x = "Time (s)")

        }

    } else { # or PPC per participant

        ppc_plot <- brms::pp_check(
            object = fit,
            ndraws = ndraws,
            type = "ribbon_grouped",
            x = "time",
            group = "participant",
            prob = 0.5,
            prob_outer = 0.5,
            alpha = 0.5,
            cores = cores
            ) +
            ggplot2::theme_bw() +
            ggplot2::labs(x = "Time (s)")

    }

    # returning the plot
    print(ppc_plot)
    invisible(ppc_plot)

}
