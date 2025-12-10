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
#'   containing a \emph{binary} predictor (e.g., group or condition). If
#'   \code{predictor_id = NA}, the function tests whether the outcome differs
#'   from \code{0 + chance_level + sesoi} over time.
#' @param family A \pkg{brms} family object describing the response
#'   distribution to be used in the model (defaults to \code{gaussian()}).
#' @param kvalue Numeric; basis dimension \code{k} passed to the smooth term
#'   \code{s(time, ..., k = kvalue)}.
#' @param bs Character; Character scalar; type of spline basis to be used by \pkg{brms}
#'   (passed to \code{s()}, e.g., \code{"tp"} for thin-plate splines).
#' @param multilevel Character; which model to fit. One of
#'   \itemize{
#'     \item \code{"full"}: Full GAMM with participant-level random/varying effects;
#'     \item \code{"summary"}: GAMM fitted to participant-level summary
#'       statistics (mean outcome and its standard deviation);
#'     \item \code{"group"}: Group-level GAM fitted to participant-averaged
#'       data (no random/varying effects).
#'   }
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
#'   smoothed time course and detected clusters.
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
        multilevel = c("summary", "full", "group"),
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

    if (multilevel == "full") {

        # construct the smooth term dynamically
        if (is.na(predictor_id) ) {

            smooth_term <- glue::glue("s({time_id}, bs = '{bs}', k = {kvalue})")

        } else {

            smooth_term <- glue::glue("s({time_id}, bs = '{bs}', k = {kvalue}, by = {predictor_id})")

        }

        # full formula
        if (is.na(predictor_id) ) {

            formula_str <- glue::glue("{outcome_id} ~ 1 + {smooth_term} + (1 | {participant_id})")

        } else {

            formula_str <- glue::glue("{outcome_id} ~ 1 + {predictor_id} + {smooth_term} + (1 + {predictor_id} | {participant_id})")

        }

        # convert to formula
        formula_obj <- brms::bf(formula_str)

        # fitting the model
        brms_gam <- brms::brm(
            formula = formula_obj,
            data = data,
            family = family,
            warmup = warmup,
            iter = iter,
            chains = chains,
            cores = cores,
            backend = backend
            )

        if (is.null(n_post_samples) ) {

            n_post_samples <- brms::ndraws(brms_gam)

        }

        # computing the posterior odds over time
        if (is.na(predictor_id) ) {

            # retrieving posterior predictions (draws)
            post_draws <- brms_gam$data |>
                tidybayes::add_epred_draws(
                    object = brms_gam,
                    ndraws = n_post_samples
                    ) |>
                data.frame()

            # # computing the posterior odds over time
            # prob_y_above <- post_draws |>
            #     dplyr::select(.data$participant, .data$time, .data$.epred, .data$.draw) |>
            #     # computing posterior probability at the group level
            #     dplyr::group_by(.data$time) |>
            #     dplyr::summarise(
            #         prob_above = mean(.data$.epred > (0 + chance_level + sesoi) )
            #         ) |>
            #     dplyr::mutate(prob_ratio = .data$prob_above / (1 - .data$prob_above) ) |>
            #     dplyr::ungroup() |>
            #     # ensuring there is no 0 or +Inf values
            #     dplyr::mutate(prob_ratio = pmin(.data$prob_ratio, n_post_samples) ) |>
            #     dplyr::mutate(prob_ratio = pmax(.data$prob_ratio, 1 / n_post_samples) ) |>
            #     data.frame()
            #
            # # retrieving posterior predictions
            # post_prob_slope <- post_draws |>
            #     dplyr::select(.data$participant, .data$time,  .data$.epred, .data$.draw) |>
            #     # computing mean posterior probability at the group level
            #     dplyr::group_by(.data$time) |>
            #     dplyr::summarise(
            #         post_prob = stats::quantile(x = .data$.epred, probs = 0.5),
            #         lower = stats::quantile(x = .data$.epred, probs = 0.025),
            #         upper = stats::quantile(x = .data$.epred, probs = 0.975)
            #         ) |>
            #     dplyr::ungroup()
            #
            # # joining with prob_y_above
            # prob_y_above <- dplyr::left_join(prob_y_above, post_prob_slope, by = "time")

        } else {

            # retrieving posterior predictions (draws)
            post_draws <- tidybayes::add_epred_draws(
                object = brms_gam,
                newdata = tidyr::crossing(
                    # time = brms_gam$data[[time_id]],
                    # predictor = brms_gam$data[[predictor_id]]
                    time = sort(unique(data[[time_id]])),
                    predictor = sort(unique(data[[predictor_id]]))
                    ),
                ndraws = n_post_samples,
                re_formula = NA
                ) |>
                data.frame()

            # # retrieving predictor labels
            # cond1 <- levels(post_draws$predictor)[1]
            # cond2 <- levels(post_draws$predictor)[2]
            #
            # # computing the posterior odds over time
            # prob_y_above <- post_draws |>
            #     dplyr::select(.data$time, .data$predictor, .data$.epred, .data$.draw) |>
            #     tidyr::pivot_wider(names_from = .data$predictor, values_from = .data$.epred) |>
            #     dplyr::mutate(epred_diff = .data[[cond2]] - .data[[cond1]]) |>
            #     # computing posterior probability at the group level
            #     dplyr::group_by(.data$time) |>
            #     dplyr::summarise(
            #         prob_above = mean(.data$epred_diff > (0 + chance_level + sesoi) )
            #         ) |>
            #     dplyr::mutate(prob_ratio = .data$prob_above / (1 - .data$prob_above) ) |>
            #     dplyr::ungroup() |>
            #     # ensuring there is no 0 or +Inf values
            #     dplyr::mutate(prob_ratio = pmin(.data$prob_ratio, n_post_samples) ) |>
            #     dplyr::mutate(prob_ratio = pmax(.data$prob_ratio, 1 / n_post_samples) ) |>
            #     data.frame()
            #
            # # retrieving posterior predictions
            # post_prob_slope <- post_draws |>
            #     dplyr::select(.data$time, .data$predictor, .data$.epred, .data$.draw) |>
            #     tidyr::pivot_wider(names_from = .data$predictor, values_from = .data$.epred) |>
            #     dplyr::mutate(epred_diff = .data[[cond2]] - .data[[cond1]]) |>
            #     # computing mean posterior probability at the group level
            #     dplyr::group_by(.data$time) |>
            #     dplyr::summarise(
            #         post_prob = stats::quantile(x = .data$epred_diff, probs = 0.5),
            #         lower = stats::quantile(x = .data$epred_diff, probs = 0.025),
            #         upper = stats::quantile(x = .data$epred_diff, probs = 0.975)
            #         ) |>
            #     dplyr::ungroup()
            #
            # # joining with prob_y_above
            # prob_y_above <- dplyr::left_join(prob_y_above, post_prob_slope, by = "time")

        }

    } else if (multilevel == "summary") {

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

            # define the smooth term
            smooth_term <- glue::glue("s(time, bs = '{bs}', k = {kvalue})")

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

            # define the smooth term
            smooth_term <- glue::glue("s(time, bs = '{bs}', k = {kvalue}, by = predictor)")

        }

        # full formula
        if (is.na(predictor_id) ) {

            formula_str <- glue::glue("outcome_mean | se(outcome_sd) ~ 1 + {smooth_term} + (1 | participant)")

        } else {

            formula_str <- glue::glue("outcome_mean | se(outcome_sd) ~ 1 + predictor + {smooth_term} + (1 + predictor | participant)")

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
            backend = backend
            )

        if (is.null(n_post_samples) ) {

            n_post_samples <- brms::ndraws(brms_gam)

        }

        # computing the posterior odds over time
        if (is.na(predictor_id) ) {

            # # retrieving posterior predictions (draws)
            # post_draws <- brms_gam$data |>
            #     tidybayes::add_epred_draws(
            #         object = brms_gam,
            #         ndraws = n_post_samples
            #         ) |>
            #     data.frame()

            # newdata grid over time
            newdata_grid <- tidyr::crossing(
                time = sort(unique(brms_gam$data$time) )
                # predictor = levels(brms_gam$data$predictor)
                ) |>
                # dummy outcome_sd to satisfy | se(outcome_sd)
                dplyr::mutate(outcome_sd = 1)

            post_draws <- tidybayes::add_epred_draws(
                object = brms_gam,
                newdata = newdata_grid,
                ndraws = n_post_samples,
                re_formula = NA
                ) |>
                data.frame()

            # # computing the posterior odds over time
            # prob_y_above <- post_draws |>
            #     dplyr::select(.data$participant, .data$time, .data$.epred, .data$.draw) |>
            #     # computing posterior probability at the group level
            #     dplyr::group_by(.data$time) |>
            #     dplyr::summarise(
            #         prob_above = mean(.data$.epred > (0 + chance_level + sesoi) )
            #         ) |>
            #     dplyr::mutate(prob_ratio = .data$prob_above / (1 - .data$prob_above) ) |>
            #     dplyr::ungroup() |>
            #     # ensuring there is no 0 or +Inf values
            #     dplyr::mutate(prob_ratio = pmin(.data$prob_ratio, n_post_samples) ) |>
            #     dplyr::mutate(prob_ratio = pmax(.data$prob_ratio, 1 / n_post_samples) ) |>
            #     data.frame()
            #
            # # retrieving posterior predictions
            # post_prob_slope <- post_draws |>
            #     dplyr::select(.data$participant, .data$time,  .data$.epred, .data$.draw) |>
            #     # computing mean posterior probability at the group level
            #     dplyr::group_by(.data$time) |>
            #     dplyr::summarise(
            #         post_prob = stats::quantile(x = .data$.epred, probs = 0.5, na.rm = TRUE),
            #         lower = stats::quantile(x = .data$.epred, probs = 0.025, na.rm = TRUE),
            #         upper = stats::quantile(x = .data$.epred, probs = 0.975, na.rm = TRUE)
            #         ) |>
            #     dplyr::ungroup()
            #
            # # joining with prob_y_above
            # prob_y_above <- dplyr::left_join(prob_y_above, post_prob_slope, by = "time")

        } else {

            # newdata grid over time × predictor
            newdata_grid <- tidyr::crossing(
                time = sort(unique(brms_gam$data$time) ),
                predictor = levels(brms_gam$data$predictor)
                ) |>
                # dummy outcome_sd to satisfy | se(outcome_sd)
                dplyr::mutate(outcome_sd = 1)

            post_draws <- tidybayes::add_epred_draws(
                object = brms_gam,
                newdata = newdata_grid,
                ndraws = n_post_samples,
                re_formula = NA
                ) |>
                data.frame()

            # # retrieving predictor labels
            # cond1 <- levels(brms_gam$data$predictor)[1]
            # cond2 <- levels(brms_gam$data$predictor)[2]
            #
            # # computing the posterior odds over time
            # prob_y_above <- post_draws |>
            #     # dplyr::select(.data$time, .data$participant, .data$predictor, .data$.epred, .data$.draw) |>
            #     dplyr::select(.data$time, .data$predictor, .data$.epred, .data$.draw) |>
            #     tidyr::pivot_wider(names_from = .data$predictor, values_from = .data$.epred) |>
            #     dplyr::mutate(epred_diff = .data[[cond2]] - .data[[cond1]]) |>
            #     # computing posterior probability at the group level
            #     dplyr::group_by(.data$time) |>
            #     dplyr::summarise(
            #         prob_above = mean(.data$epred_diff > (0 + chance_level + sesoi) )
            #         ) |>
            #     dplyr::mutate(prob_ratio = .data$prob_above / (1 - .data$prob_above) ) |>
            #     dplyr::ungroup() |>
            #     # ensuring there is no 0 or +Inf values
            #     dplyr::mutate(prob_ratio = pmin(.data$prob_ratio, n_post_samples) ) |>
            #     dplyr::mutate(prob_ratio = pmax(.data$prob_ratio, 1 / n_post_samples) ) |>
            #     data.frame()
            #
            # # retrieving group-level posterior predictions
            # post_prob_slope <- post_draws |>
            #     # dplyr::select(.data$time, .data$participant, .data$predictor, .data$.epred, .data$.draw) |>
            #     dplyr::select(.data$time, .data$predictor, .data$.epred, .data$.draw) |>
            #     tidyr::pivot_wider(names_from = .data$predictor, values_from = .data$.epred) |>
            #     dplyr::mutate(epred_diff = .data[[cond2]] - .data[[cond1]]) |>
            #     # computing mean posterior probability at the group level
            #     dplyr::group_by(.data$time) |>
            #     dplyr::summarise(
            #         post_prob = stats::quantile(x = .data$epred_diff, probs = 0.5),
            #         lower = stats::quantile(x = .data$epred_diff, probs = 0.025),
            #         upper = stats::quantile(x = .data$epred_diff, probs = 0.975)
            #         ) |>
            #     dplyr::ungroup()
            #
            # # joining with prob_y_above
            # prob_y_above <- dplyr::left_join(prob_y_above, post_prob_slope, by = "time")

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

            # define the smooth term
            smooth_term <- glue::glue("s(time, bs = '{bs}', k = {kvalue}, by = predictor)")

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
            backend = backend
            )

        if (is.null(n_post_samples) ) {

            n_post_samples <- brms::ndraws(brms_gam)

        }

        # computing the posterior odds over time
        if (is.na(predictor_id) ) {

            # retrieving posterior predictions (draws)
            post_draws <- brms_gam$data |>
                tidybayes::add_epred_draws(
                    object = brms_gam,
                    ndraws = n_post_samples
                    ) |>
                data.frame()

            # prob_y_above <- post_draws |>
            #     dplyr::select(.data$time, .data$.epred, .data$.draw) |>
            #     dplyr::group_by(.data$time) |>
            #     dplyr::summarise(
            #         prob_above = mean(.data$.epred > (0 + chance_level + sesoi) )
            #         ) |>
            #     dplyr::mutate(prob_ratio = .data$prob_above / (1 - .data$prob_above) ) |>
            #     dplyr::ungroup() |>
            #     # ensuring there is no 0 or +Inf values
            #     dplyr::mutate(prob_ratio = pmin(.data$prob_ratio, n_post_samples) ) |>
            #     dplyr::mutate(prob_ratio = pmax(.data$prob_ratio, 1 / n_post_samples) ) |>
            #     data.frame()
            #
            # # retrieving posterior predictions
            # post_prob_slope <- post_draws |>
            #     dplyr::select(.data$time, .data$.epred, .data$.draw) |>
            #     # computing mean posterior probability at the group level
            #     dplyr::group_by(.data$time) |>
            #     dplyr::summarise(
            #         post_prob = stats::quantile(x = .data$.epred, probs = 0.5),
            #         lower = stats::quantile(x = .data$.epred, probs = 0.025),
            #         upper = stats::quantile(x = .data$.epred, probs = 0.975)
            #         ) |>
            #     dplyr::ungroup()
            #
            # # joining with prob_y_above
            # prob_y_above <- dplyr::left_join(prob_y_above, post_prob_slope, by = "time")

        } else {

            # retrieving posterior predictions (draws)
            post_draws <- tidybayes::add_epred_draws(
                object = brms_gam,
                newdata = tidyr::crossing(
                    time = sort(unique(brms_gam$data$time) ),
                    predictor = unique(brms_gam$data$predictor)
                    ),
                ndraws = n_post_samples
                ) |>
                data.frame()

            # # retrieving predictor labels
            # cond1 <- levels(post_draws$predictor)[1]
            # cond2 <- levels(post_draws$predictor)[2]
            #
            # # computing the posterior odds over time
            # prob_y_above <- post_draws |>
            #     dplyr::select(.data$time, .data$predictor, .data$.epred, .data$.draw) |>
            #     tidyr::pivot_wider(names_from = .data$predictor, values_from = .data$.epred) |>
            #     # head()
            #     dplyr::mutate(epred_diff = .data[[cond2]] - .data[[cond1]]) |>
            #     # computing posterior probability at the group level
            #     dplyr::group_by(.data$time) |>
            #     dplyr::summarise(
            #         prob_above = mean(.data$epred_diff > (0 + chance_level + sesoi) )
            #         ) |>
            #     dplyr::mutate(prob_ratio = .data$prob_above / (1 - .data$prob_above) ) |>
            #     dplyr::ungroup() |>
            #     # ensuring there is no 0 or +Inf values
            #     dplyr::mutate(prob_ratio = pmin(.data$prob_ratio, n_post_samples) ) |>
            #     dplyr::mutate(prob_ratio = pmax(.data$prob_ratio, 1 / n_post_samples) ) |>
            #     data.frame()
            #
            # # retrieving posterior predictions
            # post_prob_slope <- post_draws |>
            #     dplyr::select(.data$time, .data$predictor, .data$.epred, .data$.draw) |>
            #     tidyr::pivot_wider(names_from = .data$predictor, values_from = .data$.epred) |>
            #     dplyr::mutate(epred_diff = .data[[cond2]] - .data[[cond1]]) |>
            #     # computing mean posterior probability at the group level
            #     dplyr::group_by(.data$time) |>
            #     dplyr::summarise(
            #         post_prob = stats::quantile(x = .data$epred_diff, probs = 0.5),
            #         lower = stats::quantile(x = .data$epred_diff, probs = 0.025),
            #         upper = stats::quantile(x = .data$epred_diff, probs = 0.975)
            #         ) |>
            #     dplyr::ungroup()
            #
            # # joining with prob_y_above
            # prob_y_above <- dplyr::left_join(prob_y_above, post_prob_slope, by = "time")

        }

    }

    # compute the posterior odds
    if (is.na(predictor_id) ) {

        prob_y_above <- .compute_one_sample_prob(
            post_draws = post_draws,
            threshold = chance_level + sesoi,
            n_post_samples = n_post_samples,
            credible_interval = credible_interval
            )

    } else {

        prob_y_above <- .compute_two_sample_prob(
            post_draws = post_draws,
            threshold = chance_level + sesoi,
            n_post_samples = n_post_samples,
            credible_interval = credible_interval
            )

    }

    # finding the clusters
    clusters <- find_clusters(
        data = prob_y_above |> dplyr::select(.data$time, value = .data$prob_ratio),
        threshold = threshold
        )

    # combining the results in a list
    clusters_results <- list(
        clusters = clusters,
        predictions = prob_y_above,
        model = brms_gam,
        multilevel = multilevel
        )

    # assigning a new class to the list
    class(clusters_results) <- "clusters_results"

    # returning the clusters and posterior probabilities
    return (clusters_results)

}

.compute_one_sample_prob <- function (post_draws, threshold, n_post_samples, credible_interval) {

    # validate credible interval
    if (!is.numeric(credible_interval) ||

        length(credible_interval) != 1L ||
        credible_interval <= 0 || credible_interval >= 1) {

        stop ("`credible_interval` must be a number strictly between 0 and 1.")

    }

    alpha <- (1 - credible_interval) / 2
    lower_q <- alpha
    upper_q <- 1 - alpha

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

    return (results)

}

.compute_two_sample_prob <- function (post_draws, threshold, n_post_samples, credible_interval) {

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

    return (results)

}

#' @export
# plot.clusters_results <- function (x, clusters_y = -Inf, clusters_colour = "black", lineend = "butt", ...) {
#
#     emp_data <- x$model$data
#
#     if (ncol(emp_data) > 2 & "predictor" %in% colnames(emp_data) ) {
#
#         cond1 <- levels(emp_data$predictor)[1]
#         cond2 <- levels(emp_data$predictor)[2]
#
#         reshaped_data <- emp_data |>
#             dplyr::summarise(
#                 outcome_mean = mean(.data$outcome_mean),
#                 .by = c(.data$time, .data$predictor)
#                 ) |>
#             tidyr::pivot_wider(
#                 names_from = .data$predictor,
#                 values_from = .data$outcome_mean
#                 ) |>
#             dplyr::mutate(outcome_mean = .data[[cond2]] - .data[[cond1]])
#
#     } else {
#
#         reshaped_data <- emp_data |>
#             dplyr::summarise(outcome_mean = mean(.data$outcome_mean), .by = .data$time)
#
#     }
#
#     reshaped_data |>
#         ggplot2::ggplot(
#             ggplot2::aes(x = .data$time, y = .data$outcome_mean)
#             ) +
#         ggplot2::geom_hline(yintercept = 0.0, linetype = 2) +
#         ggplot2::geom_ribbon(
#             data = x$predictions,
#             ggplot2::aes(x = .data$time, y = NULL, ymin = .data$lower, ymax = .data$upper),
#             fill = clusters_colour, alpha = 0.2
#             ) +
#         ggplot2::geom_line(
#             data = x$predictions,
#             ggplot2::aes(x = .data$time, y = .data$post_prob),
#             colour = clusters_colour, linewidth = 1
#             ) +
#         ggplot2::geom_segment(
#             data = x$clusters,
#             ggplot2::aes(
#                 x = .data$cluster_onset,
#                 xend = .data$cluster_offset,
#                 y = clusters_y,
#                 yend = clusters_y
#                 ),
#             colour = clusters_colour,
#             inherit.aes = FALSE,
#             lineend = lineend,
#             linewidth = 5
#             ) +
#         ggplot2::geom_line(linewidth = 0.5) +
#         ggplot2::theme_bw() +
#         ggplot2::labs(x = "Time (s)", y = "Observed and predicted effect")
#
# }
plot.clusters_results <- function (x, clusters_y = -Inf, clusters_colour = "black", lineend = "butt", ...) {

    # retrieve the empirical data
    emp_data <- x$model$data

    # reconstruct a raw data time course when possible
    if (!is.null(x$multilevel) && x$multilevel %in% c("summary", "group") && "outcome_mean" %in% names(emp_data) ) {

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

            reshaped_data <- emp_data |>
                dplyr::summarise(
                    outcome_mean = mean(.data$outcome_mean),
                    .by = .data$time
                    )

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

    p +
        ggplot2::geom_segment(
            data = x$clusters,
            ggplot2::aes(
                x    = .data$cluster_onset,
                xend = .data$cluster_offset,
                y    = clusters_y,
                yend = clusters_y
                ),
            colour = clusters_colour,
            inherit.aes = FALSE,
            lineend = lineend,
            linewidth = 5
            ) +
        ggplot2::theme_bw() +
        ggplot2::labs(x = "Time", y = "Observed and predicted effect")

}
