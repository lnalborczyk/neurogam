#' Time-resolved Bayesian testing
#'
#' Time-resolved testing using Bayesian generalised additive multilevel regression models (BGAMMs).
#'
#' @param data Dataframe, with ERP or decoding performance data (in long format).
#' @param participant_id Character, column in data specifying the participant ID.
#' @param outcome_id Character, column in data specifying the M/EEG values.
#' @param time_id Character, column in data specifying the time information.
#' @param predictor_id Character, column in data specifying the binary predictor (e.g., group or condition). If predictor_id is NA, then neurogam will test difference against 0.
#' @param family A description of the response distribution to be used in the model.
#' @param kvalue Numeric, GAM basis dimension.
#' @param bs Character, type of splines to be passed to brms.
#' @param multilevel Character, should we fit the full GAMM ("full"), GAMM with summary statistics per participant ("summary"), or group-level GAM ("group").
#' @param warmup Numeric, number of warm-up iterations per chain.
#' @param iter Numeric, number of total iterations per chain.
#' @param chains Numeric, number of MCMCs.
#' @param cores Numeric, number of parallel cores to use.
#' @param backend Character, package to use as the backend for fitting the Stan model. Options are "cmdstanr" (the default) or "rstan".
#' @param threshold Numeric, threshold for identifying clusters.
#' @param n_post_samples Numeric, number of posterior samples to use to compute the posterior probability ratio (the larger the better).
#' @param chance_level Numeric, chance level (for analysing time-resolved decoding performance).
#' @param sesoi Numeric, smallest effect size of interest (defaults to 0).
#'
#' @return A list containing the identified clusters (i.e., onset and offset) and time-series of posterior probabilities.
#'
#' @importFrom rlang .data
#' @importFrom stats gaussian
#'
#' @examples
#' \dontrun{
#' # importing some simulated EEG data
#' data(eeg_data)
#' head(eeg_data)
#'
#' # fitting the BGAMM to identify clusters
#' results <- testing_through_time(data = eeg_data)
#'
#' # displaying the identified clusters
#' print(results$clusters)
#'
#' # plotting the GAM-smoothed signal and identified clusters
#' plot(results)
#' }
#'
#' @author Ladislas Nalborczyk \email{ladislas.nalborczyk@@cnrs.fr}.
#'
#' @export

testing_through_time <- function (
        data,
        participant_id = "participant", outcome_id = "eeg",
        time_id = "time", predictor_id = "condition",
        family = gaussian(), kvalue = 20, bs = "tp",
        multilevel = c("summary", "group", "full"),
        warmup = 1000, iter = 2000, chains = 4, cores = 4,
        backend = "cmdstanr",
        threshold = 10, n_post_samples = NULL,
        chance_level = 0, sesoi = 0
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

    # checking required column names
    required_columns <- c(participant_id, outcome_id, time_id, predictor_id)
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
                tidybayes::add_epred_draws(object = brms_gam, ndraws = n_post_samples) |>
                data.frame()

            # computing the posterior odds over time
            prob_y_above <- post_draws |>
                dplyr::select(.data$participant, .data$time, .data$.epred, .data$.draw) |>
                # computing posterior probability at the group level
                dplyr::group_by(.data$time) |>
                dplyr::summarise(
                    prob_above = mean(.data$.epred > (0 + chance_level + sesoi) )
                    ) |>
                dplyr::mutate(prob_ratio = .data$prob_above / (1 - .data$prob_above) ) |>
                dplyr::ungroup() |>
                # ensuring there is no 0 or +Inf values
                dplyr::mutate(prob_ratio = pmin(.data$prob_ratio, n_post_samples) ) |>
                dplyr::mutate(prob_ratio = pmax(.data$prob_ratio, 1 / n_post_samples) ) |>
                data.frame()

            # retrieving posterior predictions
            post_prob_slope <- post_draws |>
                dplyr::select(.data$participant, .data$time,  .data$.epred, .data$.draw) |>
                # computing mean posterior probability at the group level
                dplyr::group_by(.data$time) |>
                dplyr::summarise(
                    post_prob = stats::quantile(x = .data$.epred, probs = 0.5),
                    lower = stats::quantile(x = .data$.epred, probs = 0.025),
                    upper = stats::quantile(x = .data$.epred, probs = 0.975)
                    ) |>
                dplyr::ungroup()

            # joining with prob_y_above
            prob_y_above <- dplyr::left_join(prob_y_above, post_prob_slope, by = "time")

        } else {

            # retrieving posterior predictions (draws)
            post_draws <- tidybayes::add_epred_draws(
                object = brms_gam,
                newdata = tidyr::crossing(
                    brms_gam$data$time, brms_gam$data$predictor
                    ) |>
                    dplyr::rename(time = 1, predictor = 2), ndraws = n_post_samples) |>
                data.frame()

            # retrieving predictor labels
            cond1 <- levels(post_draws$predictor)[1]
            cond2 <- levels(post_draws$predictor)[2]

            # computing the posterior odds over time
            prob_y_above <- post_draws |>
                dplyr::select(.data$time, .data$predictor, .data$.epred, .data$.draw) |>
                tidyr::pivot_wider(names_from = .data$predictor, values_from = .data$.epred) |>
                # head()
                dplyr::mutate(epred_diff = .data[[cond2]] - .data[[cond1]]) |>
                # computing posterior probability at the group level
                dplyr::group_by(.data$time) |>
                dplyr::summarise(
                    prob_above = mean(.data$epred_diff > (0 + chance_level + sesoi) )
                    ) |>
                dplyr::mutate(prob_ratio = .data$prob_above / (1 - .data$prob_above) ) |>
                dplyr::ungroup() |>
                # ensuring there is no 0 or +Inf values
                dplyr::mutate(prob_ratio = pmin(.data$prob_ratio, n_post_samples) ) |>
                dplyr::mutate(prob_ratio = pmax(.data$prob_ratio, 1 / n_post_samples) ) |>
                data.frame()

            # retrieving posterior predictions
            post_prob_slope <- post_draws |>
                dplyr::select(.data$time, .data$predictor, .data$.epred, .data$.draw) |>
                tidyr::pivot_wider(names_from = .data$predictor, values_from = .data$.epred) |>
                dplyr::mutate(epred_diff = .data[[cond2]] - .data[[cond1]]) |>
                # computing mean posterior probability at the group level
                dplyr::group_by(.data$time) |>
                dplyr::summarise(
                    post_prob = stats::quantile(x = .data$epred_diff, probs = 0.5),
                    lower = stats::quantile(x = .data$epred_diff, probs = 0.025),
                    upper = stats::quantile(x = .data$epred_diff, probs = 0.975)
                    ) |>
                dplyr::ungroup()

            # joining with prob_y_above
            prob_y_above <- dplyr::left_join(prob_y_above, post_prob_slope, by = "time")

        }

    } else if (multilevel == "summary") {

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

            # defining the smooth term
            smooth_term <- glue::glue("s(time, bs = '{bs}', k = {kvalue})")

        } else {

            # reshaping and summarising the data
            summary_data <- data |>
                # reshaping the original variables
                dplyr::mutate(predictor = as.factor(.data[[predictor_id]]) ) |>
                dplyr::mutate(participant = .data[[participant_id]]) |>
                dplyr::mutate(outcome = .data[[outcome_id]]) |>
                dplyr::mutate(time = .data[[time_id]]) |>
                # summarising per participant
                dplyr::summarise(
                    outcome_mean = mean(.data$outcome),
                    outcome_sd = stats::sd(.data$outcome),
                    .by = c(.data$participant, .data$predictor, .data$time)
                    )

            # defining the smooth term
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

        # fitting the model
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
                tidybayes::add_epred_draws(object = brms_gam, ndraws = n_post_samples) |>
                data.frame()

            # computing the posterior odds over time
            prob_y_above <- post_draws |>
                dplyr::select(.data$participant, .data$time, .data$.epred, .data$.draw) |>
                # computing posterior probability at the group level
                dplyr::group_by(.data$time) |>
                dplyr::summarise(
                    prob_above = mean(.data$.epred > (0 + chance_level + sesoi) )
                    ) |>
                dplyr::mutate(prob_ratio = .data$prob_above / (1 - .data$prob_above) ) |>
                dplyr::ungroup() |>
                # ensuring there is no 0 or +Inf values
                dplyr::mutate(prob_ratio = pmin(.data$prob_ratio, n_post_samples) ) |>
                dplyr::mutate(prob_ratio = pmax(.data$prob_ratio, 1 / n_post_samples) ) |>
                data.frame()

            # retrieving posterior predictions
            post_prob_slope <- post_draws |>
                dplyr::select(.data$participant, .data$time,  .data$.epred, .data$.draw) |>
                # computing mean posterior probability at the group level
                dplyr::group_by(.data$time) |>
                dplyr::summarise(
                    post_prob = stats::quantile(x = .data$.epred, probs = 0.5),
                    lower = stats::quantile(x = .data$.epred, probs = 0.025),
                    upper = stats::quantile(x = .data$.epred, probs = 0.975)
                    ) |>
                dplyr::ungroup()

            # joining with prob_y_above
            prob_y_above <- dplyr::left_join(prob_y_above, post_prob_slope, by = "time")

        } else {

            # retrieving posterior predictions (draws)
            post_draws <- tidybayes::add_epred_draws(
                object = brms_gam,
                # newdata = tidyr::crossing(
                #     brms_gam$data$time, brms_gam$data$predictor
                #     ) |> dplyr::rename(time = 1, predictor = 2),
                newdata = brms_gam$data,
                ndraws = n_post_samples
                ) |>
                data.frame()

            # retrieving predictor labels
            cond1 <- levels(post_draws$predictor)[1]
            cond2 <- levels(post_draws$predictor)[2]

            # computing the posterior odds over time
            prob_y_above <- post_draws |>
                # dplyr::select(.data$time, .data$predictor, .data$.epred, .data$.draw) |>
                dplyr::select(.data$time, .data$participant, .data$predictor, .data$.epred, .data$.draw) |>
                tidyr::pivot_wider(names_from = .data$predictor, values_from = .data$.epred) |>
                # head()
                dplyr::mutate(epred_diff = .data[[cond2]] - .data[[cond1]]) |>
                # computing posterior probability at the group level
                dplyr::group_by(.data$time) |>
                dplyr::summarise(
                    prob_above = mean(.data$epred_diff > (0 + chance_level + sesoi) )
                    ) |>
                dplyr::mutate(prob_ratio = .data$prob_above / (1 - .data$prob_above) ) |>
                dplyr::ungroup() |>
                # ensuring there is no 0 or +Inf values
                dplyr::mutate(prob_ratio = pmin(.data$prob_ratio, n_post_samples) ) |>
                dplyr::mutate(prob_ratio = pmax(.data$prob_ratio, 1 / n_post_samples) ) |>
                data.frame()

            # retrieving posterior predictions
            post_prob_slope <- post_draws |>
                dplyr::select(.data$time, .data$participant, .data$predictor, .data$.epred, .data$.draw) |>
                tidyr::pivot_wider(names_from = .data$predictor, values_from = .data$.epred) |>
                dplyr::mutate(epred_diff = .data[[cond2]] - .data[[cond1]]) |>
                # computing mean posterior probability at the group level
                dplyr::group_by(.data$time) |>
                dplyr::summarise(
                    post_prob = stats::quantile(x = .data$epred_diff, probs = 0.5),
                    lower = stats::quantile(x = .data$epred_diff, probs = 0.025),
                    upper = stats::quantile(x = .data$epred_diff, probs = 0.975)
                    ) |>
                dplyr::ungroup()

            # joining with prob_y_above
            prob_y_above <- dplyr::left_join(prob_y_above, post_prob_slope, by = "time")

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

            # defining the smooth term
            smooth_term <- glue::glue("s(time, bs = '{bs}', k = {kvalue})")

        } else {

            # reshaping and summarising the data
            summary_data <- data |>
                # reshaping the original variables
                dplyr::mutate(predictor = as.factor(.data[[predictor_id]]) ) |>
                dplyr::mutate(participant = .data[[participant_id]]) |>
                dplyr::mutate(outcome = .data[[outcome_id]]) |>
                dplyr::mutate(time = .data[[time_id]]) |>
                # summarising per participant
                dplyr::summarise(
                    outcome_mean = mean(.data$outcome),
                    outcome_sd = stats::sd(.data$outcome),
                    .by = c(.data$participant, .data$predictor, .data$time)
                    )

            # defining the smooth term
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

        # fitting the model
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
                tidybayes::add_epred_draws(object = brms_gam, ndraws = n_post_samples) |>
                data.frame()

            prob_y_above <- post_draws |>
                dplyr::select(.data$time, .data$.epred, .data$.draw) |>
                dplyr::group_by(.data$time) |>
                dplyr::summarise(
                    prob_above = mean(.data$.epred > (0 + chance_level + sesoi) )
                    ) |>
                dplyr::mutate(prob_ratio = .data$prob_above / (1 - .data$prob_above) ) |>
                dplyr::ungroup() |>
                # ensuring there is no 0 or +Inf values
                dplyr::mutate(prob_ratio = pmin(.data$prob_ratio, n_post_samples) ) |>
                dplyr::mutate(prob_ratio = pmax(.data$prob_ratio, 1 / n_post_samples) ) |>
                data.frame()

            # retrieving posterior predictions
            post_prob_slope <- post_draws |>
                dplyr::select(.data$time, .data$.epred, .data$.draw) |>
                # computing mean posterior probability at the group level
                dplyr::group_by(.data$time) |>
                dplyr::summarise(
                    post_prob = stats::quantile(x = .data$.epred, probs = 0.5),
                    lower = stats::quantile(x = .data$.epred, probs = 0.025),
                    upper = stats::quantile(x = .data$.epred, probs = 0.975)
                    ) |>
                dplyr::ungroup()

            # joining with prob_y_above
            prob_y_above <- dplyr::left_join(prob_y_above, post_prob_slope, by = "time")

        } else {

            # retrieving posterior predictions (draws)
            post_draws <- tidybayes::add_epred_draws(
                object = brms_gam,
                newdata = tidyr::crossing(
                    brms_gam$data$time, brms_gam$data$predictor
                    ) |> dplyr::rename(time = 1, predictor = 2),
                ndraws = n_post_samples
                ) |>
                data.frame()

            # retrieving predictor labels
            cond1 <- levels(post_draws$predictor)[1]
            cond2 <- levels(post_draws$predictor)[2]

            # computing the posterior odds over time
            prob_y_above <- post_draws |>
                dplyr::select(.data$time, .data$predictor, .data$.epred, .data$.draw) |>
                tidyr::pivot_wider(names_from = .data$predictor, values_from = .data$.epred) |>
                # head()
                dplyr::mutate(epred_diff = .data[[cond2]] - .data[[cond1]]) |>
                # computing posterior probability at the group level
                dplyr::group_by(.data$time) |>
                dplyr::summarise(
                    prob_above = mean(.data$epred_diff > (0 + chance_level + sesoi) )
                    ) |>
                dplyr::mutate(prob_ratio = .data$prob_above / (1 - .data$prob_above) ) |>
                dplyr::ungroup() |>
                # ensuring there is no 0 or +Inf values
                dplyr::mutate(prob_ratio = pmin(.data$prob_ratio, n_post_samples) ) |>
                dplyr::mutate(prob_ratio = pmax(.data$prob_ratio, 1 / n_post_samples) ) |>
                data.frame()

            # retrieving posterior predictions
            post_prob_slope <- post_draws |>
                dplyr::select(.data$time, .data$predictor, .data$.epred, .data$.draw) |>
                tidyr::pivot_wider(names_from = .data$predictor, values_from = .data$.epred) |>
                dplyr::mutate(epred_diff = .data[[cond2]] - .data[[cond1]]) |>
                # computing mean posterior probability at the group level
                dplyr::group_by(.data$time) |>
                dplyr::summarise(
                    post_prob = stats::quantile(x = .data$epred_diff, probs = 0.5),
                    lower = stats::quantile(x = .data$epred_diff, probs = 0.025),
                    upper = stats::quantile(x = .data$epred_diff, probs = 0.975)
                    ) |>
                dplyr::ungroup()

            # joining with prob_y_above
            prob_y_above <- dplyr::left_join(prob_y_above, post_prob_slope, by = "time")

        }

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
        data = brms_gam$data
        )

    # assigning a new class to the list
    class(clusters_results) <- "clusters_results"

    # returning the clusters and posterior probabilities
    return (clusters_results)

}

#' @export

plot.clusters_results <- function (x, clusters_y = -Inf, clusters_colour = "steelblue", lineend = "butt", ...) {

    if (ncol(x$data) > 2) {

        cond1 <- levels(x$data$predictor)[1]
        cond2 <- levels(x$data$predictor)[2]

        reshaped_data <- x$data |>
            dplyr::summarise(outcome_mean = mean(.data$outcome_mean), .by = c(.data$time, .data$predictor) ) |>
            tidyr::pivot_wider(names_from = .data$predictor, values_from = .data$outcome_mean) |>
            dplyr::mutate(outcome_mean = .data[[cond2]] - .data[[cond1]])

    } else {

        reshaped_data <- x$data |>
            dplyr::summarise(outcome_mean = mean(.data$outcome_mean), .by = .data$time)

    }

    reshaped_data |>
        ggplot2::ggplot(
            ggplot2::aes(x = .data$time, y = .data$outcome_mean)
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
            colour = clusters_colour, linewidth = 1
            ) +
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
        ggplot2::geom_line() +
        ggplot2::theme_bw() +
        ggplot2::labs(x = "Time (s)", y = "Difference in signal amplitude")

}
