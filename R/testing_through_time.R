#' Time-resolved Bayesian testing
#'
#' Time-resolved testing using Bayesian generalised additive multilevel regression models.
#'
#' @param data Dataframe, with ERP or decoding performance data (in long format).
#' @param participant_id Character, column in data specifying the participant ID.
#' @param meeg_id Character, column in data specifying the M/EEG values.
#' @param time_id Character, column in data specifying the time information.
#' @param predictor_id Character, column in data specifying the predictor(s).
#' @param kvalue Numeric, GAM basis dimension.
#' @param bs Character, type of splines to be passed to brms.
#' @param multilevel Logical, should we fit a simple GAM or the multilevel GAM with summary statistics per participant?
#' @param warmup Numeric, number of warm-up iterations per chain.
#' @param iter Numeric, number of total iterations per chain.
#' @param chains Numeric, number of MCMCs.
#' @param cores Numeric, number of parallel cores to use.
#' @param post_prob_ratio_threshold Numeric, threshold for identifying clusters.
#' @param chance_level Numeric, chance level (for analysing time-resolved decoding performance).
#' @param sesoi Numeric, smallest effect size of interest (defaults to 0).
#'
#' @return A list containing the identified clusters (i.e., onset and offset) and time-series of posterior probabilities.
#'
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' # importing some simulated EEG data
#' data(eeg_data)
#' head(eeg_data)
#'
#' # fitting the BGAMM to identify clusters
#' testing_through_time(data = eeg_data, post_prob_ratio_threshold = 2)
#' }
#'
#' @author Ladislas Nalborczyk \email{ladislas.nalborczyk@@gmail.com}.
#'
#' @export

testing_through_time <- function (
        data,
        participant_id = "participant", meeg_id = "eeg", time_id = "time", predictor_id = "condition",
        kvalue = 20, bs = "cr", multilevel = TRUE,
        warmup = 1000, iter = 2000, chains = 4, cores = 4,
        post_prob_ratio_threshold = 20,
        chance_level = 0, sesoi = 0
        ) {

    # some tests for variable types
    stopifnot("data must be a dataframe..." = is.data.frame(data) )
    stopifnot("kvalue must be a numeric..." = is.numeric(kvalue) )
    stopifnot("warmup must be a numeric..." = is.numeric(warmup) )
    stopifnot("iter must be a numeric..." = is.numeric(iter) )
    stopifnot("chains must be a numeric..." = is.numeric(chains) )
    stopifnot("cores must be a numeric..." = is.numeric(cores) )
    stopifnot("post_prob_ratio_threshold must be a numeric..." = is.numeric(post_prob_ratio_threshold) )
    stopifnot("chance_level must be a numeric..." = is.numeric(chance_level) )
    stopifnot("sesoi must be a numeric..." = is.numeric(sesoi) )
    stopifnot("bs must be a character..." = is.character(bs) )
    stopifnot("multilevel must be a logical..." = is.logical(multilevel) )

    # checking required column names
    required_columns <- c(participant_id, meeg_id, time_id, predictor_id)
    assertthat::assert_that(
        all(required_columns %in% colnames(data) ),
        msg = paste(
            "Missing columns:",
            paste(setdiff(required_columns, colnames(data) ), collapse = ", ")
            )
        )

    if (multilevel == TRUE) {

        # reshaping the data
        summary_data <- data |>
            # reshaping the original variables
            dplyr::mutate(predictor = as.factor(.data[[predictor_id]]) ) |>
            dplyr::mutate(participant = .data[[participant_id]]) |>
            dplyr::mutate(eeg = .data[[meeg_id]]) |>
            dplyr::mutate(time = .data[[time_id]]) |>
            # summarising per participant
            dplyr::summarise(
                eeg_mean = mean(.data$eeg),
                eeg_sd = stats::sd(.data$eeg),
                .by = c(.data$participant, .data$condition, .data$time)
                )

        # construct the smooth term dynamically
        smooth_term <- glue::glue("s(time, bs = '{bs}', k = {kvalue}, by = condition)")

        # full formula
        formula_str <- glue::glue("eeg_mean | se(eeg_sd) ~ condition + {smooth_term} + (1 | participant)")

        # convert to formula
        formula_obj <- brms::bf(formula_str)

        # fitting the GAMM
        # cat("Fitting the model...\n")
        brms_gam <- brms::brm(
            # eeg_mean | se(eeg_sd) ~
            #     condition + s(time, bs = get(bs), k = get(kvalue), by = condition) +
            #     (1 | participant),
            formula = formula_obj,
            data = summary_data,
            family = stats::gaussian(),
            warmup = warmup,
            iter = iter,
            chains = chains,
            cores = cores
            )

    } else {

        # reshaping the data
        summary_data <- data |>
            # reshaping the original variables
            dplyr::mutate(predictor = as.factor(.data[[predictor_id]]) ) |>
            dplyr::mutate(participant = .data[[participant_id]]) |>
            dplyr::mutate(eeg_mean = .data[[meeg_id]]) |>
            dplyr::mutate(time = .data[[time_id]])

        # construct the smooth term dynamically
        smooth_term <- glue::glue("s(time, bs = '{bs}', k = {kvalue}, by = condition)")

        # full formula
        formula_str <- glue::glue("eeg_mean ~ condition + {smooth_term}")

        # convert to formula
        formula_obj <- brms::bf(formula_str)

        # fitting the GAM
        # cat("Fitting the model...\n")
        brms_gam <- brms::brm(
            # eeg ~ brms::s(time, bs = bs, k = kvalue),
            formula = formula_obj,
            data = summary_data,
            family = stats::gaussian(),
            warmup = warmup,
            iter = iter,
            chains = chains,
            cores = cores
            )

    }

    # computing the posterior probability
    prob_y_above <- brms_gam$data |>
        # retrieving the posterior samples
        tidybayes::add_epred_draws(object = brms_gam, ndraws = 1e3) |>
        # converting to dataframe
        data.frame() |>
        # computing posterior probability at the group level
        dplyr::group_by(.data$time) |>
        dplyr::summarise(prob_above = mean(.data$.epred > (0 + chance_level + sesoi) ) ) |>
        dplyr::mutate(prob_ratio = .data$prob_above / (1 - .data$prob_above) ) |>
        # ensuring there is no 0 or +Inf values
        dplyr::mutate(prob_ratio = ifelse(is.infinite(.data$prob_ratio), brms::ndraws(brms_gam), .data$prob_ratio) ) |>
        dplyr::mutate(prob_ratio = ifelse(.data$prob_ratio == 0, 1 / brms::ndraws(brms_gam), .data$prob_ratio) ) |>
        dplyr::ungroup() |>
        data.frame()

    # finding the clusters
    clusters <- find_clusters(
        data = prob_y_above |> dplyr::select(-.data$prob_above),
        threshold = post_prob_ratio_threshold
        )

    # combining the results in a list
    clusters_results <- list(
        clusters = clusters,
        post_prob_timecourse = prob_y_above,
        summary_data = summary_data
        )

    # assigning a new class to the list
    class(clusters_results) <- "clusters_results"

    # returning the clusters and posterior probabilities
    return (clusters_results)

}

#' @export

plot.clusters_results <- function (x, clusters_y = 0, clusters_colour = "steelblue", ...) {

    x$summary_data |>
        dplyr::summarise(eeg_mean = mean(.data$eeg_mean), .by = .data$time) |>
        ggplot2::ggplot(
            ggplot2::aes(x = .data$time, y = .data$eeg_mean)
            ) +
        ggplot2::geom_hline(yintercept = 0.0, linetype = 2) +
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
            lineend = "round",
            linewidth = 2
            ) +
        ggplot2::geom_line() +
        ggplot2::theme_bw() +
        ggplot2::labs(x = "Time (s)", y = "M/EEEG signal")

}
