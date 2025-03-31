#' Time-resolved Bayesian testing
#'
#' Time-resolved testing using Bayesian generalised additive multilevel regression models.
#'
#' @param data Dataframe, with ERP or decoding performance data (in long format).
#' @param participant_id Numeric/Character/Factor, column in data specifying the participant ID.
#' @param meeg_id Numeric, column in data specifying the block ID.
#' @param time_id Numeric, column in data specifying the trial ID.
#' @param predictor_id Numeric, column in data specifying the trial ID.
#' @param kvalue Character, the kernel computation method, either "difference", "GLM", or "GLMM".
#' @param bs Numeric, maximum number of feature levels in the GLMM (keep it low to keep GLMM fitting feasible).
#' @param multilevel Logical, indicating whether the last block was repeated.
#' @param warmup Logical, indicating whether the last block was repeated.
#' @param iter Logical, indicating whether the last block was repeated.
#' @param chains Logical, indicating whether the last block was repeated.
#' @param cores Logical, indicating whether the last block was repeated.
#' @param post_prob_ratio_threshold Logical, indicating whether the last block was repeated.
#' @param chance_level Logical, indicating whether the last block was repeated.
#' @param sesoi Logical, indicating whether the last block was repeated.
#'
#' @return The original data augmented with the kernel.
#'
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' # importing some simulated EEG data
#' data(eeg_data)
#' head(eeg_data)
#'
#' # computing the kernel per participant
#' testing_through_time(data = eeg_data)
#' }
#'
#' @author Ladislas Nalborczyk \email{ladislas.nalborczyk@@gmail.com}.
#'
#' @export

testing_through_time <- function (
        data,
        participant_id = "participant", meeg_id = "eeg", time_id = "time", predictor_id = "condition",
        kvalue = 20, bs = "cr", multilevel = TRUE,
        warmup = 1000, iter = 5000, chains = 4, cores = 4,
        post_prob_ratio_threshold = 20,
        chance_level = 0, sesoi = 0
        ) {

    # ensuring that the method is one of the above
    # method <- match.arg(method)

    # some tests for variable types
    stopifnot("data must be a dataframe..." = is.data.frame(data) )
    stopifnot("kvalue must be a numeric..." = is.numeric(kvalue) )
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

    if (multilevel) {

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
                eeg_sd = sd(.data$eeg),
                .by = c(.data$participant, .data$condition, .data$time)
                )

        # fitting the GAM
        brms_gam <- brms::brm(
            eeg_mean ~ brms::s(time, bs = bs, k = kvalue),
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
            dplyr::mutate(eeg = .data[[meeg_id]]) |>
            dplyr::mutate(time = .data[[time_id]])

        # fitting the GAMM
        brms_gam <- brms::brm(
            eeg ~ brms::s(time, bs = bs, k = kvalue),
            data = summary_data,
            family = stats::gaussian(),
            warmup = warmup,
            iter = iter,
            chains = chains,
            cores = cores
            )

    }

    # computing the posterior probability
    prob_y_above_0 <- data.frame(time = unique(brms_gam$data$time) ) %>%
        # retrieving the posterior samples
        add_epred_draws(object = brms_gam) %>%
        # converting to dataframe
        data.frame() %>%
        # computing mean posterior probability at the group level
        group_by(time) %>%
        summarise(m = mean(.epred > (0 + chance_level + sesoi) ) ) %>%
        mutate(prob_ratio = m / (1 - m) ) %>%
        # ensuring there is no 0 or +Inf values
        mutate(prob_ratio = ifelse(is.infinite(prob_ratio), ndraws(meg_decoding_gam), prob_ratio) ) %>%
        mutate(prob_ratio = ifelse(prob_ratio == 0, 1 / ndraws(meg_decoding_gam), prob_ratio) ) %>%
        ungroup()

    # find the clusters
    clusters <- find_clusters(data = prob_y_above_0, threshold = post_prob_ratio_threshold)

    # returning the clusters
    return (clusters)

}
