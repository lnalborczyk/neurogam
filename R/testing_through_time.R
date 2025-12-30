#' Time-resolved testing based on BGAMMs
#'
#' Fits time-resolved Bayesian generalised additive (multilevel) models (BGAMMs)
#' using \pkg{brms}, and computes posterior odds for an effect at each time
#' point. The effect can be either i) a deviation of the outcome from a
#' reference value (e.g., zero or a chance level), ii) a difference between two
#' groups/conditions (varying within or between participants), or iii) the amplitude
#' of a continuous predictor varying either within (e.g., speech formants) or
#' between participants (e.g., age).
#'
#' @param data A data frame in long format containing time-resolved data.
#' @param participant_id Character; name of the column in \code{data}
#' specifying participant IDs.
#' @param outcome_id Character; name of the column in \code{data} containing
#' the outcome values (e.g., M/EEG amplitude, decoding accuracy).
#' @param outcome_sd Character; name of the column in \code{data} containing
#' the outcome SD, when \code{outcome_id} has already been summarised (default
#' value is NULL).
#' @param time_id Character; name of the column in \code{data}
#' containing time information (e.g., in seconds or samples).
#' @param predictor_id Character; name of the column in \code{data}
#'   containing either:
#'   \itemize{
#'     \item A \emph{binary} categorical predictor (e.g., group or condition),
#'       in which case the function tests, at each time point, whether the
#'       difference between the two levels differs from
#'       \code{chance_level};
#'     \item A \emph{continuous} numeric predictor, in which case the function
#'       tests, at each time point, whether the difference between the average
#'       value of the predictor +1 SD and the average value -1 SD differs from
#'       \code{chance_level}
#'       (typically with \code{chance_level = 0}).
#'     \item If \code{predictor_id = NA}, the function tests whether the outcome differs
#'       from \code{chance_level} over time (useful for decoding accuracies,
#'       for instance).
#'   }
#' @param trials_id Character; name of the column in \code{data}
#' containing the number of trials when using \code{family = binomial()}
#' and summary data. If NULL (default), the function internally summarise binary
#' data into "successes" and total number of "trials".
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
#' @param include_ar_term Logical; if \code{TRUE}, adds an AR(1) autocorrelation
#'   structure within participant via
#'   \code{autocor = brms::ar(time = "time", gr = "participant", p = 1, cov = FALSE)}.
#' @param varying_smooth Logical; should we include a varying smooth. Default is
#' \code{TRUE}. If \code{FALSE}, we only include a varying intercept and slope.
#' @param participant_clusters Logical; should we return clusters at the participant-level.
#' @param warmup Numeric; number of warm-up iterations per chain.
#' @param iter Numeric; total number of iterations per chain (including warmup).
#' @param chains Numeric; number of MCMCs.
#' @param cores Numeric; number of parallel cores to use.
#' @param backend Character; package to use as the backend for fitting the
#'   \code{Stan} model. One of \code{"cmdstanr"} (default) or \code{"rstan"}.
#' @param stan_control List; parameters to control the MCMC behaviour, using
#'   default parameters when NULL. See \code{?brm} for more details.
#' @param n_post_samples Numeric; number of posterior draws used to compute
#'   posterior probabilities. If \code{NULL} (default), all available draws
#'   from the fitted model are used.
#' @param threshold Numeric; threshold on the posterior odds used to define
#'   contiguous temporal clusters. Values greater than 1 favour the hypothesis
#'   that the effect exceeds \code{chance_level}.
#' @param threshold_type Character scalar controlling which clusters are
#'   detected. Must be one of \code{"above"}, \code{"below"}, or \code{"both"}
#'   (default). When \code{"above"}, clusters are formed where
#'   \code{value >= threshold}. When \code{"below"}, clusters are formed where
#'   \code{value <= 1/threshold}. When \code{"both"}, both types are detected
#'   and the returned data include a \code{sign} column.
#' @param chance_level Numeric; null value for the outcome (e.g., 0.5 for
#'   decoding accuracy).
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
#'     effect (or condition difference) exceeds \code{chance_level};
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
#' summary(results)
#'
#' # plot the model predictions and identified clusters
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
        participant_id = "participant", outcome_id = "eeg", outcome_sd = NULL,
        time_id = "time", predictor_id = "condition", trials_id = NULL,
        family = gaussian(), kvalue = 20, bs = "tp",
        multilevel = c("summary", "group"),
        include_ar_term = FALSE,
        participant_clusters = FALSE, varying_smooth = TRUE,
        warmup = 1000, iter = 2000, chains = 4, cores = 4,
        backend = c("cmdstanr", "rstan"),
        stan_control = NULL,
        n_post_samples = NULL,
        threshold = 10, threshold_type = c("both", "above", "below"),
        chance_level = NULL, credible_interval = 0.95
        ) {

    # some tests for variable types
    stopifnot("data must be a dataframe..." = is.data.frame(data) )
    stopifnot("kvalue must be a numeric..." = is.numeric(kvalue) )
    stopifnot("warmup must be a numeric..." = is.numeric(warmup) )
    stopifnot("iter must be a numeric..." = is.numeric(iter) )
    stopifnot("chains must be a numeric..." = is.numeric(chains) )
    stopifnot("cores must be a numeric..." = is.numeric(cores) )
    stopifnot("threshold must be a numeric..." = is.numeric(threshold) )
    stopifnot("bs must be a character..." = is.character(bs) )

    # multilevel should be one of above
    multilevel <- match.arg(multilevel)

    # backend should be one of above
    backend <- match.arg(backend)

    # threshold_type should be one of above
    threshold_type <- match.arg(threshold_type)

    # some more tests
    if (!backend %in% c("cmdstanr", "rstan") ) {

        stop ("`backend` must be either 'cmdstanr' or 'rstan'.", call. = FALSE)

    }

    if (warmup >= iter) {

        stop ("`iter` must be strictly larger than `warmup`.", call. = FALSE)

    }

    if (!is.null(n_post_samples) ) {

        if (!is.numeric(n_post_samples) || length(n_post_samples) != 1 || n_post_samples <= 0) {

            stop ("`n_post_samples` must be NULL or a positive numeric scalar.", call. = FALSE)

        }

    }

    if (include_ar_term && multilevel == "summary") {

        stop (
            "`include_ar_term` not yet implemented for multilevel = 'summary'... try multilevel = 'group'.",
            call. = FALSE
            )

    }

    # restrict supported response distributions
    fam_name <- tryCatch ({
        if (is.list(family) && !is.null(family$family) ) {
                as.character(family$family)
            } else {
                stop ("not a family object")
            }
        },
        error = function (e) {
            stop ("`family` must be a valid family object (see ?brm).", call. = FALSE)
            }
        )

    allowed_families <- c("gaussian", "binomial")

    if (!fam_name %in% allowed_families) {

        stop (
            "Unsupported `family`: '", fam_name, "'. ",
            "Currently supported families are: gaussian() and binomial().",
            call. = FALSE
            )

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

    # tests for binomial response
    is_binom <- identical(fam_name, "binomial")

    if (is_binom && is.null(trials_id) ) {

        cat("NB: When trials_id = NULL (default), neurogam assumes a binary outcome and summarise it internally in counts.\n")

        # retrieve outcome (should be 0/1)
        y <- data[[outcome_id]]

        # accept logical, integer, numeric, but enforce values are 0/1 (ignoring NA)
        ok <- all(is.na(y) | y %in% c(0, 1, FALSE, TRUE) )

        if (!ok) {

            stop (
                "For binomial() models, `", outcome_id, "` must contain trial-level 0/1 (or TRUE/FALSE).",
                call. = FALSE
                )

        }

    }

    if (is.null(chance_level) ) {

        if (is_binom) {

            # define chance_level to 0.5 by default for binomial responses
            chance_level <- 0.5

            # warning the user
            cat("Setting chance_level = 0.5 by default when family = binomial().\n")

        } else {

            # define chance_level to 0 by default for gaussian responses
            chance_level <- 0

            # warning the user
            cat("Setting chance_level = 0 by default when family = gaussian().\n")


        }

    }

    # checking predictor type
    predictor_type <- "none"

    if (!is.na(predictor_id) ) {

        pred_vec <- data[[predictor_id]]

        if (is.numeric(pred_vec) ) {

            predictor_type <- "continuous"

        } else {

            predictor_type <- "categorical"

            # optional: enforce 2 levels for categorical case
            if (length(unique(pred_vec) ) != 2) {

                stop ("For categorical `predictor_id`, there must be exactly 2 levels.", call. = FALSE)

            }

        }

    }

    if (multilevel == "summary") {

        # summarise the data appropriately
        ms <- make_summary_data(
            data = data,
            participant_id = participant_id,
            outcome_id = outcome_id,
            outcome_sd = outcome_sd,
            time_id = time_id,
            predictor_id = predictor_id,
            trials_id = trials_id,
            family = family,
            multilevel = multilevel
            )

        summary_data <- ms$data
        within_between <- ms$within_between

        # define the model formula
        formula_obj <- make_bgam_formula(
            family = family,
            multilevel = multilevel,
            predictor_type = predictor_type,
            within_between = if (!is.na(predictor_id) ) within_between$classification else NA,
            kvalue = kvalue,
            bs = bs,
            include_ar_term = include_ar_term,
            varying_smooth = varying_smooth
            )

        # include new predictor in data
        if (include_ar_term) {

            if (!is.na(predictor_id) ) {

                summary_data <- summary_data |>
                    dplyr::mutate(ar_series = interaction(.data$participant, .data$predictor) )

            } else {

                summary_data <- summary_data |>
                    dplyr::mutate(ar_series = .data$participant)

            }

        }

        # testing whether outcome_sd contains NAs
        if (!is_binom && any(is.na(summary_data$outcome_sd) ) ) {

            na_count <- sum(is.na(summary_data$outcome_sd) )

            stop (
                paste0("Internal data summary returned ", na_count, " NAs. If the input data is already summarised, please use the `outcome_sd` argument. Otherwise, make sure to input trial-by-trial data in long format (i.e., one observation/trial per row)."),
                call. = FALSE
                )

        }

        ####################################################
        # fit the model
        ####################################################
        brms_gam <- brms::brm(
            formula = formula_obj,
            data = summary_data,
            family = family,
            warmup = warmup,
            iter = iter,
            chains = chains,
            cores = cores,
            backend = backend,
            control = stan_control,
            stan_model_args = list(stanc_options = list("O1") )
            )

        if (is.null(n_post_samples) ) {

            n_post_samples <- brms::ndraws(brms_gam)

        }

        # computing the posterior odds over time
        if (is.na(predictor_id) ) {

            if (participant_clusters) {

                # newdata grid over time
                newdata_grid <- tidyr::crossing(
                    time = sort(unique(brms_gam$data$time) ),
                    participant = sort(unique(brms_gam$data$participant) )
                    ) |>
                    # add appropriate dummy to satisfy validate_data()
                    add_required_dummy(is_binom = is_binom)

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
                    # add appropriate dummy to satisfy validate_data()
                    add_required_dummy(is_binom = is_binom) |>
                    # NA participant to satisfy validate_data()
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

            if (participant_clusters) {

                if (within_between$classification == "within-subject") {

                    # newdata grid over time and predictor
                    newdata_grid <- tidyr::crossing(
                        time = sort(unique(brms_gam$data$time) ),
                        predictor = levels(brms_gam$data$predictor),
                        participant = sort(unique(brms_gam$data$participant) )
                        ) |>
                        # add appropriate dummy to satisfy validate_data()
                        add_required_dummy(is_binom = is_binom)

                    # retrieving posterior predictions (draws)
                    post_draws <- tidybayes::add_epred_draws(
                        object = brms_gam,
                        newdata = newdata_grid,
                        ndraws = n_post_samples,
                        re_formula = NULL
                        ) |>
                        data.frame()

                } else if (within_between$classification == "between-subject") {

                    cat(
                        "We can not estimate clusters at the participant level when the predictor varies across participants...\nSwitching to participant_clusters = FALSE\n"
                        )

                    participant_clusters <- FALSE

                }

            } else { # end participant_clusters

                # newdata grid over time and predictor
                newdata_grid <- tidyr::crossing(
                    time = sort(unique(brms_gam$data$time) ),
                    predictor = levels(brms_gam$data$predictor)
                    ) |>
                    # add appropriate dummy to satisfy validate_data()
                    add_required_dummy(is_binom = is_binom) |>
                    # NA participant to satisfy validate_data()
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

            if (participant_clusters) {

                # newdata grid over time and predictor +/-1 SD
                predictor_mean <- mean(brms_gam$data$predictor)
                predictor_sd <- sd(brms_gam$data$predictor)
                newdata_grid <- tidyr::crossing(
                    time = sort(unique(brms_gam$data$time) ),
                    predictor = c(predictor_mean - predictor_sd, predictor_mean + predictor_sd),
                    participant = sort(unique(brms_gam$data$participant) )
                    ) |>
                    # add appropriate dummy to satisfy validate_data()
                    add_required_dummy(is_binom = is_binom)

                # retrieving posterior predictions (draws)
                post_draws <- tidybayes::add_epred_draws(
                    object = brms_gam,
                    newdata = newdata_grid,
                    ndraws = n_post_samples,
                    re_formula = NULL
                    ) |>
                    data.frame()

            } else {

                # newdata grid over time and predictor +/-1 SD
                predictor_mean <- mean(brms_gam$data$predictor)
                predictor_sd <- sd(brms_gam$data$predictor)
                newdata_grid <- tidyr::crossing(
                    time = sort(unique(brms_gam$data$time) ),
                    predictor = c(predictor_mean - predictor_sd, predictor_mean + predictor_sd)
                    ) |>
                    # add appropriate dummy to satisfy validate_data()
                    add_required_dummy(is_binom = is_binom) |>
                    # NA participant to satisfy validate_data()
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

        }

    } else if (multilevel == "group") {

        # summarise the data appropriately
        ms <- make_summary_data(
            data = data,
            participant_id = participant_id,
            outcome_id = outcome_id,
            time_id = time_id,
            predictor_id = predictor_id,
            trials_id = trials_id,
            family = family,
            multilevel = multilevel
            )

        summary_data <- ms$data
        within_between <- ms$within_betwee

        # define the model formula
        formula_obj <- make_bgam_formula(
            family = family,
            multilevel = multilevel,
            predictor_type = predictor_type,
            within_between = if (!is.na(predictor_id) ) within_between$classification else NA,
            kvalue = kvalue,
            bs = bs,
            include_ar_term = include_ar_term,
            varying_smooth = varying_smooth
            )

        # include new predictor in data
        if (include_ar_term) {

            if (!is.na(predictor_id) ) {

                summary_data <- summary_data |>
                    dplyr::mutate(ar_series = interaction(.data$participant, .data$predictor) )

            } else {

                summary_data <- summary_data |>
                    dplyr::mutate(ar_series = .data$participant)

            }

        }

        ####################################################
        # fit the model
        ####################################################
        brms_gam <- brms::brm(
            formula = formula_obj,
            data = summary_data,
            family = family,
            warmup = warmup,
            iter = iter,
            chains = chains,
            cores = cores,
            backend = backend,
            control = stan_control
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

            # newdata grid over time and predictor +/-1 SD
            predictor_mean <- mean(brms_gam$data$predictor)
            predictor_sd <- sd(brms_gam$data$predictor)
            newdata_grid <- tidyr::crossing(
                time = sort(unique(brms_gam$data$time) ),
                predictor = c(predictor_mean - predictor_sd, predictor_mean + predictor_sd)
                )

            # retrieving posterior predictions (draws)
            post_draws <- tidybayes::add_epred_draws(
                object = brms_gam,
                newdata = newdata_grid,
                ndraws = n_post_samples,
                re_formula = NA
                ) |>
                data.frame()

        }

    }

    # compute the posterior odds
    if (is.na(predictor_id) ) {

        prob_y_above <- compute_one_sample_prob(
            post_draws = post_draws,
            participant_clusters = participant_clusters,
            threshold = chance_level,
            n_post_samples = n_post_samples,
            credible_interval = credible_interval
            )

    } else {

        prob_y_above <- compute_two_sample_prob(
            post_draws = post_draws,
            participant_clusters = participant_clusters,
            # when comparing two groups, null value should be 0
            threshold = 0,
            n_post_samples = n_post_samples,
            credible_interval = credible_interval,
            predictor_type = predictor_type
            )

    }

    if (participant_clusters) {

        # find the clusters
        clusters <- find_clusters(
            data = prob_y_above |> dplyr::select(.data$time, .data$participant, value = .data$prob_ratio),
            group = "participant",
            threshold = threshold,
            threshold_type = threshold_type
            )

    } else {

        # find the clusters
        clusters <- find_clusters(
            data = prob_y_above |> dplyr::select(.data$time, value = .data$prob_ratio),
            group = NULL,
            threshold = threshold,
            threshold_type
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

#' @export
plot.clusters_results <- function (
        x, null_value = 0,
        clusters_y = -Inf, clusters_colour = "black", lineend = "butt",
        theme = ggplot2::theme_bw(),
        ...
        ) {

    if (!ggplot2::is.theme(theme) ) {

        stop ("Argument 'theme' should be a 'theme' object.")

    }

    # retrieve the empirical data
    emp_data <- x$model$data

    # group-level or participant-level clusters?
    group_level <- ifelse(
        test = "participant" %in% names(x$clusters),
        yes = FALSE, no = TRUE
        )

    # reconstruct a raw data time course when possible
    if (!is.null(x$multilevel) && "outcome_mean" %in% names(emp_data) ) {

        if ("predictor" %in% names(emp_data) ) {

            if (is.numeric(emp_data$predictor) ) {

                stop ("plot.cluster_results() is not implemented for continuous predictors yet.")

            } else {

                cond1 <- levels(emp_data$predictor)[1]
                cond2 <- levels(emp_data$predictor)[2]

            }

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

    } else if (!is.null(x$multilevel) && "success" %in% names(emp_data) ) {

        if ("predictor" %in% names(emp_data) ) {

            cond1 <- levels(emp_data$predictor)[1]
            cond2 <- levels(emp_data$predictor)[2]

            if (group_level) {

                reshaped_data <- emp_data |>
                    dplyr::mutate(emp_prob = .data$success / .data$trials) |>
                    dplyr::summarise(
                        emp_prob = mean(.data$emp_prob),
                        .by = c(.data$time, .data$predictor)
                        ) |>
                    tidyr::pivot_wider(
                        names_from  = .data$predictor,
                        values_from = .data$emp_prob
                        ) |>
                    dplyr::mutate(outcome_mean = .data[[cond2]] - .data[[cond1]])

            } else {

                reshaped_data <- emp_data |>
                    dplyr::mutate(emp_prob = .data$success / .data$trials) |>
                    dplyr::summarise(
                        emp_prob = mean(.data$emp_prob),
                        .by = c(.data$time, .data$predictor, .data$participant)
                        ) |>
                    tidyr::pivot_wider(
                        names_from  = .data$predictor,
                        values_from = .data$emp_prob
                        ) |>
                    dplyr::mutate(outcome_mean = .data[[cond2]] - .data[[cond1]])

            }

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
        ggplot2::geom_hline(yintercept = null_value, linetype = 2) +
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
                x = .data$onset,
                xend = .data$offset,
                y = clusters_y,
                yend = clusters_y
                ),
            colour = clusters_colour,
            inherit.aes = FALSE,
            lineend = lineend,
            linewidth = 5
            ) +
        theme +
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
                onset = round(.data$onset,  digits),
                offset = round(.data$offset, digits),
                duration = round(.data$offset - .data$onset, digits)
                ) |>
            dplyr::select(
                .data$participant, .data$sign, .data$id,
                .data$onset, .data$offset, .data$duration
                )

    } else {

        clust_tbl <- x$clusters |>
            dplyr::mutate(
                onset = round(.data$onset,  digits),
                offset = round(.data$offset, digits),
                duration = round(.data$offset - .data$onset, digits)
                ) |>
            dplyr::select(
                .data$sign, .data$id, .data$onset,
                .data$offset, .data$duration
                )

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
        dplyr::mutate(duration = .data$offset - .data$onset)

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
                onset = round(.data$onset,  digits),
                offset = round(.data$offset, digits),
                duration = round(.data$duration, digits)
                ) |>
            dplyr::select(
                .data$participant, .data$sign, .data$id,
                .data$onset, .data$offset, .data$duration
                )

    } else {

        cl_print <- cl |>
            dplyr::mutate(
                onset = round(.data$onset,  digits),
                offset = round(.data$offset, digits),
                duration = round(.data$duration, digits)
                ) |>
            dplyr::select(
                .data$sign, .data$id,
                .data$onset, .data$offset, .data$duration
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
#' @param xlab Character; Label for the x-axis (usually time with some appropriate unit).
#' @param theme A \code{\link[ggplot2:theme]{theme}} object
#'   modifying the appearance of the plots.
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
        xlab = "Time (s)",
        theme = ggplot2::theme_bw(),
        ...
        ) {

    if (!ggplot2::is.theme(theme) ) {

        stop ("Argument 'theme' should be a 'theme' object.")

    }

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

    # binomial?
    is_binom <- ifelse(test = "success" %in% names(data_fit), yes = TRUE, no = FALSE)

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
                dplyr::summarise(
                    dplyr::across(dplyr::where(is.numeric), mean, na.rm = TRUE),
                    .by = .data$time
                    ) |>
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
                theme +
                ggplot2::labs(x = xlab)

        } else {

            # grid for group-level prediction
            if (is_binom) {

                # observed y on that same grid
                newdata <- tidyr::crossing(predictor = data_fit$predictor, time = data_fit$time) |>
                    add_required_dummy(is_binom = is_binom)

                y_obs <- data_fit |>
                    dplyr::mutate(emp_prob = .data$success / .data$trials) |>
                    dplyr::summarise(
                        emp_prob = mean(.data$emp_prob),
                        .by = c(.data$predictor, .data$time)
                        ) |>
                    dplyr::arrange(.data$predictor, .data$time) |>
                    dplyr::pull(.data$emp_prob)

                # predict probs of success at the group level
                yrep <- brms::posterior_epred(
                    object = fit,
                    newdata = newdata,
                    re_formula = NA,
                    ndraws = ndraws
                    )

            } else {

                newdata <- data_fit |>
                    dplyr::summarise(
                        dplyr::across(dplyr::where(is.numeric), mean, na.rm = TRUE),
                        .by = c(.data$predictor, .data$time)
                        ) |>
                    dplyr::mutate(participant = NA) |>
                    dplyr::arrange(.data$predictor, .data$time)

                # observed y on that same grid
                y_obs <- data_fit |>
                    dplyr::group_by(.data$time, .data$predictor) |>
                    dplyr::summarise(
                        y = mean(.data$outcome_mean, na.rm = TRUE),
                        .groups = "drop"
                        ) |>
                    dplyr::arrange(.data$predictor, .data$time) |>
                    dplyr::pull(.data$y)

                # simulate from posterior at the group level
                yrep <- brms::posterior_predict(
                    object = fit,
                    newdata = newdata,
                    re_formula = NA,
                    ndraws = ndraws
                    )

            }

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
                theme +
                ggplot2::labs(x = xlab)

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
            alpha = 0.5
            ) +
            theme +
            ggplot2::labs(x = xlab)

    }

    # returning the plot
    print(ppc_plot)
    invisible(ppc_plot)

}
