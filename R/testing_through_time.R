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
#' @param previous_model Optional. A previously fitted \code{brmsfit} object
#'   obtained from \code{testing_through_time()}. When provided, the model is
#'   not refitted; instead, posterior predictions and inference are recomputed
#'   using the supplied model. This is useful for exploring the effect of
#'   different \code{threshold} and \code{threshold_type} values without
#'   re-running model fitting.
#'
#'   The supplied model must be compatible with the current function call
#'   (i.e., same data structure, formula, family, and predictors). If
#'   \code{previous_model} is not \code{NULL}, arguments related to model
#'   estimation (e.g., \code{warmup}, \code{iter}, \code{chains},
#'   \code{cores}, \code{backend}, \code{stan_control}) are ignored.
#' @param participant_id Character; name of the column in \code{data}
#'   specifying participant IDs.
#' @param outcome_id Character; name of the column in \code{data} containing
#'   the outcome values (e.g., M/EEG amplitude, decoding accuracy).
#' @param outcome_sd Character; name of the column in \code{data} containing
#'   the outcome SD, when \code{outcome_id} has already been summarised (default
#'   value is NULL).
#' @param time_id Character; name of the column(s) in \code{data}
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
#'   containing the number of trials when using \code{family = binomial()}
#'   and summary data. If NULL (default), the function internally summarise binary
#'   data into "successes" and total number of "trials".
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
#' @param use_se Logical; whether to include known or internally computed
#'   measurement error via \code{y | se(outcome_sd)} in the model formula.
#' @param t2_full Logical; If TRUE, then there is a separate penalty for each
#'   combination of null space column and range space, see \code{\link[mgcv]{t2}}.
#'   Only use when fitting 2D temporal models (i.e., when \code{time_id} contains
#'   two temporal variables).
#' @param varying_smooth Logical; should we include a varying smooth. Default is
#' \code{TRUE}. If \code{FALSE}, we only include a varying intercept and slope.
#' @param participant_clusters Logical; should we return clusters at the participant-level.
#' @param warmup Numeric; number of warm-up iterations per chain.
#' @param iter Numeric; total number of iterations per chain (including warmup).
#' @param chains Numeric; number of MCMCs.
#' @param cores Numeric; number of parallel cores to use.
#' @param threads Numeric; number of threads to use in within-chain parallelisation.
#'   See \code{\link[brms]{brm}} documentation for more information.
#' @param backend Character; package to use as the backend for fitting the
#'   \code{Stan} model. One of \code{"cmdstanr"} (default) or \code{"rstan"}.
#' @param stan_control List; parameters to control the MCMC behaviour, using
#'   default parameters when NULL. See \code{?brm} for more details.
#' @param file Either NULL or a character string. In the latter case, the
#'   fitted \code{brms} model object is saved via saveRDS in a file named after
#'   the string supplied in file. The \code{.rds} extension is added
#'   automatically. If the file already exists, \code{brm} will load and return
#'   the saved model object instead of refitting the model.
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
#'     \item \code{summary_data}: data used to fit the \pkg{brms} model
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
#'
#' # posterior predictive check
#' ppc(results)
#' }
#'
#' @author Ladislas Nalborczyk \email{ladislas.nalborczyk@@cnrs.fr}.
#'
#' @seealso \code{\link[brms]{brm}}
#'
#' @export
testing_through_time <- function (
        data,
        previous_model = NULL,
        participant_id = "participant", outcome_id = "eeg", outcome_sd = NULL,
        time_id = "time", predictor_id = "condition", trials_id = NULL,
        family = gaussian(), kvalue = 20, bs = "tp",
        multilevel = c("summary", "group"),
        include_ar_term = FALSE,
        use_se = TRUE, t2_full = FALSE,
        participant_clusters = FALSE, varying_smooth = TRUE,
        warmup = 1000, iter = 2000, chains = 4, cores = 4, threads = NULL,
        backend = c("cmdstanr", "rstan"),
        stan_control = NULL,
        file = NULL,
        n_post_samples = NULL,
        threshold = 10, threshold_type = c("both", "above", "below"),
        chance_level = NULL, credible_interval = 0.95
        ) {

    # some tests for variable types
    stopifnot("data must be a dataframe..." = is.data.frame(data) )
    stopifnot("kvalue must be a numeric..." = is.numeric(kvalue) )
    stopifnot("bs must be a character..." = is.character(bs) )
    stopifnot("warmup must be a numeric..." = is.numeric(warmup) )
    stopifnot("iter must be a numeric..." = is.numeric(iter) )
    stopifnot("chains must be a numeric..." = is.numeric(chains) )
    stopifnot("cores must be a numeric..." = is.numeric(cores) )
    stopifnot("threshold must be a numeric..." = is.numeric(threshold) )
    stopifnot("credible_interval must be a numeric..." = is.numeric(credible_interval) )

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

    # time_id can be 1D or 2D (character vector of length 1 or 2)
    if (!is.character(time_id) ) {

        stop (
            "`time_id` must be a character vector (e.g., \"time\" or c(\"train_time\", \"test_time\")). ",
            "Passing unquoted names like time_id = c(train_time, test_time) will fail unless ",
            "`train_time` and `test_time` objects exist in the calling environment."
            )

    }

    if (!(length(time_id) %in% c(1L, 2L) ) ) {

        stop ("`time_id` must have length 1 or 2.")

    }

    if (length(time_id) == 2L && identical(time_id[[1]], time_id[[2]]) ) {

        stop ("When `time_id` has length 2, the two names must be different.")

    }

    # check columns exist
    missing_time_cols <- setdiff(time_id, names(data) )
    if (length(missing_time_cols) > 0) {

        stop (
            "The following `time_id` column(s) are missing from `data`: ",
            paste(missing_time_cols, collapse = ", ")
            )

    }

    # check further required column names
    required_columns <- c(participant_id, outcome_id)

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

    if (is_binom) use_se <- FALSE

    if (use_se && fam_name != "gaussian") {

        stop (
            "`se()` is only supported for gaussian models.",
            call. = FALSE
            )

    }

    if (is.null(chance_level) ) {

        if (is_binom) {

            # define chance_level to 0.5 by default for binomial responses
            chance_level <- 0.5

            # warning the user
            cat("Setting chance_level = 0.5 by default when family = binomial().\n")

        } else {

            # define chance_level to 0 by default for Gaussian responses
            chance_level <- 0

            # warning the user
            cat("Assuming null_value (chance_level) = 0 by default when family = gaussian().\n")


        }

    }

    # allow reusing an already fitted model
    use_previous_model <- !is.null(previous_model)

    if (use_previous_model) {

        if (!inherits(previous_model, "brmsfit") ) {

            stop ("`previous_model` must be a valid `brmsfit` object.", call. = FALSE)

        }

        brms_gam <- previous_model

        if (is.null(n_post_samples) ) {

            n_post_samples <- brms::ndraws(brms_gam)

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

            # enforce 2 levels for categorical case
            if (length(unique(pred_vec) ) != 2) {

                stop ("For categorical `predictor_id`, there must be exactly 2 levels.", call. = FALSE)

            }

        }

    }

    # retrieving temporal variable(s)
    is_2d_time <- length(time_id) == 2L
    time_vars <- time_id

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
            time_id = time_id,
            t2_full = t2_full,
            kvalue = kvalue,
            bs = bs,
            include_ar_term = include_ar_term,
            varying_smooth = varying_smooth,
            use_se = use_se
            )

        # include new predictor in data
        if (include_ar_term) {

            # if there is a predictor and if it varies within participants
            if (!is.na(predictor_id) && within_between$classification == "within-subject") {

                summary_data <- summary_data |>
                    dplyr::mutate(ar_series = interaction(.data$participant, .data$predictor) )

            } else {

                summary_data <- summary_data |>
                    dplyr::mutate(ar_series = .data$participant)

            }

        }

        # testing whether outcome_sd contains NAs
        if (!is_binom && use_se && any(is.na(summary_data$outcome_sd) ) ) {

            na_count <- sum(is.na(summary_data$outcome_sd) )

            stop (
                paste0("Internal data summary returned ", na_count, " NAs. If the input data is already summarised, please use the `outcome_sd` argument. Otherwise, make sure to input trial-by-trial data in long format (i.e., one observation/trial per row)."),
                call. = FALSE
                )

        }

        # displays the model formula
        message (
            "Fitting model with formula: ",
            paste(utils::capture.output(print(formula_obj) ), collapse = " "), "\n"
            )

        #####################################################
        # fit the model (unless previous_model is provided) #
        #####################################################

        if (!use_previous_model) {

            brms_gam <- brms::brm(
                formula = formula_obj,
                data = summary_data,
                family = family,
                warmup = warmup,
                iter = iter,
                chains = chains,
                cores = cores,
                threads = threads,
                backend = backend,
                control = stan_control,
                file = file,
                stan_model_args = list(stanc_options = list("O1") )
                )

            if (is.null(n_post_samples) ) {

                n_post_samples <- brms::ndraws(brms_gam)

            }

        } else {

            message ("Using `previous_model`: skipping model fitting.")

        }

        ###########################################
        # build prediction grid from summary_data #
        ###########################################

        newdata_grid <- make_prediction_grid_from_summary(
            summary_data = summary_data,
            brms_gam = brms_gam,
            predictor_id = predictor_id,
            predictor_type = predictor_type,
            participant_clusters = participant_clusters,
            within_between = within_between,
            is_2d_time = is_2d_time,
            is_binom = is_binom,
            include_ar_term = include_ar_term,
            continuous_mode = "pm1sd"
            )

        ##################################
        # retrieve posterior predictions #
        ##################################

        post_draws <- tidybayes::add_epred_draws(
            object = brms_gam,
            newdata = newdata_grid,
            ndraws = n_post_samples,
            re_formula = if (participant_clusters) NULL else NA,
            incl_autocor = TRUE
            ) |>
            data.frame()

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
        within_between <- ms$within_between

        # define the model formula
        formula_obj <- make_bgam_formula(
            family = family,
            multilevel = multilevel,
            predictor_type = predictor_type,
            within_between = if (!is.na(predictor_id) ) within_between$classification else NA,
            time_id = time_id,
            t2_full = t2_full,
            kvalue = kvalue,
            bs = bs,
            include_ar_term = include_ar_term,
            varying_smooth = varying_smooth
            )

        # include new predictor in data
        if (include_ar_term) {

            # if there is a predictor and if it varies within participants
            if (!is.na(predictor_id) && within_between$classification == "within-subject") {

                summary_data <- summary_data |>
                    dplyr::mutate(ar_series = interaction(.data$participant, .data$predictor) )

            } else {

                summary_data <- summary_data |>
                    dplyr::mutate(ar_series = .data$participant)

            }

        }

        # display the model formula
        message (
            "Fitting model with formula: ",
            paste(utils::capture.output(print(formula_obj) ), collapse = " "), "\n"
            )

        #####################################################
        # fit the model (unless previous_model is provided) #
        #####################################################

        if (!use_previous_model) {

            brms_gam <- brms::brm(
                formula = formula_obj,
                data = summary_data,
                family = family,
                warmup = warmup,
                iter = iter,
                chains = chains,
                cores = cores,
                threads = threads,
                backend = backend,
                control = stan_control,
                file = file,
                stan_model_args = list(stanc_options = list("O1") )
                )

            if (is.null(n_post_samples) ) {

                n_post_samples <- brms::ndraws(brms_gam)

            }

        } else {

            message ("Using `previous_model`: skipping model fitting.")

        }

        # compute the posterior odds over time
        if (is.na(predictor_id) ) {

            if (is_2d_time) {

                if (participant_clusters) {

                    stop ("`participant_clusters = TRUE` is not supported for 2D temporal models yet.", call. = FALSE)

                }

                newdata_grid <- tidyr::crossing(
                    time1 = sort(unique(brms_gam$data$time1) ),
                    time2 = sort(unique(brms_gam$data$time2) )
                    ) |>
                    # add appropriate dummy to satisfy validate_data()
                    add_required_dummy(is_binom = is_binom) |>
                    # NA participant to satisfy validate_data()
                    dplyr::mutate(participant = NA)

            } else {

                # newdata grid over time
                newdata_grid <- tidyr::crossing(time = sort(unique(brms_gam$data$time) ) )

                }

        } else if (predictor_type == "categorical") {

            if (is_2d_time) {

                if (participant_clusters) {

                    stop ("`participant_clusters = TRUE` is not supported for 2D temporal models yet.", call. = FALSE)

                }

                newdata_grid <- tidyr::crossing(
                    time1 = sort(unique(brms_gam$data$time1) ),
                    time2 = sort(unique(brms_gam$data$time2) ),
                    predictor = levels(brms_gam$data$predictor)
                    ) |>
                    # add appropriate dummy to satisfy validate_data()
                    add_required_dummy(is_binom = is_binom) |>
                    # NA participant to satisfy validate_data()
                    dplyr::mutate(participant = NA)

            } else {

                # newdata grid over time and predictor
                newdata_grid <- tidyr::crossing(
                    time = sort(unique(brms_gam$data$time) ),
                    predictor = levels(brms_gam$data$predictor)
                    )

            }

        } else if (predictor_type == "continuous") {

            if (is_2d_time) {

                if (participant_clusters) {

                    stop ("`participant_clusters = TRUE` is not supported for 2D temporal models yet.", call. = FALSE)

                }

                predictor_mean <- mean(brms_gam$data$predictor)
                predictor_sd <- stats::sd(brms_gam$data$predictor)

                newdata_grid <- tidyr::crossing(
                    time1 = sort(unique(brms_gam$data$time1) ),
                    time2 = sort(unique(brms_gam$data$time2) ),
                    predictor = c(predictor_mean - predictor_sd, predictor_mean + predictor_sd)
                    ) |>
                    # add appropriate dummy to satisfy validate_data()
                    add_required_dummy(is_binom = is_binom) |>
                    # NA participant to satisfy validate_data()
                    dplyr::mutate(participant = NA)

            } else {

                # newdata grid over time and predictor +/-1 SD
                predictor_mean <- mean(brms_gam$data$predictor)
                predictor_sd <- sd(brms_gam$data$predictor)
                newdata_grid <- tidyr::crossing(
                    time = sort(unique(brms_gam$data$time) ),
                    predictor = c(predictor_mean - predictor_sd, predictor_mean + predictor_sd)
                    )

            }

        }

        # retrieve posterior predictions (draws)
        post_draws <- tidybayes::add_epred_draws(
            object = brms_gam,
            newdata = newdata_grid,
            ndraws = n_post_samples,
            re_formula = NA,
            incl_autocor = FALSE
            ) |>
            data.frame()

    }

    ###########################################
    # compute the posterior odds and clusters #
    ###########################################

    if (!is_2d_time) {

        if (is.na(predictor_id) ) {

            prob_y_above <- compute_one_sample_prob(
                post_draws = post_draws,
                participant_clusters = participant_clusters,
                null_value = chance_level,
                n_post_samples = n_post_samples,
                credible_interval = credible_interval
                )

        } else {

            prob_y_above <- compute_two_sample_prob(
                post_draws = post_draws,
                # when comparing two groups/conditions, null value should be 0
                null_value = 0,
                participant_clusters = participant_clusters,
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

    } else { # 2D clusters

        if (participant_clusters) {

            stop (
                "`participant_clusters = TRUE` is not supported for 2D temporal models yet. ",
                "Please set `participant_clusters = FALSE`.",
                call. = FALSE
                )

        }

        if (is.na(predictor_id) ) {

            prob_y_above <- compute_one_sample_prob_2d(
                post_draws = post_draws,
                participant_clusters = FALSE,
                null_value = chance_level,
                n_post_samples = n_post_samples,
                credible_interval = credible_interval
                )

        } else {

            prob_y_above <- compute_two_sample_prob_2d(
                post_draws = post_draws,
                # when comparing two groups/conditions, null value should be 0
                null_value = 0,
                participant_clusters = FALSE,
                n_post_samples = n_post_samples,
                credible_interval = credible_interval,
                predictor_type = predictor_type
                )

        }

        time_id <- c("time1", "time2")

        clusters <- find_clusters(
            data = prob_y_above |> dplyr::select(dplyr::all_of(time_id), value = .data$prob_ratio),
            group = NULL,
            threshold = threshold,
            threshold_type = threshold_type,
            time_id = time_id
            )

    }

    # combine the results in a list
    clusters_results <- list(
        clusters = clusters,
        predictions = prob_y_above,
        model = brms_gam,
        summary_data = summary_data,
        multilevel = multilevel
        )

    # assign a new class to the list
    if (is_2d_time) {

        class(clusters_results) <- "clusters_results_2d"

    } else {

        class(clusters_results) <- "clusters_results_1d"

    }

    # return the clusters and posterior probabilities
    return (clusters_results)

}
