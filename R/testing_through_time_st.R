#' Internal high-level fitting and inference for spatiotemporal EEG
#'
#' Called from `testing_through_time()` when `space_id` is supplied.
#'
#' @param data Input long-format EEG data.
#' @param previous_model Optional compatible fitted `brmsfit`.
#' @param participant_id,outcome_id,outcome_sem,time_id,sensor_id,space_id,
#'   predictor_id Source column names.
#' @param family Gaussian likelihood.
#' @param kvalue,bs,tensor_d,t2_full Tensor product basis specification.
#' @param multilevel Only `"summary"` is supported.
#' @param use_se Whether to model supplied/estimated participant SEM.
#' @param varying_smooth Must be `FALSE` (shared tensor smooth only).
#' @param include_ar_term Whether to model AR dependence (unsupported).
#' @param include_by_smooth Whether to include condition-varying smooths.
#' @param participant_clusters Unsupported for spatiotemporal models.
#' @param warmup,iter,chains,cores,threads,backend,stan_control,file Fit options.
#' @param n_post_samples Number of posterior draws, or `NULL` for all.
#' @param threshold,threshold_type,credible_interval Clustering parameters.
#' @param chance_level Null value for one-sample inference.
#' @param rope Symmetric fixed ROPE half-width.
#' @param contrast_levels Optional A - B ordering for two conditions.
#' @param sensor_adjacency Optional known electrode graph.
#' @param adjacency_k,adjacency_distance_factor Automatic graph settings.
#' @param chunk_size Number of sensor-time locations per prediction chunk.
#' @param re_formula Population-level predictions default to `NA`.
#' @param ... Additional arguments to `brms::brm()` when fitting.
#'
#' @return A `clusters_results_st` object, plus `summary_data`, `clusters`,
#'   `multilevel`, and `contrast_levels`.
#'
#' @keywords internal
testing_through_time_spatiotemporal <- function (
    data,
    previous_model = NULL,
    participant_id = "participant",
    outcome_id = "eeg",
    outcome_sem = NULL,
    time_id = "time",
    sensor_id = "channel",
    space_id = c("xproj", "yproj"),
    predictor_id = NA_character_,
    family = stats::gaussian(),
    kvalue = c(15L, 10L),
    bs = c("ps", "ps"),
    tensor_d = c(2L, 1L),
    t2_full = FALSE,
    multilevel = "summary",
    include_ar_term = FALSE,
    use_se = TRUE,
    varying_smooth = FALSE,
    include_by_smooth = TRUE,
    participant_clusters = FALSE,
    warmup = 1000,
    iter = 2000,
    chains = 4,
    cores = 4,
    threads = NULL,
    backend = "cmdstanr",
    stan_control = NULL,
    file = NULL,
    n_post_samples = NULL,
    threshold = 10,
    threshold_type = "both",
    odds_mode = c(
        "directional_rope",
        "legacy_one_sided",
        "legacy_directional",
        "spatiotemporal"
        ),
    chance_level = 0,
    credible_interval = 0.95,
    rope = 0,
    contrast_levels = NULL,
    sensor_adjacency = NULL,
    adjacency_k = 4L,
    adjacency_distance_factor = 1.8,
    chunk_size = 256L,
    re_formula = NA,
    ...
    ) {

    if (!is.data.frame(data) || !nrow(data) ) {

        stop ("`data` must be a nonempty data frame.", call. = FALSE)

    }

    if (length(time_id) != 1L || length(space_id) != 2L ||
        !is.character(time_id) || !is.character(space_id) ) {

        stop (
            "Spatiotemporal fitting requires one time and two spatial coordinates.",
            call. = FALSE
            )

    }

    if (length(predictor_id) != 1L) {
        stop("`predictor_id` must name one column or be `NA`.", call. = FALSE)
    }
    if (isTRUE(participant_clusters)) {
        stop("`participant_clusters = TRUE` is not supported for 3D models.",
             call. = FALSE)
    }
    if (!identical(multilevel, "summary")) {
        stop("Only multilevel = 'summary' is supported for 3D models.",
             call. = FALSE)
    }
    if (isTRUE(varying_smooth)) {
        stop("Participant-specific 3D smooths are not implemented; set varying_smooth = FALSE.",
             call. = FALSE)
    }
    if (isTRUE(include_ar_term)) {
        stop("AR residual structure is not yet supported for 3D models.",
             call. = FALSE)
    }
    threshold_type <- match.arg(threshold_type, c("both", "above", "below"))
    backend <- match.arg(backend, c("rstan", "cmdstanr"))
    if (!is.null(n_post_samples) &&
        (!is.numeric(n_post_samples) || length(n_post_samples) != 1L ||
         !is.finite(n_post_samples) || n_post_samples < 1L ||
         n_post_samples != as.integer(n_post_samples))) {
        stop("`n_post_samples` must be NULL or a positive integer.", call. = FALSE)
    }
    if (!is.numeric(chance_level) || length(chance_level) != 1L ||
        !is.finite(chance_level)) {
        stop("`chance_level` must be a finite numeric scalar.", call. = FALSE)
    }

    summarized <- make_spatiotemporal_summary_data(
        data = data,
        participant_id = participant_id,
        outcome_id = outcome_id,
        outcome_sem = outcome_sem,
        time_id = time_id,
        sensor_id = sensor_id,
        space_id = space_id,
        predictor_id = predictor_id,
        use_se = use_se
    )
    predictor_type <- if ("predictor" %in% names(summarized)) {
        "categorical"
    } else {
        "none"
    }
    within_between <- NA_character_
    if (identical(predictor_type, "categorical")) {
        n_levels <- tapply(
            as.character(summarized$predictor),
            summarized$participant,
            function(x) length(unique(x))
        )
        if (all(n_levels == 2L)) {
            within_between <- "within-subject"
        } else if (all(n_levels == 1L)) {
            within_between <- "between-subject"
        } else {
            stop("Condition factor must be consistently within- or between-subject.",
                 call. = FALSE)
        }
    }
    formula_obj <- make_bgam_formula(
        family = family,
        multilevel = "summary",
        predictor_type = predictor_type,
        within_between = within_between,
        time_id = "time",
        space_id = c("xproj", "yproj"),
        tensor_d = tensor_d,
        kvalue = kvalue,
        bs = bs,
        t2_full = t2_full,
        include_ar_term = FALSE,
        varying_smooth = FALSE,
        include_by_smooth = include_by_smooth,
        use_se = use_se
    )

    if (!is.null(previous_model)) {
        if (!inherits(previous_model, "brmsfit")) {
            stop("`previous_model` must be a fitted `brmsfit` object.",
                 call. = FALSE)
        }
        needed <- c("outcome_mean", "participant", "time", "xproj", "yproj")
        if (!all(needed %in% names(previous_model$data))) {
            stop(
                "`previous_model` is not compatible with the Stage 3 summarized-data interface. Use `analyse_bgam_spatiotemporal()` for externally fitted models.",
                call. = FALSE
            )
        }
        if (identical(predictor_type, "categorical") !=
            ("predictor" %in% names(previous_model$data))) {
            stop("`previous_model` has incompatible predictor structure.",
                 call. = FALSE)
        }
    }

    if (is.null(previous_model)) {
        message(
            "Fitting spatial-temporal BGAMM with formula: ",
            paste(deparse(formula_obj$formula), collapse = " ")
        )
    }

    fit_result <- fit_or_reuse_bgam(
        formula = formula_obj,
        data = summarized,
        family = family,
        previous_model = previous_model,
        n_post_samples = n_post_samples,
        warmup = warmup,
        iter = iter,
        chains = chains,
        cores = cores,
        threads = threads,
        backend = backend,
        stan_control = stan_control,
        file = file,
        extra_args = list(...),
        resolve_draws = FALSE,
        reuse_message = paste0(
            "Using `previous_model`: skipping spatial-temporal model fitting."
        )
    )
    fitted <- fit_result$model

    grids <- make_spatiotemporal_prediction_grid(
        summarized,
        contrast_levels = contrast_levels
    )
    # Dummy response / SEM columns are useful for brms versions that validate
    # measurement-error terms in newdata even for posterior_epred().
    grids$a$outcome_mean <- 0
    grids$a$outcome_sem <- 0
    if (!is.null(grids$b)) {
        grids$b$outcome_mean <- 0
        grids$b$outcome_sem <- 0
    }

    # draws_available <- brms::ndraws(fitted)
    # if (is.null(n_post_samples)) {
    #     draw_ids <- NULL
    # } else {
    #     if (n_post_samples > draws_available) {
    #         stop("Requested more posterior draws than the model contains.",
    #              call. = FALSE)
    #     }
    #     draw_ids <- as.integer(round(seq(
    #         1, draws_available, length.out = n_post_samples
    #     )))
    # }

    draw_ids <- select_prediction_draws(
        n_available = brms::ndraws(fitted),
        n_requested = n_post_samples
        )

    analysed <- analyse_bgam_spatiotemporal(
        model = fitted,
        newdata = grids$a,
        newdata_b = grids$b,
        sensor_id = "sensor",
        time_id = "time",
        spatial_id = c("xproj", "yproj"),
        sensor_adjacency = sensor_adjacency,
        adjacency_k = adjacency_k,
        adjacency_distance_factor = adjacency_distance_factor,
        null_value = if (is.null(grids$b)) chance_level else 0,
        rope = rope,
        threshold = threshold,
        threshold_type = threshold_type,
        odds_mode = odds_mode,
        credible_interval = credible_interval,
        draw_ids = draw_ids,
        chunk_size = chunk_size,
        re_formula = re_formula
        )
    # Restore the external sensor and spatial names expected by the plotting
    # functions. Other standard output columns are unchanged.
    rename_output <- function(tab) {
        rename <- c(sensor = sensor_id, xproj = space_id[1L], yproj = space_id[2L])
        for (old in names(rename)) {
            if (old %in% names(tab)) {
                names(tab)[names(tab) == old] <- rename[[old]]
            }
        }
        if ("time" %in% names(tab)) {
            names(tab)[names(tab) == "time"] <- time_id
        }
        tab
    }
    analysed$predictions <- rename_output(analysed$predictions)
    analysed$cluster_membership <- rename_output(analysed$cluster_membership)
    analysed$sensor_positions <- rename_output(analysed$sensor_positions)
    analysed$clusters <- analysed$cluster_summary
    analysed$summary_data <- summarized
    analysed$multilevel <- multilevel
    analysed$contrast_levels <- if (identical(grids$contrast, "one_sample")) {
        NULL
    } else {
        grids$contrast
    }

    analysed$input_columns <- list(
        sensor_id = sensor_id,
        time_id = time_id,
        space_id = space_id,
        participant_id = participant_id,
        outcome_id = outcome_id
        )

    analysed <- set_neurogam_result_metadata(
        analysed,
        make_neurogam_metadata(
            analysis_type = "spatiotemporal",
            time_id = time_id,
            space_id = space_id,
            sensor_id = sensor_id,
            participant_id = participant_id,
            outcome_id = outcome_id,
            predictor_id = if (is.na(predictor_id)) NULL else predictor_id,
            multilevel = multilevel,
            threshold = threshold,
            threshold_type = threshold_type,
            null_value = if (is.null(grids$b)) chance_level else 0,
            rope = rope,
            credible_interval = credible_interval,
            odds_mode = "legacy_one_sided",
            contrast = analysed$contrast,
            contrast_levels = analysed$contrast_levels,
            adjacency_source = analysed$adjacency_source,
            participant_clusters = FALSE,
            n_post_samples = if (is.null(draw_ids)) NULL else length(draw_ids)
        )
    )

    return(analysed)

}
