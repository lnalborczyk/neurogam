#' Summarise trial-level data for \code{testing_through_time()}
#' @keywords internal
make_summary_data <- function (
        data,
        participant_id = "participant",
        outcome_id = "eeg",
        outcome_sem = NULL,
        time_id = "time",
        predictor_id = NA,
        trials_id = NULL,
        family = gaussian(),
        multilevel = c("summary", "group"),
        na_rm = TRUE,
        # how to aggregate continuous predictors within (participant, time)
        continuous_agg = c("mean", "first")
        ) {

    multilevel <- match.arg(multilevel)
    continuous_agg <- match.arg(continuous_agg)

    # validate time_id
    if (!is.character(time_id) || !(length(time_id) %in% c(1L, 2L) ) ) {

        stop ("`time_id` must be a character vector of length 1 or 2.", call. = FALSE)

    }

    is_2d_time <- length(time_id) == 2L

    if (is_2d_time && identical(time_id[[1]], time_id[[2]]) ) {

        stop ("When `time_id` has length 2, the two names must be different.", call. = FALSE)

    }

    # checking family name
    fam_name <- tryCatch ({

        if (is.list(family) && !is.null(family$family) ) as.character(family$family) else stop ("not a family object")

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

    is_binom <- identical(fam_name, "binomial")

    # column existence checks
    required <- c(participant_id, outcome_id, time_id)

    if (!is.na(predictor_id) ) required <- c(required, predictor_id)

    # if user provides an SEM column name, require it (gaussian only)
    if (!is.null(outcome_sem) ) {

        if (!is.character(outcome_sem) || length(outcome_sem) != 1) {

            stop ("`outcome_sem` must be NULL or a single character column name.", call. = FALSE)

        }

        required <- c(required, outcome_sem)

    }

    if (is_binom && !is.null(trials_id) ) required <- c(required, trials_id)

    missing_cols <- setdiff(required, names(data) )

    if (length(missing_cols) > 0) {

        stop ("Missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)

    }

    # internal columns (UPDATED for 2D time)
    if (!is_2d_time) {

        df <- data |>
            dplyr::mutate(
                participant = .data[[participant_id]],
                time = .data[[time_id[[1]]]]
                )

    } else {

        df <- data |>
            dplyr::mutate(
                participant = .data[[participant_id]],
                time1 = .data[[time_id[[1]]]],
                time2 = .data[[time_id[[2]]]]
                )

    }

    predictor_type <- "none"

    if (!is.na(predictor_id) ) {

        pred_raw <- df[[predictor_id]]

        if (is.numeric(pred_raw) ) {

            predictor_type <- "continuous"
            df <- df |> dplyr::mutate(predictor = as.numeric(.data[[predictor_id]]) )

        } else {

            predictor_type <- "categorical"
            df <- df |> dplyr::mutate(predictor = as.factor(.data[[predictor_id]]) )

        }

    }

    # helper for aggregating continuous predictors within (participant, time, ...)
    agg_continuous <- function (x) {

        if (continuous_agg == "first") {

            # "first" non-NA (or NA if all NA)
            x2 <- x[!is.na(x)]

            agg_result <- if (length(x2) == 0L) NA_real_ else x2[[1]]

        } else {

            agg_result <- mean(x = x, na.rm = na_rm)

        }

        return (agg_result)

    }

    within_between <- NULL

    # choose grouping keys for the outcome
    if (!is_2d_time) {

        time_keys_summary <- c("participant", "time")
        time_keys_group <- c("time")

    } else {

        time_keys_summary <- c("participant", "time1", "time2")
        time_keys_group <- c("time1", "time2")

    }

    grp_outcome <- time_keys_summary
    if (predictor_type == "categorical") grp_outcome <- c(grp_outcome, "predictor")

    if (is_binom) {

        # A) trial-level 0/1 outcome
        if (is.null(trials_id) ) {

            y <- df[[outcome_id]]
            ok <- all(is.na(y) | y %in% c(0, 1, FALSE, TRUE) )

            if (!ok) {

                stop (
                    "For binomial() with trials_id = NULL, `", outcome_id,
                    "` must contain trial-level 0/1 (or TRUE/FALSE).",
                    call. = FALSE
                    )

            }

            df <- df |> dplyr::mutate(outcome_bin = as.integer(.data[[outcome_id]]) )

            summary_data <- df |>
                dplyr::group_by(dplyr::across(dplyr::all_of(grp_outcome) ) ) |>
                dplyr::summarise(
                    success = sum(.data$outcome_bin, na.rm = na_rm),
                    trials  = if (na_rm) sum(!is.na(.data$outcome_bin)) else dplyr::n(),
                    # carry continuous predictor if needed
                    predictor = if (predictor_type == "continuous") agg_continuous(.data$predictor) else dplyr::first(.data$predictor),
                    .groups = "drop"
                    ) |>
                # drop the redundant predictor column for categorical case (already grouped)
                dplyr::mutate(predictor = if (predictor_type == "categorical") .data$predictor else .data$predictor)

            if (predictor_type == "categorical") {

                # predictor already exists via grouping; keep as is

            } else if (predictor_type == "none") {

                summary_data <- summary_data |> dplyr::select(-dplyr::any_of("predictor") )

            }

        } else {

            # B) pre-aggregated success + trials in the raw data
            df <- df |>
                dplyr::mutate(
                    success_in = .data[[outcome_id]],
                    trials_in = .data[[trials_id]]
                    )

            if (!is.numeric(df$success_in) || !is.numeric(df$trials_in) ) {

                stop ("For binomial() with trials_id provided, both success and trials must be numeric.", call. = FALSE)

            }

            summary_data <- df |>
                dplyr::group_by(dplyr::across(dplyr::all_of(grp_outcome) ) ) |>
                dplyr::summarise(
                    success = sum(.data$success_in, na.rm = na_rm),
                    trials = sum(.data$trials_in,  na.rm = na_rm),
                    predictor = if (predictor_type == "continuous") agg_continuous(.data$predictor) else dplyr::first(.data$predictor),
                    .groups = "drop"
                    )

            if (predictor_type == "categorical") {

                # OK

            } else if (predictor_type == "none") {

                summary_data <- summary_data |> dplyr::select(-dplyr::any_of("predictor") )

            }

        }

        # sanity checks
        if (any(summary_data$trials < 0, na.rm = TRUE) ) stop ("Computed `trials` contains negative values.", call. = FALSE)
        if (any(summary_data$success < 0, na.rm = TRUE) ) stop ("Computed `success` contains negative values.", call. = FALSE)
        if (any(summary_data$success > summary_data$trials, na.rm = TRUE) ) stop ("Some cells have success > trials after summarisation.", call. = FALSE)

    } else { # gaussian summary

        # if outcome_sem is provided, assume data is already summarised and
        # only validate + return it in the standard format.
        if (!is.null(outcome_sem) ) {

            # build "summary-like" output directly
            summary_data <- df |>
                dplyr::transmute(
                    participant = .data$participant,
                    time = .data$time,
                    predictor = if (!is.na(predictor_id)) .data$predictor else NULL,
                    outcome_mean = .data[[outcome_id]],
                    outcome_sem = .data[[outcome_sem]]
                    )

            # type checks
            if (!is.numeric(summary_data$outcome_mean) ) {

                stop ("`", outcome_id, "` must be numeric when providing `outcome_sem`.", call. = FALSE)

            }

            if (!is.numeric(summary_data$outcome_sem) ) {

                stop ("`", outcome_sem, "` must be numeric (SEM values).", call. = FALSE)

            }

            if (any(summary_data$outcome_sem < 0, na.rm = TRUE) ) {

                stop ("`outcome_sem` contains negative values; SEM must be >= 0.", call. = FALSE)

            }

            # enforce uniqueness of cells (what summarisation would have produced)
            key_cols <- grp_outcome

            # for continuous predictors, we do NOT group by predictor, but we do require
            # predictor to be single-valued per (participant, time) row (already true
            # if keys are unique)
            dup <- summary_data |>
                dplyr::count(dplyr::across(dplyr::all_of(key_cols) ) ) |>
                dplyr::filter(.data$n > 1)

            if (nrow(dup) > 0) {

                stop (
                    "When `outcome_sem` is provided, the input data must already be summarised ",
                    "with one row per cell: ",
                    paste(key_cols, collapse = " x "),
                    ". Found duplicated cells.",
                    call. = FALSE
                    )

            }

            # if predictor absent
            if (predictor_type == "none") {

                summary_data <- summary_data |> dplyr::select(-dplyr::any_of("predictor") )

            }

        } else {

            df <- df |> dplyr::mutate(outcome = .data[[outcome_id]])

            if (!is.na(predictor_id) ) {

                summary_data <- df |>
                    dplyr::group_by(dplyr::across(dplyr::all_of(grp_outcome) ) ) |>
                    dplyr::summarise(
                        outcome_mean = mean(.data$outcome, na.rm = na_rm),
                        # outcome_sd = stats::sd(.data$outcome, na.rm = na_rm),
                        # outcome_sem = stats::sd(.data$outcome, na.rm = na_rm) / sqrt(length(.data$outcome) ),
                        outcome_sem = stats::sd(.data$outcome, na.rm = na_rm) / sqrt(sum(!is.na(.data$outcome) ) ),
                        predictor = if (predictor_type == "continuous") agg_continuous(.data$predictor) else dplyr::first(.data$predictor),
                        .groups = "drop"
                        )

            } else {

                summary_data <- df |>
                    dplyr::group_by(dplyr::across(dplyr::all_of(grp_outcome) ) ) |>
                    dplyr::summarise(
                        outcome_mean = mean(.data$outcome, na.rm = na_rm),
                        # outcome_sd = stats::sd(.data$outcome, na.rm = na_rm),
                        # outcome_sem = stats::sd(.data$outcome, na.rm = na_rm) / sqrt(length(.data$outcome) ),
                        outcome_sem = stats::sd(.data$outcome, na.rm = na_rm) / sqrt(sum(!is.na(.data$outcome) ) ),
                        .groups = "drop"
                        )

            }

        }

    }

    # within/between classification (extended to continuous predictors)
    if (!is.na(predictor_id) ) {

        if (predictor_type == "categorical") {

            if (exists("check_within_between", mode = "function") ) {

                within_between <- check_within_between(
                    data = summary_data,
                    participant = "participant",
                    predictor = "predictor"
                    )

            } else {

                # within_between <- NULL
                within_between <- NA

            }

        } else if (predictor_type == "continuous") {

            # between-subject if each participant has ~1 unique value (ignoring NA)
            wb <- summary_data |>
                dplyr::group_by(.data$participant) |>
                dplyr::summarise(
                    n_unique = dplyr::n_distinct(.data$predictor[!is.na(.data$predictor)]),
                    .groups = "drop"
                    )

            within_between <- list(
                classification = if (all(wb$n_unique <= 1) ) "between-subject" else "within-subject",
                predictor_type = "continuous"
                )

        }

    } else {

        within_between <- NA

    }

    return (
        list(
            data = summary_data,
            within_between = within_between,
            predictor_type = predictor_type
            )
        )

}

#' Select approximately equally spaced times
#' @param time_vec Numeric vector.
#' @param N Integer number of points.
#' @return Subset of \code{time_vec}.
#' @export
st_take_n_times <- function (time_vec, N) {

    stopifnot(is.numeric(time_vec), length(time_vec) > 0L)
    stopifnot(is.numeric(N), length(N) == 1L, N >= 1)

    idx <- unique(round(seq.int(1L, length(time_vec), length.out = N) ) )

    return (time_vec[idx])

}

#' Test whether a predictor varies (or not) per participant
#' @keywords internal
check_within_between <- function (
        data,
        participant,
        predictor,
        tol = 0,
        min_within_prop = 0.9
        ) {

    # basic checks
    stopifnot(
        is.character(participant), length(participant) == 1,
        is.character(predictor), length(predictor) == 1
        )

    if (!participant %in% names(data) ) {

        stop (sprintf("Column '%s' not found in `data`.", participant), call. = FALSE)

    }

    if (!predictor %in% names(data) ) {

        stop (sprintf("Column '%s' not found in `data`.", predictor), call. = FALSE)

    }

    # drop missing predictor values
    df <- data |>
        dplyr::select(
            participant = .data[[participant]],
            predictor = .data[[predictor]]
            ) |>
        dplyr::filter(!is.na(.data$predictor) )

    is_cat <- is.factor(df$predictor) || is.character(df$predictor)

    # classify per participant
    per_participant <- df |>
        dplyr::group_by(.data$participant) |>
        dplyr::summarise(
            n_obs = dplyr::n(),
            n_unique = dplyr::n_distinct(.data$predictor),
            sd_within = if (!is_cat) stats::sd(.data$predictor) else NA_real_,
            min = if (!is_cat) min(.data$predictor) else NA_real_,
            max = if (!is_cat) max(.data$predictor) else NA_real_,
            .groups = "drop"
            ) |>
        dplyr::mutate(
            varies_within = if (is_cat) {
                .data$n_unique > 1
            } else {
                (.data$max - .data$min) > tol
            })

    # summary stats
    prop_within <- mean(per_participant$varies_within)

    classification <- dplyr::case_when(
        prop_within == 0 ~ "between-subject",
        prop_within >= min_within_prop ~ "within-subject",
        TRUE ~ "mixed (within + between)"
        )

    return (
        list(
            predictor_type = if (is_cat) "categorical" else "continuous",
            classification = classification,
            proportion_within = prop_within,
            per_participant = per_participant
            )
        )

}

#' @keywords internal
st_head_radius <- function (sensors, head_expand = 1.1) {

    # enforce numeric
    x <- as.numeric(sensors$xproj)
    y <- as.numeric(sensors$yproj)

    ok <- is.finite(x) & is.finite(y)
    if (!any(ok) ) stop ("No finite sensor coordinates to compute head radius.", call. = FALSE)

    rr <- sqrt(sensors$xproj^2 + sensors$yproj^2)
    radius <- head_expand * max(rr, na.rm = TRUE)

    if (!is.finite(radius) || radius <= 0) stop ("Computed head radius is not finite/positive.", call. = FALSE)

    return (radius)

}

#' @keywords internal
st_head_outline <- function (r, n = 600) {

    th <- seq(0, 2*pi, length.out = n)

    head <- data.frame(
        x = r * cos(th),
        y = r * sin(th),
        part = "head"
        )

    # Nose base points *on* the circle
    a1 <- 80 * pi / 180
    a2 <- 100 * pi / 180
    nose_base1 <- c(r * cos(a1), r * sin(a1) )
    nose_base2 <- c(r * cos(a2), r * sin(a2) )
    nose_tip <- c(0, 1.08 * r)

    nose <- data.frame(
        x = c(nose_base1[1], nose_tip[1], nose_base2[1]),
        y = c(nose_base1[2], nose_tip[2], nose_base2[2]),
        part = "nose"
        )

    ear_arc <- function (theta_start, theta_end, side = c("left", "right") ) {

        side <- match.arg(side)
        t <- seq(theta_start, theta_end, length.out = 120)
        anchor_x <- r * cos(t)
        anchor_y <- r * sin(t)
        bump <- 0.10 * r

        data.frame(
            x = anchor_x + bump * cos(t),
            y = anchor_y + bump * sin(t),
            part = paste0("ear_", side)
            )

    }

    left_ear  <- ear_arc(160*pi/180, 200*pi/180, "left")
    right_ear <- ear_arc(-20*pi/180,  20*pi/180, "right")

    return (rbind(head, nose, left_ear, right_ear) )

}

#' @keywords internal
st_make_grid <- function (x_range, y_range, grid_res = 100) {

    return (
        expand.grid(
            xproj = seq(x_range[1], x_range[2], length.out = grid_res),
            yproj = seq(y_range[1], y_range[2], length.out = grid_res)
            )
        )

}

#' @keywords internal
st_order_time_factor <- function (time_labels) {

    # order facets numerically even if labels are like "Time: 28 s"
    nums <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", time_labels) ) )

    return (factor(time_labels, levels = unique(time_labels[order(nums)]) ) )

}

#' @keywords internal
st_build_time_labels <- function (times, unit = "s", prefix = "Time: ") {

    return (paste0(prefix, format(times, trim = TRUE), " ", unit) )

}

#' @keywords internal
#' @importFrom stats fitted
st_predict_brms <- function (fit, newdata, re_formula = NULL, probs = c(0.025, 0.975), ndraws = NULL) {

    stopifnot(inherits(fit, "brmsfit") )

    fe <- fitted(
        fit,
        newdata = newdata,
        re_formula = re_formula,
        summary = TRUE,
        probs = probs,
        ndraws = ndraws
        )

    return (cbind(newdata, fe) )

}

#' @keywords internal
st_compute_limits <- function (
        df,
        z = "Estimate",
        limits = c("global", "global_quantile", "none"),
        limit_quantiles = c(0.01, 0.99)
        ) {

    limits <- match.arg(limits)

    if (limits == "none") return (NULL)
    if (limits == "global") return (c(-1, 1) * max(abs(df[[z]]), na.rm = TRUE) )

    qtls <- as.numeric(stats::quantile(df[[z]], probs = limit_quantiles, na.rm = TRUE) )

    return (c(-1, 1) * max(abs(qtls) ) )

}

#' @keywords internal
st_make_head_grid <- function (sensors, grid_res = 100, head_expand = 1.1) {

    r_head <- st_head_radius(sensors, head_expand = head_expand)

    return(
        list(
            r_head = r_head,
            grid_xy = st_make_grid(c(-r_head, r_head), c(-r_head, r_head), grid_res = grid_res)
            )
        )

}

#' @keywords internal
st_fill_head <- function (df, r_head, fill_na = c("none", "nearest") ) {

    fill_na <- match.arg(fill_na)

    df <- df |> dplyr::mutate(.inside = sqrt(.data$xproj^2 + .data$yproj^2) <= r_head)

    if (fill_na == "nearest") {

        # fill NAs *inside* the head only
        na_idx <- which(df$.inside & is.na(df$Estimate) )
        ok_idx <- which(df$.inside & !is.na(df$Estimate) )

        if (length(na_idx) > 0 && length(ok_idx) > 0) {

            nn <- FNN::get.knnx(
                data  = as.matrix(df[ok_idx, c("xproj", "yproj")]),
                query = as.matrix(df[na_idx, c("xproj", "yproj")]),
                k = 1
                )

            df$Estimate[na_idx] <- df$Estimate[ok_idx][nn$nn.index[, 1]]

        }

    }

    # mask outside head
    df$Estimate[!df$.inside] <- NA_real_
    df$.inside <- NULL

    return (df)

}

#' @keywords internal
st_interp_to_grid <- function (
        df_t, grid_xy, value = "Estimate",
        r_head = NULL,
        fill_na = c("none", "nearest"),
        extrap = FALSE
        ) {

    fill_na <- match.arg(fill_na)

    xo <- sort(unique(grid_xy$xproj) )
    yo <- sort(unique(grid_xy$yproj) )

    # z <- akima::interp(
    #     x = df_t$xproj,
    #     y = df_t$yproj,
    #     z = df_t[[value]],
    #     xo = xo,
    #     yo = yo,
    #     duplicate = "mean",
    #     linear = TRUE,
    #     extrap = extrap
    #     )

    z <- interp::interp(
        x = df_t$xproj,
        y = df_t$yproj,
        z = df_t[[value]],
        xo = xo,
        yo = yo,
        duplicate = "mean",
        method = "linear",
        extrap = extrap
        )

    out <- expand.grid(xproj = z$x, yproj = z$y)
    out$Estimate <- as.vector(z$z)

    # optional: fill interior NAs using nearest neighbour
    if (fill_na == "nearest" && anyNA(out$Estimate) ) {

        na_idx <- which(is.na(out$Estimate) )
        ok_idx <- which(!is.na(out$Estimate) )

        nn <- FNN::get.knnx(
            data  = as.matrix(out[ok_idx, c("xproj", "yproj")]),
            query = as.matrix(out[na_idx, c("xproj", "yproj")]),
            k = 1
            )

        out$Estimate[na_idx] <- out$Estimate[ok_idx][nn$nn.index[, 1]]

    }

    # apply head mask at the end (outside head = NA)
    if (!is.null(r_head) ) {

        out$Estimate[sqrt(out$xproj^2 + out$yproj^2) > r_head] <- NA_real_

    }

    return (out)

}

#' @keywords internal
get_baseline_null <- function (post_diff, baseline_window, credible_interval) {

    alpha <- (1 - credible_interval) / 2
    lower_q <- alpha
    upper_q <- 1 - alpha

    post_diff |>
        dplyr::filter(
            .data$time >= baseline_window[1],
            .data$time <= baseline_window[2]
            ) |>
        dplyr::group_by(.data$time) |>
        dplyr::summarise(
            lower = stats::quantile(.data$epred_diff, lower_q),
            median = stats::quantile(.data$epred_diff, 0.5),
            upper = stats::quantile(.data$epred_diff, upper_q),
            .groups = "drop"
            ) |>
        # dplyr::summarise(null_value = max(abs(c(.data$lower, .data$upper) ) ) ) |>
        dplyr::summarise(null_value = max(abs(.data$median) ) ) |>
        dplyr::pull(.data$null_value)

}

#' @keywords internal
resolve_temporal_rope <- function (
    post_draws,
    rope = NULL,
    baseline_window = NULL,
    credible_interval = 0.95,
    value_col = ".epred",
    null_center = 0
    ) {

    # An explicit ROPE, including zero, overrides baseline estimation.
    if (!is.null(rope)) {

        if (!is.numeric(rope) || length(rope) != 1L ||
            !is.finite(rope) || rope < 0) {

            stop(
                "`rope` must be NULL or a nonnegative finite number.",
                call. = FALSE
            )

        }

        return(list(rope = rope, source = "user"))

    }

    # No baseline and no explicit ROPE: point-null inference.
    if (is.null(baseline_window)) {

        return(list(rope = 0, source = "point_null"))

    }

    if (!is.numeric(baseline_window) ||
        length(baseline_window) != 2L ||
        any(!is.finite(baseline_window)) ||
        baseline_window[1L] >= baseline_window[2L]) {

        stop(
            "`baseline_window` must contain two increasing finite numbers.",
            call. = FALSE
        )

    }

    if (!all(c("time", value_col) %in% names(post_draws))) {

        stop(
            "Baseline estimation requires time and posterior draws.",
            call. = FALSE
        )

    }

    baseline_draws <- data.frame(
        time = post_draws$time,
        epred_diff = post_draws[[value_col]] - null_center
    )

    effective_rope <- get_baseline_null(
        post_diff = baseline_draws,
        baseline_window = baseline_window,
        credible_interval = credible_interval
    )

    if (length(effective_rope) != 1L ||
        !is.finite(effective_rope) ||
        effective_rope < 0) {

        stop(
            "Could not estimate the ROPE from the baseline window.",
            call. = FALSE
        )

    }

    return (
        list(
            rope = effective_rope,
            source = "baseline"
            )
        )

}

#' @keywords internal
compute_one_sample_prob <- function(
    post_draws,
    null_value,
    participant_clusters,
    n_post_samples,
    credible_interval,
    baseline_window = NULL,
    odds_mode = "directional_rope",
    rope = NULL
) {
    groups <- if (participant_clusters) c("time", "participant") else "time"

    # Preserve the historical one-sided baseline convention when explicitly
    # selected and no user ROPE has been supplied.
    legacy_baseline <- identical(odds_mode, "legacy_one_sided") &&
        is.null(rope) && !is.null(baseline_window)
    effective_null <- null_value

    if (legacy_baseline) {
        baseline <- post_draws[
            post_draws$time >= baseline_window[1L] &
                post_draws$time <= baseline_window[2L],
            , drop = FALSE
        ]
        if (!nrow(baseline)) {
            stop("No observations found within `baseline_window`.",
                 call. = FALSE)
        }
        alpha <- (1 - credible_interval) / 2
        limits <- baseline |>
            dplyr::group_by(dplyr::across(dplyr::all_of(groups))) |>
            dplyr::summarise(
                lower = stats::quantile(.data$.epred, probs = alpha),
                upper = stats::quantile(.data$.epred, probs = 1 - alpha),
                .groups = "drop"
            )
        effective_null <- max(abs(c(limits$lower, limits$upper)))
        effective <- list(rope = 0, source = "legacy_baseline_null")
    } else {
        effective <- resolve_temporal_rope(
            post_draws = post_draws,
            rope = rope,
            baseline_window = baseline_window,
            credible_interval = credible_interval,
            value_col = ".epred",
            null_center = null_value
        )
    }

    out <- summarise_posterior_long(
        data = post_draws,
        groups = groups,
        value = ".epred",
        null_value = effective_null,
        rope = effective$rope,
        credible_interval = credible_interval,
        odds_mode = odds_mode,
        odds_cap = n_post_samples
    )
    out$post_prob <- out$estimate
    out <- out[c(
        groups, "prob_above", "prob_below", "prob_inside",
        "odds_above", "odds_below", "prob_ratio",
        "post_prob", "lower", "upper"
    )]
    attr(out, "effective_rope") <- effective$rope
    attr(out, "rope_source") <- effective$source
    out
}

#' @keywords internal
compute_two_sample_prob <- function(
    post_draws,
    null_value,
    participant_clusters,
    n_post_samples,
    credible_interval,
    predictor_type,
    baseline_window = NULL,
    odds_mode = "directional_rope",
    rope = NULL
) {
    groups <- if (participant_clusters) c("time", "participant") else "time"
    paired <- pair_posterior_contrasts(post_draws, groups, predictor_type)
    # Preserve the old low-level interpretation of `null_value` as a
    # contrast width if no baseline and no explicit `rope` were supplied.
    selected_rope <- if (is.null(rope) && is.null(baseline_window)) {
        null_value
    } else {
        rope
    }
    effective <- resolve_temporal_rope(
        post_draws = paired,
        rope = selected_rope,
        baseline_window = baseline_window,
        credible_interval = credible_interval,
        value_col = "epred_diff",
        null_center = 0
    )

    out <- summarise_posterior_long(
        data = paired,
        groups = groups,
        value = "epred_diff",
        null_value = 0,
        rope = effective$rope,
        credible_interval = credible_interval,
        odds_mode = odds_mode,
        odds_cap = n_post_samples
    )
    out$post_prob <- out$estimate
    out <- out[c(
        groups, "prob_above", "prob_below", "prob_inside",
        "odds_above", "odds_below", "prob_ratio",
        "post_prob", "lower", "upper"
    )]
    attr(out, "effective_rope") <- effective$rope
    attr(out, "rope_source") <- effective$source
    out
}

#' @keywords internal
add_required_dummy <- function (df, is_binom) {

    if (is_binom) {

        # required by success | trials(trials)
        if (!"trials" %in% names(df) ) df <- dplyr::mutate(df, trials = 1)

    } else {

        # required by outcome_mean | se(outcome_sem)
        if (!"outcome_sem" %in% names(df) ) df <- dplyr::mutate(df, outcome_sem = 1)

    }

    return (df)

}

#' @keywords internal
compute_one_sample_prob_2d <- function (
    post_draws,
    null_value,
    participant_clusters,
    n_post_samples,
    credible_interval,
    baseline_window = NULL,
    odds_mode = "directional_rope",
    rope = 0
    ) {

    groups <- c("time1", "time2")

    if (participant_clusters) {

        groups <- c(groups, "participant")

    }

    # Existing 2D one-sample implementation does not use baseline_window.
    out <- summarise_posterior_long(
        data = post_draws,
        groups = groups,
        value = ".epred",
        null_value = null_value,
        rope = rope,
        credible_interval = credible_interval,
        odds_mode = odds_mode,
        odds_cap = n_post_samples
        )

    out$post_prob <- out$estimate

    # return (
    #     out[c(
    #         groups,
    #         "prob_above", "prob_below", "prob_inside",
    #         "odds_above", "odds_below", "prob_ratio",
    #         "post_prob", "lower", "upper"
    #         )]
    #     )

    out <- out[c(
        groups,
        "prob_above", "prob_below", "prob_inside",
        "odds_above", "odds_below", "prob_ratio",
        "post_prob", "lower", "upper"
    )]

    attr(out, "effective_rope") <- rope
    attr(out, "rope_source") <- if (rope == 0) {
        "point_null"
    } else {
        "user"
    }

    return (out)

}

#' @keywords internal
compute_two_sample_prob_2d <- function (
    post_draws,
    null_value,
    participant_clusters,
    n_post_samples,
    credible_interval,
    predictor_type,
    time_id = c("time1", "time2"),
    baseline_window = NULL,
    odds_mode = "directional_rope",
    rope = 0
    ) {

    groups <- c("time1", "time2")

    if (participant_clusters) {

        groups <- c(groups, "participant")

    }

    paired <- pair_posterior_contrasts(post_draws, groups, predictor_type)

    # Existing 2D contrast implementation uses one-sided posterior odds,
    # unlike the 1D directional contrast implementation.
    out <- summarise_posterior_long(
        data = paired,
        groups = groups,
        value = "epred_diff",
        null_value = null_value,
        rope = rope,
        credible_interval = credible_interval,
        odds_mode = odds_mode,
        odds_cap = n_post_samples
        )

    out$post_prob <- out$estimate

    # return (
    #     out[c(
    #         groups,
    #         "prob_above", "prob_below", "prob_inside",
    #         "odds_above", "odds_below", "prob_ratio",
    #         "post_prob", "lower", "upper"
    #         )]
    #     )

    out <- out[c(
        groups,
        "prob_above", "prob_below", "prob_inside",
        "odds_above", "odds_below", "prob_ratio",
        "post_prob", "lower", "upper"
    )]

    attr(out, "effective_rope") <- rope
    attr(out, "rope_source") <- if (rope == 0) {
        "point_null"
    } else {
        "user"
    }

    return (out)

}

# null-coalescing (avoid importing rlang just for %||%)
#' @keywords internal
`%||%` <- function (a, b) if (!is.null(a) ) a else b

# plotting 2D temporal clusters
#' @keywords internal
plot_clusters_2d <- function (
        field_df,
        clusters_df = NULL,
        value_col = "prob_ratio",
        palette = "vik",
        midpoint = 0.5,
        fill_name = NULL,
        axes_labels = c("Training time", "Testing time"),
        cluster_colour = "white"
        ) {

    stopifnot(length(axes_labels) == 2)

    hit_df <- clusters_df |> dplyr::transmute(.data$time1, .data$time2, hit = 1)

    plot_df <- field_df |>
        dplyr::left_join(hit_df, by = c("time1", "time2") ) |>
        dplyr::mutate(hit = tidyr::replace_na(.data$hit, 0) )

    p <- field_df |>
        ggplot2::ggplot(
            ggplot2::aes(
                x = .data$time1,
                y = .data$time2,
                fill = .data[[value_col]]
                )
            ) +
        ggplot2::geom_raster() +
        ggplot2::coord_equal(expand = FALSE) +
        ggplot2::geom_abline(intercept = 0, slope = 1, linetype = 2, color = cluster_colour) +
        ggplot2::labs(x = axes_labels[1], y = axes_labels[2], fill = value_col) +
        scico::scale_fill_scico(
            palette = palette,
            midpoint = midpoint,
            oob = scales::squish,
            name = fill_name
            )

    if (!is.null(clusters_df) && nrow(clusters_df) > 0) {

        p <- p +
            ggplot2::geom_contour(
                data = plot_df,
                ggplot2::aes(x = .data$time1, y = .data$time2, z = .data$hit),
                breaks = 0.5,
                colour = cluster_colour,
                linewidth = 0.5,
                inherit.aes = FALSE
                )

    }

    return (p)

}

# create data grid for model predictions
#' @keywords internal
make_prediction_grid_from_summary <- function (
        summary_data,
        brms_gam,
        predictor_id,
        predictor_type = c("none", "categorical", "continuous"),
        participant_clusters = FALSE,
        within_between = NULL,
        is_2d_time = FALSE,
        is_binom = FALSE,
        include_ar_term = FALSE,
        continuous_mode = c("pm1sd", "observed")
        ) {

    # argument checks (match testing_through_time() style)
    stopifnot("summary_data must be a dataframe..." = is.data.frame(summary_data) )
    stopifnot("brms_gam must be a brmsfit..." = inherits(brms_gam, "brmsfit") )
    stopifnot("predictor_id must be a character (or NA)..." = is.character(predictor_id) || is.na(predictor_id) )
    stopifnot("participant_clusters must be logical..." = is.logical(participant_clusters) )
    stopifnot("is_2d_time must be logical..." = is.logical(is_2d_time) )
    stopifnot("is_binom must be logical..." = is.logical(is_binom) )
    stopifnot("include_ar_term must be logical..." = is.logical(include_ar_term) )

    predictor_type <- match.arg(predictor_type)
    continuous_mode <- match.arg(continuous_mode)

    # expected time columns
    time_cols <- if (is_2d_time) c("time1", "time2") else "time"

    assertthat::assert_that(
        all(time_cols %in% names(summary_data) ),
        msg = paste(
            "Missing time column(s) in summary_data:",
            paste(setdiff(time_cols, names(summary_data) ), collapse = ", ")
            )
        )

    # base grid = time support
    base_grid <- summary_data |>
        dplyr::distinct(dplyr::across(dplyr::all_of(time_cols) ) ) |>
        dplyr::arrange(dplyr::across(dplyr::all_of(time_cols) ) )

    # add predictor dimension
    if (is.na(predictor_id) || predictor_type == "none") {

        newdata_grid <- base_grid

    } else if (predictor_type == "categorical") {

        assertthat::assert_that(
            "predictor" %in% names(summary_data),
            msg = "summary_data must contain a 'predictor' column for categorical predictors."
            )

        pred_levels <- if (is.factor(summary_data$predictor) ) {
            levels(summary_data$predictor)
        } else {
            sort(unique(summary_data$predictor) )
        }

        newdata_grid <- tidyr::crossing(
            base_grid,
            predictor = pred_levels
            )

    } else if (predictor_type == "continuous") {

        assertthat::assert_that(
            "predictor" %in% names(summary_data),
            msg = "summary_data must contain a 'predictor' column for continuous predictors."
            )

        if (continuous_mode == "observed") {

            pred_vals <- sort(unique(summary_data$predictor) )

        } else {

            # use model data for mean/sd
            predictor_mean <- mean(brms_gam$data$predictor, na.rm = TRUE)
            predictor_sd <- stats::sd(brms_gam$data$predictor, na.rm = TRUE)
            pred_vals <- c(predictor_mean - predictor_sd, predictor_mean + predictor_sd)

        }

        newdata_grid <- tidyr::crossing(
            base_grid,
            predictor = pred_vals
            )

    }

    # detect whether ar_series is required
    needs_ar_series <- include_ar_term || ("ar_series" %in% names(brms_gam$data) )

    # participant dimension
    if (participant_clusters) {

        assertthat::assert_that(
            "participant" %in% names(summary_data),
            msg = "summary_data must contain 'participant' when participant_clusters = TRUE."
            )

        # no participant clusters for between-subject predictors
        if (!is.null(within_between) &&
            !is.na(predictor_id) &&
            is.list(within_between) &&
            !is.null(within_between$classification) &&
            identical(within_between$classification, "between-subject") ) {

            stop (
                "participant_clusters = TRUE is not supported for between-subject predictors.",
                call. = FALSE
                )

        }

        newdata_grid <- tidyr::crossing(
            newdata_grid,
            participant = sort(unique(summary_data$participant) )
            )

    } else {

        # For population-level predictions:
        # - if AR is on and participant exists: use a reference participant so ar_series is well-defined
        # - otherwise: keep old behaviour (participant = NA)
        if ("participant" %in% names(summary_data) ) {

            if (needs_ar_series) {

                ref_participant <- brms_gam$data$participant[!is.na(brms_gam$data$participant)][1]

                newdata_grid <- newdata_grid |>
                    dplyr::mutate(participant = ref_participant)

            } else {

                newdata_grid <- newdata_grid |>
                    dplyr::mutate(participant = NA)

            }

        }

    }

    # ar_series (when needed)
    if (needs_ar_series) {

        if ("participant" %in% names(newdata_grid) ) {

            if ("predictor" %in% names(newdata_grid) ) {

                newdata_grid <- newdata_grid |>
                    dplyr::mutate(ar_series = interaction(.data$participant, .data$predictor, drop = TRUE) )

            } else {

                newdata_grid <- newdata_grid |>
                    dplyr::mutate(ar_series = .data$participant)

            }

        } else {

            # fall-back dummy
            newdata_grid <- newdata_grid |>
                dplyr::mutate(ar_series = "ar_dummy")

        }

        # match factor levels to model data if needed
        if ("ar_series" %in% names(brms_gam$data) && is.factor(brms_gam$data$ar_series) ) {

            newdata_grid$ar_series <- factor(newdata_grid$ar_series, levels = levels(brms_gam$data$ar_series) )

        }

    }

    # add required dummy columns for validate_data()
    newdata_grid <- newdata_grid |>
        add_required_dummy(is_binom = is_binom)

    return (newdata_grid)

}

# compute directional odds
#' @keywords internal
compute_directional_odds <- function (
    post_diff,
    null_value = 0,
    n_post_samples,
    groups = "time",
    odds_mode = "directional_rope"
    ) {

    if (!is.numeric(null_value) || length(null_value) != 1L ||
        !is.finite(null_value) || null_value < 0) {

        stop (
            "`null_value` must be a nonnegative ROPE half-width.",
            call. = FALSE
            )

    }

    out <- summarise_posterior_long(
        data = post_diff,
        groups = groups,
        value = "epred_diff",
        null_value = 0,
        rope = null_value,
        odds_mode = odds_mode,
        odds_cap = n_post_samples
        )

    return (
        out[c(
            groups, "prob_above", "prob_below", "prob_inside",
            "odds_above", "odds_below", "prob_ratio"
            )]
        )

}
