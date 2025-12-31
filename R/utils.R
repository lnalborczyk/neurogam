#' Summarise trial-level data for \code{testing_through_time()}
#' @keywords internal
make_summary_data <- function (
        data,
        participant_id = "participant",
        outcome_id = "eeg",
        outcome_sd = NULL,
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

    # family name (robust)
    fam_name <- tryCatch({

        if (is.list(family) && !is.null(family$family) ) as.character(family$family) else stop ("not a family object")

    }, error = function(e) {

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

    # if user provides an SD column name, require it (gaussian only)
    if (!is.null(outcome_sd) ) {

        if (!is.character(outcome_sd) || length(outcome_sd) != 1) {

            stop ("`outcome_sd` must be NULL or a single character column name.", call. = FALSE)

        }

        required <- c(required, outcome_sd)

    }

    if (is_binom && !is.null(trials_id) ) required <- c(required, trials_id)

    missing_cols <- setdiff(required, names(data) )

    if (length(missing_cols) > 0) {

        stop ("Missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)

    }

    # internal columns
    df <- data |>
        dplyr::mutate(
            participant = .data[[participant_id]],
            time = .data[[time_id]]
            )

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
    # - categorical predictor: group by predictor
    # - continuous predictor: DO NOT group by predictor; keep as a covariate column
    grp_outcome <- c("participant", "time")
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
                    trials_in  = .data[[trials_id]]
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

        # if outcome_sd is provided, assume data is already summarised and
        # only validate + return it in the standard format.
        if (!is.null(outcome_sd) ) {

            # build "summary-like" output directly
            summary_data <- df |>
                dplyr::transmute(
                    participant = .data$participant,
                    time = .data$time,
                    predictor = if (!is.na(predictor_id)) .data$predictor else NULL,
                    outcome_mean = .data[[outcome_id]],
                    outcome_sd = .data[[outcome_sd]]
                    )

            # type checks
            if (!is.numeric(summary_data$outcome_mean) ) {

                stop ("`", outcome_id, "` must be numeric when providing `outcome_sd`.", call. = FALSE)

            }

            if (!is.numeric(summary_data$outcome_sd) ) {

                stop ("`", outcome_sd, "` must be numeric (SD values).", call. = FALSE)

            }

            if (any(summary_data$outcome_sd < 0, na.rm = TRUE) ) {

                stop ("`outcome_sd` contains negative values; SD must be >= 0.", call. = FALSE)

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
                    "When `outcome_sd` is provided, the input data must already be summarised ",
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
                        outcome_sd = stats::sd(.data$outcome, na.rm = na_rm),
                        predictor = if (predictor_type == "continuous") agg_continuous(.data$predictor) else dplyr::first(.data$predictor),
                        .groups = "drop"
                        )

            } else {

                summary_data <- df |>
                    dplyr::group_by(dplyr::across(dplyr::all_of(grp_outcome) ) ) |>
                    dplyr::summarise(
                        outcome_mean = mean(.data$outcome, na.rm = na_rm),
                        outcome_sd = stats::sd(.data$outcome, na.rm = na_rm),
                        .groups = "drop"
                        )

            }

            # if (predictor_type == "categorical") {
            #
            #     # predictor is already grouped; keep
            #
            # } else if (predictor_type == "none") {
            #
            #     summary_data <- summary_data |> dplyr::select(-dplyr::any_of("predictor") )
            #
            # }

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

    z <- akima::interp(
        x = df_t$xproj,
        y = df_t$yproj,
        z = df_t[[value]],
        xo = xo,
        yo = yo,
        duplicate = "mean",
        linear = TRUE,
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
compute_one_sample_prob <- function (
        post_draws, null_value, participant_clusters, n_post_samples, credible_interval
        ) {

    # validate credible interval
    if (!is.numeric(credible_interval) ||

        length(credible_interval) != 1 ||
        credible_interval <= 0 || credible_interval >= 1) {

        stop ("`credible_interval` must be a number strictly between 0 and 1.")

    }

    alpha <- (1 - credible_interval) / 2
    lower_q <- alpha
    upper_q <- 1 - alpha

    if (participant_clusters) {

        prob_y_above <- post_draws |>
            dplyr::group_by(.data$time, .data$participant) |>
            dplyr::summarise(prob_above = mean(.data$.epred > null_value) ) |>
            dplyr::mutate(prob_ratio = .data$prob_above / (1 - .data$prob_above) ) |>
            dplyr::mutate(
                prob_ratio = pmin(.data$prob_ratio, n_post_samples),
                prob_ratio = pmax(.data$prob_ratio, 1 / n_post_samples)
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

        results <- dplyr::left_join(
            prob_y_above, post_prob_slope, by = c("time", "participant")
            )

    } else {

        prob_y_above <- post_draws |>
            dplyr::group_by(.data$time) |>
            dplyr::summarise(
                prob_above = mean(.data$.epred > null_value)
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

#' @keywords internal
compute_two_sample_prob <- function (
        post_draws, null_value, participant_clusters, n_post_samples,
        credible_interval, predictor_type
        ) {

    # validate credible interval
    if (!is.numeric(credible_interval) ||

        length(credible_interval) != 1 ||
        credible_interval <= 0 || credible_interval >= 1) {

        stop ("`credible_interval` must be a number strictly between 0 and 1.")

    }

    # convert continuous predictor to factor
    if (predictor_type == "continuous") {

        post_draws$predictor <- as.character(round(x = post_draws$predictor, digits = 3) )

    }

    alpha <- (1 - credible_interval) / 2
    lower_q <- alpha
    upper_q <- 1 - alpha

    cond1 <- unique(post_draws$predictor)[1]
    cond2 <- unique(post_draws$predictor)[2]

    if (participant_clusters) {

        post_diff <- post_draws |>
            dplyr::select(
                .data$time, .data$predictor, .data$participant,
                .data$.epred, .data$.draw
                ) |>
            tidyr::pivot_wider(names_from = .data$predictor, values_from = .data$.epred) |>
            dplyr::mutate(epred_diff = .data[[cond2]] - .data[[cond1]])

        prob_y_above <- post_diff |>
            dplyr::group_by(.data$time, .data$participant) |>
            dplyr::summarise(
                prob_above = mean(.data$epred_diff > null_value)
                ) |>
            dplyr::mutate(prob_ratio = .data$prob_above / (1 - .data$prob_above) ) |>
            dplyr::ungroup() |>
            dplyr::mutate(
                prob_ratio = pmin(.data$prob_ratio, n_post_samples),
                prob_ratio = pmax(.data$prob_ratio, 1 / n_post_samples)
                ) |>
            data.frame()

        post_prob_slope <- post_diff |>
            dplyr::group_by(.data$time, .data$participant) |>
            dplyr::summarise(
                post_prob = stats::quantile(x = .data$epred_diff, probs = 0.5),
                lower = stats::quantile(x = .data$epred_diff, probs = lower_q),
                upper = stats::quantile(x = .data$epred_diff, probs = upper_q)
                ) |>
            dplyr::ungroup()

        results <- dplyr::left_join(
            prob_y_above, post_prob_slope,
            by = c("time", "participant")
            )

    } else {

        post_diff <- post_draws |>
            dplyr::select(.data$time, .data$predictor, .data$.epred, .data$.draw) |>
            tidyr::pivot_wider(names_from = .data$predictor, values_from = .data$.epred) |>
            dplyr::mutate(epred_diff = .data[[cond2]] - .data[[cond1]])

        prob_y_above <- post_diff |>
            dplyr::group_by(.data$time) |>
            dplyr::summarise(
                prob_above = mean(.data$epred_diff > null_value)
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
                post_prob = stats::quantile(x = .data$epred_diff, probs = 0.5),
                lower = stats::quantile(x = .data$epred_diff, probs = lower_q),
                upper = stats::quantile(x = .data$epred_diff, probs = upper_q)
                ) |>
            dplyr::ungroup()

        results <- dplyr::left_join(prob_y_above, post_prob_slope, by = "time")

    }

    return (results)

}

#' @keywords internal
add_required_dummy <- function (df, is_binom) {

    if (is_binom) {

        # required by success | trials(trials)
        if (!"trials" %in% names(df) ) df <- dplyr::mutate(df, trials = 1)

    } else {

        # required by outcome_mean | se(outcome_sd)
        if (!"outcome_sd" %in% names(df) ) df <- dplyr::mutate(df, outcome_sd = 1)

    }

    return (df)

}
