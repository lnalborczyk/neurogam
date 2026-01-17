#' Check residual autocorrelation for brms models
#'
#' Compute and visualise residual autocorrelation diagnostics (ACF and partial ACF)
#' for \code{brms} models, optionally comparing raw residuals to "corrected" residuals
#' obtained by removing the contribution of autoregressive (AR) terms.
#'
#' This function is designed for models fit with \code{brms} that may include
#' autoregressive structures (e.g., \code{ar()} terms). Residual autocorrelation is
#' evaluated per time series (defined by \code{series_id} or fallback groupings) and
#' summarised using posterior draws of residuals.
#'
#' @param fit A fitted \code{brmsfit} object.
#'
#' @param data Optional \code{data.frame} used to compute diagnostics. Defaults to
#'   \code{fit$data}. Must match the data used to fit \code{fit} (same rows and order),
#'   otherwise residual draws will not align with \code{data}.
#'
#' @param time_id Name of the time column in \code{data}. Must be present. The time
#'   variable is used to order observations within each series prior to computing
#'   correlations.
#'
#' @param series_id Name of the series (time-series grouping) column in \code{data}.
#'   If present, it is used as the grouping variable. If absent, a series factor is
#'   constructed using \code{participant_id} and \code{predictor_id} if available (see
#'   Details).
#'
#' @param participant_id Name of a participant identifier column. Used to construct
#'   a fallback series grouping if \code{series_id} is not present.
#'
#' @param predictor_id Name of a predictor identifier column. Used to construct
#'   a fallback series grouping if \code{series_id} is not present.
#'
#' @param resid_method Residual method passed to \code{residuals.brmsfit()} (argument
#'   \code{method}). Defaults to \code{"posterior_epred"}.
#'
#' @param resid_type Residual type passed to \code{residuals.brmsfit()} (argument
#'   \code{resid_type}). One of \code{"ordinary"} or \code{"pearson"}.
#'
#' @param ndraws Number of posterior draws of residuals to request from
#'   \code{residuals.brmsfit(..., summary = FALSE)}. These draws are used both to compute
#'   the pointwise residual summary stored in \code{data} and to compute ACF/PACF bands.
#'
#' @param use_posterior_mean Logical; if \code{TRUE}, raw residuals stored in
#'   \code{data} are computed as the posterior mean across residual draws. If \code{FALSE},
#'   the posterior median is used.
#'
#' @param ar_summary Summary used to form a single vector of AR coefficients for the
#'   "point" corrected residuals stored in \code{data}. One of \code{"mean"} or
#'   \code{"median"}.
#'
#' @param max_lag Maximum lag for ACF/PACF computation (positive integer). Lags are
#'   computed from 1 to \code{max_lag}.
#'
#' @param n_series_plot Number of series to randomly sample (without replacement) for
#'   plotting ACF/PACF panels and residual time series panels.
#'
#' @param seed Optional random seed used for sampling series and posterior draws. If
#'   \code{NULL}, sampling is not seeded.
#'
#' @param acf_ndraws Number of posterior draws to subsample when computing ACF/PACF
#'   summaries. If \code{acf_ndraws > ndraws}, the maximum available
#'   draws are used.
#'
#' @param acf_ci Credible interval width for ACF/PACF bands (e.g., \code{0.95} for a 95%
#'   interval). Must lie in \code{(0, 1)}.
#'
#' @param theme A \code{ggplot2} theme object added to plots.
#'
#' @param verbose Logical; if \code{TRUE}, emit messages and warnings about grouping,
#'   duplicated time values, or posterior draw alignment.
#'
#' @details
#' \strong{Series grouping.} If \code{series_id} is present in \code{data}, it defines
#' the grouping. Otherwise, the function tries to construct \code{.series} as:
#' \enumerate{
#'   \item \code{interaction(participant_id, predictor_id)} if both columns exist;
#'   \item \code{participant_id} if only participant exists;
#'   \item \code{predictor_id} if only predictor exists;
#'   \item a single-level factor \code{"group"} if none exist.
#' }
#'
#' \strong{Corrected residuals.} If AR parameters are detected in the posterior draws
#' (columns matching \code{"^ar\\["}), corrected residuals are computed per series as:
#' \deqn{e_t^{(corr)} = e_t - \sum_{k=1}^{p} \phi_k e_{t-k}}
#' using either the posterior mean or median of \eqn{\phi_k} (controlled by
#' \code{ar_summary}) for the \emph{point} corrected residuals stored in \code{data}.
#'
#' @return
#' A list of class \code{"neurogam_autocor_check"} with components:
#' \describe{
#'   \item{ar_summary}{A \code{data.frame} summarising posterior AR parameters (if present).}
#'   \item{rho1_by_series_raw}{Lag-1 autocorrelation per series computed from \code{.resid_raw}.}
#'   \item{rho1_by_series_corrected}{Lag-1 autocorrelation per series computed from
#'     \code{.resid_corrected} (NA if no AR parameters).}
#'   \item{acf_raw}{A \code{data.frame} containing ACF summaries for raw residuals, with
#'     columns \code{.series}, \code{lag}, \code{q_low}, \code{q50}, \code{q_high}, \code{point}.}
#'   \item{acf_corrected}{Same as \code{acf_raw} but for corrected residuals (empty/NA if no AR).}
#'   \item{pacf_raw}{A \code{data.frame} containing PACF summaries for raw residuals, with
#'     columns \code{.series}, \code{lag}, \code{q_low}, \code{q50}, \code{q_high}, \code{point}.}
#'   \item{pacf_corrected}{Same as \code{pacf_raw} but for corrected residuals (empty/NA if no AR).}
#'   \item{plots}{A named list of \code{ggplot} objects (ACF/PACF panels, rho1 summaries,
#'     residual time-series plots, and AR density plots if applicable).}
#' }
#'
#' @examples
#' \dontrun{
#' # fit is a brmsfit with an AR(1) term
#' out <- check_residual_autocorrelation(
#'     fit = fit,
#'     time_id = "time",
#'     series_id = "ar_series"
#'     )
#'
#' # inspect ACF/PACF summaries
#' head(out$acf_raw)
#' head(out$pacf_raw)
#'
#' # plot ACF/PACF panels
#' out$plots$acf_raw
#' out$plots$acf_corrected
#' }
#'
#' @author Ladislas Nalborczyk \email{ladislas.nalborczyk@@cnrs.fr}.
#'
#' @export
check_residual_autocorrelation <- function (
        fit,
        data = NULL,
        time_id = "time",
        series_id = "ar_series",
        participant_id = "participant",
        predictor_id = "predictor",
        resid_method = "posterior_epred",
        resid_type = c("ordinary", "pearson"),
        ndraws = 1000,
        use_posterior_mean = TRUE,
        ar_summary = c("mean", "median"),
        max_lag = 20,
        n_series_plot = 9,
        seed = 666,
        acf_ndraws = 1000,
        acf_ci = 0.95,
        theme = ggplot2::theme_bw(),
        verbose = TRUE
        ) {

    if (!inherits(fit, "brmsfit") ) {

        stop ("`fit` must be a brmsfit object.", call. = FALSE)

    }

    if (!ggplot2::is_theme(theme) ) {

        stop ("Argument 'theme' should be a 'theme' object.", call. = FALSE)

    }

    if (!is.numeric(max_lag) || length(max_lag) != 1 || max_lag < 1) {

        stop ("`max_lag` must be a positive integer.", call. = FALSE)

    }

    if (!is.numeric(n_series_plot) || length(n_series_plot) != 1 || n_series_plot < 1) {

        stop ("`n_series_plot` must be a positive integer.", call. = FALSE)

    }

    if (!is.numeric(acf_ndraws) || length(acf_ndraws) != 1 || acf_ndraws < 10) {

        stop ("`acf_ndraws` must be a numeric scalar >= 10.", call. = FALSE)

    }

    if (!is.numeric(acf_ci) || length(acf_ci) != 1 || acf_ci <= 0 || acf_ci >= 1) {

        stop ("`acf_ci` must be a numeric scalar in (0, 1).", call. = FALSE)

    }

    resid_type <- match.arg(resid_type)
    ar_summary <- match.arg(ar_summary)

    if (is.null(data) ) {

        data <- fit$data

    }

    if (!is.data.frame(data) ) {

        stop ("`data` must be a data.frame.", call. = FALSE)

    }

    checks <- list()

    #############################
    # ensure time column exists #
    #############################

    if (!time_id %in% names(data) ) {

        stop ("`data` must contain the time column '", time_id, "'.", call. = FALSE)

    }

    ##############################
    # build/ensure series column #
    ##############################

    df <- data

    if (series_id %in% names(df) ) {

        df[[".series"]] <- df[[series_id]]
        checks$series_source <- series_id

    } else {

        has_participant <- participant_id %in% names(df)
        has_predictor <- predictor_id %in% names(df)

        if (has_participant && has_predictor) {

            df[[".series"]] <- interaction(df[[participant_id]], df[[predictor_id]], drop = TRUE)
            checks$series_source <- paste0("interaction(", participant_id, ", ", predictor_id, ")")

            if (isTRUE(verbose) ) {

                message ("`", series_id, "` not found; using participant x predictor as series grouping.")

            }

        } else if (has_participant) {

            df[[".series"]] <- df[[participant_id]]
            checks$series_source <- participant_id

            if (isTRUE(verbose) ) {

                message ("`", series_id, "` not found; using participant as series grouping.")

            }

        } else if (has_predictor) {

            df[[".series"]] <- df[[predictor_id]]
            checks$series_source <- predictor_id

            if (isTRUE(verbose) ) {

                message ("`", series_id, "` and `", participant_id, "` not found; using predictor as series grouping.")

            }

        } else {

            df[[".series"]] <- factor("group")
            checks$series_source <- "single-series"

            if (isTRUE(verbose) ) {

                message (
                    "No `", series_id, "`, `", participant_id, "`, or `", predictor_id,
                    "` found; treating all rows as a single time series ('group')."
                    )

            }

        }

    }

    #######################################
    # time uniqueness check within series #
    #######################################

    dup_df <- df |>
        dplyr::count(.data$.series, .data[[time_id]], name = "n") |>
        dplyr::filter(.data$n > 1L)

    checks$time_unique_within_series <- (nrow(dup_df) == 0L)
    checks$time_duplicates <- dup_df

    if (!checks$time_unique_within_series) {

        msg <- paste0(
            "Time points within series are not unique (n = ", nrow(dup_df), " duplicated series x time cells). ",
            "This violates brms::ar() requirements if you use the same grouping."
            )

        if (isTRUE(verbose) ) warning (msg, call. = FALSE)

    }

    #####################################
    # residuals via residuals.brmsfit() #
    #####################################

    res_draws <- stats::residuals(
        object = fit,
        method = resid_method,
        resid_type = resid_type,
        ndraws = ndraws,
        summary = FALSE
        )

    if (!is.matrix(res_draws) ) {

        stop ("Expected residuals.brmsfit(..., summary = FALSE) to return a matrix.", call. = FALSE)

    }

    if (ncol(res_draws) != nrow(df) ) {

        stop (
            "Row mismatch: residuals(fit) has ", ncol(res_draws),
            " observations but `data` has ", nrow(df), " rows. ",
            "Make sure `data` is the exact data used to fit `fit` (same rows and order).",
            call. = FALSE
            )

    }

    if (isTRUE(use_posterior_mean) ) {

        df[[".resid_raw"]] <- colMeans(res_draws)
        checks$resid_raw_summary <- "mean"

    } else {

        df[[".resid_raw"]] <- apply(res_draws, 2, stats::median)
        checks$resid_raw_summary <- "median"

    }

    checks$resid_method <- resid_method
    checks$resid_type <- resid_type
    checks$ndraws_residuals <- nrow(res_draws)

    ##################################
    # extract AR parameters (if any) #
    ##################################

    draws_df <- brms::as_draws_df(fit)
    draws_df <- as.data.frame(draws_df)
    ar_names <- names(draws_df)[grepl("^ar\\[", names(draws_df) )]

    if (length(ar_names) == 0L) {

        checks$has_ar_parameters <- FALSE
        ar_summary_df <- data.frame()

        if (isTRUE(verbose) ) {

            message ("No parameters named like 'ar[.]' found in the fitted model draws; cannot compute corrected residuals.")

        }

    } else {

        checks$has_ar_parameters <- TRUE

        lag_index <- as.integer(gsub("^ar\\[|\\]$", "", ar_names) )
        o <- order(lag_index)
        ar_names <- ar_names[o]
        lag_index <- lag_index[o]

        ar_summary_df <- do.call(
            rbind,
            lapply(ar_names, function (nm) {
                x <- draws_df[[nm]]
                data.frame(
                    parameter = nm,
                    mean = mean(x),
                    sd = stats::sd(x),
                    q025 = stats::quantile(x, 0.025),
                    q50  = stats::quantile(x, 0.5),
                    q975 = stats::quantile(x, 0.975),
                    row.names = NULL
                    )
            })
        )

    }

    ##############################################
    # compute corrected residuals for one series #
    ##############################################

    compute_corrected_1d <- function (r, phi) {

        p <- length(phi)
        out <- rep(NA_real_, length(r) )

        if (length(r) <= p) {

            return (out)

        }

        for (t in seq_len(length(r) ) ) {

            if (t <= p) next

            acc <- 0

            for (k in seq_len(p) ) {

                acc <- acc + phi[k] * r[t - k]

            }

            out[t] <- r[t] - acc

        }

        return (out)

    }

    ##########################################################
    # corrected residuals (point summary, for rho1/timeplot) #
    ##########################################################

    if (isTRUE(checks$has_ar_parameters) ) {

        phi_hat <- vapply(
            ar_names,
            function (nm) {
                x <- draws_df[[nm]]
                if (ar_summary == "mean") mean(x) else stats::median(x)
            },
            numeric(1)
            )

        df[[".resid_corrected"]] <- NA_real_

        df <- df |>
            dplyr::arrange(.data$.series, .data[[time_id]]) |>
            dplyr::group_by(.data$.series) |>
            dplyr::mutate(
                .resid_corrected = compute_corrected_1d(.data[[".resid_raw"]], phi_hat)
                ) |>
            dplyr::ungroup()

        checks$ar_phi_summary <- ar_summary

    } else {

        df[[".resid_corrected"]] <- NA_real_

    }

    #############################################
    # lag-1 residual autocorrelation per series #
    #############################################

    rho1_by_series_raw <- df |>
        dplyr::arrange(.data$.series, .data[[time_id]]) |>
        dplyr::group_by(.data$.series) |>
        dplyr::summarise(
            n_time = dplyr::n(),
            rho1 = {
                r <- .data[[".resid_raw"]]
                r <- r[is.finite(r)]
                if (length(r) < 3L) NA_real_
                else stats::cor(x = r[-1], y = r[-length(r)], use = "complete.obs")
            },
            .groups = "drop"
            )

    rho1_by_series_corrected <- df |>
        dplyr::arrange(.data$.series, .data[[time_id]]) |>
        dplyr::group_by(.data$.series) |>
        dplyr::summarise(
            n_time = dplyr::n(),
            rho1 = {
                r <- .data[[".resid_corrected"]]
                r <- r[is.finite(r)]
                if (length(r) < 3L) NA_real_
                else stats::cor(x = r[-1], y = r[-length(r)], use = "complete.obs")
            },
            .groups = "drop"
            )

    ##############################
    # choose series to visualise #
    ##############################

    series_levels <- unique(df[[".series"]])

    if (!is.null(seed) ) set.seed(seed)

    n_pick <- min(n_series_plot, length(series_levels) )
    pick <- if (n_pick > 0) sample(series_levels, n_pick) else character(0)

    ###########################################################
    # compute ACF + PACF point and bands from posterior draws #
    ###########################################################

    compute_acf_pacf_from_draws <- function (
            df,
            res_draws,
            phi_draws = NULL,
            compute_corrected = FALSE,
            pick,
            max_lag,
            time_id,
            acf_ndraws,
            acf_ci
            ) {

        if (length(pick) == 0) {

            return (list(acf = data.frame(), pacf = data.frame()) )

        }

        if (isTRUE(compute_corrected) ) {

            if (is.null(phi_draws) ) {

                return (list(acf = data.frame(), pacf = data.frame()) )

            }

            n_min <- min(nrow(res_draws), nrow(phi_draws) )

            if (nrow(res_draws) != nrow(phi_draws) && isTRUE(verbose) ) {

                message (
                    "Draw mismatch between residual draws (", nrow(res_draws),
                    ") and AR draws (", nrow(phi_draws),
                    "). Using the first ", n_min, " draws."
                    )

            }

            res_draws <- res_draws[seq_len(n_min), , drop = FALSE]
            phi_draws <- phi_draws[seq_len(n_min), , drop = FALSE]

        }

        n_post <- nrow(res_draws)
        draw_ids <- sample(seq_len(n_post), min(acf_ndraws, n_post) )

        q_lo <- (1 - acf_ci) / 2
        q_hi <- 1 - q_lo

        acf_out <- list()
        pacf_out <- list()

        for (sid in pick) {

            sub <- df[df[[".series"]] == sid, , drop = FALSE]
            ord <- order(sub[[time_id]])
            idx <- which(df[[".series"]] == sid)[ord]

            if (length(idx) < 3L) next

            rr <- res_draws[draw_ids, idx, drop = FALSE]

            if (isTRUE(compute_corrected) ) {

                phi_sub <- phi_draws[draw_ids, , drop = FALSE]

                corrected <- matrix(NA_real_, nrow = nrow(rr), ncol = ncol(rr))

                for (i in seq_len(nrow(rr) ) ) {

                    corrected[i, ] <- compute_corrected_1d(rr[i, ], phi_sub[i, ])

                }

                rr <- corrected

            }

            acf_mat <- matrix(NA_real_, nrow = nrow(rr), ncol = max_lag)
            pacf_mat <- matrix(NA_real_, nrow = nrow(rr), ncol = max_lag)

            for (i in seq_len(nrow(rr) ) ) {

                r_i <- as.numeric(rr[i, ])
                r_i <- r_i[is.finite(r_i)]

                if (length(r_i) < 3L) next

                ac <- stats::acf(r_i, plot = FALSE, lag.max = max_lag)
                acf_vec <- as.numeric(ac$acf)[-1]

                if (length(acf_vec) >= max_lag) {

                    acf_mat[i, ] <- acf_vec[seq_len(max_lag)]

                }

                pc <- stats::pacf(r_i, plot = FALSE, lag.max = max_lag)
                pacf_vec <- as.numeric(pc$acf)

                if (length(pacf_vec) >= max_lag) {

                    pacf_mat[i, ] <- pacf_vec[seq_len(max_lag)]

                }

            }

            lag_vec <- seq_len(max_lag)

            acf_out[[as.character(sid)]] <- data.frame(
                .series = sid,
                lag = lag_vec,
                q_low = apply(acf_mat, 2, stats::quantile, probs = q_lo, na.rm = TRUE),
                q50 = apply(acf_mat, 2, stats::quantile, probs = 0.5, na.rm = TRUE),
                q_high = apply(acf_mat, 2, stats::quantile, probs = q_hi, na.rm = TRUE),
                # point = apply(acf_mat, 2, stats::quantile, probs = 0.5, na.rm = TRUE),
                stringsAsFactors = FALSE
                )

            pacf_out[[as.character(sid)]] <- data.frame(
                .series = sid,
                lag = lag_vec,
                q_low = apply(pacf_mat, 2, stats::quantile, probs = q_lo, na.rm = TRUE),
                q50 = apply(pacf_mat, 2, stats::quantile, probs = 0.5, na.rm = TRUE),
                q_high = apply(pacf_mat, 2, stats::quantile, probs = q_hi, na.rm = TRUE),
                # point = apply(pacf_mat, 2, stats::quantile, probs = 0.5, na.rm = TRUE),
                stringsAsFactors = FALSE
                )

        }

        acf_df <- if (length(acf_out) > 0) do.call(rbind, acf_out) else data.frame()
        pacf_df <- if (length(pacf_out) > 0) do.call(rbind, pacf_out) else data.frame()

        return (list(acf = acf_df, pacf = pacf_df) )

    }

    ##############################
    # compute ACF/PACF dataframes #
    ##############################

    acf_raw <- data.frame()
    acf_corrected <- data.frame()
    pacf_raw <- data.frame()
    pacf_corrected <- data.frame()

    raw_res <- compute_acf_pacf_from_draws(
        df = df,
        res_draws = res_draws,
        phi_draws = NULL,
        compute_corrected = FALSE,
        pick = pick,
        max_lag = max_lag,
        time_id = time_id,
        acf_ndraws = acf_ndraws,
        acf_ci = acf_ci
        )

    acf_raw <- raw_res$acf
    pacf_raw <- raw_res$pacf

    if (isTRUE(checks$has_ar_parameters) ) {

        phi_draws <- as.matrix(draws_df[, ar_names, drop = FALSE])

        corr_res <- compute_acf_pacf_from_draws(
            df = df,
            res_draws = res_draws,
            phi_draws = phi_draws,
            compute_corrected = TRUE,
            pick = pick,
            max_lag = max_lag,
            time_id = time_id,
            acf_ndraws = acf_ndraws,
            acf_ci = acf_ci
            )

        acf_corrected <- corr_res$acf
        pacf_corrected <- corr_res$pacf

    }

    #########
    # plots #
    #########

    plots <- list()

    # AR parameter density (if present)
    if (nrow(ar_summary_df) > 0) {

        ar_long <- do.call(
            rbind,
            lapply(ar_names, function (nm) {
                data.frame(parameter = nm, value = draws_df[[nm]])
            })
        )

        plots$ar_density <- ar_long |>
            ggplot2::ggplot(ggplot2::aes(x = .data$value) ) +
            ggplot2::geom_density() +
            ggplot2::facet_wrap(~parameter, scales = "free") +
            theme +
            ggplot2::labs(
                title = "Posterior density of AR parameter(s)",
                x = "AR coefficient",
                y = "Density"
                )

    }

    # rho1 by series (raw)
    plots$rho1_series_raw <- rho1_by_series_raw |>
        ggplot2::ggplot(ggplot2::aes(x = .data$.series, y = .data$rho1) ) +
        ggplot2::geom_segment(ggplot2::aes(y = 0, yend = .data$rho1) ) +
        ggplot2::geom_point() +
        ggplot2::coord_flip() +
        theme +
        ggplot2::labs(
            title = "Lag-1 residual autocorrelation per series (raw residuals)",
            x = "Series",
            y = "rho1"
            )

    if (isTRUE(checks$has_ar_parameters) ) {

        plots$rho1_series_corrected <- rho1_by_series_corrected |>
            ggplot2::ggplot(ggplot2::aes(x = .data$.series, y = .data$rho1) ) +
            ggplot2::geom_segment(ggplot2::aes(y = 0, yend = .data$rho1) ) +
            ggplot2::geom_point() +
            ggplot2::coord_flip() +
            theme +
            ggplot2::labs(
                title = "Lag-1 residual autocorrelation per series (corrected residuals)",
                x = "Series",
                y = "rho1"
                )

    }

    # ACF plot: raw (bands + point)
    if (nrow(acf_raw) > 0) {

        plots$acf_raw <- acf_raw |>
            ggplot2::ggplot(ggplot2::aes(x = .data$lag) ) +
            ggplot2::geom_hline(yintercept = 0, linetype = 2) +
            ggplot2::geom_segment(
                colour = "steelblue",
                linewidth = 2,
                alpha = 0.5,
                ggplot2::aes(y = .data$q_low, yend = .data$q_high, xend = .data$lag)
                ) +
            ggplot2::geom_segment(
                # ggplot2::aes(y = .data$point, yend = 0, xend = .data$lag)
                ggplot2::aes(y = .data$q50, yend = 0, xend = .data$lag)
                ) +
            ggplot2::geom_point(
                # ggplot2::aes(y = .data$point),
                ggplot2::aes(y = .data$q50),
                size = 1
                ) +
            ggplot2::facet_wrap(~.data$.series) +
            theme +
            ggplot2::labs(
                title = paste0("Raw residuals ACF (method = '", resid_method, "', ci_width = ", acf_ci, ")"),
                x = "Lag",
                y = "ACF"
                )

    }

    # ACF plot: corrected (bands + point)
    if (nrow(acf_corrected) > 0) {

        plots$acf_corrected <- acf_corrected |>
            ggplot2::ggplot(ggplot2::aes(x = .data$lag) ) +
            ggplot2::geom_hline(yintercept = 0, linetype = 2) +
            ggplot2::geom_segment(
                colour = "steelblue",
                linewidth = 2,
                alpha = 0.5,
                ggplot2::aes(y = .data$q_low, yend = .data$q_high, xend = .data$lag)
                ) +
            ggplot2::geom_segment(
                # ggplot2::aes(y = .data$point, yend = 0, xend = .data$lag)
                ggplot2::aes(y = .data$q50, yend = 0, xend = .data$lag)
                ) +
            ggplot2::geom_point(
                # ggplot2::aes(y = .data$point),
                ggplot2::aes(y = .data$q50),
                size = 1
                ) +
            ggplot2::facet_wrap(~.data$.series) +
            theme +
            ggplot2::labs(
                title = paste0("Corrected residuals ACF (method = '", resid_method, "', ci_width = ", acf_ci, ")"),
                x = "Lag",
                y = "ACF"
                )

    }

    # PACF plot: raw (bands + point)
    if (nrow(pacf_raw) > 0) {

        plots$pacf_raw <- pacf_raw |>
            ggplot2::ggplot(ggplot2::aes(x = .data$lag) ) +
            ggplot2::geom_hline(yintercept = 0, linetype = 2) +
            ggplot2::geom_segment(
                colour = "steelblue",
                linewidth = 2,
                alpha = 0.5,
                ggplot2::aes(y = .data$q_low, yend = .data$q_high, xend = .data$lag)
                ) +
            ggplot2::geom_segment(
                # ggplot2::aes(y = .data$point, yend = 0, xend = .data$lag)
                ggplot2::aes(y = .data$q50, yend = 0, xend = .data$lag)
                ) +
            ggplot2::geom_point(
                # ggplot2::aes(y = .data$point),
                ggplot2::aes(y = .data$q50),
                size = 1
                ) +
            ggplot2::facet_wrap(~.data$.series) +
            theme +
            ggplot2::labs(
                title = paste0("Raw residuals PACF (method = '", resid_method, "', ci_width = ", acf_ci, ")"),
                x = "Lag",
                y = "PACF"
                )

    }

    # PACF plot: corrected (bands + point)
    if (nrow(pacf_corrected) > 0) {

        plots$pacf_corrected <- pacf_corrected |>
            ggplot2::ggplot(ggplot2::aes(x = .data$lag) ) +
            ggplot2::geom_hline(yintercept = 0, linetype = 2) +
            ggplot2::geom_segment(
                colour = "steelblue",
                linewidth = 2,
                alpha = 0.5,
                ggplot2::aes(y = .data$q_low, yend = .data$q_high, xend = .data$lag)
                ) +
            ggplot2::geom_segment(
                # ggplot2::aes(y = .data$point, yend = 0, xend = .data$lag)
                ggplot2::aes(y = .data$q50, yend = 0, xend = .data$lag)
                ) +
            ggplot2::geom_point(
                # ggplot2::aes(y = .data$point),
                ggplot2::aes(y = .data$q50),
                size = 1
                ) +
            ggplot2::facet_wrap(~.data$.series) +
            theme +
            ggplot2::labs(
                title = paste0("Corrected residuals PACF (method = '", resid_method, "', ci_width = ", acf_ci, ")"),
                x = "Lag",
                y = "PACF"
                )

    }

    # raw residuals vs time for sampled series
    if (length(pick) > 0) {

        subp <- df[df[[".series"]] %in% pick, , drop = FALSE]

        plots$resid_time_raw <- subp |>
            ggplot2::ggplot(ggplot2::aes(x = .data[[time_id]], y = .data[[".resid_raw"]]) ) +
            ggplot2::geom_hline(yintercept = 0, linetype = 2) +
            ggplot2::geom_line() +
            ggplot2::facet_wrap(~.data$.series, scales = "free_y") +
            theme +
            ggplot2::labs(
                title = "Raw residuals through time (sampled series)",
                x = time_id,
                y = "Residual"
                )

    }

    # corrected residuals vs time for sampled series
    if (length(pick) > 0 && isTRUE(checks$has_ar_parameters) ) {

        subp <- df[df[[".series"]] %in% pick, , drop = FALSE]

        plots$resid_time_corrected <- subp |>
            ggplot2::ggplot(ggplot2::aes(x = .data[[time_id]], y = .data[[".resid_corrected"]]) ) +
            ggplot2::geom_hline(yintercept = 0, linetype = 2) +
            ggplot2::geom_line() +
            ggplot2::facet_wrap(~.data$.series, scales = "free_y") +
            theme +
            ggplot2::labs(
                title = "Corrected residuals through time (sampled series)",
                x = time_id,
                y = "Residual"
                )

    }

    out <- list(
        ar_summary = ar_summary_df,
        rho1_by_series_raw = rho1_by_series_raw,
        rho1_by_series_corrected = rho1_by_series_corrected,
        acf_raw = acf_raw,
        acf_corrected = acf_corrected,
        pacf_raw = pacf_raw,
        pacf_corrected = pacf_corrected,
        plots = plots
        )

    class(out) <- c("neurogam_autocor_check", class(out) )

    return (out)

}
