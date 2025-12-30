#' Check and visualise residual autocorrelation in a fitted \pkg{neurogam} model
#'
#' Computes practical diagnostics for residual autocorrelation in a fitted
#' \pkg{brms} model, with special handling of time-resolved neurogam data where
#' each participant can have multiple independent time series (e.g., one per
#' condition).
#'
#' The function:
#' \itemize{
#'   \item checks that time points are unique within each AR series,
#'   \item extracts and summarises AR parameter(s) (e.g., \code{ar[1]}) if present,
#'   \item computes posterior-mean residuals \eqn{y - E[y \mid x]},
#'   \item quantifies lag-1 residual autocorrelation per series,
#'   \item visualises residual autocorrelation via ACF curves and summary plots.
#' }
#'
#' @param fit A fitted \code{\link[brms:brm]{brmsfit}} object.
#' @param data A data frame used to fit the model (or the same rows in the same order).
#'   Must contain a time column and (directly or indirectly) an AR series identifier.
#'   If \code{NULL} (default), uses data from \code{fit}.
#' @param time_id Character; name of the time column. Defaults to \code{"time"}.
#' @param series_id Character; name of the column defining independent time series
#'   for autocorrelation checks. Defaults to \code{"ar_series"}.
#'   If \code{series_id} is not present in \code{data}, it is created as
#'   \code{interaction(participant_id, predictor_id, drop = TRUE)} when possible,
#'   otherwise as \code{participant_id}.
#' @param participant_id Character; participant column name used to build \code{series}
#'   if needed. Defaults to \code{"participant"}.
#' @param predictor_id Character; predictor/condition column name used to build \code{series}
#'   if needed. Defaults to \code{"predictor"}.
#' @param outcome_mean_id Character; outcome column for Gaussian-type residuals.
#'   Defaults to \code{"outcome_mean"}.
#' @param success_id,trials_id Character; columns for binomial outcomes
#'   (\code{success} and \code{trials}). Defaults to \code{"success"} and \code{"trials"}.
#' @param max_lag Integer; maximum lag for ACF curves. Defaults to 25.
#' @param n_series_plot Integer; number of series to plot ACFs for (sampled).
#'   Defaults to 9.
#' @param seed Integer or \code{NULL}; random seed for sampling series to plot.
#'   Defaults to 42.
#' @param use_posterior_mean Logical; if \code{TRUE} (default) uses posterior mean
#'   \code{E[y|x]} to compute residuals. If \code{FALSE}, uses median.
#' @param theme A \code{\link[ggplot2:theme]{theme}} object
#'   modifying the appearance of the plots.
#' @param verbose Logical; if \code{TRUE}, emits informative messages and warnings.
#'
#' @return A named list with components:
#' \describe{
#'   \item{data}{A copy of \code{data} with \code{.series} and \code{.resid} columns added.}
#'   \item{checks}{A list of checks (duplicates, missing columns, etc.).}
#'   \item{ar_summary}{A data frame summarising AR parameter draws (if present).}
#'   \item{rho1_by_series}{A data frame with lag-1 residual correlation per series.}
#'   \item{acf_df}{A data frame of ACF values for sampled series (for plotting).}
#'   \item{plots}{A list of \pkg{ggplot2} objects.}
#' }
#'
#' @examples
#' \dontrun{
#' # after fitting a brms model
#' out <- check_residual_autocorrelation(fit)
#' out$ar_summary
#' out$plots$rho1_hist
#' out$plots$acf
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
        outcome_mean_id = "outcome_mean",
        success_id = "success",
        trials_id = "trials",
        max_lag = 25,
        n_series_plot = 9,
        seed = 42,
        use_posterior_mean = TRUE,
        theme = ggplot2::theme_bw(),
        verbose = TRUE
        ) {

    if (!inherits(fit, "brmsfit") ) {

        stop ("`fit` must be a brmsfit object.", call. = FALSE)

    }

    if (!ggplot2::is.theme(theme) ) {

        stop ("Argument 'theme' should be a 'theme' object.")

    }

    if (is.null(data) ) {

        data <- fit$data

    }

    if (!is.data.frame(data) ) {

        stop ("`data` must be a data.frame.", call. = FALSE)

    }

    checks <- list()

    # ensure time column exists
    if (!time_id %in% names(data) ) {

        stop ("`data` must contain the time column '", time_id, "'.", call. = FALSE)

    }

    # build/ensure series column
    df <- data

    if (series_id %in% names(df) ) {

        df[[".series"]] <- df[[series_id]]
        checks$series_source <- series_id

    } else {

        can_build_series <- participant_id %in% names(df)
        has_predictor <- predictor_id %in% names(df)

        if (!can_build_series) {

            stop (
                "Could not find `series_id` ('", series_id, "') and cannot build it because `participant_id` ('",
                participant_id, "') is missing in `data`.",
                call. = FALSE
                )

        }

        if (has_predictor) {

            df[[".series"]] <- interaction(df[[participant_id]], df[[predictor_id]], drop = TRUE)
            checks$series_source <- paste0("interaction(", participant_id, ", ", predictor_id, ")")

        } else {

            df[[".series"]] <- df[[participant_id]]
            checks$series_source <- participant_id

            if (isTRUE(verbose) ) {

                message ("`", series_id, "` not found; using participant as series grouping for autocorrelation checks.")

            }

        }

    }

    # time uniqueness check within series
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

        if (isTRUE(verbose) ) warning(msg, call. = FALSE)

    }

    # infer outcome type (gaussian vs binomial) by presence of columns (pragmatic)
    has_binom_cols <- all(c(success_id, trials_id) %in% names(df) )
    has_gauss_col <- outcome_mean_id %in% names(df)

    if (!has_binom_cols && !has_gauss_col) {

        stop (
            "Could not find outcome columns. Provide either '", outcome_mean_id,
            "' for Gaussian or both '", success_id, "' and '", trials_id, "' for Binomial.",
            call. = FALSE
            )

    }

    outcome_kind <- if (has_binom_cols) "binomial" else "gaussian"
    checks$outcome_kind <- outcome_kind

    # expected value E[y|x] from brms
    epred <- brms::posterior_epred(fit)

    if (ncol(epred) != nrow(df) ) {

        stop (
            "Row mismatch: posterior_epred(fit) has ", ncol(epred),
            " observations but `data` has ", nrow(df), " rows. ",
            "Make sure `data` is the exact data used to fit `fit` (same rows and order).",
            call. = FALSE
            )

    }

    if (isTRUE(use_posterior_mean) ) {

        mu_hat <- colMeans(epred)
        checks$epred_summary <- "mean"

    } else {

        mu_hat <- apply(epred, 2, stats::median)
        checks$epred_summary <- "median"

    }

    # compute residuals
    if (outcome_kind == "gaussian") {

        y <- df[[outcome_mean_id]]
        df[[".resid"]] <- y - mu_hat

    } else {

        # residuals on the response scale for proportions (success/trials)
        y <- df[[success_id]] / df[[trials_id]]
        df[[".resid"]] <- y - mu_hat

    }

    # lag-1 residual autocorrelation per series
    rho1_by_series <- df |>
        dplyr::arrange(.data$.series, .data[[time_id]]) |>
        dplyr::group_by(.data$.series) |>
        dplyr::summarise(
            n_time = dplyr::n(),
            rho1 = {
                r <- .data[[".resid"]]
                if (length(r) < 3L) NA_real_ else stats::cor(x = r[-1], y = r[-length(r)], use = "complete.obs")
            },
            .groups = "drop"
            )

    checks$rho1_n_na <- sum(is.na(rho1_by_series$rho1) )

    # extract AR parameter(s), if any
    draws <- brms::as_draws_df(fit)
    ar_names <- names(draws)[grepl("^ar\\[", names(draws) )]

    if (length(ar_names) == 0L) {

        ar_summary <- data.frame()
        checks$has_ar_parameters <- FALSE

        if (isTRUE(verbose) ) {

            message ("No parameters named like 'ar[.]' found in the fitted model draws.")

        }

    } else {

        checks$has_ar_parameters <- TRUE
        ar_summary <- do.call(
            rbind,
            lapply(ar_names, function (nm) {
                x <- draws[[nm]]
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

    # build ACF df for sampled series (posterior-mean residuals)
    series_levels <- unique(df[[".series"]])

    if (!is.null(seed) ) set.seed(seed)

    n_pick <- min(n_series_plot, length(series_levels) )
    pick <- if (n_pick > 0) sample(series_levels, n_pick) else character(0)

    acf_df <- data.frame()

    if (length(pick) > 0) {

        acf_df <- do.call(
            rbind,
            lapply(pick, function (sid) {
                sub <- df[df[[".series"]] == sid, , drop = FALSE]
                sub <- sub[order(sub[[time_id]]), , drop = FALSE]
                r <- sub[[".resid"]]
                # stats::acf returns lag 0..max_lag; keep 1..max_lag for plotting
                ac <- stats::acf(r, plot = FALSE, lag.max = max_lag, na.action = stats::na.pass)
                data.frame(
                    .series = sid,
                    lag = as.numeric(ac$lag),
                    acf = as.numeric(ac$acf),
                    stringsAsFactors = FALSE
                    )
            })
        )

        acf_df <- acf_df[acf_df$lag > 0, , drop = FALSE]

    }

    # plots
    plots <- list()

    # AR parameter density (if present)
    if (nrow(ar_summary) > 0) {

        ar_long <- do.call(
            rbind,
            lapply(ar_names, function (nm) {
                data.frame(parameter = nm, value = draws[[nm]])
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

    # rho1 histogram
    plots$rho1_hist <- rho1_by_series |>
        ggplot2::ggplot(ggplot2::aes(x = .data$rho1) ) +
        ggplot2::geom_histogram(bins = 30) +
        theme +
        ggplot2::labs(
            title = "Lag-1 residual autocorrelation by series",
            x = "rho1 = cor(resid[t], resid[t-1])",
            y = "Count"
            )

    # rho1 by series (dot plot)
    plots$rho1_series <- rho1_by_series |>
        ggplot2::ggplot(ggplot2::aes(x = .data$.series, y = .data$rho1) ) +
        ggplot2::geom_point() +
        ggplot2::coord_flip() +
        theme +
        ggplot2::labs(
            title = "Lag-1 residual autocorrelation per series",
            x = "Series",
            y = "rho1"
            )

    # ACF plot for sampled series
    if (nrow(acf_df) > 0) {

        plots$acf <- acf_df |>
            ggplot2::ggplot(ggplot2::aes(x = .data$lag, y = .data$acf) ) +
            ggplot2::geom_hline(yintercept = 0) +
            ggplot2::geom_segment(ggplot2::aes(xend = .data$lag, yend = 0) ) +
            ggplot2::facet_wrap(~.data$.series) +
            theme +
            ggplot2::labs(
                title = paste0("Residual ACF (posterior-", checks$epred_summary, " residuals)"),
                x = "Lag",
                y = "ACF"
                )

    }

    # helpful check plot: residuals vs time for sampled series
    if (length(pick) > 0) {

        subp <- df[df[[".series"]] %in% pick, , drop = FALSE]

        plots$resid_time <- subp |>
            ggplot2::ggplot(ggplot2::aes(x = .data[[time_id]], y = .data[[".resid"]]) ) +
            ggplot2::geom_line() +
            ggplot2::facet_wrap(~.data$.series, scales = "free_y") +
            theme +
            ggplot2::labs(
                title = "Residuals through time (sampled series)",
                x = time_id,
                y = "Residual"
                )

    }

    out <- list(
        data = df,
        checks = checks,
        ar_summary = ar_summary,
        rho1_by_series = rho1_by_series,
        acf_df = acf_df,
        plots = plots
        )

    class(out) <- c("neurogam_autocor_check", class(out) )

    return (out)

}
