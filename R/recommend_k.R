#' Recommend a smooth basis dimension k via effective complexity "knee" detection
#'
#' This function takes a \code{"clusters_results"} object produced by
#' \code{\link{testing_through_time}}, refits the underlying \pkg{brms} model
#' for a grid of \code{k} values, computes model comparison criteria
#' (\code{loo}, \code{waic}), and identifies a recommended
#' smooth basis dimension via an automatic "knee" detection procedure on the
#' effective number of parameters (either \code{p_loo} or \code{p_waic}).
#'
#' @param object An object of class \code{"clusters_results"} as returned by
#'   \code{\link{testing_through_time}}. The object must contain a fitted
#'   \pkg{brms} model in its \code{$model} slot.
#' @param k_min Numeric; minimum value of the smooth basis dimension
#'   \code{k} to consider.
#' @param k_max Numeric; maximum value of the smooth basis dimension
#'   \code{k} to consider.
#' @param k_step Numeric; step size between successive \code{k} values.
#'   The sequence of tested values is \code{seq(k_min, k_max, by = k_step)}.
#' @param criterion Character vector of model comparison criteria to add via
#'   \code{\link[brms]{add_criterion}}. Defaults to
#'   \code{c("waic", "loo")}.
#' @param knee_method Character; method for knee detection. One of:
#'   \itemize{
#'     \item \code{"geometric"}: standard geometric elbow on the raw
#'           \code{p_*} values;
#'     \item \code{"geometric_smooth"}: geometric elbow on a loess-smoothed
#'           curve of \code{p_*} vs \code{k}, with weights
#'           \code{1 / SE(p_*)^2} when available (recommended when
#'           \code{p_*} estimates are noisy).
#'   }
#' @param loess_span Numeric; smoothing parameter passed to
#'   \code{\link[stats]{loess}} when \code{knee_method = "geometric_smooth"}.
#'   Defaults to \code{0.75}.
#' @param verbose Logical; if \code{TRUE} (default), prints progress messages
#'   while refitting models for different \code{k} values.
#'
#' @details
#' For each \code{k} in \code{seq(k_min, k_max, by = k_step)}, the function:
#' \enumerate{
#'   \item updates all \code{s(..., k = ...)} terms in the original \pkg{brms}
#'     formula to use the new value \code{k}, and refits the model via
#'     \code{update()};
#'   \item calls \code{\link[brms]{add_criterion}} to compute the requested
#'     criteria, and extracts \code{p_loo} or \code{p_waic} and their standard
#'     errors (where available);
#'   \item builds a comparison table of criteria values across \code{k}.
#' }
#'
#' To recommend a basis dimension, the function treats the chosen complexity
#' measure (\code{p_loo} or \code{p_waic}) as a function of \code{k} and uses
#' an elbow (knee) method:
#' \itemize{
#'   \item \code{knee_method = "geometric"} applies the geometric knee
#'         procedure directly to \code{p_*} vs \code{k};
#'   \item \code{knee_method = "geometric_smooth"} first fits a loess curve
#'         \code{p_* ~ k} (with weights \code{1 / SE(p_*)^2} when available),
#'         then applies the geometric knee to the smoothed curve. This
#'         reduces the influence of noisy \code{p_*} estimates.
#' }
#'
#' The resulting plot shows \code{p_*} vs \code{k} with:
#' \itemize{
#'   \item points and a line for the raw \code{p_*} estimates;
#'   \item vertical error bars \code{± SE(p_*)} (when SEs are available);
#'   \item (optionally) a dashed line for the smoothed \code{p_*} curve;
#'   \item a vertical dashed line at the recommended \code{k} (knee).
#' }
#'
#' @return An object of class \code{"recommend_k_results"}, which is a list
#'   with elements:
#'   \itemize{
#'     \item \code{models}: named list of refitted \code{brmsfit} objects,
#'       one per \code{k} value;
#'     \item \code{comparison}: data frame summarising model criteria for each
#'       \code{k} (one row per \code{k}), including \code{p_*} and its SE;
#'     \item \code{k_values}: numeric vector of \code{k} values that were
#'       evaluated;
#'     \item \code{recommended_k}: the \code{k} value identified as the knee
#'       based on the chosen effective complexity measure and method;
#'     \item \code{plot}: a \code{ggplot2} object displaying \code{p_*} as a
#'       function of \code{k} with error bars and the knee highlighted;
#'     \item \code{knee_on}: the complexity measure used (\code{"p_loo"} or
#'       \code{"p_waic"});
#'     \item \code{knee_method}: the knee detection method used.
#'   }
#'
#' @examples
#' \dontrun{
#' # import some simulated EEG data
#' data(eeg_data)
#' head(eeg_data)
#'
#' # fit a time-resolved GAMM
#' res <- testing_through_time(
#'   data = eeg_data,
#'   participant_id = "participant",
#'   outcome_id = "eeg",
#'   time_id = "time",
#'   predictor_id = NA,
#'   kvalue = 20,
#'   multilevel = "summary"
#'   )
#'
#' # recommend an optimal smooth basis dimension k
#' k_res <- recommend_k(
#'   object = res,
#'   k_min = 10,
#'   k_max = 40,
#'   k_step = 5,
#'   criterion = "waic"
#'   )
#'
#' # print summary in the console
#' summary(k_res)
#'
#' # extract the recommended k value
#' k_res$recommended_k
#'
#' # access the comparison table
#' k_res$comparison
#'
#' # visualise effective complexity vs. k
#' k_res$plot
#' }
#'
#' @author Ladislas Nalborczyk \email{ladislas.nalborczyk@@cnrs.fr}.
#'
#' @seealso
#'   \code{\link{testing_through_time}},
#'   \code{\link[brms]{brm}},
#'   \code{\link[brms]{add_criterion}}
#'   \code{\link[stats]{loess}}
#'
#' @importFrom stats sd formula as.formula update loess predict
#'
#' @export
recommend_k <- function (
        object,
        k_min = 10,
        k_max = 40,
        k_step = 5,
        criterion = c("waic", "loo"),
        knee_method = c("geometric_smooth", "geometric"),
        loess_span = 0.5,
        verbose = TRUE
        ) {

    # basic checks
    if (!inherits(object, "clusters_results") ) {

        stop ("`object` must be of class 'clusters_results'.", call. = FALSE)

    }

    if (is.null(object$model) || !inherits(object$model, "brmsfit") ) {

        stop ("`object$model` must be a valid 'brmsfit' object.", call. = FALSE)

    }

    if (!is.numeric(k_min) || !is.numeric(k_max) || !is.numeric(k_step) ) {

        stop ("`k_min`, `k_max`, and `k_step` must be numeric.", call. = FALSE)

    }

    criterion <- match.arg(criterion)
    knee_method <- match.arg(knee_method)
    k_values <- seq(k_min, k_max, by = k_step)

    if (length(k_values) < 2L) {

        stop ("The sequence of k values must contain at least two points.", call. = FALSE)

    }

    # helper: update k in all s(...) terms of a brmsfit
    .update_s_k <- function (expr, new_k) {

        if (!is.call(expr) ) return(expr)

        # if this is an s(...) call, update/add k argument
        if (identical(expr[[1L]], as.name("s") ) ) {

            args <- as.list(expr)
            arg_names <- names(args)

            if ("k" %in% arg_names) {

                args[["k"]] <- new_k

            } else {

                args[["k"]] <- new_k

            }

            return (as.call(args) )

        }

        # otherwise, recurse on all components
        for (i in seq_along(expr) ) {

            expr[[i]] <- .update_s_k(expr[[i]], new_k)

        }

        return (expr)

    }

    .update_brms_k <- function (fit, new_k) {

        bf_old <- formula(fit)
        f_core <- bf_old$formula

        lhs <- f_core[[2L]]
        rhs <- f_core[[3L]]

        rhs_new <- .update_s_k(rhs, new_k)
        f_new  <- as.formula(call("~", lhs, rhs_new) )

        # update the brms model
        update(fit, formula. = f_new)

    }

    # storage objects
    model_list <- vector("list", length(k_values) )
    names(model_list) <- paste0("k_", k_values)
    comp_rows <- vector("list", length(k_values) )

    # loop over k values
    for (i in seq_along(k_values) ) {

        k_val <- k_values[i]

        if (verbose) {

            message ("Refitting model for k = ", k_val, " ...\n")

        }

        # refit with updated k
        fit_k <- .update_brms_k(object$model, new_k = k_val)

        # add criteria
        fit_k <- brms::add_criterion(
            x = fit_k,
            criterion = criterion,
            model_name = paste0("gam_k", k_val)
            )

        model_list[[i]] <- fit_k

        # build comparison row
        row <- list(
            k = k_val,
            model = paste0("gam_k", k_val)
            )

        # LOO
        if ("loo" %in% criterion && !is.null(fit_k$criteria$loo) ) {

            knee_on <- "p_loo"

            loo_obj <- fit_k$criteria$loo
            est <- loo_obj$estimates
            row$loo_elpd <- est["elpd_loo", "Estimate"]
            row$loo_elpd_se <- est["elpd_loo", "SE"]
            row$p_loo <- est["p_loo", "Estimate"]
            row$p_loo_see <- est["p_loo", "SE"]
            row$looic <- est["looic", "Estimate"]

        }

        # WAIC
        if ("waic" %in% criterion && !is.null(fit_k$criteria$waic) ) {

            knee_on <- "p_waic"

            waic_obj <- fit_k$criteria$waic
            est <- waic_obj$estimates
            row$waic_elpd <- est["elpd_waic", "Estimate"]
            row$waic_elpd_se <- est["elpd_waic", "SE"]
            row$p_waic <- est["p_waic", "Estimate"]
            row$p_waic_se <- est["p_waic", "SE"]
            row$waic <- est["waic", "Estimate"]

        }

        comp_rows[[i]] <- as.data.frame(row, stringsAsFactors = FALSE)

    }

    # assemble comparison table
    comparison <- do.call(rbind, comp_rows)
    rownames(comparison) <- NULL

    # choose which p_* column to use for knee detection
    if (!knee_on %in% names(comparison) ) {

        warning (
            "Column '", knee_on,
            "' is missing from the comparison table; knee detection ",
            "cannot be performed. `recommended_k` set to NA."
            )
        recommended_k <- NA_real_

    } else if (all(is.na(comparison[[knee_on]]) ) ) {

        warning (
            "All values of '", knee_on,
            "' are NA; knee detection cannot be performed. `recommended_k` ",
            "set to NA."
            )

        recommended_k <- NA_real_

    } else {

        # aata for knee detection (remove rows with NA in p_*):
        df_knee <- comparison[!is.na(comparison[[knee_on]]), , drop = FALSE]
        x  <- df_knee$k
        y  <- df_knee[[knee_on]]

        se_col <- paste0(knee_on, "_se")
        y_se <- if (se_col %in% names(df_knee) ) df_knee[[se_col]] else rep(NA_real_, length(y) )

        if (knee_method == "geometric_smooth") {

            # weights = 1 / SE^2 when SE available & > 0, otherwise 1
            w <- rep(1, length(y) )
            if (!all(is.na(y_se) ) ) {

                ok_se <- !is.na(y_se) & y_se > 0
                w[ok_se] <- 1 / (y_se[ok_se]^2)

            }

            lo_obj <- stats::loess(y ~ x, weights = w, span = loess_span)
            y_smooth <- stats::predict(lo_obj, newdata = data.frame(x = x) )

            # store smoothed curve back in comparison for plotting
            comparison[[paste0(knee_on, "_smooth")]] <- NA_real_
            comparison[[paste0(knee_on, "_smooth")]][match(x, comparison$k)] <- y_smooth

            y_for_knee <- y_smooth

        } else {

            # plain geometric knee on raw p_* values
            y_for_knee <- y

        }

        # rescale to [0, 1] for numerical stability
        x_scaled <- (x - min(x) ) / (max(x) - min(x) )
        y_scaled <- (y_for_knee - min(y_for_knee) ) / (max(y_for_knee) - min(y_for_knee) )

        # line between first and last point
        x1 <- x_scaled[1L]
        y1 <- y_scaled[1L]
        x2 <- x_scaled[length(x_scaled)]
        y2 <- y_scaled[length(y_scaled)]

        denom <- sqrt((y2 - y1)^2 + (x2 - x1)^2)

        if (denom == 0) {

            # degenerate case: all points identical in scaled space
            knee_idx <- which.max(y_for_knee)

        } else {

            distances <- abs(
                (y2 - y1) * x_scaled -
                    (x2 - x1) * y_scaled +
                    x2 * y1   - y2 * x1
                ) / denom

            knee_idx <- which.max(distances)

        }

        recommended_k <- x[knee_idx]

    }

    # plot: p_* vs k with SE and optional smooth
    if (knee_on %in% names(comparison) ) {

        p_col <- knee_on
        se_col <- paste0(knee_on, "_se")
        sm_col <- paste0(knee_on, "_smooth")
        has_se <- se_col %in% names(comparison)
        has_sm <- sm_col %in% names(comparison) && any(!is.na(comparison[[sm_col]]))

        p <- comparison |>
            ggplot2::ggplot(ggplot2::aes(x = .data$k, y = .data[[p_col]]) ) +
            ggplot2::geom_line() +
            ggplot2::geom_point(size = 2)

        if (has_se) {

            p <- p +
                ggplot2::geom_errorbar(
                    ggplot2::aes(
                        ymin = .data[[p_col]] - .data[[se_col]],
                        ymax = .data[[p_col]] + .data[[se_col]]
                        ),
                    width = diff(range(comparison$k, na.rm = TRUE) ) * 0.01
                    )

        }

        if (has_sm) {

            p <- p +
                ggplot2::geom_line(
                    ggplot2::aes(y = .data[[sm_col]]),
                    linetype = "dashed"
                    )

        }

        p <- p +
            ggplot2::geom_vline(
                xintercept = recommended_k,
                linetype = 2
                ) +
            ggplot2::theme_bw() +
            ggplot2::labs(
                x = "Smooth basis dimension (k)",
                y = paste0("Model complexity (", knee_on, ")"),
                title = "Tuning k based on effective number of parameters",
                subtitle = paste0(
                    "Recommended k (", knee_method, ", ", knee_on, "): ",
                    recommended_k
                    )
                )

        print(p)

    } else {

        p <- NULL

    }

    # build output object
    out <- list(
        models = model_list,
        comparison = comparison,
        k_values = k_values,
        recommended_k = recommended_k,
        plot = p,
        knee_on = knee_on,
        knee_method   = knee_method
        )

    class(out) <- "recommend_k_results"

    return (out)

}

#' Print method for recommend_k_results objects
#'
#' Pretty printer for objects of class \code{"recommend_k_results"} returned
#' by the k-recommendation procedure (e.g., based on the effective number of
#' parameters \code{p_loo} or \code{p_waic}).
#'
#' @param x An object of class \code{"recommend_k_results"}.
#' @param digits Integer; number of digits to display for numeric summaries.
#' @param ... Further arguments passed to or from other methods (currently
#'   ignored).
#'
#' @return The input object \code{x}, invisibly.
#'
#' @export
print.recommend_k_results <- function (x, digits = 3, ...) {

    if (!inherits(x, "recommend_k_results") ) {

        stop ("`x` must be of class 'recommend_k_results'.", call. = FALSE)

    }

    cat("\n==== k recommendation results ==================================\n\n")

    # basic info
    n_models <- if (!is.null(x$models) ) length(x$models) else NA_integer_
    k_range  <- if (!is.null(x$k_values) ) range(x$k_values) else c(NA, NA)

    cat("Number of models fitted  : ", n_models, "\n", sep = "")
    cat("Range of k values tested : [",
        paste(round(k_range, digits), collapse = ", "),
        "]\n", sep = "")

    if (!is.null(x$knee_on) ) {

        cat("Knee based on            : ", x$knee_on, "\n", sep = "")

    }

    # recommended k
    rec_k <- x$recommended_k
    cat("Recommended k (knee)     : ", rec_k, "\n", sep = "")
    cat("\n=================================================================\n")
    invisible(x)

}

#' Summary method for recommend_k_results objects
#'
#' Provides a concise summary of k-tuning results, including the range of
#' tested \code{k} values, the effective complexity measures
#' (e.g., \code{p_loo} or \code{p_waic}), and the recommended \code{k}
#' obtained via knee detection.
#'
#' @param object An object of class \code{"recommend_k_results"}.
#' @param digits Integer; number of digits to display for numeric summaries.
#' @param ... Further arguments passed to or from other methods (currently
#'   ignored).
#'
#' @return The input object \code{object}, invisibly.
#'
#' @export
summary.recommend_k_results <- function (object, digits = 3, ...) {

    if (!inherits(object, "recommend_k_results") ) {

        stop ("`object` must be of class 'recommend_k_results'.", call. = FALSE)

    }

    cat("\n==== Summary of k recommendation ===================================\n\n")

    n_models <- if (!is.null(object$models) ) length(object$models) else NA_integer_
    k_vals <- object$k_values %||% NA_real_
    k_range <- range(k_vals, na.rm = TRUE)

    cat("Number of models fitted  : ", n_models, "\n", sep = "")
    cat("k values tested          : ",
        paste(sort(unique(round(k_vals, digits) ) ), collapse = ", "),
        "\n", sep = "")

    cat("Range of k values        : [",
        paste(round(k_range, digits), collapse = ", "),
        "]\n", sep = "")

    if (!is.null(object$knee_on) ) {

        cat("Knee based on            : ", object$knee_on, "\n", sep = "")

    }

    rec_k <- object$recommended_k
    cat("Recommended k (knee)     : ", rec_k, "\n\n", sep = "")

    # basic complexity stats
    comp <- object$comparison

    if (!is.null(comp) && is.data.frame(comp) ) {

        # determine which complexity column was used
        p_col <- NULL

        if (!is.null(object$knee_on) && object$knee_on %in% names(comp) ) {

            p_col <- object$knee_on

        } else if ("p_loo" %in% names(comp) ) {

            p_col <- "p_loo"

        } else if ("p_waic" %in% names(comp) ) {

            p_col <- "p_waic"

        }

        # small comparison table (rounded)
        cat("Comparison table (rounded):\n\n")
        comp_print <- comp
        num_cols <- vapply(comp_print, is.numeric, logical(1L) )
        comp_print[num_cols] <- lapply(comp_print[num_cols], round, digits = digits)
        comp_print <- comp_print |> dplyr::arrange(.data$waic)

        print(comp_print, row.names = FALSE)

    } else {

        cat("No comparison table available.\n")

    }

    cat("\n====================================================================\n")
    invisible(object)

}
