#' Plot time-resolved GAM results (1D clusters)
#'
#' Visualises the fitted posterior effect through time for a
#' \code{clusters_results_1d} object. The plot shows:
#' \itemize{
#'   \item the posterior mean (\code{x$predictions$post_prob}) as a line,
#'   \item its uncertainty interval (\code{x$predictions$lower}–\code{x$predictions$upper}) as a ribbon,
#'   \item and detected significant clusters (\code{x$clusters}) as thick horizontal segments.
#' }
#'
#' When possible (depending on \code{x$multilevel} and the structure of \code{x$model$data}),
#' the method also reconstructs and overlays an empirical time course (e.g., averaged outcome
#' or difference between two predictor levels).
#'
#' If clusters are available at the participant level (i.e., a \code{participant} column exists
#' in \code{x$clusters}), the plot is facetted by participant.
#'
#' @param x A \code{clusters_results_1d} object (typically returned by
#'   \code{\link[neurogam]{testing_through_time}} for 1D time).
#' @param null_value Numeric. Reference value shown as a horizontal dashed line.
#'   Defaults to \code{0}.
#' @param clusters_y Numeric. Vertical position used to draw the cluster segments.
#'   Defaults to \code{-Inf} (draws at the bottom of the panel).
#' @param clusters_colour Character. Colour used for the prediction line/ribbon and
#'   the cluster segments. Defaults to \code{"black"}.
#' @param lineend Character. Line end style for cluster segments, passed to
#'   \code{\link[ggplot2]{geom_segment}} (e.g., \code{"butt"}, \code{"round"}).
#' @param theme A ggplot2 theme object. Defaults to \code{\link[ggplot2]{theme_bw}()}.
#' @param ... Currently unused. Included for S3 compatibility.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object.
#'
#' @examples
#' \dontrun{
#' res <- testing_through_time(...)
#' plot(res)
#' }
#'
#' @export
plot.clusters_results_1d <- function (
        x,
        null_value = 0,
        clusters_y = -Inf,
        clusters_colour = "black",
        lineend = "butt",
        theme = ggplot2::theme_bw(),
        ...
        ) {

    if (!ggplot2::is.theme(theme) ) {

        stop ("Argument 'theme' should be a 'theme' object.")

    }

    emp_data <- x$model$data
    clusters <- x$clusters
    preds <- x$predictions

    # group-level or participant-level clusters?
    group_level <- ifelse(
        test = "participant" %in% names(clusters),
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

#' Plot time-generalisation GAM results (2D clusters)
#'
#' Visualises 2D (time-by-time) posterior predictions from a
#' \code{clusters_results_2d} object as a two-panel figure:
#' \enumerate{
#'   \item summary empirical data
#'   \item summary posterior predictions
#' }
#'
#' Significant clusters stored in \code{x$clusters} are assumed to be
#' \strong{pointwise} (one row per time-by-time cell). Cluster locations are
#' overlaid as contour (which includes interpolation).
#'
#' This method relies on \pkg{patchwork} to combine the two panels using \code{+}.
#'
#' @param x A \code{clusters_results_2d} object (typically returned by
#'   \code{\link[neurogam]{testing_through_time}} for 2D time / time-generalisation).
#' @param theme A ggplot2 theme object added to each panel. Defaults to
#'   \code{\link[ggplot2]{theme_bw}()}.
#' @param palette Character, a scico palette, see scico::scico_palette_show().
#' @param midpoint Numeric, the midpoint for diverging color palettes.
#' @param axes_labels Character vector, the x- and y-axis labels.
#' @param cluster_outline_colour Character. Colour used to outline cluster cells.
#'   Defaults to \code{"white"}.
#' @param ... Currently unused. Included for S3 compatibility.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object (a patchwork composition of two ggplots).
#'
#' @examples
#' \dontrun{
#' res2d <- testing_through_time(...)
#' plot(res2d)
#' }
#'
#' @export
plot.clusters_results_2d <- function (
        x,
        theme = ggplot2::theme_bw(),
        palette = "vik",
        midpoint = 0.5,
        axes_labels = c("Training time (s)", "Testing time (s)"),
        cluster_outline_colour = "black",
        ...
        ) {

    if (!palette %in% scico::scico_palette_names() ) {

        stop ("Argument 'palette' should be a scio palette. See scico::scico_palette_show().")

    }

    if (!ggplot2::is.theme(theme) ) {

        stop ("Argument 'theme' should be a 'theme' object.")

    }

    emp_data <- x$model$data
    clusters <- x$clusters
    preds <- x$predictions

    if (!is.null(clusters) && nrow(clusters) > 0) {

        bookkeeping <- c("id", "sign", "value", "n_points", "participant")
        candidate_cols <- setdiff(names(clusters), bookkeeping)

        is_num <- vapply(clusters[, candidate_cols, drop = FALSE], is.numeric, logical(1) )
        time_cols <- candidate_cols[is_num]

        if (length(time_cols) == 2) {

            is_2d <- TRUE
            time_id <- time_cols

        }

    }

    stopifnot(!is.null(preds) )
    stopifnot(length(time_id) == 2)

    if (!all(time_id %in% names(preds) ) ) {

        stop ("x$predictions must contain the 2D time columns...")

    }

    if (!all(c("post_prob", "prob_ratio") %in% names(preds) ) ) {

        stop ("For the 2D two-panel plot, x$predictions must contain columns: 'post_prob' and 'prob_ratio'.")

    }

    # plotting summary data
    emp_data_summary <- emp_data |>
        dplyr::summarise(
            outcome_mean = mean(.data$outcome_mean),
            .by = c(.data$time1, .data$time2)
            )

    p1 <- plot_clusters_2d(
        field_df = emp_data_summary,
        clusters_df = clusters,
        value_col = "outcome_mean",
        palette = palette,
        midpoint = midpoint,
        fill_name = "Observed",
        axes_labels = axes_labels,
        cluster_colour = cluster_outline_colour
        ) + theme

    # plotting model predictions
    p2 <- plot_clusters_2d(
        field_df = preds,
        clusters_df = clusters,
        value_col = "post_prob",
        palette = palette,
        midpoint = midpoint,
        axes_labels = axes_labels,
        fill_name = "Predicted",
        cluster_colour = cluster_outline_colour
        ) + theme

    return (patchwork::wrap_plots(p1, p2, nrow = 1) )

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

    if (inherits(object, "clusters_results_2d" ) ) {

        stop ("PPCs are not supported for 2D cluster results.")

    }

    fit <- object$model
    ppc_type <- match.arg(ppc_type)

    if (!inherits(object, "clusters_results_1d") ) {

        stop ("`object` must be of class 'clusters_results_1d'.", call. = FALSE)

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
                # incl_autocor = FALSE,
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
                    # incl_autocor = FALSE,
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
