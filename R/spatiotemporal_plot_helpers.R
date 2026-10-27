# Internal helpers for plotting sensor-by-time BGAMM results.

validate_st_plot_columns <- function(data, columns, object_name) {
    if (!is.data.frame(data)) {
        stop("`", object_name, "` must be a data frame.", call. = FALSE)
    }

    missing <- setdiff(columns, names(data))
    if (length(missing)) {
        stop(
            "Missing columns in `", object_name, "`: ",
            paste(missing, collapse = ", "), call. = FALSE
        )
    }
}

summarise_observed_st <- function(
    data, outcome_id, sensor_id, time_id, participant_id = NULL
) {
    validate_st_plot_columns(
        data, c(sensor_id, time_id, outcome_id, participant_id),
        "observed_data"
    )

    if (!is.numeric(data[[outcome_id]])) {
        stop("The observed outcome must be numeric.", call. = FALSE)
    }
    if (!is.numeric(data[[time_id]]) || any(!is.finite(data[[time_id]]))) {
        stop("Observed times must be finite and numeric.", call. = FALSE)
    }

    participant <- if (is.null(participant_id)) {
        NULL
    } else {
        as.character(data[[participant_id]])
    }
    data <- data.frame(
        sensor = as.character(data[[sensor_id]]),
        time = data[[time_id]],
        value = data[[outcome_id]],
        stringsAsFactors = FALSE
    )
    if (!is.null(participant)) {
        if (anyNA(participant)) {
            stop("Observed participant IDs cannot be missing.", call. = FALSE)
        }
        data$participant <- participant
    }
    if (anyNA(data[c("sensor", "time")])) {
        stop("Observed sensor/time identifiers cannot be missing.", call. = FALSE)
    }

    if (!is.null(participant_id)) {
        # Average trials within participants before averaging participants.
        first <- stats::aggregate(
            value ~ sensor + time + participant,
            data = data,
            FUN = function(z) {
                if (all(is.na(z))) NA_real_ else mean(z, na.rm = TRUE)
            },
            na.action = stats::na.pass
        )
        data <- first
    }
    out <- stats::aggregate(
        value ~ sensor + time,
        data = data,
        FUN = function(z) {
            if (all(is.na(z))) {
                NA_real_
            } else {
                mean(z, na.rm = TRUE)
            }
        },
        na.action = stats::na.pass
    )
    names(out)[3L] <- "observed"
    out
}

#' Prepare predicted and observed EEG trajectories for plotting
#'
#' Observations are averaged within each sensor-time location. For a model
#' contrast A - B, supply `observed_data` for A and `observed_data_b` for B.
#' The displayed residual is the difference between observed *mean* and
#' posterior expected value, not a trial-level residual.
#'
#' @param results Result from [analyse_bgam_spatiotemporal()].
#' @param observed_data Observed data in long format (condition A if relevant).
#' @param observed_data_b Optional observed condition B for a model contrast.
#' @param outcome_id Outcome column of observed data.
#' @param sensor_id Sensor column name in both data and `results$predictions`.
#' @param time_id Numeric time column.
#' @param spatial_id Two projected coordinate column names.
#' @param participant_id Optional participant column for an equal-participant-
#'   weighted grand average (average within participant first).
#'
#' @return A data frame with standardized sensor, time, x, y, observed,
#'   predicted, lower, upper, residual, cluster_id and cluster_sign columns.
#' @keywords internal
prepare_spatiotemporal_plot_data <- function(
    results,
    observed_data,
    observed_data_b = NULL,
    outcome_id = "voltage",
    sensor_id = "channel",
    time_id = "time",
    spatial_id = c("xproj", "yproj"),
    participant_id = NULL
) {
    if (!is.list(results) || is.null(results$predictions) ||
        is.null(results$cluster_membership)) {
        stop("`results` needs `predictions` and `cluster_membership`.",
             call. = FALSE)
    }
    names_needed <- c(sensor_id, time_id, spatial_id)
    if (!is.character(names_needed) || anyNA(names_needed) ||
        any(!nzchar(names_needed)) || length(spatial_id) != 2L ||
        anyDuplicated(names_needed)) {
        stop("Specify distinct sensor, time, and two spatial columns.",
             call. = FALSE)
    }
    validate_st_plot_columns(
        results$predictions,
        c(names_needed, "estimate", "lower", "upper"),
        "results$predictions"
    )

    pred <- results$predictions
    if (!is.numeric(pred[[time_id]]) || any(!is.finite(pred[[time_id]])) ||
        !all(vapply(pred[spatial_id], is.numeric, logical(1))) ||
        any(!is.finite(as.matrix(pred[spatial_id])))) {
        stop("Prediction coordinates must be finite numeric values.",
             call. = FALSE)
    }
    if (anyNA(pred[c(sensor_id, time_id)])) {
        stop("Prediction sensor/time IDs must not be missing.", call. = FALSE)
    }
    if (anyDuplicated(pred[c(sensor_id, time_id)])) {
        stop("Predictions must have one row per sensor and time.",
             call. = FALSE)
    }

    pred <- data.frame(
        sensor = as.character(pred[[sensor_id]]),
        time = pred[[time_id]],
        x = pred[[spatial_id[1L]]],
        y = pred[[spatial_id[2L]]],
        predicted = pred$estimate,
        lower = pred$lower,
        upper = pred$upper,
        stringsAsFactors = FALSE
    )
    if (any(!is.finite(as.matrix(pred[c("predicted", "lower", "upper")])))) {
        stop("Predictions and credible intervals must be finite.",
             call. = FALSE)
    }
    layout <- unique(pred[c("sensor", "x", "y")])
    if (anyDuplicated(layout$sensor)) {
        stop("Each sensor must have a single pair of coordinates.",
             call. = FALSE)
    }
    if (anyDuplicated(layout[c("x", "y")])) {
        stop("Different sensors have identical projected coordinates.",
             call. = FALSE)
    }

    is_contrast <- identical(results$contrast, "A_minus_B")
    if (is_contrast && is.null(observed_data_b)) {
        stop(
            "This result is an A - B contrast; supply `observed_data_b` ",
            "to obtain the corresponding observed contrast.",
            call. = FALSE
        )
    }
    if (!is_contrast && !is.null(observed_data_b)) {
        stop("`observed_data_b` is only valid for A - B contrasts.",
             call. = FALSE)
    }

    obs <- summarise_observed_st(
        observed_data, outcome_id, sensor_id, time_id, participant_id
    )
    if (is_contrast) {
        obs_b <- summarise_observed_st(
            observed_data_b, outcome_id, sensor_id, time_id, participant_id
        )
        names(obs_b)[3L] <- "observed_b"
        obs <- merge(obs, obs_b, by = c("sensor", "time"), all = TRUE)
        obs$observed <- obs$observed - obs$observed_b
        obs$observed_b <- NULL
    }

    # all.x preserves the prediction grid, including locations with no
    # corresponding observations.
    dat <- merge(pred, obs, by = c("sensor", "time"), all.x = TRUE, sort = FALSE)
    dat <- dat[order(dat$sensor, dat$time), , drop = FALSE]
    rownames(dat) <- NULL
    if (all(is.na(dat$observed))) {
        stop("No observed sensor-time samples match the prediction grid.",
             call. = FALSE)
    }
    if (anyNA(dat$observed)) {
        warning(
            "Some predictions lack matching observed means; these are ",
            "omitted from observed and residual layers.", call. = FALSE
        )
    }
    dat$residual <- dat$observed - dat$predicted

    members <- results$cluster_membership
    if (nrow(members)) {
        validate_st_plot_columns(
            members, c(sensor_id, time_id, "id", "sign"),
            "results$cluster_membership"
        )
        members <- data.frame(
            sensor = as.character(members[[sensor_id]]),
            time = members[[time_id]],
            cluster_id = members$id,
            cluster_sign = as.character(members$sign),
            stringsAsFactors = FALSE
        )
        if (anyDuplicated(members[c("sensor", "time")])) {
            stop("Each sensor-time sample can belong to only one cluster.",
                 call. = FALSE)
        }
        dat <- merge(
            dat, members, by = c("sensor", "time"), all.x = TRUE,
            sort = FALSE
        )
        dat <- dat[order(dat$sensor, dat$time), , drop = FALSE]
        rownames(dat) <- NULL
    } else {
        dat$cluster_id <- NA_integer_
        dat$cluster_sign <- NA_character_
    }
    dat
}

st_time_settings <- function(dat, max_gap = NULL) {
    times <- sort(unique(dat$time))
    spacings <- diff(times)
    step <- if (length(spacings)) stats::median(spacings) else 1
    if (!is.finite(step) || step <= 0) {
        stop("Distinct numeric time points are required.", call. = FALSE)
    }
    if (is.null(max_gap)) {
        max_gap <- 1.5 * step
    }
    if (!is.numeric(max_gap) || length(max_gap) != 1L ||
        !is.finite(max_gap) || max_gap <= 0) {
        stop("`max_gap` must be a positive finite scalar.", call. = FALSE)
    }
    list(times = times, step = step, max_gap = max_gap)
}

#' Construct nonbridging sensor-specific cluster intervals
#'
#' Splits cluster memberships at changes of sign, cluster ID, absent
#' time samples, or time gaps exceeding `max_gap`.
#'
#' @param data Standardized output of `prepare_spatiotemporal_plot_data()`.
#' @param max_gap Largest allowed temporal gap. By default, 1.5 times the
#'   median spacing in the prediction grid.
#'
#' @return Data frame containing sensor, sign, cluster_id, xmin and xmax.
#' @keywords internal
make_cluster_intervals_per_sensor <- function(data, max_gap = NULL) {
    settings <- st_time_settings(data, max_gap)
    hits <- data[!is.na(data$cluster_id), , drop = FALSE]
    empty <- data.frame(
        sensor = character(), cluster_id = integer(), sign = character(),
        xmin = numeric(), xmax = numeric(), stringsAsFactors = FALSE
    )
    if (!nrow(hits)) {
        return(empty)
    }

    hits$time_index <- match(hits$time, settings$times)
    splits <- split(
        seq_len(nrow(hits)),
        interaction(hits$sensor, hits$cluster_sign, hits$cluster_id,
                    drop = TRUE)
    )
    pieces <- list()
    for (rows in splits) {
        z <- hits[rows, , drop = FALSE]
        z <- z[order(z$time), , drop = FALSE]
        breaks <- c(
            TRUE,
            diff(z$time_index) != 1L |
                diff(z$time) > settings$max_gap
        )
        run <- cumsum(breaks)
        for (idx in split(seq_len(nrow(z)), run)) {
            a <- z[idx, , drop = FALSE]
            pieces[[length(pieces) + 1L]] <- data.frame(
                sensor = a$sensor[1L],
                cluster_id = a$cluster_id[1L],
                sign = a$cluster_sign[1L],
                xmin = min(a$time) - settings$step / 2,
                xmax = max(a$time) + settings$step / 2,
                stringsAsFactors = FALSE
            )
        }
    }
    out <- do.call(rbind, pieces)
    rownames(out) <- NULL
    out
}

#' Construct highlighted prediction-line segments within clusters
#'
#' Each segment connects only neighboring points sharing the same sensor,
#' sign and cluster ID; no lines cross missing samples or long time gaps.
#'
#' @param data Standardized posterior/observation table.
#' @param max_gap Largest allowed temporal gap.
#'
#' @return Data frame with segment endpoints and sign.
#' @keywords internal
make_cluster_segments_per_sensor <- function(data, max_gap = NULL) {
    settings <- st_time_settings(data, max_gap)
    dat <- data[order(data$sensor, data$time), , drop = FALSE]
    dat$time_index <- match(dat$time, settings$times)
    if (nrow(dat) < 2L) {
        return(dat[FALSE, , drop = FALSE])
    }
    next_idx <- seq_len(nrow(dat) - 1L) + 1L
    start_idx <- next_idx - 1L
    same <- dat$sensor[start_idx] == dat$sensor[next_idx] &
        !is.na(dat$cluster_id[start_idx]) &
        !is.na(dat$cluster_id[next_idx]) &
        dat$cluster_id[start_idx] == dat$cluster_id[next_idx] &
        dat$cluster_sign[start_idx] == dat$cluster_sign[next_idx] &
        dat$time_index[next_idx] - dat$time_index[start_idx] == 1L &
        dat$time[next_idx] - dat$time[start_idx] <= settings$max_gap
    same[is.na(same)] <- FALSE
    starts <- start_idx[same]
    ends <- next_idx[same]
    if (!length(starts)) {
        out <- dat[FALSE, , drop = FALSE]
        out$timeend <- numeric()
        out$predictedend <- numeric()
        return(out)
    }
    out <- dat[starts, , drop = FALSE]
    out$timeend <- dat$time[ends]
    out$predictedend <- dat$predicted[ends]
    rownames(out) <- NULL
    out
}

st_head_geometry <- function(data, n = 180L) {
    layout <- unique(data[c("sensor", "x", "y")])
    center <- c(mean(range(layout$x)), mean(range(layout$y)))
    radius <- max(sqrt((layout$x - center[1L])^2 +
                       (layout$y - center[2L])^2)) * 1.04
    theta <- seq(0, 2 * pi, length.out = n)
    outline <- data.frame(
        x = center[1L] + radius * cos(theta),
        y = center[2L] + radius * sin(theta)
    )
    nose <- data.frame(
        x = center[1L] + radius * c(-0.13, 0, 0.13),
        y = center[2L] + radius * c(0.99, 1.13, 0.99)
    )
    list(center = center, radius = radius, outline = outline, nose = nose)
}

# Inverse-distance weighting (IDW) avoids a dependency on akima. It is a
# visualization method, not an inferential or physiological forward model.
interpolate_st_topography <- function(
    data, x_grid, y_grid, power = 2
) {
    usable <- is.finite(data$value)
    data <- data[usable, , drop = FALSE]
    if (!nrow(data)) {
        return(rep(NA_real_, length(x_grid)))
    }
    delta_x <- outer(x_grid, data$x, "-")
    delta_y <- outer(y_grid, data$y, "-")
    d2 <- delta_x^2 + delta_y^2
    exact <- d2 < 1e-20
    weights <- 1 / pmax(d2, 1e-20)^(power / 2)
    out <- as.vector(weights %*% data$value / rowSums(weights))
    exact_rows <- which(rowSums(exact) > 0)
    if (length(exact_rows)) {
        out[exact_rows] <- vapply(
            exact_rows,
            function(i) mean(data$value[which(exact[i, ])]),
            numeric(1)
        )
    }
    out
}
