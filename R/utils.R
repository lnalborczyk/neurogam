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
        is.character(predictor), length(predictor)   == 1
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

    # per_participant <- df |>
    #     dplyr::group_by(!!participant) |>
    #     dplyr::summarise(
    #         n_obs = dplyr::n(),
    #         n_unique = dplyr::n_distinct(!!predictor),
    #         sd_within = if (!is_cat) sd(!!predictor) else NA_real_,
    #         min = if (!is_cat) min(!!predictor) else NA_real_,
    #         max = if (!is_cat) max(!!predictor) else NA_real_,
    #         .groups = "drop"
    #         )

    # classify per participant
    # per_participant <- per_participant |>
    #     dplyr::mutate(
    #         varies_within = if (is_cat) {
    #             .data$n_unique > 1
    #         } else {
    #             (max - min) > tol
    #         }
    #     )

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
    if (limits == "global") return (range(df[[z]], na.rm = TRUE))

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
