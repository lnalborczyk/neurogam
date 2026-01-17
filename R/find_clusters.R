#' Find contiguous clusters in a time series
#'
#' Identify contiguous clusters of time points in a time series where a
#' variable exceeds a positive threshold, falls below a negative threshold,
#' or both. Clusters are defined as consecutive time points satisfying the
#' thresholding condition.
#'
#' By default, the function detects \strong{both positive and negative clusters}
#' in a single call and returns a column indicating the cluster sign.
#'
#' If a grouping variable is provided (e.g., \code{"participant"}), clusters
#' are detected independently within each group level.
#'
#' @param data A data frame containing at least two columns named
#'   \code{"time"} and \code{"value"}. If \code{group} is not \code{NULL},
#'   \code{data} must also contain that grouping column.
#' @param threshold Numeric scalar specifying the (positive) threshold used to
#'   define clusters. Positive clusters are defined where
#'   \code{value >= threshold}, and negative clusters where
#'   \code{value <= 1/threshold}. Must be non-negative when
#'   \code{threshold_type = "both"}.
#' @param group Optional grouping column name (character scalar) used to find
#'   clusters independently within each group level (e.g.,
#'   \code{"participant"}). Set to \code{NULL} (default) to ignore grouping.
#' @param threshold_type Character scalar controlling which clusters are
#'   detected. Must be one of \code{"above"}, \code{"below"}, or \code{"both"}
#'   (default). When \code{"above"}, clusters are formed where
#'   \code{value >= threshold}. When \code{"below"}, clusters are formed where
#'   \code{value <= 1/threshold}. When \code{"both"}, both types are detected
#'   and the returned data include a \code{sign} column.
#' @param time_id Character; name of the column(s) in \code{data}
#' containing time information (e.g., in seconds or samples).
#'
#' @return A data frame with one row per detected cluster and columns:
#'   \itemize{
#'     \item \code{id}: integer cluster index (starting at 1). When
#'       \code{group} is provided, \code{id} restarts at 1 within each
#'       group level;
#'     \item \code{onset}: time of the first point in the cluster;
#'     \item \code{offset}: time of the last point in the cluster;
#'     \item \code{n_points}: number of time points in the cluster;
#'     \item \code{sign}: character indicating cluster sign
#'       (\code{"positive"} or \code{"negative"}).
#'   }
#'
#'   If \code{group} is not \code{NULL}, the returned data frame also contains
#'   the grouping column (named as in \code{group}).
#'
#'   If no clusters are found, an empty data frame with the same column structure
#'   is returned.
#'
#' @details
#' The function assumes that the \code{time} variable is numeric and that
#' consecutive rows correspond to consecutive time points (within each group
#' if grouping is used). Internally, the data are:
#' \enumerate{
#'   \item filtered to remove rows with missing values;
#'   \item arranged by \code{time} (and by \code{group} then \code{time}, if used);
#'   \item thresholded to identify positive and/or negative excursions;
#'   \item segmented into runs of consecutive threshold-exceeding values,
#'         which define clusters.
#' }
#'
#' @importFrom rlang .data
#' @importFrom dplyr all_of
#'
#' @author Ladislas Nalborczyk \email{ladislas.nalborczyk@@cnrs.fr}
#'
#' @examples
#' \dontrun{
#' set.seed(666)
#' df <- data.frame(
#'   time = seq(0, 1, length.out = 100),
#'   value = c(
#'     rnorm(30, 0, 1),
#'     rnorm(20,  4, 1), # positive cluster
#'     rnorm(20, -4, 1), # negative cluster
#'     rnorm(30, 0, 1)
#'     )
#'   )
#'
#' # Detect both positive and negative clusters
#' find_clusters(data = df, threshold = 3, threshold_type = "both")
#'
#' # One-sided detection (positive only)
#' find_clusters(data = df, threshold = 3, threshold_type = "above")
#'
#' # One-sided detection (negative only)
#' find_clusters(data = df, threshold = 3, threshold_type = "below")
#'
#' # Grouped example (e.g., per participant)
#' df_g <- rbind(
#'   transform(df, participant = "P01"),
#'   transform(df, participant = "P02")
#'   )
#'
#' find_clusters(
#'   data = df_g,
#'   threshold = 3,
#'   group = "participant",
#'   threshold_type = "both"
#'   )
#' }
#'
#' @author Ladislas Nalborczyk \email{ladislas.nalborczyk@@cnrs.fr}.
#'
#' @export
find_clusters <- function (
        data,
        threshold = 10,
        group = NULL,
        threshold_type = c("both", "above", "below"),
        time_id = "time"
        ) {

    stopifnot("`data` must be a data frame." = is.data.frame(data) )
    threshold_type <- match.arg(threshold_type)

    if (!is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold) ) {

        stop ("`threshold` must be a finite numeric scalar.", call. = FALSE)

    }

    if (threshold_type == "both" && threshold < 0) {

        stop ("`threshold` must be >= 0 when `threshold_type = \"both\"`.", call. = FALSE)

    }

    # time_id validation
    if (!is.character(time_id) || !(length(time_id) %in% c(1L, 2L) ) ) {

        stop ("`time_id` must be a character vector of length 1 or 2.", call. = FALSE)

    }

    if (length(time_id) == 2L && identical(time_id[[1]], time_id[[2]]) ) {

        stop ("When `time_id` has length 2, the two names must be different.", call. = FALSE)

    }

    if (!is.null(group) ) {

        if (!is.character(group) || length(group) != 1L || group == "") {

            stop ("`group` must be NULL or a single non-empty character string.", call. = FALSE)

        }

        if (!group %in% names(data) ) {

            stop ("Grouping column '", group, "' not found in `data`.", call. = FALSE)

        }

    }

    # required columns
    required_columns <- c(time_id, "value")
    if (!is.null(group) ) required_columns <- c(required_columns, group)

    missing_cols <- setdiff(required_columns, names(data) )

    if (length(missing_cols) > 0L) {

        stop (
            "Data is missing required column(s): ",
            paste(missing_cols, collapse = ", "),
            call. = FALSE
            )

    }

    # subset columns and filter missing
    data <- data |>
        dplyr::select(dplyr::all_of(required_columns) ) |>
        # dplyr::filter(
        #     # !is.na(.data$time),
        #     # if (length(time_id) == 1L) !is.na(.data[[time_id[[1]]]]),
        #     !is.na(.data$value)
        #     # if (is.null(group) ) TRUE else !is.na(.data[[group]])
        #     )
        dplyr::filter(!is.na(.data$value) )

    ####################
    # 1D temporal case #
    ####################

    if (length(time_id) == 1L) {

        time_col <- time_id[[1]]
        data_1d <- data |> dplyr::rename(time = dplyr::all_of(time_col) )

        data_1d <- if (is.null(group) ) {
            dplyr::arrange(data_1d, .data$time)
        } else {
            dplyr::arrange(data_1d, .data[[group]], .data$time)
        }

        compute_clusters_1d <- function (dat, direction = c("positive", "negative"), sign_label) {

            direction <- match.arg(direction)

            if (is.null(group) ) {
                out <- dat |>
                    dplyr::mutate(
                        hit = if (direction == "positive") .data$value >= threshold else .data$value <= 1 / threshold,
                        change = dplyr::lag(.data$hit, default = FALSE) != .data$hit,
                        id = cumsum(.data$change & .data$hit)
                        ) |>
                    dplyr::filter(.data$hit) |>
                    dplyr::group_by(.data$id) |>
                    dplyr::summarise(
                        onset = dplyr::first(.data$time),
                        offset = dplyr::last(.data$time),
                        n_points = dplyr::n(),
                        .groups = "drop"
                        ) |>
                    dplyr::mutate(sign = sign_label) |>
                    data.frame()

            } else {

                out <- dat |>
                    dplyr::group_by(.data[[group]]) |>
                    dplyr::mutate(
                        hit = if (direction == "positive") .data$value >= threshold else .data$value <= 1 / threshold,
                        change = dplyr::lag(.data$hit, default = FALSE) != .data$hit,
                        id = cumsum(.data$change & .data$hit)
                        ) |>
                    dplyr::filter(.data$hit) |>
                    dplyr::group_by(.data[[group]], .data$id) |>
                    dplyr::summarise(
                        onset = dplyr::first(.data$time),
                        offset = dplyr::last(.data$time),
                        n_points = dplyr::n(),
                        .groups = "drop"
                        ) |>
                    dplyr::mutate(sign = sign_label) |>
                    data.frame()

            }

            return (out)

        }

        clusters <- switch (
            threshold_type,
            both  = dplyr::bind_rows(
                compute_clusters_1d(data_1d, "positive", "positive"),
                compute_clusters_1d(data_1d, "negative", "negative")
                ),
            above = compute_clusters_1d(data_1d, "positive", "positive"),
            below = compute_clusters_1d(data_1d, "negative", "negative")
            )

        if (!is.null(group) && nrow(clusters) > 0) clusters <- dplyr::arrange(clusters, .data[[group]])

        return (clusters)

    }

    ######################################################
    # 2D temporal case                                   #
    # Return per-cell memberships (t1,t2) with id + sign #
    ######################################################

    t1 <- time_id[[1]]
    t2 <- time_id[[2]]

    data <- if (is.null(group) ) {

        dplyr::arrange(data, .data[[t1]], .data[[t2]])

    } else {

        dplyr::arrange(data, .data[[group]], .data[[t1]], .data[[t2]])

    }

    # 4-neighbour connected components on a logical matrix
    label_components <- function (mask) {

        nr <- nrow(mask)
        nc <- ncol(mask)
        lab <- matrix(0L, nr, nc)
        cur <- 0L
        qx <- integer(nr * nc); qy <- integer(nr * nc)

        for (r in seq_len(nr) ) for (c in seq_len(nc) ) {

            if (!mask[r, c] || lab[r, c] != 0L) next

            cur <- cur + 1L
            head <- 1L
            tail <- 1L
            qx[1L] <- r
            qy[1L] <- c
            lab[r, c] <- cur

            while (head <= tail) {

                x <- qx[head]; y <- qy[head]; head <- head + 1L
                if (x > 1L  && mask[x-1L, y] && lab[x-1L, y] == 0L) { tail <- tail + 1L; qx[tail] <- x-1L; qy[tail] <- y;   lab[x-1L, y] <- cur }
                if (x < nr  && mask[x+1L, y] && lab[x+1L, y] == 0L) { tail <- tail + 1L; qx[tail] <- x+1L; qy[tail] <- y;   lab[x+1L, y] <- cur }
                if (y > 1L  && mask[x, y-1L] && lab[x, y-1L] == 0L) { tail <- tail + 1L; qx[tail] <- x;   qy[tail] <- y-1L; lab[x, y-1L] <- cur }
                if (y < nc  && mask[x, y+1L] && lab[x, y+1L] == 0L) { tail <- tail + 1L; qx[tail] <- x;   qy[tail] <- y+1L; lab[x, y+1L] <- cur }

            }

        }

        return (lab)

    }

    compute_membership_2d <- function (dat, direction = c("positive", "negative"), sign_label) {

        direction <- match.arg(direction)

        hit_fun <- if (direction == "positive") {

            function (v) v >= threshold

        } else {

            function (v) v <= 1 / threshold

        }

        build_for_slice <- function (slice_df) {

            u1 <- sort(unique(slice_df[[t1]]) )
            u2 <- sort(unique(slice_df[[t2]]) )

            i1 <- match(slice_df[[t1]], u1)
            i2 <- match(slice_df[[t2]], u2)

            mat_hit <- matrix(FALSE, nrow = length(u1), ncol = length(u2) )
            mat_hit[cbind(i1, i2)] <- hit_fun(slice_df$value)

            labs <- label_components(mat_hit)
            nlab <- max(labs)

            if (nlab == 0L) return (data.frame() )

            # map each hit cell back to (t1,t2) with its component id
            idx <- which(labs > 0L, arr.ind = TRUE)

            out <- tibble::tibble(
                id = labs[idx],
                sign = sign_label,
                .name_repair = "minimal"
                )

            out <- dplyr::bind_cols(
                out,
                tibble::as_tibble(stats::setNames(list(u1[idx[, 1]]), t1) ),
                tibble::as_tibble(stats::setNames(list(u2[idx[, 2]]), t2) )
                )

            # attach value if you want it available for debugging/plotting
            # (we join from slice_df to avoid matrix reconstruction issues)
            out <- out |>
                dplyr::left_join(
                    slice_df |> dplyr::select(dplyr::all_of(c(t1, t2) ), .data$value),
                    by = c(t1, t2)
                    )

            # n_points per component (repeated per row, 1D-like metadata)
            out <- out |>
                dplyr::group_by(.data$id) |>
                dplyr::mutate(n_points = dplyr::n()) |>
                dplyr::ungroup()

            return (data.frame(out) )

        }

        if (is.null(group) ) {

            build_for_slice(dat)

        } else {

            dat |>
                dplyr::group_by(.data[[group]]) |>
                dplyr::group_modify(\(d, ...) build_for_slice(d)) |>
                dplyr::ungroup() |>
                data.frame()

        }

    }

    clusters <- switch (
        threshold_type,
        both = {
            out_pos <- compute_membership_2d(data, "positive", "positive")
            out_neg <- compute_membership_2d(data, "negative", "negative")
            # make ids unique across sign within each group (optional, but avoids collisions)
            if (!is.null(group) ) {
                out_neg <- out_neg |>
                    dplyr::group_by(.data[[group]]) |>
                    dplyr::mutate(id = .data$id + dplyr::coalesce(max(out_pos$id[out_pos[[group]] == .data[[group]][1]]), 0L)) |>
                    dplyr::ungroup()
            } else {
                if (nrow(out_pos) > 0 && nrow(out_neg) > 0) out_neg$id <- out_neg$id + max(out_pos$id)
            }

            dplyr::bind_rows(out_pos, out_neg)
        },
        above = compute_membership_2d(data, "positive", "positive"),
        below = compute_membership_2d(data, "negative", "negative")
        )

    if (!is.null(group) && nrow(clusters) > 0) clusters <- dplyr::arrange(clusters, .data[[group]])

    return (clusters)

}
