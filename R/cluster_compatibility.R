#' Convert multidimensional clustering results to legacy output
#'
#' Internal adapter retaining the public `clusters` data-frame contract.
#' The standardised cluster summaries and memberships remain available
#' separately from the output of `find_clusters_nd()`.
#'
#' @param result A list returned by `find_clusters_nd()`.
#' @param time_id One or two temporal coordinate column names.
#' @param group An optional grouping column name.
#'
#' @return A legacy-compatible data frame: one row per cluster for 1D,
#'   or one row per cluster member for 2D.
#'
#' @keywords internal
cluster_nd_to_legacy <- function (result, time_id = "time", group = NULL) {

    stopifnot(
        is.list(result),
        is.data.frame(result$clusters),
        is.data.frame(result$membership),
        is.character(time_id),
        length(time_id) %in% c(1L, 2L)
        )

    group_cols <- if (is.null(group) ) character() else group

    if (length(time_id) == 1L) {

        out <- result$clusters
        min_col <- paste0(time_id, "_min")
        max_col <- paste0(time_id, "_max")

        # `find_clusters_nd()` supplies onset/offset automatically only
        # when its coordinate is literally named `time`.
        out$onset <- out[[min_col]]
        out$offset <- out[[max_col]]

        if (nrow(out) > 0L) {

            order_cols <- c(group_cols, "sign", "onset", "offset")
            ord <- do.call(order, unname(out[order_cols]) )
            out <- out[ord, , drop = FALSE]

            # Legacy 1D cluster IDs start at 1 separately by sign/group.
            strata <- interaction(
                out[c(group_cols, "sign")],
                drop = TRUE,
                lex.order = TRUE
                )

            out$id <- as.integer(stats::ave(seq_len(nrow(out) ), strata, FUN = seq_along) )

        }

        columns <- c(group_cols, "id", "onset", "offset", "n_points", "sign")
        out <- out[, columns, drop = FALSE]

    } else {

        out <- result$membership
        key_cols <- c(group_cols, "id")
        sizes <- result$clusters[, c(key_cols, "n_points"), drop = FALSE]

        if (nrow(out) > 0L) {

            # Match on the complete group + cluster key; IDs may be
            # repeated in different groups.
            key <- function (dat) {

                do.call(paste, c(unname(dat[key_cols]), sep = "\r") )

            }

            positions <- match(key(out), key(sizes) )

            if (anyNA(positions) ) {

                stop ("Cluster membership contains an unknown cluster ID.", call. = FALSE)

            }

            out$n_points <- sizes$n_points[positions]
            order_cols <- c(group_cols, "sign", time_id)
            out <- out[do.call(order, unname(out[order_cols]) ), , drop = FALSE]

        } else {

            out$n_points <- integer()

        }

        columns <- c(group_cols, "id", "sign", time_id, "value", "n_points")
        out <- out[, columns, drop = FALSE]

    }

    rownames(out) <- NULL

    return (as.data.frame(out) )

}
