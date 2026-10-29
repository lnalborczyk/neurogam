#' Internal dispatch classification for neurogam analysis dimensions
#'
#' This function routes analyses without altering any data or model settings.
#' Validation of 1D/2D and spatiotemporal data stays in the respective
#' analysis implementations.
#'
#' @param space_id Optional spatial-coordinate column names.
#' @return Either `"temporal"` or `"spatiotemporal"`.
#' @keywords internal
classify_bgam_analysis <- function (space_id = NULL) {

    if (is.null(space_id) ) {

        return ("temporal")

    }

    return ("spatiotemporal")

}

#' Normalise the options supplied to the spatiotemporal dispatcher
#'
#' @param sensor_id Name of the EEG sensor identifier column.
#' @param baseline_window Baseline-ROPE window or NULL.
#' @param varying_smooth Whether participant-specific smooths are requested.
#' @param varying_smooth_missing Whether the public caller omitted that option.
#' @param kvalue,bs Tensor marginal settings supplied by the user.
#' @param kvalue_missing,bs_missing Whether those options were omitted.
#' @return A list containing normalized `kvalue` and `bs`.
#' @keywords internal
normalize_spatiotemporal_dispatch <- function (
    sensor_id,
    baseline_window,
    varying_smooth,
    varying_smooth_missing,
    kvalue,
    kvalue_missing,
    bs,
    bs_missing
    ) {

    if (is.null(sensor_id) || !is.character(sensor_id) ||
        length(sensor_id) != 1L || is.na(sensor_id) ) {

        stop ("Provide `sensor_id` when using `space_id`.", call. = FALSE)

    }

    if (!identical(baseline_window, c(0, 0.1) ) &&
        !is.null(baseline_window) ) {

        stop (
            "Baseline-derived ROPEs are not implemented for 3D models; use `rope` instead.",
            call. = FALSE
            )

    }

    if (isTRUE(varying_smooth) ) {

        if (!isTRUE(varying_smooth_missing) ) {

            stop (
                "Participant-specific 3D smooths are not implemented; set varying_smooth = FALSE.",
                call. = FALSE
                )

        }

        message (
            "Spatiotemporal model: using participant-level intercepts without 3D varying smooths."
            )

    }

    return (
        list(
            kvalue = if (isTRUE(kvalue_missing) ) c(15L, 10L) else kvalue,
            bs = if (isTRUE(bs_missing) ) c("ps", "ps") else bs
            )
        )

}
