#' Print a temporal neurogam analysis
#'
#' @param x A result of class `clusters_results_1d`.
#' @param digits Number of decimal places to print.
#' @param ... Reserved for future use.
#' @return The input object, invisibly.
#' @export
print.clusters_results_1d <- function(x, digits = 3L, ...) {
    neurogam_print_report(x, digits = digits, detailed = FALSE)
}

#' Print a temporal-generalization neurogam analysis
#'
#' @param x A result of class `clusters_results_2d`.
#' @param digits Number of decimal places to print.
#' @param ... Reserved for future use.
#' @return The input object, invisibly.
#' @export
print.clusters_results_2d <- function(x, digits = 3L, ...) {
    neurogam_print_report(x, digits = digits, detailed = FALSE)
}

#' Print a spatiotemporal neurogam analysis
#'
#' @param x A result of class `clusters_results_st`.
#' @param digits Number of decimal places to print.
#' @param ... Reserved for future use.
#' @return The input object, invisibly.
#' @export
print.clusters_results_st <- function(x, digits = 3L, ...) {
    neurogam_print_report(x, digits = digits, detailed = FALSE)
}

#' Summarize a temporal neurogam analysis
#'
#' @param object A result of class `clusters_results_1d`.
#' @param digits Number of decimal places to print.
#' @param ... Reserved for future use.
#' @return The input object, invisibly.
#' @export
summary.clusters_results_1d <- function(object, digits = 3L, ...) {
    neurogam_print_report(object, digits = digits, detailed = TRUE)
}

#' Summarize a temporal-generalization neurogam analysis
#'
#' @param object A result of class `clusters_results_2d`.
#' @param digits Number of decimal places to print.
#' @param ... Reserved for future use.
#' @return The input object, invisibly.
#' @export
summary.clusters_results_2d <- function(object, digits = 3L, ...) {
    neurogam_print_report(object, digits = digits, detailed = TRUE)
}

#' Summarize a spatiotemporal neurogam analysis
#'
#' @param object A result of class `clusters_results_st`.
#' @param digits Number of decimal places to print.
#' @param ... Reserved for future use.
#' @return The input object, invisibly.
#' @export
summary.clusters_results_st <- function(object, digits = 3L, ...) {
    neurogam_print_report(object, digits = digits, detailed = TRUE)
}
