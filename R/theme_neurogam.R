#' neurogam ggplot2 theme
#'
#' Default ggplot2 theme used by neurogam figures.
#'
#' @param base_size Base font size.
#' @param base_family Base font family.
#'
#' @return A ggplot2 theme.
#'
#' @export

theme_neurogam <- function (base_size = 12, base_family = "") {

    ggplot2::theme_bw(
        base_size = base_size,
        base_family = base_family
        ) +
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major.x = ggplot2::element_blank(),
            panel.grid.major.y = ggplot2::element_blank()
            )

}
