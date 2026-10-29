precompile_vignette <- function (name) {

    stopifnot(
        file.exists(
            file.path("vignettes-raw", paste0(name, ".qmd") )
            )
        )

    dir.create(
        "vignettes/figures",
        recursive = TRUE,
        showWarnings = FALSE
        )

    old_dir <- getwd()
    old_fig_path <- knitr::opts_chunk$get("fig.path")

    on.exit({
        setwd(old_dir)
        knitr::opts_chunk$set(fig.path = old_fig_path)
        }, add = TRUE)

    setwd("vignettes-raw")

    knitr::opts_chunk$set(
        fig.path = paste0(
            "../vignettes/figures/",
            name,
            "-"
            )
        )

    knitr::knit(
        input = paste0(name, ".qmd"),
        output = file.path(
            "..",
            "vignettes",
            paste0(name, ".qmd")
            ),
        envir = new.env(parent = globalenv() ),
        quiet = FALSE
        )

    # Adjust figure references for the new document location
    staged <- file.path(
        "..",
        "vignettes",
        paste0(name, ".qmd")
        )

    content <- readLines(staged, warn = FALSE)

    content <- gsub(
        "../vignettes/figures/",
        "figures/",
        content,
        fixed = TRUE
        )

    writeLines(content, staged, useBytes = TRUE)
    invisible(staged)

}

precompile_vignette("eeg_temporal")
precompile_vignette("temporal_2d")
precompile_vignette("eeg_spatio_temporal")
precompile_vignette("residual_autocorrelation")
