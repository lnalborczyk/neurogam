test_that("spatiotemporal formula uses a 2D spatial and 1D temporal marginal", {
    f <- make_bgam_formula(
        family = gaussian(),
        multilevel = "summary",
        predictor_type = "none",
        time_id = "time",
        space_id = c("xproj", "yproj"),
        tensor_d = c(2, 1),
        kvalue = c(15, 10),
        bs = c("ps", "ps"),
        varying_smooth = FALSE,
        use_se = FALSE
    )
    text <- paste(deparse(f$formula), collapse = " ")
    expect_match(text, "t2\\(xproj, yproj, time", fixed = FALSE)
    expect_match(text, "d = c\\(2, 1\\)")
    expect_match(text, "k = c\\(15, 10\\)")
    expect_match(text, "\\(1 \\| participant\\)")
})

test_that("categorical spatiotemporal formula contains condition smooth", {
    f <- make_bgam_formula(
        family = gaussian(),
        multilevel = "summary",
        predictor_type = "categorical",
        within_between = "within-subject",
        time_id = "time",
        space_id = c("xproj", "yproj"),
        kvalue = c(15, 10),
        bs = c("ps", "ps"),
        varying_smooth = FALSE,
        use_se = TRUE
    )
    text <- paste(deparse(f$formula), collapse = " ")
    expect_match(text, "by = predictor", fixed = TRUE)
    #expect_match(text, "1 + predictor | participant", fixed = TRUE)
    expect_match(
        gsub("\\s+", " ", text),
        "\\(1 \\+ predictor \\| participant\\)"
    )
    expect_match(text, "se(outcome_sem", fixed = TRUE)
})

test_that("unsupported spatiotemporal configurations error explicitly", {
    basic <- function(...) make_bgam_formula(
        family = gaussian(),
        multilevel = "summary", predictor_type = "none",
        space_id = c("xproj", "yproj"),
        varying_smooth = FALSE,
        ...
    )
    expect_error(basic(tensor_d = c(1, 1, 1)), "tensor_d")
    expect_error(basic(kvalue = c(15, 10, 5)), "kvalue")
    expect_error(basic(include_ar_term = TRUE), "AR")
    expect_error(make_bgam_formula(
        family = gaussian(),
        multilevel = "summary", predictor_type = "none",
        space_id = c("xproj", "yproj"),
        varying_smooth = TRUE
    ), "varying_smooth")
})

test_that("spatiotemporal feature is opt-in", {
    f <- make_bgam_formula(
        family = gaussian(), multilevel = "summary",
        predictor_type = "none", time_id = "time",
        varying_smooth = FALSE
    )
    expect_match(paste(deparse(f$formula), collapse = " "), "s(time", fixed = TRUE)
})
