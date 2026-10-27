test_that("temporal formula parts retain response conventions", {
    expect_identical(
        bgam_formula_response_lhs("gaussian", "summary", TRUE, FALSE),
        "outcome_mean | se(outcome_sem, sigma = TRUE)"
    )
    expect_identical(
        bgam_formula_response_lhs("gaussian", "summary", FALSE, FALSE),
        "outcome_mean"
    )
    expect_identical(
        bgam_formula_response_lhs("gaussian", "summary", FALSE, TRUE),
        "outcome_mean | se(outcome_sem, sigma = TRUE)"
    )
    expect_identical(
        bgam_formula_response_lhs("gaussian", "group", TRUE, FALSE),
        "outcome_mean"
    )
    expect_identical(
        bgam_formula_response_lhs("binomial", "summary", FALSE, FALSE),
        "success | trials(trials)"
    )
})

test_that("temporal smooth construction retains original basis syntax", {
    one_dim <- bgam_formula_temporal_smooth_terms(
        is_2d_time = FALSE,
        predictor_type = "categorical",
        within_between = "within-subject",
        kvalue = 20,
        bs = "tp",
        t2_full = FALSE,
        include_by_smooth = TRUE,
        include_by_varying_smooth = TRUE
    )

    expect_identical(
        one_dim$smooth,
        "s(time, bs = 'tp', k = 20, by = predictor)"
    )
    expect_identical(
        one_dim$varying,
        "s(time, participant, bs = 'fs', m = 1, k = 10, by = predictor)"
    )

    two_dim <- bgam_formula_temporal_smooth_terms(
        is_2d_time = TRUE,
        predictor_type = "categorical",
        within_between = "between-subject",
        kvalue = 16,
        bs = "tp",
        t2_full = TRUE,
        include_by_smooth = TRUE,
        include_by_varying_smooth = FALSE
    )

    expect_identical(
        two_dim$smooth,
        "t2(time1, time2, by = predictor, bs = 'tp', k = c(16, 16), full = TRUE)"
    )
    expect_identical(
        two_dim$varying,
        "t2(time1, time2, participant, bs = 'fs', k = 8)"
    )
})

test_that("temporal RHS preserves existing ordering and random effects", {
    expect_identical(
        bgam_formula_temporal_rhs(
            multilevel = "summary",
            predictor_type = "none",
            within_between = NA_character_,
            varying_smooth = TRUE,
            smooth = "s(time)",
            varying = "s(time, participant)"
        ),
        "1 + (1 | participant) + s(time) + s(time, participant)"
    )

    expect_identical(
        bgam_formula_temporal_rhs(
            multilevel = "summary",
            predictor_type = "categorical",
            within_between = "within-subject",
            varying_smooth = TRUE,
            smooth = "s(time, by = predictor)",
            varying = "s(time, participant)"
        ),
        "1 + predictor + s(time, by = predictor) + s(time, participant)"
    )

    expect_identical(
        bgam_formula_temporal_rhs(
            multilevel = "summary",
            predictor_type = "categorical",
            within_between = "within-subject",
            varying_smooth = FALSE,
            smooth = "s(time, by = predictor)",
            varying = "s(time, participant)"
        ),
        "1 + predictor + (1 + predictor | participant) + s(time, by = predictor)"
    )

    expect_identical(
        bgam_formula_temporal_rhs(
            multilevel = "group",
            predictor_type = "categorical",
            within_between = "between-subject",
            varying_smooth = FALSE,
            smooth = "s(time)",
            varying = "unused"
        ),
        "1 + predictor + s(time)"
    )
})

test_that("public formula builder still constructs 1D/2D/3D models", {
    one <- make_bgam_formula(
        family = stats::gaussian(),
        multilevel = "summary",
        predictor_type = "none",
        varying_smooth = FALSE,
        use_se = FALSE
    )
    txt_one <- paste(deparse(one$formula), collapse = " ")
    expect_match(txt_one, "s\\(time,", fixed = FALSE)
    expect_match(txt_one, "participant", fixed = TRUE)

    two <- make_bgam_formula(
        family = stats::gaussian(),
        multilevel = "group",
        predictor_type = "categorical",
        within_between = "between-subject",
        time_id = c("time1", "time2"),
        varying_smooth = FALSE,
        use_se = FALSE
    )
    txt_two <- paste(deparse(two$formula), collapse = " ")
    expect_match(txt_two, "t2\\(time1, time2", fixed = FALSE)

    three <- make_bgam_formula(
        family = stats::gaussian(),
        multilevel = "summary",
        predictor_type = "none",
        time_id = "time",
        space_id = c("xproj", "yproj"),
        kvalue = c(15, 10),
        bs = c("ps", "ps"),
        varying_smooth = FALSE,
        use_se = FALSE
    )
    txt_three <- paste(deparse(three$formula), collapse = " ")
    expect_match(txt_three, "t2\\(xproj, yproj, time", fixed = FALSE)
})

test_that("dimension dispatch and ST default normalization are stable", {
    expect_identical(classify_bgam_analysis(), "temporal")
    expect_identical(classify_bgam_analysis(c("xproj", "yproj")), "spatiotemporal")

    defaults <- suppressMessages(normalize_spatiotemporal_dispatch(
        sensor_id = "channel",
        baseline_window = c(0, 0.1),
        varying_smooth = TRUE,
        varying_smooth_missing = TRUE,
        kvalue = 20,
        kvalue_missing = TRUE,
        bs = "tp",
        bs_missing = TRUE
    ))
    expect_identical(defaults$kvalue, c(15L, 10L))
    expect_identical(defaults$bs, c("ps", "ps"))

    explicit <- normalize_spatiotemporal_dispatch(
        sensor_id = "channel",
        baseline_window = NULL,
        varying_smooth = FALSE,
        varying_smooth_missing = FALSE,
        kvalue = c(12, 8),
        kvalue_missing = FALSE,
        bs = c("tp", "cr"),
        bs_missing = FALSE
    )
    expect_identical(explicit$kvalue, c(12, 8))
    expect_identical(explicit$bs, c("tp", "cr"))

    expect_error(
        normalize_spatiotemporal_dispatch(
            sensor_id = NULL,
            baseline_window = NULL,
            varying_smooth = FALSE,
            varying_smooth_missing = FALSE,
            kvalue = 20,
            kvalue_missing = FALSE,
            bs = "tp",
            bs_missing = FALSE
        ),
        "sensor_id"
    )
    expect_error(
        normalize_spatiotemporal_dispatch(
            sensor_id = "channel",
            baseline_window = NULL,
            varying_smooth = TRUE,
            varying_smooth_missing = FALSE,
            kvalue = 20,
            kvalue_missing = FALSE,
            bs = "tp",
            bs_missing = FALSE
        ),
        "varying_smooth"
    )
})
