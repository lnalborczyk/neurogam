test_that("the inference core computes finite-sample capped legacy odds", {
    draws <- cbind(
        positive = c(1, 2, 3, 4),
        negative = c(-1, -2, -3, -4),
        mixed = c(-1, -1, 1, 1)
    )

    out <- summarise_posterior_matrix(
        draws,
        odds_mode = "legacy_one_sided",
        odds_cap = 4
    )

    expect_equal(out$prob_above, c(1, 0, 0.5))
    expect_equal(out$prob_ratio, c(4, 0.25, 1))
    expect_equal(out$estimate, c(2.5, -2.5, 0))
    expect_equal(nrow(out), 3L)
})

test_that("directional ROPE odds retain a neutral region", {
    draws <- cbind(
        positive = rep(2, 8),
        negative = rep(-2, 8),
        neutral = rep(0, 8),
        uncertain = c(rep(-2, 4), rep(2, 4))
    )

    out <- summarise_posterior_matrix(
        draws,
        rope = 0.5,
        odds_mode = "legacy_directional",
        odds_cap = 8
    )

    expect_equal(out$prob_ratio, c(8, 1 / 8, 1, 1))
    expect_equal(out$prob_inside[3L], 1)
    expect_equal(out$prob_above[4L], 0.5)
    expect_equal(out$prob_below[4L], 0.5)
})

test_that("spatiotemporal mode preserves uncapped odds", {

    draws <- cbind(
        rep(2, 4),
        rep(-2, 4),
        c(2, -2, 2, -2)
    )

    # Explicitly request the original uncapped strategy.
    out <- summarise_spatiotemporal_draws(
        draws,
        rope = 0,
        odds_mode = "spatiotemporal"
    )

    expect_true(is.infinite(out$prob_ratio[1L]))
    expect_equal(out$prob_ratio[2L], 0)
    expect_equal(out$prob_ratio[3L], 1)

    expect_named(out, c(
        "estimate", "lower", "upper",
        "prob_above", "prob_below", "prob_inside",
        "odds_above", "odds_below", "prob_ratio"
    ))
})

test_that("one-sample adapters preserve temporal columns and grouping", {
    df <- expand.grid(time = c(0, 0.1), participant = c("s1", "s2"), .draw = 1:6)
    df$.epred <- ifelse(df$time == 0, -2, 2)

    # one <- compute_one_sample_prob(
    #     df, 0, FALSE, 6, 0.95, baseline_window = NULL
    # )

    one <- compute_one_sample_prob(
        df,
        null_value = 0,
        participant_clusters = FALSE,
        n_post_samples = 6,
        credible_interval = 0.95,
        baseline_window = NULL,
        odds_mode = "legacy_one_sided"
    )

    expect_equal(one$prob_ratio, c(1 / 6, 6))

    # expect_named(one, c("time", "prob_above", "prob_ratio",
    #                     "post_prob", "lower", "upper"))
    expect_named(one, c(
        "time",
        "prob_above", "prob_below", "prob_inside",
        "odds_above", "odds_below", "prob_ratio",
        "post_prob", "lower", "upper"
    ))
    # expect_equal(one$prob_ratio, c(1 / 6, 6))

    by_subject <- compute_one_sample_prob(
        df, 0, TRUE, 6, 0.95, baseline_window = NULL
    )
    expect_equal(nrow(by_subject), 4L)
    expect_true("participant" %in% names(by_subject))
})

test_that("paired contrasts respect draw matching and predictor order", {
    data <- expand.grid(
        time = c(0.1, 0.2),
        predictor = c("A", "B"),
        .draw = 1:8,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    data$.epred <- ifelse(data$predictor == "B", data$.draw + 1, data$.draw)

    paired <- pair_posterior_contrasts(data, coordinates = "time")
    expect_true(all(paired$epred_diff == 1))

    # out <- compute_two_sample_prob(
    #     data,
    #     null_value = 0,
    #     participant_clusters = FALSE,
    #     n_post_samples = 8,
    #     credible_interval = 0.95,
    #     predictor_type = "categorical",
    #     baseline_window = NULL
    # )

    out <- compute_two_sample_prob(
        data,
        null_value = 0,
        participant_clusters = FALSE,
        n_post_samples = 8,
        credible_interval = 0.95,
        predictor_type = "categorical",
        baseline_window = NULL,
        odds_mode = "legacy_one_sided"
    )

    expect_equal(out$prob_ratio, c(8, 8))
    expect_equal(out$post_prob, c(1, 1))

    set.seed(123)
    mixed <- data[sample(nrow(data)), , drop = FALSE]
    paired_mixed <- pair_posterior_contrasts(mixed, coordinates = "time")
    expect_true(all(paired_mixed$epred_diff == 1))
})

test_that("2D adapters preserve one-sided odds and coordinate columns", {
    data <- expand.grid(
        time1 = c(1, 2),
        time2 = c(1, 2),
        predictor = c("A", "B"),
        .draw = 1:5,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    data$.epred <- ifelse(data$predictor == "B", 2, 0)

    # out <- compute_two_sample_prob_2d(
    #     data,
    #     null_value = 0,
    #     participant_clusters = FALSE,
    #     n_post_samples = 5,
    #     credible_interval = 0.95,
    #     predictor_type = "categorical"
    # )

    out <- compute_two_sample_prob_2d(
        data,
        null_value = 0,
        participant_clusters = FALSE,
        n_post_samples = 5,
        credible_interval = 0.95,
        predictor_type = "categorical",
        odds_mode = "legacy_one_sided"
    )

    expect_equal(nrow(out), 4L)
    expect_true(all(out$prob_ratio == 5))
    # expect_named(out, c("time1", "time2", "prob_above", "prob_ratio",
    #                     "post_prob", "lower", "upper"))
    expect_named(out, c(
        "time1", "time2",
        "prob_above", "prob_below", "prob_inside",
        "odds_above", "odds_below", "prob_ratio",
        "post_prob", "lower", "upper"
    ))
})

test_that("core rejects invalid draws and unsupported settings", {

    expect_error(
        summarise_posterior_matrix(matrix(NA_real_, 2, 2)),
        "finite"
    )

    expect_error(
        summarise_posterior_matrix(matrix(1, 2, 2), rope = -1),
        "rope"
    )

    expect_error(
        summarise_posterior_matrix(
            matrix(1, 2, 2),
            odds_mode = "unsupported"
        )
    )

    # The finite-sample odds cap defaults to the number of draws.
    # out <- summarise_posterior_matrix(
    #     matrix(1, 2, 2)
    # )
    #
    # expect_equal(out$prob_ratio, c(2, 2))
    # expect_true(all(is.finite(out$prob_ratio)))
    out <- summarise_posterior_matrix(
        matrix(1, 2, 2),
        odds_mode = "legacy_one_sided"
    )

    expect_equal(out$prob_ratio, c(2, 2))
    expect_true(all(is.finite(out$prob_ratio)))
})
