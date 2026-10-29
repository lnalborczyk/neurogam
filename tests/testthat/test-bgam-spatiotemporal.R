test_that("posterior summaries recover one-sided odds and ROPE behavior", {

    draws <- cbind(
        rep(2, 40),
        rep(-2, 40),
        rep(c(-1, 1), each = 20),
        rep(0, 40)
    )

    actual <- summarise_spatiotemporal_draws(
        draws,
        null_value = 0,
        rope = 0.5
    )

    # Posterior probabilities are unchanged.
    expect_equal(actual$prob_above, c(1, 0, 0.5, 0))
    expect_equal(actual$prob_below, c(0, 1, 0.5, 0))
    expect_equal(actual$prob_inside, c(0, 0, 0, 1))

    # Default: finite-sample-capped one-sided odds.
    # expect_equal(
    #     actual$prob_ratio,
    #     c(40, 1 / 40, 1, 1 / 40)
    # )

    # New default: continuous directional posterior odds.
    expect_equal(actual$odds_above, c(Inf, 0, 1, 0))
    expect_equal(actual$odds_below, c(0, Inf, 1, 0))
    expect_equal(actual$prob_ratio, c(Inf, 0, 1, 0))

    # The historical method remains available explicitly.
    legacy <- summarise_spatiotemporal_draws(
        draws,
        null_value = 0,
        rope = 0.5,
        odds_mode = "legacy_one_sided"
    )

    expect_equal(
        legacy$prob_ratio,
        c(40, 1 / 40, 1, 1 / 40)
    )

    expect_equal(actual$estimate[1:2], c(2, -2))

    # Historical spatiotemporal behavior remains available.
    historical <- summarise_spatiotemporal_draws(
        draws,
        null_value = 0,
        rope = 0.5,
        odds_mode = "spatiotemporal"
    )

    expect_true(is.infinite(historical$prob_ratio[1]))
    expect_equal(historical$prob_ratio[2], 0)
    expect_equal(historical$prob_ratio[3:4], c(1, 1))
})

test_that("nonzero null values are applied to posterior locations", {
    draws <- cbind(rep(5, 20), rep(3, 20))
    actual <- summarise_spatiotemporal_draws(
        draws,
        null_value = 4,
        rope = 0.25
    )

    expect_equal(actual$prob_above, c(1, 0))
    expect_equal(actual$prob_below, c(0, 1))
})

test_that("paired differences preserve posterior draw dependence", {
    set.seed(12)
    shared <- rnorm(100)
    condition_a <- cbind(shared + 2, shared - 3)
    condition_b <- cbind(shared, shared)

    difference <- condition_a - condition_b
    actual <- summarise_spatiotemporal_draws(difference)

    expect_equal(actual$estimate, c(2, -3), tolerance = 1e-12)
    expect_equal(actual$prob_above, c(1, 0))
    expect_equal(actual$prob_below, c(0, 1))
    expect_equal(actual$lower, c(2, -3), tolerance = 1e-12)
    expect_equal(actual$upper, c(2, -3), tolerance = 1e-12)
})

test_that("posterior summary helper rejects invalid input", {
    expect_error(summarise_spatiotemporal_draws(1:10), "matrix")
    expect_error(
        summarise_spatiotemporal_draws(matrix(1, 2, 2), rope = -0.5),
        "nonnegative"
    )
    expect_error(
        summarise_spatiotemporal_draws(
            matrix(1, 2, 2), credible_interval = 1
        ),
        "strictly between"
    )
})
