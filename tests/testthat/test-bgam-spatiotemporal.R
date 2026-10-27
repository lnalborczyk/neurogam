test_that("posterior draw summaries recover directional odds and ROPE behavior", {
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

    expect_equal(actual$prob_above, c(1, 0, 0.5, 0))
    expect_equal(actual$prob_below, c(0, 1, 0.5, 0))
    expect_equal(actual$prob_inside, c(0, 0, 0, 1))
    expect_true(is.infinite(actual$prob_ratio[1]))
    expect_equal(actual$prob_ratio[2], 0)
    expect_equal(actual$prob_ratio[3:4], c(1, 1))
    expect_equal(actual$estimate[1:2], c(2, -2))
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
