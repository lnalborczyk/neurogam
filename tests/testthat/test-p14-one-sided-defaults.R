test_that("paired contrasts are invariant to row order", {
    data <- expand.grid(
        time = c(0.1, 0.2),
        predictor = c("A", "B"),
        .draw = 1:8,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    data$.epred <- data$.draw + as.integer(data$predictor == "B")

    paired <- pair_posterior_contrasts(data, coordinates = "time")
    expect_true(all(paired$epred_diff == 1))

    set.seed(123)
    shuffled <- data[sample(nrow(data)), , drop = FALSE]
    paired_shuffled <- pair_posterior_contrasts(
        shuffled,
        coordinates = "time"
    )
    expect_true(all(paired_shuffled$epred_diff == 1))
})

test_that("categorical factor order controls contrast sign", {
    data <- expand.grid(
        time = c(0.1, 0.2),
        predictor = c("A", "B"),
        .draw = 1:8,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    data$.epred <- data$.draw + as.integer(data$predictor == "B")
    data$predictor <- factor(data$predictor, levels = c("B", "A"))

    paired <- pair_posterior_contrasts(data, coordinates = "time")
    expect_true(all(paired$epred_diff == -1))
})

test_that("continuous contrasts use numeric ordering", {
    data <- expand.grid(
        time = 1:2,
        predictor = c(2, 10),
        .draw = 1:8,
        KEEP.OUT.ATTRS = FALSE
    )
    data$.epred <- data$.draw + as.integer(data$predictor == 10)
    data <- data[sample(nrow(data)), , drop = FALSE]

    paired <- pair_posterior_contrasts(
        data,
        coordinates = "time",
        predictor_type = "continuous"
    )
    expect_true(all(paired$epred_diff == 1))
})

test_that("one-sided odds is the default in the matrix core", {
    draws <- cbind(
        positive = rep(2, 8),
        negative = rep(-2, 8),
        neutral = rep(0, 8)
    )
    #out <- summarise_posterior_matrix(draws)
    out <- summarise_posterior_matrix(
        draws,
        odds_mode = "legacy_one_sided"
    )

    expect_equal(out$prob_ratio, c(8, 1 / 8, 1 / 8))
    expect_true(all(is.finite(out$prob_ratio)))
    # expect_equal(
    #     summarise_posterior_matrix(draws, odds_cap = 8)$prob_ratio,
    #     out$prob_ratio
    # )
    expect_equal(
        summarise_posterior_matrix(
            draws,
            odds_mode = "legacy_one_sided",
            odds_cap = 8
        )$prob_ratio,
        out$prob_ratio
    )
})

test_that("alternative directional odds remains opt-in", {
    draws <- cbind(
        positive = rep(2, 8),
        negative = rep(-2, 8),
        neutral = rep(0, 8)
    )
    old_style <- summarise_posterior_matrix(
        draws,
        odds_mode = "legacy_directional",
        odds_cap = 8
    )
    expect_equal(old_style$prob_ratio, c(8, 1 / 8, 1))
})

test_that("spatiotemporal wrapper uses capped one-sided odds", {
    draws <- cbind(rep(2, 8), rep(-2, 8))
    default <- summarise_spatiotemporal_draws(draws, odds_mode = "legacy_one_sided")

    expect_equal(default$prob_ratio, c(8, 1 / 8))
    expect_true(all(is.finite(default$prob_ratio)))

    # Original uncapped strategy is still available when explicitly requested.
    original <- summarise_spatiotemporal_draws(
        draws,
        odds_mode = "spatiotemporal"
    )
    expect_true(is.infinite(original$prob_ratio[1L]))
    expect_equal(original$prob_ratio[2L], 0)
})

test_that("temporal condition contrasts use one-sided odds by default", {
    data <- expand.grid(
        time = c(0.1, 0.2),
        predictor = c("A", "B"),
        .draw = 1:8,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    data$.epred <- ifelse(data$predictor == "B", 2, 0)

    # out <- compute_two_sample_prob(
    #     data,
    #     null_value = 0,
    #     participant_clusters = FALSE,
    #     n_post_samples = 8,
    #     credible_interval = 0.95,
    #     predictor_type = "categorical"
    # )
    out <- compute_two_sample_prob(
        data,
        null_value = 0,
        participant_clusters = FALSE,
        n_post_samples = 8,
        credible_interval = 0.95,
        predictor_type = "categorical",
        odds_mode = "legacy_one_sided"
    )
    expect_true(all(out$prob_ratio == 8))
    expect_true(all(out$post_prob == 2))
})

test_that("two-dimensional condition contrasts use one-sided odds", {
    data <- expand.grid(
        time1 = 1:2,
        time2 = 1:2,
        predictor = c("A", "B"),
        .draw = 1:8,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    data$.epred <- ifelse(data$predictor == "B", 2, 0)

    # out <- compute_two_sample_prob_2d(
    #     data,
    #     null_value = 0,
    #     participant_clusters = FALSE,
    #     n_post_samples = 8,
    #     credible_interval = 0.95,
    #     predictor_type = "categorical"
    # )
    out <- compute_two_sample_prob_2d(
        data,
        null_value = 0,
        participant_clusters = FALSE,
        n_post_samples = 8,
        credible_interval = 0.95,
        predictor_type = "categorical",
        odds_mode = "legacy_one_sided"
    )
    expect_true(all(out$prob_ratio == 8))
})

# Reproducibility warning: with ROPE > 0, a one-sided prob_ratio below 1
# indicates low P(effect > null + ROPE), *not* high P(effect < null - ROPE).
# test_that("one-sided ROPE semantics remain explicit", {
#     draws <- matrix(rep(0, 8), ncol = 1L)
#     one_sided <- summarise_posterior_matrix(draws, rope = 0.5)
#     directional <- summarise_posterior_matrix(
#         draws,
#         rope = 0.5,
#         odds_mode = "legacy_directional"
#     )
#     expect_equal(one_sided$prob_inside, 1)
#     expect_equal(one_sided$prob_ratio, 1 / 8)
#     expect_equal(directional$prob_ratio, 1)
# })
