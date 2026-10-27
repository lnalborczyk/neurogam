test_that("directional_rope retains continuous probabilities and odds", {
    draws <- cbind(
        positive = rep(2, 8),
        negative = rep(-2, 8),
        inside = rep(0, 8),
        split = c(rep(-2, 4), rep(2, 4))
    )

    result <- summarise_posterior_matrix(
        draws,
        rope = 0.5,
        odds_mode = "directional_rope"
    )

    expect_equal(result$prob_above, c(1, 0, 0, 0.5))
    expect_equal(result$prob_below, c(0, 1, 0, 0.5))
    expect_equal(result$prob_inside, c(0, 0, 1, 0))
    expect_equal(result$odds_above, c(Inf, 0, 0, 1))
    expect_equal(result$odds_below, c(0, Inf, 0, 1))

    score <- directional_rope_cluster_score(
        result$odds_above,
        result$odds_below,
        threshold = 10
    )

    expect_equal(score, c(Inf, 0, 1, 1))
})

test_that("directional_rope clustering distinguishes negative from ROPE", {
    predictions <- data.frame(
        time = 1:6,
        odds_above = c(Inf, Inf, 0, 0, 1, 1),
        odds_below = c(0, 0, 0, Inf, 1, 1),
        prob_ratio = c(Inf, Inf, 0, 0, 1, 1)
    )

    result <- assemble_temporal_clusters(
        predictions = predictions,
        is_2d_time = FALSE,
        participant_clusters = FALSE,
        threshold = 10,
        threshold_type = "both",
        odds_mode = "directional_rope"
    )

    expect_equal(nrow(result$cluster_summary), 2L)
    expect_equal(nrow(result$cluster_membership), 3L)
    expect_setequal(result$cluster_summary$sign, c("positive", "negative"))
    expect_false(3L %in% result$cluster_membership$time)
    expect_false(5L %in% result$cluster_membership$time)
    expect_false(6L %in% result$cluster_membership$time)
})

test_that("directional_rope requires both directional odds columns", {
    data <- data.frame(time = 1:3, prob_ratio = c(20, 1, 0.05))
    expect_error(
        assemble_temporal_clusters(
            data, FALSE, FALSE, 10, "both", odds_mode = "directional_rope"
        ),
        "(directional|odds|required|columns)"
    )

    legacy <- assemble_temporal_clusters(
        data, FALSE, FALSE, 10, "both", odds_mode = "legacy_one_sided"
    )
    expect_equal(nrow(legacy$cluster_summary), 2L)
})
