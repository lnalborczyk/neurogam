test_that("1D cluster assembly preserves legacy and standardized outputs", {
    predictions <- data.frame(
        time = seq_len(8),
        prob_ratio = c(1, 20, 20, 1, 0.02, 0.02, 1, 1)
    )

    clustering <- assemble_temporal_clusters(
        predictions = predictions,
        is_2d_time = FALSE,
        participant_clusters = FALSE,
        threshold = 10,
        threshold_type = "both"
    )

    expect_named(
        clustering,
        c("clusters", "cluster_summary", "cluster_membership")
    )
    expect_named(
        clustering$clusters,
        c("id", "onset", "offset", "n_points", "sign")
    )
    expect_equal(nrow(clustering$cluster_summary), 2L)
    expect_equal(nrow(clustering$cluster_membership), 4L)
    expect_setequal(clustering$clusters$sign, c("positive", "negative"))
    expect_equal(sum(clustering$cluster_summary$n_points), 4L)

    result <- assemble_temporal_results(
        clustering = clustering,
        predictions = predictions,
        model = NULL,
        summary_data = predictions,
        multilevel = "summary",
        is_2d_time = FALSE
    )

    expect_s3_class(result, "clusters_results_1d")
    expect_identical(result$predictions, predictions)
    expect_named(
        result,
        c(
            "clusters", "cluster_summary", "cluster_membership",
            "predictions", "model", "summary_data", "multilevel"
        )
    )
})

test_that("2D cluster assembly does not join diagonal neighbors", {
    predictions <- expand.grid(time1 = 1:2, time2 = 1:2)
    predictions$prob_ratio <- ifelse(
        predictions$time1 == predictions$time2, 20, 1
    )

    clustering <- assemble_temporal_clusters(
        predictions = predictions,
        is_2d_time = TRUE,
        participant_clusters = FALSE,
        threshold = 10,
        threshold_type = "both"
    )

    expect_equal(nrow(clustering$cluster_summary), 2L)
    expect_equal(nrow(clustering$cluster_membership), 2L)
    expect_named(
        clustering$clusters,
        c("id", "sign", "time1", "time2", "value", "n_points")
    )

    result <- assemble_temporal_results(
        clustering = clustering,
        predictions = predictions,
        model = NULL,
        summary_data = predictions,
        multilevel = "group",
        is_2d_time = TRUE
    )
    expect_s3_class(result, "clusters_results_2d")
})

test_that("participant clusters remain separated by group", {
    predictions <- expand.grid(
        time = 1:3,
        participant = c("A", "B"),
        KEEP.OUT.ATTRS = FALSE
    )
    predictions$prob_ratio <- 20

    clustering <- assemble_temporal_clusters(
        predictions = predictions,
        is_2d_time = FALSE,
        participant_clusters = TRUE,
        threshold = 10,
        threshold_type = "above"
    )

    expect_equal(nrow(clustering$cluster_summary), 2L)
    expect_equal(nrow(clustering$cluster_membership), 6L)
    expect_setequal(clustering$clusters$participant, c("A", "B"))
})

test_that("empty hits are allowed but malformed input is rejected", {
    predictions <- data.frame(time = 1:4, prob_ratio = rep(1, 4))
    clustering <- assemble_temporal_clusters(
        predictions, FALSE, FALSE, 10, "both"
    )
    expect_equal(nrow(clustering$cluster_summary), 0L)
    expect_equal(nrow(clustering$cluster_membership), 0L)

    expect_error(
        assemble_temporal_clusters(
            data.frame(time = 1:3), FALSE, FALSE, 10, "both"
        ),
        "required clustering columns"
    )
    expect_error(
        assemble_temporal_clusters(
            data.frame(time = 1:3, prob_ratio = NA_real_),
            FALSE, FALSE, 10, "both"
        ),
        "No complete posterior-odds observations"
    )
})
