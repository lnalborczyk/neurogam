test_that("P2.1 read-only views preserve 1D legacy results", {
    x <- list(
        clusters = data.frame(id = 1L, sign = "positive"),
        cluster_summary = data.frame(
            id = 1L, sign = "positive", n_points = 2L
        ),
        cluster_membership = data.frame(
            time = c(0.1, 0.2), id = c(1L, 1L),
            sign = c("positive", "positive")
        ),
        predictions = data.frame(time = c(0.1, 0.2), prob_ratio = c(20, 20)),
        model = structure(list(), class = "mock_brmsfit"),
        summary_data = data.frame(time = c(0.1, 0.2)),
        multilevel = "summary"
    )
    class(x) <- "clusters_results_1d"

    original_names <- names(x)
    original_class <- class(x)
    expect_true(validate_neurogam_result(x, check_membership = TRUE))

    view <- as_neurogam_result_view(x, check_membership = TRUE)
    expect_s3_class(view, "neurogam_result_view")
    expect_identical(view$predictions, x$predictions)
    expect_identical(view$cluster_summary, x$cluster_summary)
    expect_identical(view$cluster_membership, x$cluster_membership)
    expect_identical(view$model, x$model)
    expect_identical(view$metadata$analysis_type, "temporal")
    expect_identical(view$metadata$model_coordinates, "time")
    expect_identical(view$metadata$model_ndim, 1L)
    expect_null(view$metadata$odds_mode)

    expect_identical(names(x), original_names)
    expect_identical(class(x), original_class)
    expect_identical(x$clusters, data.frame(id = 1L, sign = "positive"))
})

test_that("P2.1 correctly identifies 2D generalized temporal results", {
    x <- list(
        predictions = data.frame(
            time1 = c(0.1, 0.1), time2 = c(0.2, 0.3),
            prob_ratio = c(1, 20)
        ),
        cluster_summary = data.frame(
            id = 7L, sign = "positive", n_points = 1L
        ),
        cluster_membership = data.frame(
            time1 = 0.1, time2 = 0.3, id = 7L, sign = "positive"
        ),
        model = NULL
    )
    class(x) <- "clusters_results_2d"

    expect_true(validate_neurogam_result(x, check_membership = TRUE))
    meta <- neurogam_result_metadata(x)
    expect_identical(meta$analysis_type, "temporal_generalization")
    expect_identical(meta$model_coordinates, c("time1", "time2"))
    expect_identical(meta$cluster_coordinates, c("time1", "time2"))
    expect_identical(meta$n_clusters, 1L)
})

test_that("P2.1 identifies spatial model dimensions and observed graph geometry", {
    x <- list(
        predictions = data.frame(
            channel = c("Fz", "Cz"), latency = c(0.1, 0.1),
            xx = c(0, 1), yy = c(0, 1), prob_ratio = c(20, 1)
        ),
        cluster_summary = data.frame(
            id = 1L, sign = "positive", n_points = 1L
        ),
        cluster_membership = data.frame(
            channel = "Fz", latency = 0.1,
            id = 1L, sign = "positive"
        ),
        model = NULL,
        input_columns = list(
            sensor_id = "channel", time_id = "latency", space_id = c("xx", "yy")
        ),
        threshold = 10,
        null_value = 0,
        rope = 0,
        adjacency_source = "auto_knn"
    )
    class(x) <- "clusters_results_st"

    meta <- neurogam_result_metadata(x)
    expect_identical(meta$analysis_type, "spatiotemporal")
    expect_identical(meta$model_coordinates, c("xx", "yy", "latency"))
    expect_identical(meta$cluster_coordinates, c("channel", "latency"))
    expect_identical(meta$model_ndim, 3L)
    expect_identical(meta$adjacency_source, "auto_knn")
    expect_identical(meta$threshold, 10)
    expect_true(validate_neurogam_result(x, check_membership = TRUE))
})

test_that("P2.1 handles empty clusters without inventing inference metadata", {
    x <- list(
        predictions = data.frame(time = 1:3, prob_ratio = rep(1, 3)),
        cluster_summary = data.frame(),
        cluster_membership = data.frame(),
        model = NULL
    )
    class(x) <- "clusters_results_1d"

    view <- as_neurogam_result_view(x, check_membership = TRUE)
    expect_equal(nrow(view$cluster_summary), 0L)
    expect_equal(nrow(view$cluster_membership), 0L)
    expect_identical(view$metadata$n_clusters, 0L)
    expect_null(view$metadata$threshold)
    expect_null(view$metadata$odds_mode)
})

test_that("P2.1 validation rejects malformed results", {
    x <- list(
        predictions = data.frame(time = 1, prob_ratio = 20),
        cluster_summary = data.frame(id = 1L, sign = "positive", n_points = 2L),
        cluster_membership = data.frame(time = 1, id = 1L, sign = "positive"),
        model = NULL
    )
    class(x) <- "clusters_results_1d"

    expect_error(
        validate_neurogam_result(x, check_membership = TRUE),
        "membership counts"
    )
    expect_error(validate_neurogam_result(unclass(x)), "must inherit")
    x$predictions$prob_ratio <- NULL
    expect_error(validate_neurogam_result(x), "prob_ratio")
})
