# P2.2: fast, deterministic result-metadata contracts; no Stan required.

make_p22_temporal_fixture <- function(two_dimensional = FALSE) {
    if (two_dimensional) {
        predictions <- data.frame(
            time1 = c(0.1, 0.1), time2 = c(0.2, 0.3),
            prob_ratio = c(20, 1)
        )
        membership <- data.frame(
            time1 = 0.1, time2 = 0.2, id = 1L, sign = "positive"
        )
        type <- "clusters_results_2d"
    } else {
        predictions <- data.frame(time = c(0.1, 0.2), prob_ratio = c(20, 1))
        membership <- data.frame(time = 0.1, id = 1L, sign = "positive")
        type <- "clusters_results_1d"
    }
    x <- list(
        clusters = data.frame(id = 1L),
        cluster_summary = data.frame(
            id = 1L, sign = "positive", n_points = 1L
        ),
        cluster_membership = membership,
        predictions = predictions,
        model = NULL,
        summary_data = predictions,
        multilevel = "summary"
    )
    class(x) <- type
    x
}

test_that("P2.2 records 1D metadata without changing fields or class", {
    x <- make_p22_temporal_fixture()
    previous_names <- names(x)
    previous_class <- class(x)
    previous_predictions <- x$predictions

    meta <- make_neurogam_metadata(
        analysis_type = "temporal",
        time_id = "time",
        participant_id = "subject",
        outcome_id = "voltage",
        predictor_id = "condition",
        multilevel = "summary",
        threshold = 10,
        threshold_type = "both",
        null_value = 0,
        baseline_window = c(-0.2, 0),
        credible_interval = 0.95,
        odds_mode = "legacy_one_sided",
        contrast = "condition_difference",
        n_post_samples = 500
    )

    annotated <- set_neurogam_result_metadata(x, meta)
    expect_identical(names(annotated), previous_names)
    expect_identical(class(annotated), previous_class)
    expect_identical(annotated$predictions, previous_predictions)
    expect_null(attr(x, "neurogam_metadata", exact = TRUE))
    expect_true(validate_neurogam_result(annotated, check_membership = TRUE))

    view <- as_neurogam_result_view(annotated)
    expect_identical(view$metadata$analysis_type, "temporal")
    expect_identical(view$metadata$odds_mode, "legacy_one_sided")
    expect_identical(view$metadata$baseline_window, c(-0.2, 0))
    expect_identical(view$metadata$participant_id, "subject")
    expect_identical(view$metadata$threshold, 10)
    expect_null(view$metadata$contrast_levels)
})

test_that("P2.2 annotates 2D constructor without modifying its output fields", {
    source <- make_p22_temporal_fixture(two_dimensional = TRUE)
    cluster <- list(
        clusters = source$clusters,
        cluster_summary = source$cluster_summary,
        cluster_membership = source$cluster_membership
    )
    metadata <- make_neurogam_metadata(
        analysis_type = "temporal_generalization",
        time_id = c("time1", "time2"),
        threshold = 20,
        odds_mode = "legacy_one_sided"
    )
    result <- assemble_temporal_results(
        clustering = cluster,
        predictions = source$predictions,
        model = NULL,
        summary_data = source$summary_data,
        multilevel = "summary",
        is_2d_time = TRUE,
        metadata = metadata
    )

    expect_identical(names(result), names(source))
    expect_identical(class(result), class(source))
    expect_identical(neurogam_result_metadata(result)$model_ndim, 2L)
    expect_identical(
        neurogam_result_metadata(result)$cluster_coordinates,
        c("time1", "time2")
    )
    expect_identical(neurogam_result_metadata(result)$threshold, 20)
})

test_that("P2.2 records separate model and cluster geometry in 3D", {
    x <- list(
        predictions = data.frame(
            channel = "Fz", latency = 0.1, xx = 0, yy = 1,
            prob_ratio = 20
        ),
        cluster_summary = data.frame(
            id = 1L, sign = "positive", n_points = 1L
        ),
        cluster_membership = data.frame(
            channel = "Fz", latency = 0.1, id = 1L, sign = "positive"
        ),
        model = NULL
    )
    class(x) <- "clusters_results_st"
    metadata <- make_neurogam_metadata(
        analysis_type = "spatiotemporal",
        time_id = "latency",
        space_id = c("xx", "yy"),
        sensor_id = "channel",
        threshold = 10,
        null_value = 0,
        rope = 0,
        odds_mode = "legacy_one_sided",
        contrast = "A_minus_B",
        contrast_levels = c("A", "B"),
        adjacency_source = "auto_knn"
    )
    out <- set_neurogam_result_metadata(x, metadata)
    observed <- neurogam_result_metadata(out)
    expect_identical(observed$model_coordinates, c("xx", "yy", "latency"))
    expect_identical(observed$cluster_coordinates, c("channel", "latency"))
    expect_identical(observed$contrast_levels, c("A", "B"))
    expect_identical(observed$adjacency_source, "auto_knn")
    expect_identical(observed$model_ndim, 3L)
    expect_identical(class(out), "clusters_results_st")
    expect_true(validate_neurogam_result(out, check_membership = TRUE))
})

test_that("P2.2 metadata survive saveRDS and keep counts live", {
    x <- make_p22_temporal_fixture()
    x <- set_neurogam_result_metadata(x, make_neurogam_metadata(
        analysis_type = "temporal", time_id = "time", threshold = 10
    ))
    path <- tempfile(fileext = ".rds")
    on.exit(unlink(path), add = TRUE)
    saveRDS(x, path)
    restored <- readRDS(path)
    expect_identical(
        neurogam_result_metadata(restored)$threshold,
        10
    )
    restored$cluster_summary <- restored$cluster_summary[FALSE, ]
    restored$cluster_membership <- restored$cluster_membership[FALSE, ]
    expect_identical(neurogam_result_metadata(restored)$n_clusters, 0L)
    expect_identical(neurogam_result_metadata(restored)$n_cluster_members, 0L)
})

test_that("P2.2 does not invent settings on older result objects", {
    old <- make_p22_temporal_fixture()
    meta <- neurogam_result_metadata(old)
    expect_null(meta$odds_mode)
    expect_null(meta$threshold)
    expect_null(meta$baseline_window)
    expect_identical(meta$analysis_type, "temporal")
})

test_that("P2.2 rejects mismatched analysis metadata", {
    x <- make_p22_temporal_fixture()
    meta <- make_neurogam_metadata(
        analysis_type = "temporal_generalization",
        time_id = c("time1", "time2")
    )
    expect_error(
        set_neurogam_result_metadata(x, meta),
        "conflicts"
    )
    expect_error(
        make_neurogam_metadata("spatiotemporal", time_id = "time"),
        "two spatial"
    )
    expect_error(
        make_neurogam_metadata("temporal", time_id = "time", threshold = 1),
        "greater than 1"
    )
})
