make_p23_result <- function(type = c("temporal", "temporal_generalization", "spatiotemporal"),
                            empty = FALSE, recorded = TRUE) {
    type <- match.arg(type)
    if (identical(type, "temporal")) {
        predictions <- data.frame(time = c(0.1, 0.2, 0.3), prob_ratio = c(20, 20, 1))
        clusters <- data.frame(
            id = c(1L, 1L), sign = c("positive", "negative"),
            n_points = c(2L, 1L), onset = c(0.1, 0.3),
            offset = c(0.2, 0.3)
        )
        members <- data.frame(
            id = c(1L, 1L, 1L),
            sign = c("positive", "positive", "negative"),
            time = c(0.1, 0.2, 0.3)
        )
        legacy <- clusters
        class_name <- "clusters_results_1d"
    } else if (identical(type, "temporal_generalization")) {
        predictions <- data.frame(
            time1 = c(0.1, 0.2), time2 = c(0.1, 0.2),
            prob_ratio = c(20, 0.05)
        )
        clusters <- data.frame(
            id = c(1L, 1L), sign = c("positive", "negative"),
            n_points = c(1L, 1L), time1_min = c(0.1, 0.2),
            time1_max = c(0.1, 0.2), time2_min = c(0.1, 0.2),
            time2_max = c(0.1, 0.2)
        )
        members <- data.frame(
            id = c(1L, 1L), sign = c("positive", "negative"),
            time1 = c(0.1, 0.2), time2 = c(0.1, 0.2)
        )
        legacy <- members
        class_name <- "clusters_results_2d"
    } else {
        predictions <- data.frame(
            channel = c("Fz", "Cz"), time = c(0.1, 0.2),
            xproj = c(0, 0.5), yproj = c(1, 0),
            prob_ratio = c(20, 0.05)
        )
        clusters <- data.frame(
            id = c(1L, 2L), sign = c("positive", "negative"),
            n_points = c(1L, 1L), time_min = c(0.1, 0.2),
            time_max = c(0.1, 0.2), n_sensors = c(1L, 1L)
        )
        members <- data.frame(
            channel = c("Fz", "Cz"), time = c(0.1, 0.2),
            id = c(1L, 2L), sign = c("positive", "negative")
        )
        legacy <- NULL
        class_name <- "clusters_results_st"
    }

    if (empty) {
        clusters <- clusters[FALSE, , drop = FALSE]
        members <- members[FALSE, , drop = FALSE]
        if (!is.null(legacy)) {
            legacy <- legacy[FALSE, , drop = FALSE]
        }
    }

    x <- list(
        clusters = legacy,
        predictions = predictions,
        cluster_summary = clusters,
        cluster_membership = members,
        model = structure(list(), class = "mock_brmsfit"),
        multilevel = "summary"
    )
    class(x) <- class_name

    if (recorded) {
        metadata <- if (identical(type, "spatiotemporal")) {
            make_neurogam_metadata(
                type, time_id = "time", space_id = c("xproj", "yproj"),
                sensor_id = "channel", threshold = 10,
                threshold_type = "both", odds_mode = "legacy_one_sided",
                adjacency_source = "auto_knn", n_post_samples = 500
            )
        } else {
            make_neurogam_metadata(
                type,
                time_id = if (identical(type, "temporal")) "time" else c("time1", "time2"),
                threshold = 10, threshold_type = "both",
                odds_mode = "legacy_one_sided", n_post_samples = 500
            )
        }
        x <- set_neurogam_result_metadata(x, metadata)
    }

    x
}

test_that("P2.3 print and summary handle all three classes", {
    for (type in c("temporal", "temporal_generalization", "spatiotemporal")) {
        x <- make_p23_result(type)
        old_names <- names(x)
        old_class <- class(x)
        old_meta <- attr(x, "neurogam_metadata", exact = TRUE)

        printed <- capture.output(returned <- print(x))
        detailed <- capture.output(summarized <- summary(x))

        expect_identical(returned, x)
        expect_identical(summarized, x)
        expect_identical(names(x), old_names)
        expect_identical(class(x), old_class)
        expect_identical(attr(x, "neurogam_metadata", exact = TRUE), old_meta)
        expect_true(any(grepl("Clusters found: 2", printed, fixed = TRUE)))
        expect_true(any(grepl("legacy_one_sided", detailed, fixed = TRUE)))
        expect_true(any(grepl("Cluster statistics", detailed, fixed = TRUE)))
        expect_true(any(grepl("Prediction locations", detailed, fixed = TRUE)))
    }
})

test_that("P2.3 uses cluster summary rows, not arbitrary component IDs", {
    x <- make_p23_result("temporal_generalization")
    x$cluster_summary$id <- c(1L, 1L)
    x$cluster_summary$sign <- c("positive", "negative")
    printed <- capture.output(print(x))
    expect_true(any(grepl("Clusters found: 2", printed, fixed = TRUE)))
})

test_that("P2.3 is robust to no clusters in each dimension", {
    for (type in c("temporal", "temporal_generalization", "spatiotemporal")) {
        x <- make_p23_result(type, empty = TRUE)
        expect_output(print(x), "Clusters found: 0")
        expect_output(summary(x), "No clusters exceeded")
    }
})

test_that("P2.3 does not infer unknown historical settings", {
    x <- make_p23_result("temporal", recorded = FALSE)
    printed <- capture.output(print(x))
    expect_true(any(grepl("Posterior odds mode: not recorded", printed, fixed = TRUE)))
    expect_false(any(grepl("legacy_one_sided", printed, fixed = TRUE)))
})

test_that("P2.3 report survives malformed or absent model objects", {
    x <- make_p23_result("spatiotemporal")
    x$model <- NULL
    expect_output(summary(x), "Clusters found: 2")

    x$model <- structure(list(), class = "not_a_brmsfit")
    expect_output(summary(x), "Backend: not_a_brmsfit")
})

test_that("P2.3 validates report options", {
    x <- make_p23_result("temporal")
    expect_error(neurogam_print_report(x, digits = -1), "digits")
    expect_error(neurogam_print_report(x, max_clusters = 0), "max_clusters")
    expect_error(neurogam_print_report(x, max_clusters = NA_real_), "max_clusters")
    expect_output(neurogam_print_report(x, max_clusters = 1), "additional cluster")
    expect_output(neurogam_print_report(x, max_clusters = Inf), "Clusters found: 2")
})

test_that("P2.3 S3 methods are registered", {
    for (suffix in c("1d", "2d", "st")) {
        cls <- paste0("clusters_results_", suffix)
        expect_true(is.function(utils::getS3method("print", cls, optional = TRUE)))
        expect_true(is.function(utils::getS3method("summary", cls, optional = TRUE)))
    }
})
