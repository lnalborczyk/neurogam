make_p24_result <- function(type = c("temporal", "temporal_generalization"),
                            empty = FALSE, with_model = FALSE) {
    type <- match.arg(type)
    if (identical(type, "temporal")) {
        predictions <- data.frame(
            time = c(0.1, 0.2, 0.3), post_prob = c(0, 0.3, 0.5),
            lower = c(-0.1, 0.1, 0.2), upper = c(0.1, 0.5, 0.8),
            prob_ratio = c(1, 10, 20)
        )
        clusters <- data.frame(
            id = 1L, sign = "positive", n_points = 2L,
            onset = 0.2, offset = 0.3
        )
        members <- data.frame(
            id = rep(1L, 2), sign = rep("positive", 2),
            time = c(0.2, 0.3)
        )
        class_name <- "clusters_results_1d"
    } else {
        predictions <- expand.grid(
            time1 = c(0.1, 0.2), time2 = c(0.1, 0.2)
        )
        predictions$post_prob <- c(0.1, 0.2, 0.3, 0.4)
        predictions$prob_ratio <- c(1, 20, 20, 1)
        clusters <- data.frame(
            id = 1L, sign = "positive", n_points = 2L,
            time1_min = 0.1, time1_max = 0.2,
            time2_min = 0.1, time2_max = 0.1
        )
        members <- data.frame(
            id = rep(1L, 2), sign = rep("positive", 2),
            time1 = c(0.1, 0.2), time2 = c(0.1, 0.1)
        )
        class_name <- "clusters_results_2d"
    }
    if (empty) {
        clusters <- clusters[FALSE, , drop = FALSE]
        members <- members[FALSE, , drop = FALSE]
    }
    result <- list(
        predictions = predictions,
        cluster_summary = clusters,
        cluster_membership = members,
        clusters = members
    )
    if (with_model) {
        result$model <- structure(list(), class = "mock_fit")
    }
    class(result) <- class_name
    result
}

test_that("P2.4 typed temporal plots work without a model", {
    x <- make_p24_result("temporal")
    before <- x

    predictions <- plot(x, type = "predictions")
    clusters <- plot(x, type = "clusters")
    expect_s3_class(predictions, "ggplot")
    expect_s3_class(clusters, "ggplot")
    expect_equal(nrow(predictions$data), nrow(x$predictions))
    expect_equal(nrow(clusters$data), nrow(x$cluster_summary))
    expect_identical(x, before)

    expect_error(plot(x, type = "unexpected"), "Unsupported `type`")
    expect_error(
        plot(x, type = "predictions", surprising_argument = 1),
        "Unused plotting argument"
    )
})

test_that("P2.4 typed cross-temporal plots work with empty clusters", {
    x <- make_p24_result("temporal_generalization")
    expect_s3_class(plot(x, type = "predictions"), "ggplot")
    expect_s3_class(plot(x, type = "clusters"), "ggplot")

    empty <- make_p24_result("temporal_generalization", empty = TRUE)
    p <- plot(empty, type = "clusters")
    expect_s3_class(p, "ggplot")
    expect_equal(nrow(p$data), 0L)
    expect_error(plot(x, type = "predictions", palette = "not_a_palette"),
                 "scico palette")
})

test_that("P2.4 handles an empty 1D cluster table", {
    x <- make_p24_result("temporal", empty = TRUE)
    p <- plot(x, type = "clusters")
    expect_s3_class(p, "ggplot")
    expect_equal(nrow(p$data), 0L)
})

test_that("P2.4 requires saved estimates for predictions plots", {
    x <- make_p24_result("temporal")
    x$predictions$post_prob <- NULL
    expect_error(plot(x, type = "predictions"), "post_prob.*estimate")
})

test_that("P2.4 spatiotemporal dispatch distinguishes prediction and clusters", {
    prediction <- neurogam_st_plot_request("predictions")
    expect_identical(prediction$type, "trajectories")
    expect_identical(prediction$dots$cluster_style, "none")

    clusters <- neurogam_st_plot_request("clusters", list(layout = "facets"))
    expect_identical(clusters$type, "trajectories")
    expect_identical(clusters$dots$cluster_style, "segments")
    expect_identical(clusters$dots$layout, "facets")

    topography <- neurogam_st_plot_request("topography", list(n_times = 4L))
    expect_identical(topography$type, "topography")
    expect_identical(topography$dots$n_times, 4L)

    expect_error(neurogam_st_plot_request("predictions", list(cluster_style = "segments")),
                 "cluster_style")
    expect_error(neurogam_st_plot_request("topography", list(results = 1)),
                 "Do not repeat")
    expect_error(neurogam_st_plot_request("topography", list(1)),
                 "unique names")
})

test_that("P2.4 registers all three plot methods", {
    for (suffix in c("1d", "2d", "st")) {
        method <- utils::getS3method(
            "plot", paste0("clusters_results_", suffix), optional = TRUE
        )
        expect_true(is.function(method))
        expect_true("type" %in% names(formals(method)))
    }
})

test_that("P2.4 reports missing spatiotemporal observations helpfully", {
    x <- structure(list(), class = "clusters_results_st")
    expect_error(plot(x, type = "topography"), "observed_data")
    expect_error(plot(x, type = "trajectories"), "observed_data")
    expect_error(plot(x, type = "unknown"), "Unsupported `type`")
})

test_that("P2.4 delegates all 3D views to existing plotters", {
    sensors <- data.frame(
        channel = c("Fz", "Cz", "Pz"),
        xproj = c(-0.5, 0.5, 0),
        yproj = c(0, 0, 1)
    )
    predictions <- merge(sensors, data.frame(time = c(0.1, 0.2)))
    predictions$estimate <- seq_len(nrow(predictions)) / 10
    predictions$lower <- predictions$estimate - 0.1
    predictions$upper <- predictions$estimate + 0.1
    predictions$prob_ratio <- 1
    observed <- predictions[c("channel", "time")]
    observed$voltage <- predictions$estimate

    members <- data.frame(
        channel = "Fz", time = 0.1, id = 1L, sign = "positive"
    )
    x <- list(
        predictions = predictions,
        cluster_summary = data.frame(
            id = 1L, sign = "positive", n_points = 1L,
            n_sensors = 1L, time_min = 0.1, time_max = 0.1
        ),
        cluster_membership = members,
        model = NULL
    )
    class(x) <- "clusters_results_st"

    for (type in c("trajectories", "predictions", "clusters")) {
        p <- plot(x, type = type, observed_data = observed, layout = "facets")
        expect_s3_class(p, "ggplot")
    }
    p_topo <- plot(
        x, type = "topography", observed_data = observed,
        n_times = 2L, grid_n = 15L, contour = FALSE
    )
    expect_s3_class(p_topo, "ggplot")
})
