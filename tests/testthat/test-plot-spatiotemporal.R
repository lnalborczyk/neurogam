make_plot_test_data <- function(contrast = FALSE, empty_clusters = FALSE) {
    locations <- data.frame(
        channel = c("F3", "F4", "P3", "P4"),
        xproj = c(-1, 1, -1, 1),
        yproj = c(1, 1, -1, -1)
    )
    grid <- merge(locations, data.frame(time = seq(0, 0.07, by = 0.01)))
    grid <- grid[order(grid$channel, grid$time), , drop = FALSE]
    grid$estimate <- sin(grid$time * 25) + grid$xproj * 0.1
    grid$lower <- grid$estimate - 0.15
    grid$upper <- grid$estimate + 0.15

    observed <- grid[c("channel", "time")]
    observed$voltage <- grid$estimate + 0.08
    if (contrast) {
        observed_b <- observed
        observed_b$voltage <- 0.08
    } else {
        observed_b <- NULL
    }
    members <- data.frame(
        channel = c("F3", "F3", "F3", "F3", "P4", "P4"),
        time = c(0, 0.01, 0.03, 0.04, 0.05, 0.06),
        id = c(1L, 1L, 1L, 1L, 2L, 2L),
        sign = c("positive", "positive", "positive", "positive",
                 "negative", "negative"),
        value = c(15, 20, 20, 15, 0.03, 0.05)
    )
    if (empty_clusters) {
        members <- members[FALSE, , drop = FALSE]
    }
    results <- list(
        predictions = grid,
        cluster_membership = members,
        contrast = if (contrast) "A_minus_B" else "one_sample"
    )
    list(results = results, observed = observed, observed_b = observed_b)
}

test_that("observed and predicted arrays are matched correctly", {
    ex <- make_plot_test_data()
    d <- prepare_spatiotemporal_plot_data(
        ex$results, ex$observed,
        outcome_id = "voltage", sensor_id = "channel"
    )
    expect_equal(nrow(d), nrow(ex$results$predictions))
    expect_equal(d$residual, rep(0.08, nrow(d)), tolerance = 1e-9)
    expect_equal(sum(!is.na(d$cluster_id)), 6L)
})

test_that("subject-weighted observed means differ from trial-weighted means", {
    ex <- make_plot_test_data()
    first <- ex$observed
    first$subject <- "A"
    second <- ex$observed
    second$subject <- "B"
    second$voltage <- second$voltage + 2
    first <- rbind(first, first, second)
    d <- prepare_spatiotemporal_plot_data(
        ex$results, first, outcome_id = "voltage", sensor_id = "channel",
        participant_id = "subject"
    )
    expect_equal(d$residual, rep(1.08, nrow(d)), tolerance = 1e-9)
})

test_that("condition contrasts require condition B observations", {
    ex <- make_plot_test_data(contrast = TRUE)
    expect_error(
        prepare_spatiotemporal_plot_data(
            ex$results, ex$observed, sensor_id = "channel"
        ),
        "observed_data_b"
    )
    d <- prepare_spatiotemporal_plot_data(
        ex$results, ex$observed, ex$observed_b,
        sensor_id = "channel"
    )
    expect_equal(d$residual, rep(0, nrow(d)), tolerance = 1e-9)
})

test_that("sensor intervals do not cross omitted samples", {
    ex <- make_plot_test_data()
    d <- prepare_spatiotemporal_plot_data(
        ex$results, ex$observed, sensor_id = "channel"
    )
    intervals <- make_cluster_intervals_per_sensor(d)
    expect_equal(nrow(intervals), 3L)
    expect_equal(sum(intervals$sensor == "F3"), 2L)
    expect_true(all(intervals$xmin < intervals$xmax))

    segments <- make_cluster_segments_per_sensor(d)
    expect_equal(nrow(segments), 3L)
    expect_true(all(segments$timeend > segments$time))
    expect_false(any(segments$time == 0.01 & segments$timeend == 0.03))
})

test_that("topography displays selected time panels and residuals", {
    ex <- make_plot_test_data()
    p <- plot_spatiotemporal_topography(
        ex$results, ex$observed, sensor_id = "channel",
        times = c(0.011, 0.041), grid_n = 15L, contour = FALSE,
        cluster_overlay = "all"
    )
    expect_s3_class(p, "ggplot")
    expect_equal(length(unique(p$data$time_label)), 2L)
    expect_equal(length(unique(p$data$view)), 3L)
    expect_silent(ggplot2::ggplot_build(p))
})

test_that("both cluster styles render in scalp and facet layouts", {
    ex <- make_plot_test_data()
    for (plot_layout in c("scalp", "facets")) {
        for (style in c("rectangles", "segments", "none")) {
            p <- plot_sensor_trajectories(
                ex$results, ex$observed,
                sensor_id = "channel", layout = plot_layout,
                cluster_style = style
            )
            expect_s3_class(p, "ggplot")
            expect_silent(ggplot2::ggplot_build(p))
        }
    }
})

test_that("empty clusters do not break either visualization", {
    ex <- make_plot_test_data(empty_clusters = TRUE)
    p1 <- plot_spatiotemporal_topography(
        ex$results, ex$observed, sensor_id = "channel",
        grid_n = 15L, n_times = 2L, contour = FALSE
    )
    p2 <- plot_sensor_trajectories(
        ex$results, ex$observed,
        sensor_id = "channel", cluster_style = "segments"
    )
    expect_silent(ggplot2::ggplot_build(p1))
    expect_silent(ggplot2::ggplot_build(p2))
})

test_that("missing sensors and duplicate predictions fail explicitly", {
    ex <- make_plot_test_data()
    expect_error(
        plot_sensor_trajectories(
            ex$results, ex$observed,
            sensor_id = "channel", sensors = "unknown"
        ),
        "Unknown sensor"
    )
    ex$results$predictions <- rbind(
        ex$results$predictions, ex$results$predictions[1L, ]
    )
    expect_error(
        prepare_spatiotemporal_plot_data(
            ex$results, ex$observed, sensor_id = "channel"
        ),
        "one row per sensor and time"
    )
})
