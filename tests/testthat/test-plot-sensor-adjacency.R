test_that("plot_sensor_adjacency returns an editable ggplot", {
    sensors <- data.frame(
        sensor = c("Pz", "Fz", "Cz"),
        xproj = c(0, 0, 0),
        yproj = c(-1, 1, 0)
    )
    adjacency <- make_sensor_adjacency(
        sensors,
        k = 1,
        distance_factor = 1.5
    )

    plot <- plot_sensor_adjacency(sensors, adjacency)

    expect_s3_class(plot, "ggplot")
    expect_length(plot$layers, 3L)  # Edges, points, labels

    edge_data <- ggplot2::layer_data(plot, 1L)
    point_data <- ggplot2::layer_data(plot, 2L)

    expect_equal(nrow(edge_data), 2L)
    expect_equal(nrow(point_data), 3L)
    expect_true(all(edge_data$x == edge_data$xend))
    expect_equal(sort(point_data$y), c(-1, 0, 1))

    # Regular ggplot2 customisation should work.
    customised <- plot + ggplot2::labs(subtitle = "Custom subtitle")
    expect_s3_class(customised, "ggplot")
    expect_identical(customised$labels$subtitle, "Custom subtitle")
})

test_that("sensor order is matched by names, not row order", {
    sensors <- data.frame(
        channel = c("Pz", "Fz", "Cz"),
        xproj = c(0, 0, 0),
        yproj = c(-1, 1, 0)
    )
    adjacency <- make_sensor_adjacency(
        sensors,
        sensor_id = "channel",
        k = 1,
        distance_factor = 1.5
    )

    plot1 <- plot_sensor_adjacency(
        sensors,
        adjacency,
        sensor_id = "channel"
    )
    plot2 <- plot_sensor_adjacency(
        sensors[3:1, , drop = FALSE],
        adjacency[3:1, 3:1, drop = FALSE],
        sensor_id = "channel"
    )

    expect_equal(
        ggplot2::layer_data(plot1, 1L),
        ggplot2::layer_data(plot2, 1L)
    )
})

test_that("isolated electrodes and empty edge lists can be plotted", {
    sensors <- data.frame(
        sensor = c("A", "B"),
        xproj = c(-1, 1),
        yproj = c(0, 0)
    )
    adjacency <- matrix(
        FALSE,
        2, 2,
        dimnames = list(sensors$sensor, sensors$sensor)
    )

    plot <- plot_sensor_adjacency(
        sensors,
        adjacency,
        show_labels = FALSE,
        show_axes = FALSE
    )

    expect_s3_class(plot, "ggplot")
    expect_length(plot$layers, 1L)  # Only the sensor markers
    expect_equal(nrow(ggplot2::layer_data(plot, 1L)), 2L)
})

test_that("invalid adjacency matrices are rejected", {
    sensors <- data.frame(
        sensor = c("A", "B"),
        xproj = c(0, 1),
        yproj = c(0, 0)
    )
    adjacency <- matrix(
        c(FALSE, TRUE, TRUE, FALSE),
        2, 2,
        dimnames = list(sensors$sensor, sensors$sensor)
    )

    expect_error(
        plot_sensor_adjacency(sensors, unname(adjacency)),
        "named logical or numeric matrix"
    )

    asymmetric <- adjacency
    asymmetric["A", "B"] <- FALSE
    expect_error(
        plot_sensor_adjacency(sensors, asymmetric),
        "symmetric"
    )

    bad_names <- adjacency
    rownames(bad_names) <- c("A", "Z")
    expect_error(
        plot_sensor_adjacency(sensors, bad_names),
        "covering all sensors"
    )
})
