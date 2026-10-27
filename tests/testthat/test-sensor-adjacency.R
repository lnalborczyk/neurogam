test_that("automatic adjacency is named, symmetric, and has no self-links", {
    sensors <- data.frame(
        sensor = c("Pz", "Fz", "Cz"),
        xproj = c(3, 0, 1),
        yproj = c(0, 0, 0)
    )

    adjacency <- make_sensor_adjacency(
        sensors,
        k = 2,
        distance_factor = 2
    )

    expect_type(adjacency, "logical")
    expect_equal(dim(adjacency), c(3L, 3L))
    expect_identical(rownames(adjacency), c("Cz", "Fz", "Pz"))
    expect_identical(rownames(adjacency), colnames(adjacency))
    expect_identical(adjacency, t(adjacency))
    expect_false(any(diag(adjacency)))

    expect_true(adjacency["Fz", "Cz"])
    expect_true(adjacency["Cz", "Pz"])
    expect_false(adjacency["Fz", "Pz"])
})

test_that("sensor-time repetition and input ordering do not affect topology", {
    sensors <- expand.grid(
        xproj = 0:2,
        yproj = 0:2
    )
    sensors$sensor <- sprintf("S%02d", seq_len(nrow(sensors)))

    repeated <- merge(
        sensors,
        data.frame(time = c(0.1, 0.2, 0.3)),
        by = NULL
    )

    baseline <- make_sensor_adjacency(sensors)
    repeated_adj <- make_sensor_adjacency(repeated)

    set.seed(20261027)
    shuffled <- make_sensor_adjacency(
        repeated[sample(nrow(repeated)), , drop = FALSE]
    )

    expect_identical(repeated_adj, baseline)
    expect_identical(shuffled, baseline)
    expect_true(any(baseline))
})

test_that("distance pruning prevents unnecessary long-distance links", {
    positions <- data.frame(
        sensor = c("A", "B", "C"),
        xproj = c(0, 1, 3),
        yproj = c(0, 0, 0)
    )

    filtered <- make_sensor_adjacency(
        positions,
        k = 2,
        distance_factor = 2
    )
    unfiltered <- make_sensor_adjacency(
        positions,
        k = 2,
        distance_factor = NULL
    )

    expect_false(filtered["A", "C"])
    expect_true(unfiltered["A", "C"])
})

test_that("one-sensor montages are well defined", {
    one <- data.frame(sensor = "Cz", xproj = 0, yproj = 0)
    adjacency <- make_sensor_adjacency(one)

    expect_identical(dim(adjacency), c(1L, 1L))
    expect_false(adjacency[1, 1])
    expect_identical(rownames(adjacency), "Cz")
})

test_that("invalid montage geometry is rejected", {
    same_position <- data.frame(
        sensor = c("A", "B"),
        xproj = c(0, 0),
        yproj = c(0, 0)
    )
    expect_error(
        make_sensor_adjacency(same_position),
        "identical projected positions"
    )

    conflicting <- data.frame(
        sensor = c("A", "A"),
        xproj = c(0, 1),
        yproj = c(0, 0)
    )
    expect_error(
        make_sensor_adjacency(conflicting),
        "consistent coordinate"
    )

    incomplete <- data.frame(
        sensor = c("A", "B"),
        xproj = c(0, NA_real_),
        yproj = c(0, 1)
    )
    expect_error(
        make_sensor_adjacency(incomplete),
        "finite numeric"
    )

    valid <- data.frame(
        sensor = c("A", "B"),
        xproj = c(0, 1),
        yproj = c(0, 0)
    )
    expect_error(make_sensor_adjacency(valid, k = 0), "positive integer")
    expect_error(
        make_sensor_adjacency(valid, distance_factor = -1),
        "positive finite"
    )
})

test_that("disconnected graphs raise a warning", {
    positions <- data.frame(
        sensor = c("A", "B", "C"),
        xproj = c(0, 1, 100),
        yproj = c(0, 0, 0)
    )

    expect_warning(
        make_sensor_adjacency(positions, k = 2),
        "disconnected components"
    )
})

test_that("auto-generated adjacency works with multidimensional clustering", {
    skip_if_not_installed("igraph")

    sensors <- data.frame(
        sensor = c("Fz", "Cz", "Pz"),
        xproj = c(0, 1, 3),
        yproj = c(0, 0, 0)
    )
    adjacency <- make_sensor_adjacency(
        sensors,
        k = 2,
        distance_factor = 2
    )
    values <- expand.grid(sensor = sensors$sensor, time = 1:3)
    values$value <- 1
    values$value[
        values$time == 2 & values$sensor %in% c("Fz", "Cz")
    ] <- 20
    values$value[values$time == 3 & values$sensor == "Pz"] <- 0.02

    result <- find_clusters_nd(
        values,
        coordinates = "time",
        sensor_id = "sensor",
        sensor_adjacency = adjacency,
        threshold = 10
    )

    expect_equal(nrow(result$clusters), 2L)
    expect_equal(sort(result$clusters$n_points), c(1L, 2L))
    expect_setequal(result$clusters$sign, c("positive", "negative"))
})

test_that("analysis function defaults to automatic adjacency", {
    expect_identical(
        formals(analyse_bgam_spatiotemporal)$sensor_adjacency,
        NULL
    )
})
