test_that(
    desc = "one-dimensional intervals are connected and signs stay separate",
    code = {
        d <- data.frame(time = 1:8, value = c(1, 20, 30, 1, 0.01, 0.02, 1, 20) )
        r <- find_clusters_nd(d, threshold = 10)
        expect_equal(nrow(r$clusters), 3L)
        expect_equal(sort(r$clusters$n_points), c(1L, 2L, 2L) )
        expect_setequal(r$clusters$sign, c("positive", "negative") )
        expect_equal(nrow(r$membership), 5L)
        expect_equal(sort(r$membership$time), c(2, 3, 5, 6, 8) )
        })

test_that(
    desc = "four-neighbor grid connectivity excludes diagonal links",
    code = {
        d <- expand.grid(train_time = 1:3, test_time = 1:3)
        d$value <- 1
        d$value[(d$train_time == 1 & d$test_time == 1) | (d$train_time == 2 & d$test_time == 2)] <- 20
        r <- find_clusters_nd(d, coordinates = c("train_time", "test_time"), threshold = 10)
        expect_equal(nrow(r$clusters), 2L)
        expect_true(all(r$clusters$n_points == 1L) )
        })

test_that(
    desc = "sensor-time topology connects spatial and temporal neighbors",
    code = {
        adj <- matrix(FALSE, 3, 3, dimnames = list(c("Fz", "Cz", "Pz"), c("Fz", "Cz", "Pz") ) )
        adj["Fz", "Cz"] <- adj["Cz", "Fz"] <- TRUE
        adj["Cz", "Pz"] <- adj["Pz", "Cz"] <- TRUE
        d <- expand.grid(time = 1:4, sensor = c("Fz", "Cz", "Pz"), stringsAsFactors = FALSE)
        d$value <- 1
        d$value[(d$sensor == "Fz" & d$time == 1) |
              (d$sensor == "Cz" & d$time %in% c(1, 2)) |
              (d$sensor == "Pz" & d$time == 2)] <- 20
        d <- d[sample(nrow(d)), ]
        r <- find_clusters_nd(d, coordinates = "time", sensor_id = "sensor", sensor_adjacency = adj, threshold = 10)
        expect_equal(nrow(r$clusters), 1L)
        expect_equal(r$clusters$n_points, 4L)
        expect_equal(r$clusters$n_sensors, 3L)
        expect_equal(r$clusters$onset, 1)
        expect_equal(r$clusters$offset, 2)
        })

test_that(
    desc = "missing cells cannot create phantom connectivity",
    code = {
    d <- data.frame(train_time = c(1, 3), test_time = c(1, 1), value = 20)
    r <- find_clusters_nd(d, coordinates = c("train_time", "test_time"), max_gap = list(train_time = 1) )
    expect_equal(nrow(r$clusters), 2L)
    })

test_that(
    desc = "group boundaries are respected",
    code = {
    d <- expand.grid(time = 1:3, participant = c("P1", "P2"), stringsAsFactors = FALSE)
    d$value <- 20
    r <- find_clusters_nd(d, group = "participant")
    expect_equal(nrow(r$clusters), 2L)
    expect_setequal(r$clusters$participant, c("P1", "P2") )
    })

test_that(
    desc = "no hits return typed empty tables",
    code = {
    d <- data.frame(time = 1:5, value = rep(1, 5) )
    r <- find_clusters_nd(d)
    expect_equal(nrow(r$clusters), 0L)
    expect_equal(nrow(r$membership), 0L)
    expect_true(all(c("id", "sign", "onset", "offset") %in% names(r$clusters) ) )
    })

test_that(
    desc = "duplicated locations and asymmetric adjacency fail",
    code = {
    expect_error(find_clusters_nd(data.frame(time = c(1, 1), value = c(20, 20) ) ), "Duplicate")
    a <- matrix(c(0, 1, 0, 0), 2, dimnames = list(c("A", "B"), c("A", "B") ) )
    d <- data.frame(time = c(1, 1), sensor = c("A", "B"), value = 20)
    expect_error(find_clusters_nd(d, sensor_id = "sensor", sensor_adjacency = a), "symmetric")
    })
