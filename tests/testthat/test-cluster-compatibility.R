test_that("1D legacy schema is preserved", {
    dat <- data.frame(
        time = seq_len(8),
        value = c(20, 20, 1, 0.01, 0.01, 1, 20, 1)
    )
    result <- find_clusters_nd(dat, coordinates = "time", threshold = 10)
    legacy <- cluster_nd_to_legacy(result)

    expect_named(legacy, c("id", "onset", "offset", "n_points", "sign"))
    expect_equal(nrow(legacy), 3L)
    expect_equal(sum(legacy$n_points), nrow(result$membership))
})

test_that("2D legacy memberships retain the required columns", {
    dat <- expand.grid(time1 = 1:4, time2 = 1:4)
    dat$value <- ifelse(dat$time1 == 2 & dat$time2 %in% 2:3, 20, 1)
    result <- find_clusters_nd(
        dat,
        coordinates = c("time1", "time2"),
        threshold = 10
    )
    legacy <- cluster_nd_to_legacy(result, time_id = c("time1", "time2"))

    expect_named(
        legacy,
        c("id", "sign", "time1", "time2", "value", "n_points")
    )
    expect_equal(nrow(legacy), 2L)
    expect_true(all(legacy$n_points == 2L))
})

test_that("The adapter handles no detected clusters", {
    dat <- data.frame(time = seq_len(10), value = rep(1, 10))
    result <- find_clusters_nd(dat, threshold = 10)
    legacy <- cluster_nd_to_legacy(result)
    expect_equal(nrow(legacy), 0L)
    expect_named(legacy, c("id", "onset", "offset", "n_points", "sign"))
})
