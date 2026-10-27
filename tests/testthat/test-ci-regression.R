# Fast regression checks. Do not fit brms models here.

# canonical_components <- function(membership, coordinates) {
#     if (nrow(membership) == 0L) {
#         return(character())
#     }
#
#     pieces <- split(
#         membership,
#         interaction(membership[c("id", "sign")], drop = TRUE)
#     )
#
#     sort(vapply(
#         pieces,
#         function(component) {
#             coordinates_text <- do.call(
#                 paste,
#                 c(unname(component[coordinates]), sep = ":")
#             )
#             paste(sort(coordinates_text), collapse = "|")
#         },
#         character(1)
#     ))
# }

canonical_components <- function(membership, coordinates) {

    if (nrow(membership) == 0L) {
        return(character())
    }

    pieces <- split(
        membership,
        interaction(
            membership[c("id", "sign")],
            drop = TRUE
        )
    )

    signatures <- vapply(
        pieces,
        function(component) {

            coordinates_text <- do.call(
                paste,
                c(
                    unname(component[coordinates]),
                    sep = ":"
                )
            )

            paste0(
                component$sign[1L],
                ":",
                paste(sort(coordinates_text), collapse = "|")
            )
        },
        character(1)
    )

    # Ignore arbitrary component IDs and ordering.
    unname(sort(signatures))
}

test_that("1D clusters and legacy summaries remain compatible", {
    d <- data.frame(
        time = 1:8,
        value = c(1, 20, 20, 1, 0.05, 0.05, 1, 1)
    )

    nd <- find_clusters_nd(d, threshold = 10)
    old <- find_clusters(d, threshold = 10)
    legacy <- cluster_nd_to_legacy(nd, time_id = "time")

    expect_equal(nrow(nd$clusters), 2L)
    expect_equal(nrow(nd$membership), 4L)
    expect_setequal(nd$clusters$sign, c("positive", "negative"))

    columns <- c("sign", "onset", "offset", "n_points")
    ordered <- function(x) {
        x <- as.data.frame(x[, columns, drop = FALSE])
        x <- x[order(x$sign, x$onset), , drop = FALSE]
        rownames(x) <- NULL
        x
    }
    expect_equal(ordered(legacy), ordered(old))
})

test_that("2D grid components do not connect diagonally", {
    d <- expand.grid(
        time1 = 1:4,
        time2 = 1:4,
        KEEP.OUT.ATTRS = FALSE
    )
    d$value <- 1
    d$value[d$time1 %in% 1:2 & d$time2 == 1] <- 20
    d$value[d$time1 == 4 & d$time2 == 4] <- 20

    nd <- find_clusters_nd(
        d,
        coordinates = c("time1", "time2"),
        threshold = 10
    )

    shuffled <- find_clusters_nd(
        d[nrow(d):1, , drop = FALSE],
        coordinates = c("time1", "time2"),
        threshold = 10
    )

    expect_equal(nrow(nd$clusters), 2L)
    expect_equal(sort(nd$clusters$n_points), c(1L, 2L))
    expect_identical(
        canonical_components(nd$membership, c("time1", "time2")),
        canonical_components(shuffled$membership, c("time1", "time2"))
    )
})

test_that("sensor adjacency and spatiotemporal clustering are deterministic", {
    sensors <- data.frame(
        channel = c("Fz", "Cz", "Pz"),
        xproj = c(0, 1, 2),
        yproj = c(0, 0, 0)
    )

    adjacency <- make_sensor_adjacency(
        sensors,
        sensor_id = "channel",
        k = 1,
        distance_factor = NULL
    )
    expect_true(isSymmetric(adjacency))
    expect_false(any(diag(adjacency)))
    expect_true(adjacency["Fz", "Cz"])
    expect_true(adjacency["Cz", "Pz"])
    expect_false(adjacency["Fz", "Pz"])

    eeg <- expand.grid(
        channel = c("Fz", "Cz", "Pz"),
        time = 1:6,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    eeg$value <- 1
    eeg$value[(eeg$channel == "Fz" & eeg$time %in% c(2, 3)) |
              (eeg$channel == "Cz" & eeg$time == 3)] <- 20
    eeg$value[eeg$channel == "Pz" & eeg$time == 6] <- 20
    eeg$value[eeg$channel %in% c("Cz", "Pz") & eeg$time == 5] <- 0.02

    run <- function(dat) {
        find_clusters_nd(
            dat,
            coordinates = "time",
            sensor_id = "channel",
            sensor_adjacency = adjacency,
            threshold = 10
        )
    }

    nd <- run(eeg)
    reordered <- run(eeg[nrow(eeg):1, , drop = FALSE])

    expect_equal(nrow(nd$clusters), 3L)
    expect_equal(nrow(nd$membership), 6L)
    expect_equal(sort(nd$clusters$n_points), c(1L, 2L, 3L))
    expect_identical(
        canonical_components(nd$membership, c("channel", "time")),
        canonical_components(reordered$membership, c("channel", "time"))
    )
})

test_that("posterior summaries use one-sided odds by default", {

    draws <- cbind(
        positive = rep(3, 40),
        negative = rep(-2, 40),
        inside = rep(0, 40)
    )

    x <- summarise_spatiotemporal_draws(
        draws,
        null_value = 0,
        rope = 0.1,
        credible_interval = 0.95
    )

    expect_equal(x$estimate, c(3, -2, 0))
    expect_equal(x$prob_above, c(1, 0, 0))
    expect_equal(x$prob_below, c(0, 1, 0))
    expect_equal(x$prob_inside, c(0, 0, 1))

    # One-sided odds, capped using 40 posterior draws.
    # expect_equal(
    #     x$prob_ratio,
    #     c(40, 1 / 40, 1 / 40)
    # )

    expect_equal(x$prob_ratio, c(Inf, 0, 0))
    expect_equal(x$odds_above, c(Inf, 0, 0))
    expect_equal(x$odds_below, c(0, Inf, 0))

    # Confirm the historical mode is still available.
    historical <- summarise_spatiotemporal_draws(
        draws,
        null_value = 0,
        rope = 0.1,
        credible_interval = 0.95,
        odds_mode = "spatiotemporal"
    )

    expect_true(is.infinite(historical$prob_ratio[1]))
    expect_equal(historical$prob_ratio[2:3], c(0, 1))
})


test_that("cluster comparison ignores IDs but preserves signs", {

    original <- data.frame(
        time1 = c(1, 2, 4),
        time2 = c(1, 1, 4),
        id = c(1L, 1L, 2L),
        sign = c("positive", "positive", "negative")
    )

    reordered <- original[c(3, 2, 1), ]

    # Swap arbitrary cluster IDs.
    reordered$id <- c(1L, 2L, 2L)

    expect_identical(
        canonical_components(
            original,
            c("time1", "time2")
        ),
        canonical_components(
            reordered,
            c("time1", "time2")
        )
    )

    # Changing cluster polarity must be detected.
    reordered$sign[1] <- "positive"

    expect_false(identical(
        canonical_components(
            original,
            c("time1", "time2")
        ),
        canonical_components(
            reordered,
            c("time1", "time2")
        )
    ))
})
