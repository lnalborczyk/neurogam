test_that("participant-level summary keeps sensors and conditions separate", {
    d <- expand.grid(
        subject = c("s1", "s2"),
        channel = c("Fz", "Cz"),
        time = c(0.1, 0.2),
        condition = c("A", "B"),
        trial = 1:3,
        KEEP.OUT.ATTRS = FALSE
    )
    d$xproj <- ifelse(d$channel == "Fz", 0, 1)
    d$yproj <- ifelse(d$channel == "Fz", 1, 0)
    d$voltage <- ifelse(d$condition == "B", 1, 0) + d$trial * 0.1
    s <- make_spatiotemporal_summary_data(
        d,
        participant_id = "subject",
        outcome_id = "voltage",
        sensor_id = "channel",
        predictor_id = "condition"
    )
    expect_equal(nrow(s), 2L * 2L * 2L * 2L)
    expect_true(all(s$outcome_sem > 0))
    expect_equal(levels(s$predictor), c("A", "B"))

    g <- make_spatiotemporal_prediction_grid(s)
    expect_equal(nrow(g$a), 4L)
    expect_equal(as.character(unique(g$a$predictor)), "B")
    expect_equal(as.character(unique(g$b$predictor)), "A")
    expect_identical(g$a$time, g$b$time)
    expect_identical(g$a$sensor, g$b$sensor)
})

test_that("single-observation bins require use_se = FALSE or supplied SEM", {
    d <- data.frame(
        participant = "s1", channel = "Cz", time = 0.1,
        xproj = 0, yproj = 0, voltage = 1
    )
    expect_error(make_spatiotemporal_summary_data(
        d, "participant", "voltage", sensor_id = "channel"
    ), "at least two trials")
    expect_silent(make_spatiotemporal_summary_data(
        d, "participant", "voltage", sensor_id = "channel", use_se = FALSE
    ))
})

test_that("coordinate conflicts are detected", {
    d <- data.frame(
        participant = c("a", "a"), channel = c("Fz", "Fz"),
        time = c(0.1, 0.2), xproj = c(0, 1),
        yproj = c(0, 0), voltage = c(1, 2)
    )
    expect_error(make_spatiotemporal_summary_data(
        d, "participant", "voltage", sensor_id = "channel", use_se = FALSE
    ), "one pair")
})
