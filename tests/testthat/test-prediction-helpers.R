test_that("posterior draw selection preserves endpoints and count", {
    expect_null(select_prediction_draws(100, NULL))
    expect_identical(select_prediction_draws(100, 1), 1L)
    expect_identical(select_prediction_draws(100, 100), seq_len(100))
    x <- select_prediction_draws(100, 7)
    expect_length(x, 7)
    expect_identical(x[c(1, 7)], c(1L, 100L))
    expect_true(all(diff(x) > 0))
    expect_error(select_prediction_draws(5, 6), "n_requested")
    expect_error(select_prediction_draws(5, 1.5), "n_requested")
})

test_that("chunks cover rows exactly once and in order", {
    chunks <- prediction_chunks(9, 4)
    expect_identical(unlist(chunks, use.names = FALSE), 1:9)
    #expect_identical(lengths(chunks), c(4L, 4L, 1L))
    expect_identical(
        unname(lengths(chunks)),
        c(4L, 4L, 1L)
    )
    expect_error(prediction_chunks(0, 4), "positive integers")
})

test_that("paired contrasts use the same draw ids and row ordering", {
    a <- data.frame(location = 1:5, value = 11:15)
    b <- data.frame(location = 1:5, value = 1:5)
    calls <- list()
    mock_epred <- function(model, newdata, re_formula, draw_ids) {
        calls[[length(calls) + 1L]] <<- list(
            locations = newdata$location,
            draw_ids = draw_ids,
            re_formula = re_formula
        )
        matrix(rep(newdata$value, each = length(draw_ids)),
               nrow = length(draw_ids)) + matrix(draw_ids, length(draw_ids),
                                                nrow(newdata))
    }
    out <- chunked_epred(
        NULL, a, newdata_b = b, draw_ids = c(2L, 5L),
        chunk_size = 2L, predict_fun = mock_epred
    )
    expect_length(out, 3)
    expect_equal(do.call(cbind, out), matrix(10, 2, 5))
    expect_identical(lapply(calls, `[[`, "draw_ids"),
                     rep(list(c(2L, 5L)), 6))
    expect_identical(lapply(calls, `[[`, "locations"),
                     list(1:2, 1:2, 3:4, 3:4, 5L, 5L))
})

test_that("chunked summaries preserve prediction location order", {
    mock_epred <- function(model, newdata, re_formula, draw_ids) {
        rbind(newdata$value, newdata$value + 2)
    }
    data <- data.frame(value = c(4, 8, 2, 6, 10))
    result <- chunked_epred(
        NULL, data, chunk_size = 2,
        predict_fun = mock_epred,
        summarise_fun = function(x) data.frame(estimate = colMeans(x))
    )
    expect_equal(do.call(rbind, result)$estimate, data$value + 1)
})
