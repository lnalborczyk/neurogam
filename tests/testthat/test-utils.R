test_that("st_take_n_times selects endpoints", {

    x <- seq(0, 1, length.out = 101)
    out <- st_take_n_times(x, 8)
    expect_true(out[1] == x[1])
    expect_true(out[length(out)] == x[length(x)])

    })

test_that("st_take_n_times handles n=1", {

    x <- 1:10
    out <- st_take_n_times(x, 1)
    expect_length(out, 1)

    })
