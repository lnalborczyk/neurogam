test_that("fit_or_reuse_bgam forwards fitting arguments without sampling", {
    recorded <- NULL
    fake_model <- structure(list(source = "fit"), class = "brmsfit")

    fake_fit <- function(...) {
        recorded <<- list(...)
        fake_model
    }

    ans <- fit_or_reuse_bgam(
        formula = stats::as.formula("y ~ x"),
        data = data.frame(y = 1:3, x = 4:6),
        family = stats::gaussian(),
        warmup = 500,
        iter = 1500,
        chains = 2,
        cores = 1,
        threads = NULL,
        backend = "cmdstanr",
        stan_control = list(adapt_delta = 0.95),
        file = "cache/example",
        extra_args = list(seed = 123L, refresh = 0L),
        fit_fun = fake_fit,
        ndraws_fun = function(model) 500L
    )

    expect_identical(ans$model, fake_model)
    expect_identical(ans$n_post_samples, 500L)
    expect_false(ans$reused)
    expect_identical(recorded$warmup, 500)
    expect_identical(recorded$iter, 1500)
    expect_identical(recorded$backend, "cmdstanr")
    expect_identical(recorded$control, list(adapt_delta = 0.95))
    expect_identical(recorded$file, "cache/example")
    expect_identical(recorded$seed, 123L)
    expect_identical(recorded$refresh, 0L)
    expect_identical(
        recorded$stan_model_args,
        list(stanc_options = list("O1"))
    )
    expect_identical(names(recorded), c(
        "formula", "data", "family", "warmup", "iter", "chains",
        "cores", "threads", "backend", "control", "file",
        "stan_model_args", "seed", "refresh"
    ))
})

test_that("reusing a model never calls the fitter", {
    previous <- structure(list(source = "previous"), class = "brmsfit")
    fail_fit <- function(...) stop("Fitting was called unexpectedly")
    draw_calls <- 0L
    fake_ndraws <- function(model) {
        draw_calls <<- draw_calls + 1L
        1200L
    }

    expect_message(
        ans <- fit_or_reuse_bgam(
            formula = stats::as.formula("y ~ 1"),
            data = data.frame(y = 1),
            family = stats::gaussian(),
            previous_model = previous,
            fit_fun = fail_fit,
            ndraws_fun = fake_ndraws
        ),
        "skipping model fitting"
    )
    expect_identical(ans$model, previous)
    expect_true(ans$reused)
    expect_identical(ans$n_post_samples, 1200L)
    expect_identical(draw_calls, 1L)

    expect_message(
        ans2 <- fit_or_reuse_bgam(
            formula = stats::as.formula("y ~ 1"),
            data = data.frame(y = 1),
            family = stats::gaussian(),
            previous_model = previous,
            n_post_samples = 100L,
            fit_fun = fail_fit,
            ndraws_fun = fake_ndraws
        ),
        "skipping model fitting"
    )
    expect_identical(ans2$n_post_samples, 100L)
    expect_identical(draw_calls, 1L)
})

test_that("draw resolution can be left to the caller", {
    fake_fit <- function(...) structure(list(), class = "brmsfit")

    ans <- fit_or_reuse_bgam(
        formula = stats::as.formula("y ~ 1"),
        data = data.frame(y = 1),
        family = stats::gaussian(),
        resolve_draws = FALSE,
        fit_fun = fake_fit,
        ndraws_fun = function(model) stop("ndraws must not be called")
    )

    expect_null(ans$n_post_samples)
})

test_that("rstan backend receives its original default compilation args", {
    recorded <- NULL
    fit_fn <- function(...) {
        recorded <<- list(...)
        structure(list(), class = "brmsfit")
    }

    fit_or_reuse_bgam(
        formula = stats::as.formula("y ~ 1"),
        data = data.frame(y = 1),
        family = stats::gaussian(),
        backend = "rstan",
        n_post_samples = 50L,
        fit_fun = fit_fn
    )

    expect_identical(recorded$stan_model_args, list())
})

test_that("invalid previous models and duplicate fitting args are rejected", {
    args <- list(
        formula = stats::as.formula("y ~ 1"),
        data = data.frame(y = 1),
        family = stats::gaussian()
    )

    expect_error(
        do.call(fit_or_reuse_bgam, c(args, list(previous_model = list()))),
        "brmsfit"
    )
    expect_error(
        do.call(fit_or_reuse_bgam, c(args, list(extra_args = list(cores = 7)))),
        "cores"
    )
    expect_error(
        do.call(fit_or_reuse_bgam, c(args, list(warmup = 100, iter = 100))),
        "larger than"
    )
})
