#' Fit a BGAM or reuse an existing brms model
#'
#' Centralizes the brms fitting arguments and the `previous_model` path.
#' This helper does not validate model/data compatibility, construct prediction
#' grids, or compute any posterior quantities; callers retain those duties.
#'
#' @param formula A brms model formula.
#' @param data The model data frame.
#' @param family A brms family.
#' @param previous_model Optional existing `brmsfit`.
#' @param n_post_samples Optional requested number of posterior draws.
#' @param warmup,iter,chains,cores,threads Sampling settings.
#' @param backend Backend used by brms.
#' @param stan_control Sampling control options forwarded as `control`.
#' @param file Optional model cache filename.
#' @param stan_model_args Optional Stan compilation options. If NULL, uses
#'   `stanc_options = list("O1")` for cmdstanr, and an empty list for rstan.
#' @param extra_args Named or unnamed arguments forwarded to `brms::brm()`.
#' @param resolve_draws Whether to populate `n_post_samples` from the model
#'   if the caller supplied NULL. Set FALSE when draw selection is handled
#'   separately by the caller.
#' @param reuse_message Message printed if an existing fit is reused.
#' @param fit_fun,ndraws_fun Callable dependencies; injectable for fast tests.
#'
#' @return A list containing `model`, `n_post_samples`, and `reused`.
#' @keywords internal
fit_or_reuse_bgam <- function (
    formula,
    data,
    family,
    previous_model = NULL,
    n_post_samples = NULL,
    warmup = 1000,
    iter = 2000,
    chains = 4,
    cores = 4,
    threads = NULL,
    backend = "cmdstanr",
    stan_control = NULL,
    file = NULL,
    stan_model_args = NULL,
    extra_args = list(),
    resolve_draws = TRUE,
    reuse_message = "Using `previous_model`: skipping model fitting.",
    fit_fun = brms::brm,
    ndraws_fun = brms::ndraws
    ) {

    if (!is.list(extra_args) ) {

        stop ("`extra_args` must be a list.", call. = FALSE)

    }

    backend <- match.arg(backend, c("rstan", "cmdstanr") )

    reserved <- c(
        "formula", "data", "family", "warmup", "iter", "chains",
        "cores", "threads", "backend", "control", "file",
        "stan_model_args"
        )


    duplicates <- intersect(names(extra_args), reserved)

    if (length(duplicates) > 0L) {

        stop (
            "The following arguments cannot be passed through `...`: ",
            paste(duplicates, collapse = ", "),
            call. = FALSE
            )

    }

    reused <- !is.null(previous_model)

    if (reused) {

        if (!inherits(previous_model, "brmsfit") ) {

            stop (
                "`previous_model` must be a valid `brmsfit` object.",
                call. = FALSE
                )

        }

        model <- previous_model

        message (reuse_message)

    } else {
        if (!is.numeric(warmup) || length(warmup) != 1L ||
            !is.numeric(iter) || length(iter) != 1L ||
            !is.finite(warmup) || !is.finite(iter) || warmup >= iter) {

            stop ("`iter` must be strictly larger than `warmup`.", call. = FALSE)
        }

        if (is.null(stan_model_args) ) {

            stan_model_args <- if (identical(backend, "cmdstanr") ) {

                list(stanc_options = list("O1") )

            } else {

                list()

            }

        }

        brm_args <- c(
            list(
                formula = formula,
                data = data,
                family = family,
                warmup = warmup,
                iter = iter,
                chains = chains,
                cores = cores,
                threads = threads,
                backend = backend,
                control = stan_control,
                file = file,
                stan_model_args = stan_model_args
                ),
            extra_args
            )

        model <- do.call(fit_fun, brm_args)

    }

    if (isTRUE(resolve_draws) && is.null(n_post_samples) ) {

        n_post_samples <- ndraws_fun(model)

    }

    return (
        list(
            model = model,
            n_post_samples = n_post_samples,
            reused = reused
            )
        )

}
