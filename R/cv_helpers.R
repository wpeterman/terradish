.terradish_fit_converged <- function(fit)
{
  # Outer stationarity does not establish convergence of the profiled model.
  # Missing legacy diagnostics cannot certify a successful fit either.
  isTRUE(fit$convergence$code == 0L) &&
    isTRUE(fit$fit$subproblem$convergence == 0L)
}

.terradish_measurement_model <- function(model, subset = NULL)
{
  if (inherits(model, c("terradish_measurement_model",
                        "radish_measurement_model")))
  {
    subsetter <- attr(model, "subsetter", exact = TRUE)
    if (!is.null(subset) && is.function(subsetter))
      return(subsetter(subset))
    return(model)
  }

  if (!is.character(model) || length(model) != 1L)
    stop("`model` must be a measurement-model function or one of 'mlpe', 'wishart', or 'ls'")

  switch(tolower(model),
         mlpe = mlpe,
         wishart = generalized_wishart,
         generalized_wishart = generalized_wishart,
         ls = leastsquares,
         leastsquares = leastsquares,
         stop("Unknown measurement model: ", model))
}

.fit_terradish_with_fallback <- function(formula, data, measurement_model, nu = NULL,
                                         conductance_model = loglinear_conductance,
                                         theta = NULL, cores = 1L, dots = list(),
                                         response_matrix = NULL)
{
  # An explicit optimizer must be forwarded once. The default still tries
  # Newton first and falls back to BFGS only when that call raises an error.
  requested_optimizer <- dots$optimizer
  dots$optimizer <- NULL
  fit_once <- function(optimizer)
  {
    args <- c(list(formula = formula,
                   data = data,
                   conductance_model = conductance_model,
                   measurement_model = measurement_model,
                   nu = nu,
                   optimizer = optimizer,
                   cores = cores),
              dots)
    if (!is.null(theta))
      args$theta <- theta
    eval_env <- list2env(list(gd_mat = response_matrix), parent = parent.frame())
    call <- as.call(c(list(as.name("terradish")), args))
    tryCatch(eval(call, envir = eval_env), error = identity)
  }

  fit <- fit_once(if (is.null(requested_optimizer)) "newton" else requested_optimizer)
  if (is.null(requested_optimizer) && inherits(fit, "error"))
    fit <- fit_once("bfgs")
  if (inherits(fit, "error"))
    stop("Could not optimize terradish model: ", conditionMessage(fit), call. = FALSE)

  fit
}
