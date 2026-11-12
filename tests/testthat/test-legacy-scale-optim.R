# terra_radish_scale_optim() was a silent alias. It now behaves like the other
# legacy wrappers: it signals a deprecation warning and forwards its arguments
# unchanged to terradish_scale_optim().

test_that("terra_radish_scale_optim() warns and forwards to terradish_scale_optim()", {
  dat <- melip_fixture(1:6, aggregate_factor = 4)
  melip.Fst <- dat$melip.Fst
  covariates <- dat$covariates[["altitude"]]
  names(covariates) <- "altitude"
  # A trivial scale function keeps the search to two cheap refits.
  identity_scale <- function(x, scale) x
  arguments <- list(
    melip.Fst ~ altitude,
    covariates = covariates,
    coords = dat$coords,
    scale_grid = list(altitude = c(1, 2)),
    search = "grid",
    scale_fun = identity_scale,
    measurement_model = leastsquares,
    control = NewtonRaphsonControl(maxit = 20, verbose = FALSE)
  )

  legacy <- expect_fixture_warnings(
    do.call(terra_radish_scale_optim, arguments),
    "deprecated"
  )
  modern <- do.call(terradish_scale_optim, arguments)

  expect_s3_class(legacy, "terradish_scale_optim")
  expect_identical(legacy$par, modern$par)
  expect_identical(legacy$value, modern$value)
  expect_identical(coef(legacy$fit), coef(modern$fit))
})

test_that("the deprecation warning names the replacement", {
  expect_warning(
    tryCatch(terra_radish_scale_optim(), error = function(e) NULL),
    "terradish_scale_optim"
  )
})
