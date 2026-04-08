toy_scale_fun <- function(r, scale)
{
  if (scale <= 1.5)
    return(r)
  terra::focal(r, w = matrix(1, 3, 3), fun = mean, na.rm = TRUE)
}

test_that("terradish_scale_optim evaluates a full grid of candidate scales", {
  dat <- melip_fixture(1:6)
  melip.Fst <- dat$melip.Fst

  opt <- terradish_scale_optim(
    melip.Fst ~ altitude + forestcover,
    covariates = dat$covariates,
    coords = dat$coords,
    lower = c(altitude = 1, forestcover = 1),
    upper = c(altitude = 2, forestcover = 2),
    search = "grid",
    grid_points = 2,
    scale_fun = toy_scale_fun,
    postprocess = scale_covariates,
    measurement_model = leastsquares,
    optimizer = "bfgs",
    control = NewtonRaphsonControl(maxit = 2, verbose = FALSE)
  )

  expect_s3_class(opt, "terradish_scale_optim")
  expect_equal(nrow(opt$evaluations), 4)
  expect_true(all(c("altitude", "forestcover", "objective", "aic", "loglik", "status") %in%
                    names(opt$evaluations)))
  expect_true(all(opt$par >= 1 & opt$par <= 2))
  expect_s3_class(opt$fit, "terradish")
  expect_equal(opt$search, "grid")
})

test_that("terra_radish_scale_optim alias supports coordinate search", {
  dat <- melip_fixture(1:6)
  melip.Fst <- dat$melip.Fst
  covariates <- dat$covariates[["altitude"]]
  names(covariates) <- "altitude"

  opt <- terra_radish_scale_optim(
    melip.Fst ~ altitude,
    covariates = covariates,
    coords = dat$coords,
    scales = c(altitude = 1.5),
    lower = c(altitude = 1),
    upper = c(altitude = 2),
    search = "coordinate",
    maxit = 1,
    tol = 0,
    scale_fun = toy_scale_fun,
    postprocess = scale_covariates,
    measurement_model = leastsquares,
    optimizer = "bfgs",
    control = NewtonRaphsonControl(maxit = 2, verbose = FALSE),
    verbose = FALSE
  )

  expect_true(opt$par[["altitude"]] >= 1)
  expect_true(opt$par[["altitude"]] <= 2)
  expect_gte(nrow(opt$evaluations), 2)
  expect_equal(opt$search, "coordinate")
  expect_true(is.finite(opt$value))
})
