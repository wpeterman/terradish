melip_fixture <- function(keep = NULL)
{
  data(melip, package = "terradish", envir = environment())
  melip.altitude <- terra::unwrap(melip.altitude)
  melip.forestcover <- terra::unwrap(melip.forestcover)
  melip.coords <- terra::unwrap(melip.coords)

  if (!is.null(keep))
  {
    melip.Fst <- melip.Fst[keep, keep, drop = FALSE]
    melip.coords <- melip.coords[keep]
  }

  covariates <- c(melip.altitude, melip.forestcover)
  names(covariates) <- c("altitude", "forestcover")
  covariates <- scale_covariates(covariates)

  list(melip.Fst = melip.Fst,
       covariates = covariates,
       coords = melip.coords)
}

fit_fixture <- function(keep = 1:12,
                        formula = melip.Fst ~ altitude + forestcover,
                        measurement_model = leastsquares,
                        ...)
{
  dat <- melip_fixture(keep)
  melip.Fst <- dat$melip.Fst
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  fit <- suppressWarnings(
    terradish(formula,
              data = surface,
              conductance_model = loglinear_conductance,
              measurement_model = measurement_model,
              ...)
  )

  list(data = dat, surface = surface, fit = fit)
}
