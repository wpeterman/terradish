# "Removed unused spatial covariates" is raised once per user-level call, not
# once for every internal refit. A direct call still sees it.

collect_warnings <- function(expr) {
  seen <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    seen <<- c(seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = seen)
}

unused_surface <- function() {
  dat <- melip_fixture(1:12, aggregate_factor = 3)
  keep <- !duplicated(terra::cellFromXY(dat$covariates,
                                        terra::crds(dat$coords)))
  list(Fst = dat$melip.Fst[keep, keep],
       surface = conductance_surface(dat$covariates, dat$coords[keep]))
}

test_that("a direct terradish() call warns once about unused covariates", {
  fx <- unused_surface()
  Fst <- fx$Fst
  out <- collect_warnings(
    terradish(Fst ~ altitude, fx$surface, measurement_model = leastsquares))
  unused <- grep("^Removed unused spatial covariates", out$warnings,
                 value = TRUE)
  expect_identical(unused, "Removed unused spatial covariates: forestcover")
})

test_that("a conductance-model factory called directly always warns", {
  fx <- unused_surface()
  for (i in 1:2)
    expect_warning(loglinear_conductance(~ altitude, fx$surface$x),
                   "^Removed unused spatial covariates: forestcover$")
})

test_that("cross-validation warns once per distinct set of unused covariates", {
  skip_on_cran()
  fx <- unused_surface()
  Fst <- fx$Fst
  folds <- rep(1:3, length.out = nrow(Fst))
  out <- collect_warnings(terradish_cv_folds(
    fx$surface,
    list(altitude = Fst ~ altitude, forest = Fst ~ forestcover),
    folds = folds, model = "ls"))
  unused <- grep("^Removed unused spatial covariates", out$warnings,
                 value = TRUE)
  # Three folds, two formulas, and baseline refits used to repeat each
  # message many times.
  expect_setequal(unused, c("Removed unused spatial covariates: forestcover",
                            "Removed unused spatial covariates: altitude"))
  expect_length(unused, 2L)
  expect_s3_class(out$value, "terradish_cv_folds")
})

test_that("the warning scope closes after a call, including after an error", {
  fx <- unused_surface()
  Fst <- fx$Fst
  state <- .unused_covariate_state

  suppressWarnings(
    terradish(Fst ~ altitude, fx$surface, measurement_model = leastsquares))
  expect_identical(state$depth, 0L)
  expect_length(state$seen, 0L)

  # An error inside the scoped call must not leave the scope open.
  expect_error(
    terradish(Fst ~ altitude, fx$surface, measurement_model = leastsquares,
              conductance = FALSE))
  expect_identical(state$depth, 0L)

  # A second user-level call warns again.
  expect_warning(
    terradish(Fst ~ altitude, fx$surface, measurement_model = leastsquares),
    "^Removed unused spatial covariates: forestcover$")
})
