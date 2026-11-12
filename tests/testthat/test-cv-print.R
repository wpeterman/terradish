# print() for a cross-validation result names each measurement model once,
# however many candidate formulas share it.

test_that("print.terradish_cv_folds names each measurement model once", {
  skip_on_cran()
  dat <- melip_fixture(1:9)
  melip.Fst <- dat$melip.Fst
  surface <- conductance_surface(dat$covariates[["altitude"]], dat$coords)
  folds <- terradish_folds(dat$coords, k = 3, method = "random", seed = 1)

  cv <- suppressWarnings(terradish_cv_folds(
    surface,
    list(uniform = melip.Fst ~ 1, altitude = melip.Fst ~ altitude),
    folds = folds, model = "ls"))

  printed <- capture.output(returned <- print(cv))
  model_line <- grep("^Measurement model:", printed, value = TRUE)
  expect_length(model_line, 1L)
  expect_identical(trimws(model_line), "Measurement model: leastsquares")
  expect_identical(returned, cv)

  # The summary table follows, one row per candidate formula.
  expect_true(any(grepl("uniform", printed)))
  expect_true(any(grepl("altitude", printed)))
})

test_that("the printed model line collapses repeated names", {
  # A hand-built object isolates the print method from fitting.
  cv <- structure(
    list(measurement_model = c("mlpe", "mlpe", "mlpe"), repeats = 1L,
         nuisance = "reprofile",
         summary = data.frame(model = c("a", "b", "c"), folds = 3L)),
    class = "terradish_cv_folds")
  printed <- capture.output(print(cv))
  expect_identical(trimws(grep("^Measurement model:", printed, value = TRUE)),
                   "Measurement model: mlpe")
})
