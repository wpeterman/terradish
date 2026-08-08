# The `radish_*` wrappers exist so code written against the old package keeps
# working. Nothing else in the suite calls them, so a mistake in
# `.terradish_forward_call()` would silently break every legacy script without
# failing a test. These check both halves of the contract: the deprecation
# warning fires, and the call really does reach the new function with the
# arguments intact.

test_that("radish() warns and forwards to terradish()", {
  res <- fit_fixture(keep = 1:8)
  surface <- res$surface
  melip.Fst <- res$data$melip.Fst

  expect_warning(
    legacy <- radish(melip.Fst ~ altitude + forestcover,
                     data = surface,
                     conductance_model = loglinear_conductance,
                     measurement_model = leastsquares,
                     control = NewtonRaphsonControl(maxit = 3, verbose = FALSE)),
    "deprecated"
  )
  modern <- suppressWarnings(
    terradish(melip.Fst ~ altitude + forestcover,
              data = surface,
              conductance_model = loglinear_conductance,
              measurement_model = leastsquares,
              control = NewtonRaphsonControl(maxit = 3, verbose = FALSE))
  )
  expect_equal(unname(coef(legacy)), unname(coef(modern)))
  expect_equal(unname(legacy$loglik), unname(modern$loglik))
  expect_s3_class(legacy, "terradish")
})

test_that("radish_algorithm() warns and forwards to terradish_algorithm()", {
  res <- fit_fixture(keep = 1:8)
  surface <- res$surface
  S <- ifelse(res$data$melip.Fst < 0, 0, res$data$melip.Fst)
  f <- loglinear_conductance(~ altitude + forestcover, surface$x)

  expect_warning(
    legacy <- radish_algorithm(f, leastsquares, surface, S, c(-0.3, 0.3),
                               partial = FALSE),
    "deprecated"
  )
  modern <- terradish_algorithm(f, leastsquares, surface, S, c(-0.3, 0.3),
                                partial = FALSE)
  expect_equal(legacy$objective, modern$objective)
  expect_equal(legacy$gradient, modern$gradient)
  expect_equal(legacy$hessian, modern$hessian)
})

test_that("radish_distance() and radish_grid() warn and forward", {
  res <- fit_fixture(keep = 1:8)
  surface <- res$surface
  melip.Fst <- res$data$melip.Fst
  theta <- matrix(c(-0.3, 0.3), nrow = 1,
                  dimnames = list(NULL, c("altitude", "forestcover")))

  expect_warning(
    legacy_dist <- radish_distance(theta, ~ altitude + forestcover,
                                   data = surface,
                                   conductance_model = loglinear_conductance),
    "deprecated"
  )
  modern_dist <- terradish_distance(theta, ~ altitude + forestcover,
                                    data = surface,
                                    conductance_model = loglinear_conductance)
  expect_equal(legacy_dist$distance, modern_dist$distance)

  expect_warning(
    legacy_grid <- radish_grid(theta, melip.Fst ~ altitude + forestcover,
                               data = surface,
                               conductance_model = loglinear_conductance,
                               measurement_model = leastsquares),
    "deprecated"
  )
  modern_grid <- terradish_grid(theta, melip.Fst ~ altitude + forestcover,
                                data = surface,
                                conductance_model = loglinear_conductance,
                                measurement_model = leastsquares)
  expect_equal(legacy_grid$loglik, modern_grid$loglik)
})

test_that(".terradish_forward_call rewrites the function but keeps the arguments", {
  target <- function(a, b = 2) list(a = a, b = b)
  caller <- function(a, b) terradish:::.terradish_forward_call(match.call(), "target")
  # force the calls outside expect_*() so the frame the rewritten call is
  # evaluated in is this test block, not a promise inside the expectation
  both <- caller(a = 5, b = 7)
  defaulted <- caller(a = 5)
  expect_equal(both, list(a = 5, b = 7))
  expect_equal(defaulted, list(a = 5, b = 2))
})

test_that("every legacy wrapper is exported and deprecated", {
  legacy <- c("radish", "radish_algorithm", "radish_cv", "radish_distance",
              "radish_grid", "radish_multiscale", "radish_parameters")
  exported <- getNamespaceExports(asNamespace("terradish"))
  expect_true(all(legacy %in% exported))
  # each body must route through the shared deprecation helper, so none can
  # drift into a silent alias
  for (fn in legacy) {
    body_text <- paste(deparse(body(get(fn, envir = asNamespace("terradish")))),
                       collapse = " ")
    expect_true(grepl(".terradish_deprecate", body_text, fixed = TRUE),
                info = fn)
  }
})
