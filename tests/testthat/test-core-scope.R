test_that("screening stages require full-graph refinement", {
  unrefined <- list(exact_refine = FALSE)
  expect_error(.normalize_coarse_raster_control(unrefined),
               "not supported in the terradish core; available on the experimental branch",
               fixed = TRUE)
  expect_error(.normalize_landmark_control(unrefined, n_focal = 12L),
               "not supported in the terradish core; available on the experimental branch",
               fixed = TRUE)
  expect_true(.normalize_coarse_raster_control(NULL)$exact_refine)
  expect_true(.normalize_landmark_control(NULL, 12L)$exact_refine)
})

test_that("grid evaluation retains the full graph and response", {
  graph <- list(marker = "full graph")
  response <- diag(4)
  full <- .terradish_grid_approximation(graph, response)
  expect_identical(full$data, graph)
  expect_identical(full$S, response)
  expect_false(full$info$used)
  for (approximation in c("landmark", "coarse_raster")) {
    expect_error(.terradish_grid_approximation(graph, response, approximation),
                 "arg.*none")
  }
})
