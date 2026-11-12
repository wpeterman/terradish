# plot(fit) with neither `type` nor `data` used to stop, because the default
# type needs the graph. It now falls back to the fit plot with a message. A
# type that needs the graph and is requested explicitly must still stop.

test_that("plot(fit) without data falls back to the fit plot with a message", {
  fx <- fit_fixture(keep = 1:8, aggregate_factor = 4)

  expect_message(p <- plot(fx$fit), "type = \"fit\"")
  expect_s3_class(p, "ggplot")

  # Same object as an explicit request for the fit plot, which is silent.
  expect_silent(explicit <- plot(fx$fit, type = "fit"))
  expect_equal(p$data, explicit$data)
})

test_that("an explicitly requested type that needs data still stops", {
  fx <- fit_fixture(keep = 1:8, aggregate_factor = 4)

  expect_error(plot(fx$fit, type = "marginal_response"), "require `data`")
  expect_error(plot(fx$fit, type = "marginal"), "require `data`")
  expect_error(plot(fx$fit, type = "surface"), "requires `data`")
})

test_that("supplying data keeps the default marginal-response plot", {
  skip_on_cran()
  fx <- fit_fixture(keep = 1:8, aggregate_factor = 4)

  # No fallback message when the graph is supplied.
  expect_no_message(
    p <- plot(fx$fit, data = fx$surface, n = 5),
    message = "type = \"fit\""
  )
  expect_true(inherits(p, "ggplot") || inherits(p, "terradish_plot_list"))
})
