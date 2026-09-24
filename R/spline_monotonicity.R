.spline_monotonicity <- function(model, theta, focal_x) {
  info <- attr(model, "smooth_loglinear_info", exact = TRUE)
  if (is.null(info) || !length(info$smooth_specs)) return(NULL)
  rows <- lapply(info$smooth_specs, function(spec) {
    limits <- range(focal_x[[spec$variable]], finite = TRUE)
    cuts <- sort(unique(c(limits, spec$knots[spec$knots > limits[1] & spec$knots < limits[2]],
      spec$Boundary.knots[spec$Boundary.knots > limits[1] & spec$Boundary.knots < limits[2]])))
    degree <- if (spec$basis == "ns") 3L else spec$degree
    slopes <- numeric()
    for (j in seq_len(max(0L, length(cuts) - 1L))) {
      span <- diff(cuts[c(j, j + 1L)])
      t <- seq(0, 1, length.out = degree + 1L)
      grid <- setNames(data.frame(cuts[j] + span * t), spec$variable)
      B <- .smooth_loglinear_basis(spec$label, grid, spec$df, spec$basis,
        spec$degree, spec$intercept, parent.frame(), spec)
      values <- as.vector(B %*% theta[spec$columns])
      coefficients <- qr.solve(outer(t, 0:degree, `^`), values)
      derivative <- coefficients[-1] * seq_len(degree)
      roots <- if (any(abs(derivative) > 1e-12)) polyroot(derivative) else complex()
      roots <- Re(roots[abs(Im(roots)) < 1e-8 & Re(roots) > 0 & Re(roots) < 1])
      bounds <- sort(unique(c(0, roots, 1)))
      mid <- (utils::head(bounds, -1L) + utils::tail(bounds, -1L)) / 2
      slopes <- c(slopes, as.vector(outer(mid, 0:(degree - 1L), `^`) %*% derivative) / span)
    }
    tolerance <- 1e-8 * max(1, abs(slopes))
    signs <- sign(slopes[abs(slopes) > tolerance])
    changes <- if (length(signs) > 1L) sum(diff(signs) != 0) else 0L
    data.frame(term = spec$label, variable = spec$variable,
      lower = limits[1], upper = limits[2], monotone = changes == 0,
      direction = if (!length(signs)) "constant" else if (changes) "nonmonotone" else
        if (signs[1] > 0) "increasing" else "decreasing",
      derivative_sign_changes = changes)
  })
  do.call(rbind, rows)
}
