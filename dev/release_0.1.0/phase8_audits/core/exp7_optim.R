source("C:/Users/peterman.73/OneDrive - The Ohio State University/R/Packages/terradish/dev/release_0.1.0/phase8_audits/core/helper.R")
L <- make_land(25, seed = 4, n_sites = 15)
surf <- L$surf
sim <- simulate_covariance_response(c(x1 = 0.7, x2 = -0.3), ~x1 + x2, surf, tau = 1, sigma = 0.2, nu = 200, seed = 1)
D <- dist_from_cov(sim$covariance)
starts <- list(c(0, 0), c(2, 2), c(-2, 2), c(-2, -2), c(3, -3))
for (st in starts) {
  f <- tryCatch(suppressWarnings(terradish(D ~ x1 + x2, surf, loglinear_conductance, mlpe, theta = st)), error = function(e) NULL)
  if (is.null(f)) { cat("start", st, ": ERROR\n"); next }
  cat(sprintf("start (%g,%g): theta=(%s) loglik=%.4f |grad|=%.2e steps=%d boundary=%s\n", st[1], st[2],
      if (is.null(coef(f))) "none" else paste(round(coef(f), 4), collapse = ","), f$loglik,
      if (is.null(f$mle$gradient)) NA else sqrt(sum(f$mle$gradient^2)), f$cost[1], f$fit$boundary))
}
for (opt in c("newton", "bfgs")) {
  f <- suppressWarnings(terradish(D ~ s(x1, df = 3) + x2, surf, smooth_loglinear_conductance, mlpe, optimizer = opt))
  cat(opt, "spline: loglik", round(f$loglik, 4), "|grad|", signif(sqrt(sum(f$mle$gradient^2)), 3), "steps", f$cost[1], "coef", round(coef(f), 3), "\n")
}
# Wishart with large nu: gradient scale vs ctol
f <- terradish(sim$covariance ~ x1 + x2, surf, loglinear_conductance, wishart_covariance, nu = 20000)
cat("WC nu=2e4: |grad| =", signif(sqrt(sum(f$mle$gradient^2)), 3), " steps", f$cost[1], "\n")
f <- terradish(sim$covariance ~ x1 + x2, surf, loglinear_conductance, wishart_covariance, nu = 20)
cat("WC nu=20: |grad| =", signif(sqrt(sum(f$mle$gradient^2)), 3), " steps", f$cost[1], "coef", coef(f), "\n")
