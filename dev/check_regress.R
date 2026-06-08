suppressMessages(devtools::load_all(".", quiet = TRUE))
set.seed(42)
fail <- 0
chk <- function(label, cond) {
  cat(sprintf("%-55s %s\n", label, ifelse(isTRUE(cond), "OK", "**FAIL**")))
  if (!isTRUE(cond)) fail <<- fail + 1
}
approx <- function(a, b, tol = 1e-8) max(abs(a - b)) < tol

E <- matrix(c(1.2, 0.3, 0.2,
              0.3, 1.5, 0.4,
              0.2, 0.4, 1.1), nrow = 3, byrow = TRUE)

## --- regression: wishart_covariates covariance interface + base match ---
site_env <- data.frame(altitude = c(-1, 0, 1))
g <- wishart_covariates(site_env, model = "wishart_covariance")
K <- attr(g, "kernel_covariates")[, , 1]
Sigma <- 0.8 * E + 0.3 * K + 0.2 * diag(3)
S <- rWishart(1, df = 25, Sigma = Sigma)[, , 1] / 25
start <- g(E, S, nu = 25)
fit <- g(E, S, phi = start$phi, nu = 25)
chk("wc: phi names", identical(names(start$phi), c("tau", "lambda_altitude", "sigma")))
chk("wc: finite objective", is.finite(fit$objective))
chk("wc: hessian dim", all(dim(fit$hessian) == c(3, 3)))

# zero kernel weight == base wishart_covariance
wrapped <- g(E, S, phi = c(tau = 0.8, lambda_altitude = 0, sigma = log(0.2)), nu = 25)
base <- wishart_covariance(E, S, phi = c(tau = 0.8, sigma = log(0.2)), nu = 25)
chk("wc: objective == base", approx(wrapped$objective, base$objective))
chk("wc: gradient_E == base", approx(wrapped$gradient_E, base$gradient_E))
chk("wc: hessian sane (no NaN)", all(is.finite(wrapped$hessian)))

## --- regression: wishart_covariates generalized ---
Sd <- diag(Sigma) %*% t(rep(1, 3)) + rep(1, 3) %*% t(diag(Sigma)) - 2 * Sigma
gg <- wishart_covariates(site_env, model = "generalized_wishart")
wrapped_g <- gg(E, Sd, phi = c(tau = 0.8, lambda_altitude = 0, sigma = log(0.2)), nu = 25)
base_g <- generalized_wishart(E, Sd, phi = c(tau = 0.8, sigma = log(0.2)), nu = 25)
chk("wc gen: objective == base", approx(wrapped_g$objective, base_g$objective))
chk("wc gen: fitted == base", approx(wrapped_g$fitted, base_g$fitted))

## --- regression: wishart_covariance numDeriv (still correct after refactor) ---
suppressMessages(library(numDeriv))
phi <- c(tau = 0.7, sigma = log(0.25))
fo <- function(p) wishart_covariance(E, S, phi = p, nu = 25,
                                     gradient = FALSE, hessian = FALSE, partial = FALSE)$objective
wc <- wishart_covariance(E, S, phi = phi, nu = 25)
chk("base wc: grad vs numDeriv", approx(c(wc$gradient), numDeriv::grad(fo, phi), 1e-6))
chk("base wc: hess vs numDeriv", approx(wc$hessian, numDeriv::hessian(fo, phi), 1e-4))

cat(sprintf("\n%d failure(s)\n", fail))
