# Precompute for the directional vignette: model-based recovery on a small,
# well-conditioned lattice (reliable; the coalescent-data difficulty is discussed
# separately and honestly in the vignette text).
suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressMessages(library(terra)); suppressMessages(library(numDeriv))
set.seed(1)

DIM <- 7L
r <- rast(nrows = DIM, ncols = DIM, xmin = 0, xmax = DIM, ymin = 0, ymax = DIM)
gx <- xFromCell(r, seq_len(ncell(r))); gy <- yFromCell(r, seq_len(ncell(r)))
v1_values <- scale(gx + 0.5 * gy)[, 1]
elev_values <- scale(gx + 0.35 * gy +
                       0.35 * sin(pi * gx / DIM) -
                       0.25 * cos(pi * gy / DIM))[, 1]
covs <- c(setValues(r, v1_values), setValues(r, elev_values))
names(covs) <- c("v1", "elev")
fc <- c(1L, DIM, DIM * DIM, DIM * (DIM - 1L) + 1L, (DIM * DIM) %/% 2L, DIM * 3L + 2L)
coords <- xyFromCell(r, fc)
surface <- conductance_surface(covs, coords, directions = 8)
dir_cov <- edge_gradient(covs[["elev"]], surface)
dir_cov$d <- cbind(elev = dir_cov$d)
gen <- terradish:::.directed_generator(~ v1, surface, dir_cov)
nf <- length(surface$demes)

theta_true <- 0.5; gamma_true <- 0.6; nu <- 5000
E <- terradish_directed_algorithm(gen, NULL, surface, NULL, par = c(theta_true, gamma_true))$covariance
set.seed(7); S <- rWishart(1, df = nu, Sigma = (E + 0.2 * diag(nf)) / nu)[, , 1]

# numDeriv gradient check (engine)
par <- c(0.3, 0.4)
an <- terradish_directed_algorithm(gen, wishart_covariance, surface, S, par, nu = nu, gradient = TRUE)$gradient
nd <- numDeriv::grad(function(pp) terradish_directed_algorithm(gen, wishart_covariance, surface, S, pp, nu = nu, gradient = FALSE)$objective, par)
grad_err <- max(abs(an - nd))

fit <- terradish_directed(S ~ v1, data = surface, directional = dir_cov,
                          measurement_model = wishart_covariance, nu = nu)
fit_reversible <- terradish_directed(S ~ v1, data = surface,
                                     directional = dir_cov,
                                     measurement_model = wishart_covariance,
                                     nu = nu, gamma_bound = 0)
est <- c(fit$theta, fit$gamma)
coef_table <- data.frame(
  parameter = c("v1", "gamma_elev"),
  symbol = c("theta", "gamma"),
  truth = c(theta_true, gamma_true),
  estimate = unname(est[c("v1", "gamma_elev")]),
  std_error = unname(fit$se[c("v1", "gamma_elev")]),
  row.names = NULL
)
coef_table$z_value <- coef_table$estimate / coef_table$std_error
coef_table$p_value <- pmin(2 * (1 - pnorm(abs(coef_table$z_value))), 1)

lrt_stat <- 2 * (as.numeric(logLik(fit)) - as.numeric(logLik(fit_reversible)))
lrt_df <- length(fit$gamma)
lrt_table <- data.frame(
  comparison = "directed vs. reversible",
  logLik_directed = as.numeric(logLik(fit)),
  logLik_reversible = as.numeric(logLik(fit_reversible)),
  statistic = lrt_stat,
  df = lrt_df,
  p_value = stats::pchisq(lrt_stat, df = lrt_df, lower.tail = FALSE)
)

edge_rates <- directed_rates(fit, data = surface, directional = dir_cov)
edge_summary <- data.frame(
  quantity = c(
    "Number of graph edges",
    "Median absolute log rate ratio",
    "Maximum absolute log rate ratio",
    "Median rate multiplier in favored direction",
    "Maximum rate multiplier in favored direction"
  ),
  value = c(
    nrow(edge_rates),
    median(edge_rates$abs_log_rate_ratio),
    max(edge_rates$abs_log_rate_ratio),
    median(exp(edge_rates$abs_log_rate_ratio)),
    max(exp(edge_rates$abs_log_rate_ratio))
  )
)

res <- list(
  fit = fit,
  example = list(DIM = DIM, focal_cells = fc,
                 v1_values = v1_values, elev_values = elev_values),
  coef_table = coef_table,
  lrt_table = lrt_table,
  edge_rates = edge_rates,
  edge_summary = edge_summary,
  grad_err = grad_err,
  tau = unname(fit$phi["tau"]),
  loglik = fit$loglik,
  truth = c(theta = theta_true, gamma = gamma_true),
  nu = nu,
  dim = c(vertices = nrow(surface$x), focal = length(surface$demes),
          edges = nrow(surface$edge_pairs))
)
saveRDS(res, "vignettes/vignette-directional.rds")
cat(sprintf("grad_err=%.2e  tau=%.3f\n", grad_err, res$tau))
print(coef_table, digits = 3)
print(lrt_table, digits = 3)
print(edge_summary, digits = 3)
