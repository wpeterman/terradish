# Precompute for the directional vignette: model-based recovery on a small,
# well-conditioned lattice (reliable; the coalescent-data difficulty is discussed
# separately and honestly in the vignette text).
suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressMessages(library(terra)); suppressMessages(library(numDeriv))
set.seed(1)

DIM <- 7L
r <- rast(nrows = DIM, ncols = DIM, xmin = 0, xmax = DIM, ymin = 0, ymax = DIM)
gx <- xFromCell(r, seq_len(ncell(r))); gy <- yFromCell(r, seq_len(ncell(r)))
covs <- c(setValues(r, scale(gx + 0.5 * gy)[, 1]), setValues(r, scale(gx)[, 1]))
names(covs) <- c("v1", "elev")
fc <- c(1L, DIM, DIM * DIM, DIM * (DIM - 1L) + 1L, (DIM * DIM) %/% 2L, DIM * 3L + 2L)
coords <- xyFromCell(r, fc)
surface <- conductance_surface(covs, coords, directions = 8)
dir_cov <- edge_gradient(covs[["elev"]], surface)
gen <- terradish:::.directed_generator(~ v1, surface, dir_cov)
nf <- length(surface$demes)

theta_true <- 0.5; gamma_true <- 0.6; nu <- 1500
E <- terradish_directed_algorithm(gen, NULL, surface, NULL, par = c(theta_true, gamma_true))$covariance
set.seed(7); S <- rWishart(1, df = nu, Sigma = (E + 0.2 * diag(nf)) / nu)[, , 1]

# numDeriv gradient check (engine)
par <- c(0.3, 0.4)
an <- terradish_directed_algorithm(gen, wishart_covariance, surface, S, par, nu = nu, gradient = TRUE)$gradient
nd <- numDeriv::grad(function(pp) terradish_directed_algorithm(gen, wishart_covariance, surface, S, pp, nu = nu, gradient = FALSE)$objective, par)
grad_err <- max(abs(an - nd))

fit <- terradish_directed(S ~ v1, data = surface, directional = dir_cov,
                          measurement_model = wishart_covariance, nu = nu)
coef_table <- rbind(true = c(v1 = theta_true, gamma_elev = gamma_true),
                    estimated = c(fit$theta[1], fit$gamma[1]))

res <- list(coef_table = coef_table, grad_err = grad_err,
            tau = unname(fit$phi["tau"]), loglik = fit$loglik)
saveRDS(res, "vignettes/vignette-directional.rds")
cat(sprintf("grad_err=%.2e  tau=%.3f\n", grad_err, res$tau))
print(round(coef_table, 3))
