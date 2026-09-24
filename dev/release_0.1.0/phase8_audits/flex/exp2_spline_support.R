source("C:/Users/peterman.73/OneDrive - The Ohio State University/R/Packages/terradish/dev/release_0.1.0/phase8_audits/flex/common.R")
L <- make_land(25, seed = 3, nsite = 16)
s <- L$surface
# non-monotone truth in x1
logc <- -0.8 * s$x$x1^2 + 0.3 * s$x$x2
S <- sim_S(s, logc - mean(logc), nu = 500)
fit <- suppressWarnings(terradish(S ~ x2 + s(x1, df = 4), data = s,
          conductance_model = smooth_loglinear_conductance,
          measurement_model = wishart_covariance, nu = 500, leverage = FALSE,
          control = NewtonRaphsonControl(maxit = 50, verbose = FALSE)))
print(round(coef(fit), 3))
c_none  <- conductance(s, fit)
# clamp x1 to the FULL graph range: clamping is then a no-op on covariate values
rng_focal <- range(s$x$x1[s$demes]); rng_all <- range(s$x$x1)
cat("focal x1 range:", round(rng_focal,3), " graph x1 range:", round(rng_all,3), "\n")
c_focal <- conductance(s, fit, support = "focal", clamp_covariates = "x1")
v_none <- values(c_none$est)[,1]; v_focal <- values(c_focal$est)[,1]
inside <- s$x$x1 >= rng_focal[1] & s$x$x1 <= rng_focal[2]
cat("cells whose x1 is unchanged by clamping:", sum(inside), "of", length(inside), "\n")
cat("max |log ratio| of est on UNCHANGED cells (should be 0):",
    signif(max(abs(log(v_focal[!is.na(v_focal)][inside] / v_none[!is.na(v_none)][inside]))), 4), "\n")
cat("cor(log est none, log est focal) on unchanged cells:",
    round(cor(log(v_none[!is.na(v_none)][inside]), log(v_focal[!is.na(v_focal)][inside])), 4), "\n")
# the same check for a loglinear fit (control: should be exactly 0)
fitl <- suppressWarnings(terradish(S ~ x1 + x2, data = s, conductance_model = loglinear_conductance,
          measurement_model = wishart_covariance, nu = 500, leverage = FALSE,
          control = NewtonRaphsonControl(maxit = 50, verbose = FALSE)))
a <- values(conductance(s, fitl)$est)[,1]; b <- values(conductance(s, fitl, support="focal", clamp_covariates="x1")$est)[,1]
a <- a[!is.na(a)]; b <- b[!is.na(b)]
cat("loglinear control max |log ratio| on unchanged cells:", signif(max(abs(log(b[inside]/a[inside]))), 4), "\n")

# plot(): automatic vs user-supplied factory (as documented in ?plot.terradish)
p_auto <- plot(fit, type = "marginal", data = s, n = 25, marginal_covariates = "x1")
p_user <- plot(fit, type = "marginal", data = s, n = 25, marginal_covariates = "x1",
               conductance_model = smooth_loglinear_conductance)
d1 <- p_auto$data; d2 <- p_user$data
cat("marginal curve, auto   :", round(d1$est[c(1,7,13,19,25)], 3), "\n")
cat("marginal curve, user f :", round(d2$est[c(1,7,13,19,25)], 3), "\n")
# truth-shaped reference at same x: exp(-0.8 x^2) up to scale
