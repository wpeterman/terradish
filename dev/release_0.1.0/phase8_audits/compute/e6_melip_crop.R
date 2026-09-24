source("C:/Users/peterman.73/OneDrive - The Ohio State University/R/Packages/terradish/dev/release_0.1.0/phase8_audits/compute/common.R")
data(melip)
alt <- unwrap(melip.altitude); fc <- unwrap(melip.forestcover); co <- unwrap(melip.coords)
cv <- c(alt, fc); names(cv) <- c("altitude", "forestcover"); cv <- scale_covariates(cv)
cat("melip raster:", dim(cv)[1:2], "res", res(cv), " extent", as.vector(ext(cv)), "\n")
cat("site bbox:", range(crds(co)[,1]), range(crds(co)[,2]), "\n")
sf <- conductance_surface(cv, co, directions = 8)
sc <- conductance_surface(cv, co, directions = 8, crop_buffer = 0.05)
cat("vertices full", nrow(sf$x), " cropped(0.05)", nrow(sc$x), "\n")
th <- matrix(c(0, 0), 1, dimnames = list(NULL, c("altitude", "forestcover")))
Rf <- terradish_distance(th, ~ altitude + forestcover, sf, loglinear_conductance)$distance[,,1]
Rc <- terradish_distance(th, ~ altitude + forestcover, sc, loglinear_conductance)$distance[,,1]
off <- upper.tri(Rf) & Rf > 0
cat("Zero-resistance duplicate-cell pairs excluded from relative changes:",
    sum(upper.tri(Rf) & Rf == 0), "\n")
cat(sprintf("IBD resistance: mean rel increase %.3f, max %.3f, cor %.4f\n", mean(Rc[off]/Rf[off]-1), max(Rc[off]/Rf[off]-1), cor(Rc[off], Rf[off])))
ff <- suppressWarnings(terradish(melip.Fst ~ altitude + forestcover, sf, measurement_model = mlpe))
fcr <- suppressWarnings(terradish(melip.Fst ~ altitude + forestcover, sc, measurement_model = mlpe))
cat("mlpe coef full   :", round(coef(ff), 3), " SE", round(sqrt(diag(vcov(ff))), 3), "\n")
cat("mlpe coef cropped:", round(coef(fcr), 3), " SE", round(sqrt(diag(vcov(fcr))), 3), "\n")

