library(terradish)
library(terra)
data(melip)

# Unwrap the portable example data and standardize each raster layer.
altitude <- terra::unwrap(melip.altitude)
forestcover <- terra::unwrap(melip.forestcover)
sites <- terra::unwrap(melip.coords)
covariates <- c(altitude, forestcover)
names(covariates) <- c("altitude", "forestcover")
covariates <- scale_covariates(covariates)
graph <- conductance_surface(covariates, sites)  # eight neighbors by default

# MLPE accounts for shared sampling sites among pairwise genetic distances.
fit <- terradish(melip.Fst ~ forestcover + altitude, graph,
                 measurement_model = mlpe)
fit$convergence
summary(fit)
coef(fit)
confint(fit)
vcov(fit)
plot(fit, type = "surface", data = graph)
plot(fit, type = "fit")


# Raster covariates come first; coordinates come second.
environment <- pairwise_endpoint_covariates(covariates, sites)
lonlat <- terra::crds(terra::project(sites, "EPSG:4326"))
local_crs <- sprintf("+proj=aeqd +lat_0=%f +lon_0=%f +datum=WGS84 +units=m",
                     mean(lonlat[, 2]), mean(lonlat[, 1]))
projected_sites <- terra::project(sites, local_crs)
pairwise <- pairwise_covariates(environment,
  geographic_100km = dist(terra::crds(projected_sites)) / 100000)

joint <- terradish(melip.Fst ~ forestcover + altitude, graph,
                   measurement_model = mlpe_covariates(pairwise))
terradish_ibe_ratio(joint)

# The same object can enter a model for a coherent covariance response.
wishart_measurement <- wishart_covariates(pairwise, model = "wishart_covariance")

