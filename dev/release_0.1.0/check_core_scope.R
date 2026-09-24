# Phase 1b scope check. Exact exported names avoid matching the retained
# terradish_cv_folds through the legacy radish_cv substring.
pkgload::load_all()
removed <- c("terradish_directed", "directed_rates", "edge_gradient", "edge_flow",
  "terradish_hierarchical", "conductance_field", "wishart_drift_covariates",
  "linear_conductance", "pair_subset_measurement_model", "terradish_multiscale",
  "radish_multiscale", "terradish_kron_reduce", "terradish_kron_reduce_tiled",
  "terradish_cv", "terradish_cv_replicates", "radish_cv", "cv_model_selection",
  "terradish_results", "terradish_parameters", "radish_parameters")
stopifnot(!any(removed %in% getNamespaceExports("terradish")))
ns <- asNamespace("terradish")
for (name in c("terradish", "terradish_algorithm", "covariance_response_power")) {
  choices <- eval(formals(get(name, ns))$solver)
  stopifnot(setequal(choices, c("direct", "auto", "amg")))
}
stopifnot(identical(eval(formals(terradish_grid)$approximation), "none"))
stopifnot(!any(c("block_cg_reduced_laplacian", "pcg_reduced_laplacian",
                 "pcg_reduced_laplacian_ic") %in% getNamespaceExports("terradish")))
cat("Removed exports absent; retained solver choices and exact grids verified.\n")
