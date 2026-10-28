# Local teaching utilities, not additions to the public package API.
comparison_table <- function(comparisons) {
  out <- lapply(comparisons, function(x) data.frame(
    Comparison = paste(if (x$arm == "compensatory_500") "Compensatory" else "Control",
      if (x$panel == "SNP2500") "SNPs" else "microsatellites"),
    "Error reduction (%)" = 100 * x$RMSE_reduction,
    "Covariance score gain" = x$total_score_gain,
    "Positive folds" = paste(x$positive_folds, "of 4"), check.names = FALSE))
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

curve_from_record <- function(record, D, metadata) {
  theta <- unlist(record$theta)
  exp(theta["drying1"] * (D - metadata$center) / metadata$sd_linear +
      theta["drying2"] * (D - metadata$center)^2 / metadata$sd_quadratic)
}

predict_from_record <- function(record, graph) {
  theta <- unlist(record$theta)
  # Both zero coefficients give uniform conductance on the same graph.
  if (record$model == "uniform") theta <- c(drying1 = 0, drying2 = 0)
  columns <- names(theta)
  E <- terradish_distance(
    matrix(unname(theta), 1, dimnames = list(NULL, columns)),
    reformulate(columns), graph,
    conductance_model = loglinear_conductance, covariance = TRUE
  )$covariance[, , 1]
  phi <- unlist(record$phi)
  phi["tau"] * E + exp(phi["sigma"]) * diag(length(graph$demes))
}
