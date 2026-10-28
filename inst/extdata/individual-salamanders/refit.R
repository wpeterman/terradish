# Extended reproduction. Copy this file to a writable working directory.
# It fits existing data only, writes local checkpoints, and never runs SLiM.
if (Sys.getenv("TERRADISH_RUN_EXTENDED_EXAMPLE") != "true")
  stop("Set TERRADISH_RUN_EXTENDED_EXAMPLE=true to opt into the serial refits.")
library(terradish)
data(individual_salamanders)
files <- system.file("extdata", "individual-salamanders", package = "terradish",
                     mustWork = TRUE)
terms <- terra::rast(file.path(files, "drying_terms_10m.tif"))
dir.create("refit-results", showWarnings = FALSE)

for (arm in names(individual_salamanders$data)) {
  animals <- individual_salamanders$data[[arm]]
  samples <- animals$samples
  graph <- conductance_surface(terms,
    as.matrix(samples[c("easting", "northing")]), directions = 8)
  for (panel in names(animals$covariance)) {
    for (model in c("uniform", "IBR")) {
      for (fold in 0:4) {
        train <- if (fold == 0) seq_len(nrow(samples)) else which(samples$fold != fold)
        # Subset focal locations, retaining the entire fixed landscape graph.
        training_graph <- graph
        training_graph$demes <- graph$demes[train]
        training_graph$rhs <- graph$rhs[, train, drop = FALSE]
        S <- animals$covariance[[panel]][train, train]
        S <- S - rowMeans(S) - rep(colMeans(S), each = nrow(S)) + mean(S)
        columns <- if (model == "IBR") c("drying1", "drying2") else character()
        width <- unlist(individual_salamanders$metadata$start_halfwidth)
        starts <- if (length(columns)) list(zero = setNames(c(0, 0), columns),
          negative = -width[columns], positive = width[columns]) else list(uniform = NULL)
        for (label in names(starts)) {
          path <- file.path("refit-results",
            paste0(paste(arm, panel, model, fold, label, sep = "__"), ".rds"))
          # Resume from completed attempts, never replace a different start.
          if (file.exists(path)) next
          fit <- terradish(reformulate(columns, response = "S"), training_graph,
            conductance_model = loglinear_conductance,
            measurement_model = wishart_covariance, nu = 1,
            theta = starts[[label]], optimizer = "bfgs", cores = 1,
            leverage = FALSE,
            control = NewtonRaphsonControl(maxit = 100, ctol = 1e-6, ftol = 1e-8),
            measurement_control = NewtonRaphsonControl(
              maxit = 100, ctol = 1e-7, ftol = 1e-9))
          saveRDS(slim_terradish(fit), path)
          cat(arm, panel, model, fold, label, fit$convergence$message, "\n")
          flush.console()
        }
      }
    }
  }
}
