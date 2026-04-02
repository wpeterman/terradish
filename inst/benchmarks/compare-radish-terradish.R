#!/usr/bin/env Rscript

# Compare original `radish` with `terradish` on 1 and 6 cores using
# microbenchmark. Each package is benchmarked in its own fresh R session to
# avoid namespace and S3-method conflicts.
#
# Usage:
#   source("inst/benchmarks/compare-radish-terradish.R")
# or
#   Rscript inst/benchmarks/compare-radish-terradish.R

CONFIG <- list(
  radish_lib = NULL,
  terradish_lib = NULL,
  times = 5L,
  terradish_cores = c(1L, 6L),
  measurement_model = "mlpe",
  optimizer = "newton",
  maxit = 20L,
  verbose = FALSE
)

run_one_benchmark <- function(package_name,
                              lib_loc = NULL,
                              cores = NULL,
                              times = 5L,
                              measurement_model = "leastsquares",
                              optimizer = "newton",
                              maxit = 20L,
                              verbose = FALSE) {
  if (!requireNamespace("callr", quietly = TRUE)) {
    stop("Package `callr` is required.")
  }
  
  callr::r(
    function(package_name,
             lib_loc,
             cores,
             times,
             measurement_model,
             optimizer,
             maxit,
             verbose) {
      load_pkg_ns <- function(pkg, lib = NULL) {
        if (is.null(lib)) {
          return(loadNamespace(pkg))
        }
        loadNamespace(pkg, lib.loc = lib)
      }
      
      pkg_fun <- function(ns, name) {
        getExportedValue(getNamespaceName(ns), name)
      }
      
      load_melip_data <- function(pkg, lib = NULL) {
        envir <- new.env(parent = emptyenv())
        utils::data("melip", package = pkg, lib.loc = lib, envir = envir)
        as.list(envir)
      }
      
      if (!requireNamespace("microbenchmark", quietly = TRUE)) {
        stop("Package `microbenchmark` is required in child session.")
      }
      
      ns <- load_pkg_ns(package_name, lib_loc)
      
      if (identical(package_name, "radish")) {
        if (!requireNamespace("raster", quietly = TRUE)) {
          stop("Package `raster` is required for original `radish`.")
        }
        
        melip <- load_melip_data("radish", lib_loc)
        
        covariates <- raster::stack(list(
          altitude = raster::scale(melip[["melip.altitude"]]),
          forestcover = raster::scale(melip[["melip.forestcover"]])
        ))
        
        surface <- pkg_fun(ns, "conductance_surface")(
          covariates,
          melip[["melip.coords"]],
          directions = 8
        )
        
        S <- melip[["melip.Fst"]]
        
        control <- pkg_fun(ns, "NewtonRaphsonControl")(
          maxit = maxit,
          verbose = verbose
        )
        
        bench <- microbenchmark::microbenchmark(
          fit = pkg_fun(ns, "terradish")(
            S ~ altitude * forestcover,
            data = surface,
            conductance_model = pkg_fun(ns, "loglinear_conductance"),
            measurement_model = pkg_fun(ns, measurement_model),
            optimizer = optimizer,
            control = control
          ),
          times = times,
          unit = "s"
        )
        
        fit <- pkg_fun(ns, "terradish")(
          S ~ altitude * forestcover,
          data = surface,
          conductance_model = pkg_fun(ns, "loglinear_conductance"),
          measurement_model = pkg_fun(ns, measurement_model),
          optimizer = optimizer,
          control = control
        )
        
        implementation <- "radish"
      } else if (identical(package_name, "terradish")) {
        if (!requireNamespace("terra", quietly = TRUE)) {
          stop("Package `terra` is required for `terradish`.")
        }
        
        melip <- load_melip_data("terradish", lib_loc)
        
        covariates <- c(
          terra::scale(terra::unwrap(melip[["melip.altitude"]])),
          terra::scale(terra::unwrap(melip[["melip.forestcover"]]))
        )
        names(covariates) <- c("altitude", "forestcover")
        
        surface <- pkg_fun(ns, "conductance_surface")(
          covariates,
          terra::unwrap(melip[["melip.coords"]]),
          directions = 8
        )
        
        S <- melip[["melip.Fst"]]
        
        control <- pkg_fun(ns, "NewtonRaphsonControl")(
          maxit = maxit,
          verbose = verbose
        )
        
        bench <- microbenchmark::microbenchmark(
          fit = pkg_fun(ns, "terradish")(
            S ~ altitude * forestcover,
            data = surface,
            conductance_model = pkg_fun(ns, "loglinear_conductance"),
            measurement_model = pkg_fun(ns, measurement_model),
            optimizer = optimizer,
            control = control,
            cores = as.integer(cores)
          ),
          times = times,
          unit = "s"
        )
        
        fit <- pkg_fun(ns, "terradish")(
          S ~ altitude * forestcover,
          data = surface,
          conductance_model = pkg_fun(ns, "loglinear_conductance"),
          measurement_model = pkg_fun(ns, measurement_model),
          optimizer = optimizer,
          control = control,
          cores = as.integer(cores)
        )
        
        implementation <- paste0("terradish_", as.integer(cores), "_core")
      } else {
        stop("Unsupported package: ", package_name)
      }
      
      data.frame(
        implementation = implementation,
        time_sec = bench$time / 1e9,
        logLik = unname(fit$loglik),
        df = fit$df,
        stringsAsFactors = FALSE
      )
    },
    args = list(
      package_name,
      lib_loc,
      cores,
      as.integer(times),
      measurement_model,
      optimizer,
      as.integer(maxit),
      verbose
    )
  )
}

summarize_results <- function(results) {
  blocks <- split(results, results$implementation)
  
  summary <- do.call(
    rbind,
    lapply(blocks, function(x) {
      data.frame(
        implementation = x$implementation[[1]],
        n = nrow(x),
        median_sec = stats::median(x$time_sec),
        mean_sec = mean(x$time_sec),
        min_sec = min(x$time_sec),
        max_sec = max(x$time_sec),
        logLik = x$logLik[[1]],
        df = x$df[[1]],
        stringsAsFactors = FALSE
      )
    })
  )
  
  rownames(summary) <- NULL
  
  if ("radish" %in% summary$implementation) {
    baseline <- summary$median_sec[summary$implementation == "radish"][1]
    summary$speedup_vs_radish <- baseline / summary$median_sec
  } else {
    summary$speedup_vs_radish <- NA_real_
  }
  
  summary
}

benchmark_analysis <- function(config = CONFIG) {
  res_radish <- run_one_benchmark(
    package_name = "radish",
    lib_loc = config$radish_lib,
    times = config$times,
    measurement_model = config$measurement_model,
    optimizer = config$optimizer,
    maxit = config$maxit,
    verbose = config$verbose
  )
  
  res_terra_1 <- run_one_benchmark(
    package_name = "terradish",
    lib_loc = config$terradish_lib,
    cores = config$terradish_cores[[1]],
    times = config$times,
    measurement_model = config$measurement_model,
    optimizer = config$optimizer,
    maxit = config$maxit,
    verbose = config$verbose
  )
  
  res_terra_6 <- run_one_benchmark(
    package_name = "terradish",
    lib_loc = config$terradish_lib,
    cores = config$terradish_cores[[2]],
    times = config$times,
    measurement_model = config$measurement_model,
    optimizer = config$optimizer,
    maxit = config$maxit,
    verbose = config$verbose
  )
  
  raw <- rbind(res_radish, res_terra_1, res_terra_6)
  summary <- summarize_results(raw)
  
  list(raw = raw, summary = summary)
}

main <- function(config = CONFIG) {
  results <- benchmark_analysis(config)
  print(results$raw)
  cat("\nSummary:\n")
  print(results$summary, row.names = FALSE)
  invisible(results)
}

if (sys.nframe() == 0L) {
  main()
}
