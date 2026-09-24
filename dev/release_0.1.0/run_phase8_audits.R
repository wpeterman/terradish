root <- getwd()
.libPaths(c(file.path(root, "dev/check/phase2-library"), .libPaths()))
pkgload::load_all(file.path(root, "dev/check/phase7-package"), quiet = TRUE)
batch <- commandArgs(trailingOnly = TRUE)
stopifnot(length(batch) == 1L)
files <- switch(batch,
  basic = c("core/exp1b_wishart_clustered.R", "core/exp7_optim.R", "core/exp9_rho0.R",
    "scale/shift_check.R", "scale/ksr_units.R", "flex/exp2_spline_support.R",
    "compute/e1_graph.R", "compute/e4_misc.R", "compute/e6_melip_crop.R"),
  extensions = c("ext/check1_bruteforce.R", "ext/check2_behaviour.R", "ext/check5_kernel_diag.R"),
  cv = c("cv/e3_converted.R", "cv/e5.R", "cv/e6.R", "cv/e2_modes.R"),
  cv_v2 = c("cv/e6.R", "cv/e2_modes.R"),
  corrections = c("compute/e6_melip_crop.R", "ext/check1_bruteforce.R", "ext/check2_behaviour.R"),
  rho = "core/exp9_rho0.R",
  stop("Unknown batch"))
out <- file.path(root, "dev/check", paste0("phase8-audit-", batch))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
results <- list()
for (file in files) {
  cat("\nAUDIT", file, "\n")
  started <- Sys.time()
  wd <- getwd()
  setwd(out)
  env <- new.env(parent = globalenv())
  result <- tryCatch({
    sys.source(file.path(root, "dev/release_0.1.0/phase8_audits", file), envir = env)
    "completed"
  }, error = function(e) list(error = conditionMessage(e)))
  setwd(wd)
  results[[file]] <- list(started = started, finished = Sys.time(), result = result)
  saveRDS(results, file.path(root, "dev/release_0.1.0", paste0("phase8_", batch, ".rds")))
  print(results[[file]])
}
