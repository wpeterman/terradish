old <- readRDS("dev/release_0.1.0/core_baseline_0.0.47.rds")$models
fixed <- readRDS("dev/release_0.1.0/core_phase2.rds")$models
final <- readRDS("dev/release_0.1.0/core_phase8.rds")$models
rows <- lapply(names(final), function(name) {
  x <- final[[name]]; y <- fixed[[name]]; z <- old[[name]]
  data.frame(model = name, original_logLik = z$logLik, final_logLik = x$logLik,
    delta_from_corrected_logLik = x$logLik - y$logLik,
    max_delta_theta = max(abs(x$theta - y$theta)),
    max_delta_phi = max(abs(x$phi - y$phi)),
    max_delta_se = max(abs(x$se - y$se)), convergence = x$convergence$code)
})
table <- do.call(rbind, rows)
write.csv(table, "dev/release_0.1.0/phase8_baseline_comparison.csv", row.names = FALSE)
print(table, digits = 9)
