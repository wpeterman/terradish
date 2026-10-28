devtools::document(quiet = TRUE)
results <- devtools::test(stop_on_failure = TRUE)
saveRDS(results, 'dev/pond_vignette_20261013/full_test_results.rds')
writeLines('complete', 'dev/pond_vignette_20261013/full_tests.complete')

# Only numerical records are exported. No fitted closures or raster pointers
# are needed to inspect the independent vignette execution receipts.
receipts <- list.files('dev/pond_vignette_20261013/receipts', '[.]rds$', full.names = TRUE)
rows <- list()
summary <- list()
fits <- list()
for (path in receipts) {
  x <- readRDS(path)
  for (key in c('cv', 'cv_cov')) {
    if (is.null(x[[key]])) next
    rows[[length(rows) + 1L]] <- cbind(guide = x$guide, x[[key]]$results)
    summary[[length(summary) + 1L]] <- cbind(guide = x$guide, x[[key]]$summary)
  }
  for (nm in names(x$fits)) {
    z <- x$fits[[nm]]
    fits[[length(fits) + 1L]] <- data.frame(guide = x$guide, fit = nm,
      outer = z$outer, nuisance = z$nuisance, logLik = z$logLik, AIC = z$AIC,
      nu = if (is.null(z$nu)) NA_real_ else z$nu)
  }
}
scores <- do.call(rbind, rows)
stopifnot(nrow(scores) == 42L, all(scores$error == ''),
          all(scores$nuisance_iterations == 0), all(is.finite(scores$loglik)))
write.csv(scores, 'dev/pond_vignette_20261013/vignette_cv_scores.csv', row.names = FALSE)
write.csv(do.call(rbind, summary), 'dev/pond_vignette_20261013/vignette_cv_summary.csv', row.names = FALSE)
selected <- do.call(rbind, fits)
stopifnot(all(selected$outer == 0), all(selected$nuisance == 0))
write.csv(selected, 'dev/pond_vignette_20261013/vignette_fit_validation.csv', row.names = FALSE)
writeLines('complete', 'dev/pond_vignette_20261013/receipt_validation.complete')
