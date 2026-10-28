pkgload::load_all('.', quiet = TRUE, export_all = FALSE, compile = FALSE)
knitr::purl('inst/vignette-source/categorical-conductance.Rmd',
             output = 'dev/pond_vignette_20261013/categorical_benchmark.R', quiet = TRUE)
lines <- readLines('dev/pond_vignette_20261013/categorical_benchmark.R')
last <- grep('quick-reference', lines, fixed = TRUE)
if (length(last)) lines <- lines[seq_len(last[1] - 1L)]
tryCatch(source(textConnection(lines)), error = function(e) {
  if (exists('cv')) print(cv$results)
  stop(conditionMessage(e))
})
write.csv(cv$results, 'dev/pond_vignette_20261013/categorical_benchmark_cv.csv', row.names = FALSE)
writeLines('complete', 'dev/pond_vignette_20261013/categorical_benchmark.complete')
