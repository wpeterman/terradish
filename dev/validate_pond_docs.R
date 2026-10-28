devtools::document(quiet = TRUE)
devtools::test(filter = 'ohio-ponds|amphibian-cv|nuisance-variance',
               stop_on_failure = TRUE)
writeLines('complete', 'dev/pond_vignette_20261013/focused_tests.complete')
