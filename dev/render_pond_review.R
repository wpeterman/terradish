review <- 'C:/Users/peterman.73/Downloads/terradish_pond_vignettes_review_20261013'
dir.create(review, recursive = TRUE, showWarnings = FALSE)
guides <- c('v01-getting-started', 'v07-gaussian-scale-optimization', 'v04-model-comparison',
            'v06-spline-conductance', 'v08-ibe-ibr-workflow', 'v05-categorical-conductance')
for (guide in guides) {
  rmarkdown::render(file.path('vignettes', paste0(guide, '.Rmd')),
    output_dir = review, output_file = paste0(guide, '.html'),
    output_options = list(self_contained = TRUE), quiet = TRUE,
    envir = new.env(parent = globalenv()))
  message('HTML: ', guide)
}
receipts <- setNames(lapply(guides, function(guide)
  readRDS(file.path('dev/pond_vignette_20261013/receipts', paste0(guide, '.rds')))), guides)
receipts <- lapply(receipts, function(x) {
  x$fits <- lapply(x$fits, function(z) {
    z$theta <- as.list(z$theta)
    z$phi <- as.list(setNames(as.numeric(z$phi), rownames(z$phi)))
    z
  })
  if (!is.null(x$profile)) x$profile$interval <- as.list(x$profile$interval)
  x
})
jsonlite::write_json(receipts, file.path(review, 'EXECUTION_RECEIPTS.json'),
                     auto_unbox = TRUE, pretty = TRUE, digits = 16, null = 'null')
writeLines('complete', 'dev/pond_vignette_20261013/html.complete')
