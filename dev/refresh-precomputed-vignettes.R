# Regenerate computational guides from their executable source files.
# Run from the package root after installing the package being documented.
# The displayed code, output, and figures stay in the CRAN vignettes, while
# the fully executable R Markdown sources ship in inst/vignette-source.

guides <- c(
  "categorical-conductance",
  "gaussian-scale-optimization",
  "getting-started",
  "model-comparison",
  "ibe-ibr-workflow",
  "mlpe-pairwise-distance",
  "simulation-design",
  "spline-conductance"
)
requested <- commandArgs(trailingOnly = TRUE)
if (length(requested)) {
  unknown <- setdiff(requested, guides)
  if (length(unknown)) stop("Unknown guide: ", paste(unknown, collapse = ", "))
  guides <- requested
}

if (!file.exists("DESCRIPTION") || !dir.exists("vignettes"))
  stop("Run this script from the terradish package root.")

dir.create("vignettes/precomputed", recursive = TRUE, showWarnings = FALSE)

for (guide in guides) {
  source_file <- file.path("inst", "vignette-source", paste0(guide, ".Rmd"))
  target_file <- file.path("vignettes", paste0(guide, ".Rmd"))
  if (!file.exists(source_file)) stop("Missing source: ", source_file)

  # Prefix each figure so separate guides cannot overwrite one another.
  knitr::opts_chunk$set(
    dev = "png", dpi = 96, error = FALSE,
    fig.path = paste0("precomputed/", guide, "-")
  )
  knitted <- tempfile(fileext = ".md")
  evaluation <- new.env(parent = globalenv())
  # Each guide must be self-contained; objects from earlier guides cannot
  # supply an accidentally omitted setup step. Stop on any chunk error.
  knitr::knit(source_file, output = knitted, quiet = TRUE,
              envir = evaluation)
  lines <- readLines(knitted, warn = FALSE, encoding = "UTF-8")
  yaml_end <- which(lines == "---")[2L]
  if (is.na(yaml_end)) stop("Could not find vignette metadata in ", source_file)

  notice <- c(
    "",
    "The displayed results and figures in this guide were computed from the",
    "executable R Markdown source shipped with the package. Routine package",
    "checks use these precomputed results to keep the CRAN check time modest.",
    paste0("To rerun every code chunk, copy `system.file(\"vignette-source\",",
           " \"", guide, ".Rmd\", package = \"terradish\")` to a writable",
           " directory and render that copy with `rmarkdown::render()`."),
    "",
    "---",
    ""
  )
  lines <- append(lines, notice, after = yaml_end)
  # Printed R output often has line-end padding; keep the committed sources
  # whitespace-clean without changing the values displayed to readers.
  lines <- sub("[ \t]+$", "", lines)

  figures <- list.files("precomputed",
                        pattern = paste0("^", guide, "-.*[.]png$"),
                        full.names = TRUE)
  for (figure in figures) {
    destination <- file.path("vignettes", "precomputed", basename(figure))
    if (!file.copy(figure, destination, overwrite = TRUE))
      stop("Could not copy figure: ", figure)
    unlink(figure)
  }
  writeLines(lines, target_file, useBytes = TRUE)
  receipt_dir <- getOption('terradish.vignette_receipts_dir')
  if (!is.null(receipt_dir)) {
    dir.create(receipt_dir, recursive = TRUE, showWarnings = FALSE)
    receipt <- list(guide = guide, source_md5 = unname(tools::md5sum(source_file)),
                    session = capture.output(sessionInfo()))
    for (name in c('cv', 'cv_cov')) {
      score <- get0(name, envir = evaluation, inherits = FALSE)
      if (inherits(score, 'terradish_cv_folds'))
        receipt[[name]] <- list(results = score$results, summary = score$summary)
    }
    profile <- get0('profile', envir = evaluation, inherits = FALSE)
    if (is.list(profile) && !is.null(profile$interval))
      receipt$profile <- profile[c('estimate', 'interval', 'grid')]
    objects <- mget(ls(evaluation), envir = evaluation, inherits = FALSE)
    selected <- Filter(function(x) inherits(x, 'terradish'), objects)
    model_set <- get0('fits', envir = evaluation, inherits = FALSE)
    if (is.list(model_set) && all(vapply(model_set, inherits, logical(1), 'terradish')))
      selected <- c(selected, setNames(model_set, paste0('candidate_', names(model_set))))
    receipt$fits <- lapply(selected, function(x) list(theta = coef(x), phi = x$fit$phi,
      logLik = as.numeric(logLik(x)), AIC = AIC(x), outer = x$convergence$code,
      nuisance = x$fit$subproblem$convergence, nu = x$comparison$nu))
    saveRDS(receipt, file.path(receipt_dir, paste0(guide, '.rds')))
  }
  message(guide, ": ", length(figures), " precomputed figures")
}
