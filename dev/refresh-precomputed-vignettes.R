# Regenerate the five computational guides from their executable source files.
# Run from the package root after installing the package being documented.
# The displayed code, output, and figures stay in the CRAN vignettes, while
# the fully executable R Markdown sources ship in inst/vignette-source.

guides <- c(
  "gaussian-scale-optimization",
  "getting-started",
  "ibe-ibr-workflow",
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
    dev = "png", dpi = 96,
    fig.path = paste0("precomputed/", guide, "-")
  )
  knitted <- tempfile(fileext = ".md")
  knitr::knit(source_file, output = knitted, quiet = TRUE)
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
  message(guide, ": ", length(figures), " precomputed figures")
}
