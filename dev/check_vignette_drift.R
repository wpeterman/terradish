# Check precomputed vignettes for prose and code drift from their sources.
#
# Each guide in inst/vignette-source/ is the executable source of a static copy
# in vignettes/ (see dev/refresh-precomputed-vignettes.R). Hand edits must be
# made in both files. This script compares the two without running any model:
# it reports every prose or code line present in one file and absent from the
# other. Run from the package root:
#
#   Rscript dev/check_vignette_drift.R
#
# Not compared: YAML headers, printed output (lines starting with "#>"), the
# precomputed-results notice, and figure captions (a chunk option in the source,
# an image line in the copy). Check captions by eye when you change one.

if (!file.exists("DESCRIPTION") || !dir.exists("vignettes"))
  stop("Run this script from the terradish package root.")

sources <- list.files(file.path("inst", "vignette-source"), pattern = "[.]Rmd$")

# Drop the YAML header and trailing whitespace; return the body lines.
read_body <- function(path) {
  lines <- sub("[ \t]+$", "", readLines(path, warn = FALSE, encoding = "UTF-8"))
  yaml_end <- which(lines == "---")[2L]
  lines[seq_along(lines) > yaml_end]
}

# Remove the knitr setup chunk and every hidden chunk (include = FALSE).
# knitr leaves hidden chunks out of a precomputed copy, so they are not drift.
drop_setup <- function(lines) {
  repeat {
    hidden <- grep("^```\\{r (setup|[^}]*include *= *FALSE)", lines)
    if (!length(hidden)) return(lines)
    start <- hidden[1L]
    end <- start + which(lines[-seq_len(start)] == "```")[1L]
    lines <- lines[-(start:end)]
  }
}

# Escape regular-expression metacharacters in literal text.
escape <- function(x) gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)

# Turn a source line with inline R (`r ...`) into a pattern that matches the
# same line with each expression replaced by its knitted value.
inline_pattern <- function(line) {
  parts <- strsplit(line, "`r [^`]+`")[[1]]
  n_inline <- lengths(regmatches(line, gregexpr("`r [^`]+`", line)))
  parts <- c(parts, rep("", n_inline + 1L - length(parts)))
  paste0("^", paste(escape(parts), collapse = ".*"), "$")
}

drift <- 0L
for (guide in sources) {
  src <- read_body(file.path("inst", "vignette-source", guide))
  copy <- read_body(file.path("vignettes", guide))

  # Source: drop the hidden setup chunk and all chunk fences.
  src <- drop_setup(src)
  src <- src[!grepl("^```", src) & nzchar(src)]

  # Copy: a precomputed copy starts with an inserted notice that ends at a
  # horizontal rule; drop it. A copy that still runs at build time (no
  # notice) keeps its own setup chunk; drop that instead. Then drop printed
  # output, image lines, and fences.
  if (any(grepl("^The displayed results and figures in this guide", copy))) {
    notice_end <- which(copy == "---")[1L]
    copy <- copy[seq_along(copy) > notice_end]
  } else {
    copy <- drop_setup(copy)
  }
  copy <- copy[!grepl("^(#>|!\\[|```)", copy) & nzchar(copy)]

  has_inline <- grepl("`r [^`]+`", src)
  literal <- src[!has_inline]
  patterns <- vapply(src[has_inline], inline_pattern, character(1))

  # Source lines missing from the copy.
  missing_in_copy <- setdiff(literal, copy)
  unmatched_pattern <- names(patterns)[!vapply(patterns, function(p)
    any(grepl(p, copy)), logical(1))]
  # Copy lines that match neither a literal source line nor an inline pattern.
  extra <- setdiff(copy, literal)
  if (length(extra) && length(patterns))
    extra <- extra[!vapply(extra, function(line)
      any(vapply(patterns, grepl, logical(1), x = line)), logical(1))]

  problems <- c(
    if (length(missing_in_copy)) paste("  source only:", missing_in_copy),
    if (length(unmatched_pattern)) paste("  inline line unmatched:", unmatched_pattern),
    if (length(extra)) paste("  copy only:  ", extra))
  if (length(problems)) {
    drift <- drift + length(problems)
    cat("\n", guide, "\n", paste(problems, collapse = "\n"), "\n", sep = "")
  } else {
    cat(guide, ": no drift\n", sep = "")
  }
}

if (drift > 0L) stop(drift, " drifted line(s) found.", call. = FALSE)
cat("\nAll precomputed vignettes match their sources.\n")
