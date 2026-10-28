# Consume only named fixture diagnostics; unrelated warnings remain test warnings.
expect_fixture_warnings <- function(expr, patterns, required = patterns) {
  seen <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    message <- conditionMessage(w)
    if (any(vapply(patterns, grepl, logical(1), x = message))) {
      seen <<- c(seen, message)
      invokeRestart("muffleWarning")
    }
  })
  for (pattern in required)
    expect_true(any(grepl(pattern, seen)), info = paste("Expected fixture warning:", pattern))
  value
}
