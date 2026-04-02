.terradish_deprecate <- function(old, new)
{
  .Deprecated(new = new, package = "terradish", old = old)
}

.terradish_forward_call <- function(call, new)
{
  call[[1L]] <- as.name(new)
  eval(call, parent.frame(2L))
}

.terradish_set_class <- function(x, ...)
{
  class(x) <- unique(c(..., class(x)))
  x
}
