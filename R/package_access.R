.require_suggested_package <- function(pkg, context = NULL)
{
  if (!requireNamespace(pkg, quietly = TRUE))
  {
    if (is.null(context))
      stop("Package '", pkg, "' is required for this functionality.", call. = FALSE)
    stop("Package '", pkg, "' is required for ", context, ".", call. = FALSE)
  }
}

.suggested_export <- function(pkg, fun, context = NULL)
{
  .require_suggested_package(pkg, context = context)
  getExportedValue(pkg, fun)
}

.adegenet_tab <- function(...)
{
  .suggested_export("adegenet", "tab", context = "`pca_dist()`")(...)
}

.numderiv_grad <- function(...)
{
  .suggested_export("numDeriv", "grad", context = "numerical validation")(...)
}

.numderiv_hessian <- function(...)
{
  .suggested_export("numDeriv", "hessian", context = "numerical validation")(...)
}

.numderiv_jacobian <- function(...)
{
  .suggested_export("numDeriv", "jacobian", context = "numerical validation")(...)
}

# Accessors for packages used only by the benchmark harness in
# inst/benchmarks/generalized-wishart-timing.R live in that file, not here.
