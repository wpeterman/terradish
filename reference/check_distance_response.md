# Check that a squared-distance response is admissible

Tests whether a matrix supplied as the response to
[`generalized_wishart`](https://wpeterman.github.io/terradish/reference/generalized_wishart.md)
has squared-Euclidean geometry. The check evaluates the positive
semidefiniteness of the centered Gram matrix implied by the observed
response. It does not assess whether the chosen genetic summary is
biologically appropriate for a particular study.

## Usage

``` r
check_distance_response(S, tol = 1e-07)
```

## Arguments

- S:

  Numeric square matrix. Entries are interpreted as the
  **squared-distance response** \\D\\ that would be supplied to
  `generalized_wishart`, not as ordinary distances \\d\\. If you have
  ordinary distances, check `S^2` and supply that squared matrix to the
  model.

- tol:

  Positive numeric scalar. Relative tolerance used for structural checks
  and for classifying negative eigenvalues. An eigenvalue is treated as
  numerical zero when it is no smaller than
  \\-\mathrm{tol}\\\lambda\_{\max}\\, where \\\lambda\_{\max}\\ is the
  largest eigenvalue of the centered Gram matrix.

## Value

A list with the following elements:

- admissible:

  Logical. `TRUE` when the centered Gram matrix is positive semidefinite
  within `tol`.

- min_eigenvalue:

  Smallest eigenvalue of the centered Gram matrix, on the
  squared-distance scale of `S`.

- relative_min_eigenvalue:

  The smallest eigenvalue divided by the largest eigenvalue. Values
  close to zero from below usually reflect numerical rounding;
  materially negative values indicate incompatibility.

- n_negative:

  Number of eigenvalues below the scaled negative tolerance.

- tolerance:

  The relative tolerance supplied in `tol`.

- eigenvalue_threshold:

  The absolute negative-eigenvalue magnitude used for the decision,
  equal to `tol * max_eigenvalue`.

- max_eigenvalue:

  Largest eigenvalue of the centered Gram matrix.

- n_sites:

  Number of rows and columns in `S`.

## Details

For \\n\\ sites, the function forms

\$\$B = -\frac{1}{2} H D H, \qquad H = I - \frac{1}{n}
\mathbf{1}\mathbf{1}^{\mathsf T}.\$\$

The response is admissible when \\B\\ is positive semidefinite within
the requested tolerance. Symmetry, nonnegative entries, and a zero
diagonal are necessary but are not sufficient. Passing this geometry
check does not create a Wishart sampling model. Use covariance-derived
squared distances with a defensible effective information parameter.
\\F\_{ST}\\ ratio estimators have no Wishart degrees of freedom merely
because their distance matrix passes the check; consider
[`mlpe`](https://wpeterman.github.io/terradish/reference/mlpe.md) for
those responses.

This function checks the matrix that you actually intend to analyze
because scaling, transformations, locus aggregation, and missing-data
handling can change its geometry. It never applies Cailliez, Lingoes,
nearest-Euclidean, or other corrections. Those procedures alter the
observed response and may alter its biological interpretation.

For an independent check with ade4, note that
[`ade4::is.euclid()`](https://adeverse.github.io/ade4/reference/is.euclid.html)
expects ordinary distances and squares them internally. For an already
squared response `S`, the corresponding call is
`ade4::is.euclid(as.dist(sqrt(S)))`, not `ade4::is.euclid(as.dist(S))`.

## See also

[`generalized_wishart`](https://wpeterman.github.io/terradish/reference/generalized_wishart.md),
[`wishart_covariance`](https://wpeterman.github.io/terradish/reference/wishart_covariance.md),
[`mlpe`](https://wpeterman.github.io/terradish/reference/mlpe.md)

Other genetic responses:
[`lower()`](https://wpeterman.github.io/terradish/reference/lower.md),
[`pca_dist()`](https://wpeterman.github.io/terradish/reference/pca_dist.md)

## Examples

``` r
# Squared Euclidean distances from known coordinates are admissible.
xy <- matrix(c(0, 0,
               1, 0,
               0, 1,
               1, 1), ncol = 2, byrow = TRUE)
D <- as.matrix(dist(xy))^2
check_distance_response(D)
#> $admissible
#> [1] TRUE
#> 
#> $min_eigenvalue
#> [1] -3.855283e-16
#> 
#> $relative_min_eigenvalue
#> [1] -3.855283e-16
#> 
#> $n_negative
#> [1] 0
#> 
#> $tolerance
#> [1] 1e-07
#> 
#> $eigenvalue_threshold
#> [1] 1e-07
#> 
#> $max_eigenvalue
#> [1] 1
#> 
#> $n_sites
#> [1] 4
#> 

# Symmetry and a zero diagonal alone are not sufficient.
bad <- matrix(c(0, 1, 1,
                1, 0, 9,
                1, 9, 0), nrow = 3, byrow = TRUE)
check_distance_response(bad)$admissible
#> [1] FALSE
```
