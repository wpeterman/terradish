# Extract the strict lower triangle of a matrix

Returns the entries strictly below the main diagonal of a square matrix,
scanned column-by-column (standard R column-major order).

## Usage

``` r
lower(x)
```

## Arguments

- x:

  A square numeric matrix.

## Value

A numeric vector of length \\n(n-1)/2\\ containing the strict
lower-triangular entries of `x`, where \\n\\ is the dimension of `x`.

## Details

This is a thin wrapper around `x[lower.tri(x)]`. It is provided as a
convenience for extracting vectorized pairwise measurements from
symmetric matrices such as distance or covariance matrices. The order
matches the rows of
[`pairwise_covariates`](https://wpeterman.github.io/terradish/reference/pairwise_covariates.md)
and of R's [`dist()`](https://rdrr.io/r/stats/dist.html).

## See also

[`lower.tri`](https://rdrr.io/r/base/lower.tri.html)

Other genetic responses:
[`check_distance_response()`](https://wpeterman.github.io/terradish/reference/check_distance_response.md),
[`pca_dist()`](https://wpeterman.github.io/terradish/reference/pca_dist.md)

## Examples

``` r
m <- matrix(1:9, nrow = 3)
lower(m)  # entries (2,1), (3,1), (3,2)
#> [1] 2 3 6

# Useful for extracting pairwise distances:
D <- matrix(c(0, 0.1, 0.3,
              0.1, 0, 0.2,
              0.3, 0.2, 0), nrow = 3)
lower(D)
#> [1] 0.1 0.3 0.2
```
