# Control settings for Armijo backtracking line search

Tuning parameters for an objective-only Armijo line search. This line
search is currently intended for `optimizer = "bfgs"` fits where trial
objective evaluations are cheaper than repeated full
objective-plus-gradient evaluations during Hager-Zhang line search.

## Usage

``` r
ArmijoControl(
  initial = 1,
  contraction = 0.5,
  sufficient_decrease = 1e-04,
  maxit = 25L,
  min_alpha = sqrt(.Machine$double.eps),
  verbose = FALSE
)
```

## Arguments

- initial:

  Initial trial step length.

- contraction:

  Multiplicative shrinkage applied after rejected trial steps. Must be
  strictly between 0 and 1.

- sufficient_decrease:

  Armijo sufficient-decrease constant.

- maxit:

  Maximum number of backtracking steps.

- min_alpha:

  Smallest step length to try before stopping.

- verbose:

  Should line-search progress be printed?

## Value

A list of class `"terradish_armijo_control"` holding the validated
Armijo line-search settings (`type`, `initial`, `contraction`,
`sufficient_decrease`, `maxit`, `min_alpha`, and `verbose`), suitable
for the `ls.control` element of a
[`NewtonRaphsonControl`](https://wpeterman.github.io/terradish/reference/NewtonRaphsonControl.md)
object.

## Details

`ArmijoControl()` is opt-in. The default optimizer control continues to
use
[`HagerZhangControl`](https://wpeterman.github.io/terradish/reference/HagerZhangControl.md).
Armijo trial steps only evaluate the objective, so rejected steps avoid
the gradient/backpropagation work required by Hager-Zhang. After a step
is accepted, BFGS still evaluates the full gradient once at the accepted
point so the quasi-Newton update remains well-defined.

## See also

[`NewtonRaphsonControl`](https://wpeterman.github.io/terradish/reference/NewtonRaphsonControl.md),
[`HagerZhangControl`](https://wpeterman.github.io/terradish/reference/HagerZhangControl.md)

Other optimizer controls:
[`HagerZhangControl()`](https://wpeterman.github.io/terradish/reference/HagerZhangControl.md),
[`NewtonRaphsonControl()`](https://wpeterman.github.io/terradish/reference/NewtonRaphsonControl.md)

## Examples

``` r
ctrl <- NewtonRaphsonControl(
  maxit = 25,
  verbose = FALSE,
  ls.control = ArmijoControl()
)
str(ctrl$ls.control)
#> List of 7
#>  $ type               : chr "armijo"
#>  $ initial            : num 1
#>  $ contraction        : num 0.5
#>  $ sufficient_decrease: num 1e-04
#>  $ maxit              : int 25
#>  $ min_alpha          : num 1.49e-08
#>  $ verbose            : logi FALSE
#>  - attr(*, "class")= chr [1:2] "terradish_armijo_control" "list"
```
