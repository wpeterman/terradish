# Control settings for Hager-Zhang line search

Tuning parameters for the Hager-Zhang line-search implementation.

## Usage

``` r
HagerZhangControl(
  delta = 0.1,
  sigma = 0.9,
  alphamax = Inf,
  rho = 5,
  epsilon = 1e-06,
  gamma = 0.66,
  linesearchmax = 50,
  psi3 = 0.1,
  c = 1,
  verbose = FALSE
)
```

## Arguments

- delta:

  Sufficient decrease constant (Wolfe condition 1). Must satisfy
  `0 < delta < sigma < 1`. Smaller values accept shorter steps more
  readily.

- sigma:

  Curvature condition constant (Wolfe condition 2). This line search
  constant is unrelated to the `sigma` parameters of the measurement and
  conductance models. Must satisfy `delta < sigma < 1`. Values close to
  1 relax the curvature requirement and allow fewer function evaluations
  per step.

- alphamax:

  Upper bound on the trial step length. The default `Inf` lets the
  bracketing phase grow the step without limit.

- rho:

  Multiplicative expansion factor applied during the initial bracketing
  phase. Must be `> 1`.

- epsilon:

  Non-negative tolerance defining an "approximate Wolfe" regime near a
  local minimum: a small positive value (default `1e-6`) allows the
  algorithm to accept a slightly increased objective near convergence.

- gamma:

  Minimum fractional reduction of the bracket required after each secant
  step. Steps that shrink by less than this factor trigger bisection
  instead. Must be in `(0, 1)`.

- linesearchmax:

  Maximum number of line-search function evaluations. Raise this if the
  optimizer reports repeated line-search failures on difficult problems.

- psi3:

  Backtracking multiplier applied when the initial trial step produces a
  non-finite objective. Must be in `(0, 1)`.

- c:

  Initial trial step length passed to the line search. The default `1.0`
  is appropriate for Newton and quasi-Newton methods that produce
  unit-scale search directions.

- verbose:

  Should line-search progress (trial step lengths and objective values)
  be printed to the console?

## Value

A named list of class `"terradish_hager_zhang_control"` holding the line
search settings. The elements are the arguments above, stored unchanged:
`delta`, `sigma`, `alphamax`, `rho`, `epsilon`, `gamma`,
`linesearchmax`, `psi3`, `c`, and `verbose`. Pass the object to the
`ls.control` argument of
[`NewtonRaphsonControl`](https://wpeterman.github.io/terradish/reference/NewtonRaphsonControl.md)
rather than to
[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md)
directly; inspect it with [`str()`](https://rdrr.io/r/utils/str.html).
The list carries no fitted quantities, so nothing in it needs
back-transforming.

## Details

`HagerZhangControl()` is the default line search used by
[`NewtonRaphsonControl`](https://wpeterman.github.io/terradish/reference/NewtonRaphsonControl.md).
It implements the approximate Wolfe conditions from Hager and Zhang
(2006), which combine sufficient-decrease and curvature checks with a
bracketing-plus-secant strategy.

The default parameters work well for most resistance-surface
optimization problems. Common reasons to change them:

- **Line-search failures**: increase `linesearchmax` or loosen `sigma`
  toward 1.

- **Slow convergence with many small steps**: tighten `sigma` (e.g. 0.5;
  it must stay above `delta`) to enforce a stronger curvature condition.

- **Initial non-finite objectives**: increase `psi3` to take larger
  backtracking steps when `c` is too large.

## References

Hager WW, Zhang H. 2006. Algorithm 851: `CG_DESCENT`, a conjugate
gradient method with guaranteed descent. ACM Transactions on
Mathematical Software 32(1):113-137.
[doi:10.1145/1132973.1132979](https://doi.org/10.1145/1132973.1132979) .

## See also

[`NewtonRaphsonControl`](https://wpeterman.github.io/terradish/reference/NewtonRaphsonControl.md),
[`ArmijoControl`](https://wpeterman.github.io/terradish/reference/ArmijoControl.md)

Other optimizer controls:
[`ArmijoControl()`](https://wpeterman.github.io/terradish/reference/ArmijoControl.md),
[`NewtonRaphsonControl()`](https://wpeterman.github.io/terradish/reference/NewtonRaphsonControl.md)

## Examples

``` r
ctrl <- HagerZhangControl(verbose = FALSE, linesearchmax = 10)
str(ctrl)
#> List of 10
#>  $ delta        : num 0.1
#>  $ sigma        : num 0.9
#>  $ alphamax     : num Inf
#>  $ rho          : num 5
#>  $ epsilon      : num 1e-06
#>  $ gamma        : num 0.66
#>  $ linesearchmax: num 10
#>  $ psi3         : num 0.1
#>  $ c            : num 1
#>  $ verbose      : logi FALSE
#>  - attr(*, "class")= chr [1:2] "terradish_hager_zhang_control" "list"

# Use inside NewtonRaphsonControl:
opt_ctrl <- NewtonRaphsonControl(
  maxit = 50,
  ls.control = HagerZhangControl(sigma = 0.5, linesearchmax = 30)
)
```
