# Tests, R CMD check, block-CG, and tiled-Schur

## Unit tests (added)
- `tests/testthat/test-curvature.R` (12 checks, all pass): exact vs
  `gauss_newton` share objective + gradient; GN curvature symmetric everywhere
  and positive semidefinite at the optimum (for `leastsquares` and
  `generalized_wishart`); GN -> exact Hessian at a well-fitting optimum
  (rel < 0.05); `terradish(curvature="gauss_newton")` yields a usable `vcov`;
  invalid `curvature` is rejected.
- `tests/testthat/test-block-cg.R` (8 checks, all pass): `block_cg_reduced_laplacian`
  matches the direct Cholesky solve (<1e-6); warm start from a nearby parameter
  value converges; block-CG converges in fewer iterations than the summed
  single-RHS count; the `solver="block_cg"` dispatch matches the direct solve.

All run against the compiled package in-sandbox. The `curvature` change is
regression-clean: existing suites still pass (test-s3-methods 84/84,
test-wishart-covariance 54/54) once ggplot2 is attached.

## R CMD check
A full `R CMD check` needs a normal install (the GitHub-only `multiScaleR`/
`corMLPE`, vignette build, and a full amgcl recompile), which the sandbox can't
complete. Verified instead, on the actual package code: every edited R file and
both regenerated `.Rd` files parse; the new `curvature` argument is documented
in `man/terradish_algorithm.Rd` and `man/terradish.Rd` (usage + \item) so the
"undocumented arguments" check passes; trailing-newline lint fixed. NOTE: the
working tree on the synced folder is CRLF while HEAD is LF (pre-existing); the
edited files are LF (matching HEAD). Run `devtools::document()` + `R CMD check`
on a normal checkout to confirm a clean check.

## A (continued): true block conjugate gradient
- `src/radish.cpp`: `block_cg_reduced_laplacian()` — O'Leary block-CG, Jacobi
  preconditioned, one shared Krylov space for all focal RHS, pseudo-inverse on
  the small s x s coefficient systems to survive block-CG breakdown as columns
  converge. Wired as `solver = "block_cg"` in `terradish_algorithm()`/`terradish()`.
- Verified: matches the direct solve to ~1e-13; on a 40x40 grid it converged in
  64 shared iterations vs up to 202 per column (1583 summed) for independent CG.
  Like the existing Jacobi-PCG it is for small/well-conditioned graphs; large
  heterogeneous Laplacians remain the AMG solver's regime.

## C: tiled Schur-complement prototype
- `dev/large_landscape/schur_tiled.R` + `test_schur.R`: Schur/Kron elimination of
  non-focal cells, single-shot and tile-by-tile (sequential Schur composes).
- Validated EXACT: focal-site resistance recovered to rel error ~2e-14 for both
  single-shot and 4-tile elimination, with the largest per-tile interior solve
  224 cells vs 899 full. This is the contrast B's results pointed to: node
  lumping (B) biased E by ~25%; Schur elimination is exact and decomposes into
  small, parallel, out-of-core tile solves. (`terradish::terradish_kron_reduce()`
  already exposes Kron reduction; this prototype is the tiled/DD form.)
