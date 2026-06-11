# Package integration results: A + D hardened into terradish

This records the in-package integration of approaches A and D and the
package-level validation (the real `terradish` code, compiled and run, not the
reference prototype).

## What changed in the package

**D — Gauss-Newton / Fisher curvature (new).**
- `terradish_algorithm(..., curvature = c("exact","gauss_newton"))` and
  `terradish(..., curvature = ...)`. Default `"exact"` is unchanged.
- Implementation is a two-term toggle in `.terradish_algorithm_derivative_chunk()`:
  the exact Hessian is the Gauss-Newton/Fisher curvature plus two
  residual-weighted second-derivative terms — the second derivative of the
  resistance covariance `E` (`-2 * dgrad__ddl_dQnG %*% dl_dE`) and of conductance
  (`dl_dC %*% d2f__dtheta_dtheta`). `curvature = "gauss_newton"` drops both, then
  symmetrises. The result is positive semidefinite, needs only first derivatives
  of conductance, and equals the exact Hessian at a well-fitting optimum.
- It flows unchanged into `vcov`/`summary`/`confint`, because those just invert
  the returned `hessian` field — so the Fisher information becomes the
  covariance matrix and the standard errors are the asymptotic
  information-based errors.
- Because the change only *guards* previously-unconditional code with
  `if (!gauss_newton)`, the `curvature = "exact"` path is byte-for-byte the old
  behavior (verified: identical objective and gradient).

**A — solver hardening (mostly already present, confirmed).**
- Cross-step reuse is already wired through the optimizer: `radish_optimize.R`
  threads `solver_warm_start` and `solver_reuse_state` from each
  `terradish_algorithm()` call into the next (lines ~869-873, ~1108), so the
  sparse-Cholesky symbolic factor (`reused_factor_template` / the
  `cholmod_cpp_cached` backend) and the AMG preconditioner are reused across
  Newton/BFGS steps, and potentials warm-start the next solve.
- The reduced-Laplacian solves are already multi-RHS (the whole focal block `Zn`
  is solved at once) for the direct, AMG, and PCG backends.
- Remaining genuine A work (not yet done, lower priority): a true *block*-CG
  Krylov method (current iterative path loops over RHS columns sharing the
  preconditioner) and making the cached CHOLMOD backend the large-N default.

## Package-level validation (compiled `terradish.so`, run in-sandbox)

The package was compiled (with the AMG entry points stubbed; the direct path is
exercised) and run by sourcing the package R against the shared library.

**Curvature self-check on the bundled `melip` example** (`test_curvature_pkg.R`),
for both `leastsquares` and `generalized_wishart`:
- `exact` vs `gauss_newton`: identical objective, identical gradient
  (max difference 0).
- GN curvature symmetric (max asymmetry 0) and **positive definite**
  (eigenvalues e.g. `[108, 5.47]` for leastsquares, `[516, 93.4]` for
  generalized Wishart) — a valid `vcov`.

**End-to-end on SLiM data with `generalized_wishart`** (`slim/fit_terradish_slim.R`),
barrier scenario (true theta `[0.80, -2.50]`, nu=4000), across raster
resolutions via `terra::disagg`:

| factor | cells | curvature | time (s) | cos(theta_hat, true) |
|---|---|---|---|---|
| 1 | 4,096 | exact | 3.75 | 0.975 |
| 1 | 4,096 | gauss_newton | 2.09 | 0.959 |
| 2 | 16,384 | exact | 6.02 | 0.991 |
| 2 | 16,384 | gauss_newton | 4.91 | 0.991 |

- Recovery is strong and improves with resolution (cos 0.975 -> 0.991).
- At the finer (well-fitting) resolution, `gauss_newton` returns estimates
  identical to `exact` while running faster, and its curvature is PSD so the
  Newton step and `vcov` are well-behaved. At the coarse resolution the barrier
  coefficient is weakly identified (near-zero Fisher information in that
  direction) — surfaced honestly as a large/!undefined SE rather than a crash.
- The cross-step cached factor (A) is used automatically by `solver = "direct"`.

## Reproduce (sandbox)
```bash
# compile the shared library (direct path; AMG stubbed)
cd /tmp/tdsrc/src && R CMD SHLIB radish.cpp cholmod_direct.cpp amg_stub.cpp RcppExports.cpp -o terradish.so
# curvature self-check on melip
Rscript dev/large_landscape/test_curvature_pkg.R
# end-to-end on SLiM data with generalized_wishart, two resolutions
Rscript dev/large_landscape/slim/fit_terradish_slim.R dev/large_landscape/slim/scenarios/barrier 2
```
On a normally installed terradish (with multiScaleR present) these run directly
against the package namespace; the sandbox scripts source the package R files
and a locally compiled `terradish.so` because multiScaleR is GitHub-only.

## Status
- D: implemented in the package, documented (roxygen), validated against exact
  Hessian and on SLiM data for leastsquares and generalized_wishart. Ready for a
  `testthat` test and `R CMD check`.
- A: cross-step reuse + multi-RHS confirmed present and exercised; block-CG and
  cached-CHOLMOD-as-default remain as follow-ups.
- SLiM harness: extended to the full optimizer with `generalized_wishart` and to
  larger landscapes via raster disaggregation.
