# large_landscape — scaling terradish beyond ~1M cells

Design document and runnable prototypes for estimating resistance/conductance
over rasters much larger than the current ~1M-cell limit, while maximally
leveraging the information content of the data.

## Read first
- **`SCALING_DESIGN.md`** — the cost model (three walls: factorisation fill-in,
  dense-potential memory, parameter-linear Hessian), the two senses of "maximal
  information content," and a **ranked menu** of approaches with a suggested
  build order. This is the main deliverable.
- **`RESULTS.md`** — D/A/B validated end-to-end on SLiM-simulated genetic data:
  what works, the numbers, and the one finding that changed the ranking (naive
  coarsening B biases resistance; route exact inference through A or tiled Schur).
- **`slim/`** — the SLiM + tree-sequence validation harness (recipe, landscape
  generator, genotype processor, three scenarios). Toolchain stood up and run.
- **`validate_DAB.R`** — runs A/D/B against a scenario's SLiM data in R.

## Prototypes (research code; math verified, not yet wired into the package)
| File | Approach | Attacks | Verified |
|---|---|---|---|
| `verify.py` | core identities | — | adjoint grad ≈ FD 3e-9; block-PCG ≈ direct 1e-12; Hutchinson ≈ truth ~1% |
| `01_block_solver_adjoint.R` | A: matrix-free block solve + streamed adjoint gradient | W1, W2 | self-test vs dense |
| `02_adaptive_coarsen.R` | B: current-density adaptive multiresolution | N (all walls) | E resolution-convergence on a barrier landscape |
| `03_gauss_newton_fisher.R` | D: Gauss-Newton / Fisher curvature + Hutchinson | W3 + inference | GN ≈ Hessian at optimum |
| `hutchinson_diag.cpp` | D/F: compiled stochastic trace/diagonal | — | mirrors verify.py |

## Run
```bash
python3 verify.py                      # no R needed; proves the math
```
```r
# in R 4.6 with Matrix installed, from this folder:
source("01_block_solver_adjoint.R")    # prints gradient & solver errors
source("02_adaptive_coarsen.R")        # prints vertex reduction & E error
source("03_gauss_newton_fisher.R")     # prints curvature, SEs, Hutchinson trace
# hutchinson_diag.cpp: Rcpp::sourceCpp() inside the package tree, or move to src/
```

## Headline recommendations (see the design doc for the full ranking)
1. **D — Gauss-Newton / Fisher curvature first.** Cheapest; removes the
   parameter-linear Hessian cost, makes Newton steps robust, and *is* the
   information matrix behind `vcov`/SEs.
2. **A — finish the iterative engine.** Block Krylov + cross-step warm starts +
   streamed adjoint + mandatory solver-error diagnostics. Mostly hardening of
   `src/amg_solver.cpp` / `pcg_reduced_laplacian()` you already have.
3. **B — adaptive current-density multiresolution.** The structural prize;
   shrinks `N` itself. Gate on `E`-resolution-convergence before integrating.
4. **C — Schur / tiling** only when focal sites are sparse in extents too large
   to hold even coarsened (out-of-core).

On the statistics side, "maximal information" also means preferring the
`generalized_wishart` / `wishart_covariance` likelihood (full covariance + `nu`)
over MLPE-on-distances, and using the largest defensible number of focal points
(individual-level, not just populations) — which the scalable solver is what
makes affordable.

## Suggested validation
Port `verify.py`'s gradient check to `testthat`; use your local **SLiM** install
for mechanism-based parameter-recovery and CI-coverage tests on large simulated
landscapes (the decisive experiment for any approximation). Full plan in
`SCALING_DESIGN.md` §5.
