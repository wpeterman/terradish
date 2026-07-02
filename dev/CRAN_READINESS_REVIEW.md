# terradish CRAN-readiness and performance Red Team review

**Date:** 2026-07-02. **Package version:** 0.0.44. **Reviewer:** multi-agent
Red Team (CRAN compliance, code safety, correctness, performance) plus a live
`R CMD check --as-cran` on R 4.6.0 (Windows). This document is the durable
action list. Items are grouped by severity and each carries a concrete fix.

## Ground-truth check result

`R CMD check --as-cran --no-vignettes` (built from a clean `git archive` of HEAD;
vignettes skipped because they need the non-CRAN `landgraph` and precomputed
`.rds`): **Status: 1 ERROR, 4 WARNINGs, 5 NOTEs**.

**The single most useful finding: the test suite passes.** `checking tests ...
Running 'testthat.R' [10m] OK`. The full suite that appears to crash under
`devtools::test()` is a Windows `devtools`/`load_all` DLL-reload artifact, not a
package defect. Under a fresh R process with the installed package, all tests
pass. Validate with `R CMD check`, not `devtools::test()`, on this machine.

**Separate real defects from local-environment noise.** Several of the reported
items come from this machine's LaTeX/pandoc install, not from the package:

| Reported | Real defect? | Note |
|---|---|---|
| ERROR: PDF manual "without index" | No (environment) | `xkeyval Error: 'style' undefined in families 'Gin'` is a MiKTeX graphics/xkeyval mismatch. Confirm on a clean TeX. |
| WARNING: PDF manual LaTeX errors | No (environment) | Same MiKTeX issue. |
| NOTE: README/NEWS "cannot be checked without pandoc" | No (environment) | pandoc not installed here. |
| NOTE: `terradish-manual.tex` in check dir | No (environment) | Detritus from the failed PDF build. |
| NOTE: `lastMiKTeXException` in temp | No (environment) | MiKTeX artifact. |
| WARNING: vignettes present but no `inst/doc` | Partly artifact | Caused by `--no-vignettes`; a full build populates `inst/doc`. The `vignettes/*.rds` being listed IS real (see R2). |
| WARNING: incoming feasibility (landgraph) | **Yes** | The blocker (B1). |
| NOTE: undefined globals `as`, `formula`, `optim` | **Yes** | Trivial fix (M3). |
| NOTE: examples > 5s | **Yes** | Wrap slow examples (N3). |

Before submission, run an authoritative multi-platform check on a clean toolchain:
`devtools::check_win_devel()` and `rhub::rhub_check()`. The PDF-manual ERROR is
almost certainly local, but confirm it clears on CRAN's TeX.

---

## Part 1: CRAN readiness

### BLOCKERS (submission is impossible until resolved)

**B1. `landgraph` is not on CRAN.** Confirmed by the check's incoming-feasibility
WARNING: *"Strong dependencies not in the CRAN or BioC software repositories:
landgraph."* It is in `Imports` (DESCRIPTION) and re-exported from
`R/reexports.R` (7 functions: `cov_from_biallelic`, `cov_from_genetic_data`,
`fst_from_biallelic`, `dist_from_cov`, `dist_from_biallelic`, `edge_gradient`,
`edge_flow`). CRAN requires every `Imports`/`Depends` to be available on CRAN or
a listed repository. `Additional_repositories` covers Suggests, not hard
dependencies, so it does not rescue this. **Options:** (a) submit `landgraph` to
CRAN first (and any of its own non-CRAN dependencies); or (b) vendor the 7 helper
functions into terradish and drop the dependency and the re-exports. Nothing else
can be submitted until one of these is done. (`multiScaleR`, by contrast, IS on
CRAN as of v0.7.0, so it is not a blocker.)

**B2. `set.seed()` in exported functions without restoring the RNG (CRAN policy).**
This is a policy violation (`R CMD check` does not always flag it, but reviewers
reject it): package code must not change the user's global random state.
- `R/radish_cv.R:584` (`terradish_cv`) calls `set.seed(seed)` unconditionally
  (the seed is auto-generated via `sample.int` when `NULL`). Verified.
- `R/genetic_distance.R:113` (`simulate_covariance_response`) calls `set.seed`
  when a seed is supplied, no restore. Verified.
- `R/radish_optimize.R:1445` (`simulate.radish`) also reseeds; lower risk since
  it mirrors `stats:::simulate.lm`, but that idiom restores `.Random.seed`.
**Fix:** use the save/restore pattern already in this codebase at
`R/covariance_response_power.R:232-241` (save `.Random.seed`, restore with
`on.exit`), or a local RNG stream.

### MUST-FIX (WARNING-level or clear policy)

**M1. `LICENSE` file is malformed.** `License: BSD_3_clause + file LICENSE`
requires exactly two lines. The file has three (`YEAR: 2020`, `COPYRIGHT HOLDER:
Nathaniel S. Pope`, `ORGANIZATION: Nathaniel S. Pope`); the `ORGANIZATION` line is
non-standard, the year is stale, and the holder is Pope while the maintainer is
Peterman. **Fix:** reduce to two lines with the correct current holder, e.g.
`YEAR: 2026` / `COPYRIGHT HOLDER: Bill Peterman` (confirm the intended holder).

**M2. `data/melip.RData` is bzip2-compressed, not xz.** Verified (magic `BZh9`).
CRAN wants best compression and NOTEs otherwise. **Fix:**
`tools::resaveRdaFiles("data", compress = "xz")`; add `LazyData: true` and
`LazyDataCompression: xz` to DESCRIPTION (or keep LazyData off but still xz).

**M3. Undeclared imported functions (NOTE, but trivially real).** The check
reports undefined globals: `as` (methods), `formula` and `optim` (stats), used in
`.hierarchical_fit_fixed`, `terradish_directed`, `conductance_field`, and the
`kron_*` reducers. **Fix:** add roxygen `@importFrom methods as` and
`@importFrom stats formula optim` (regenerate NAMESPACE). Also add `grDevices` to
DESCRIPTION `Imports` (NAMESPACE already has `importFrom(grDevices, terrain.colors)`
but DESCRIPTION does not list it).

**M4. `.Rbuildignore` gaps ship cruft in the tarball.** Not matched by current
patterns: `.Rhistory`, `.agents/`, `.codex/`, `figure/`, `cleanup`, and the
`vignettes/vignette-*.rds` precompute files (only `_cache` is ignored). **Fix:**
add `^\.Rhistory$`, `^\.agents$`, `^\.codex$`, `^figure$`, `^cleanup$`,
`^vignettes/.*\.rds$` (keep `configure`/`configure.ac`, which are legitimate).

**M5. Committed build artifacts in `src/`.** `src/*.o` and `src/terradish.dll`
are tracked and will be packaged (non-portable, NOTE). **Fix:** add them to
`.gitignore` and `.Rbuildignore` (`^src/.*\.(o|so|dll)$`) and `git rm --cached`.

### SHOULD-FIX (NOTE-level, clear before submission)

**N1. Slow examples.** `checking examples [149s] NOTE`: `terradish` 43.6s,
`mlpe_covariates` 16.7s, `terradish_grid` 11.9s, `conductance` 9.3s,
`terradish_multiscale` 6.0s, `aic_table` 5.2s, `cv_model_selection` 5.3s. CRAN
wants each example to run in a few seconds. **Fix:** shrink the example data or
wrap the slow bodies in `\donttest{}` (still run under `--as-cran`, not in
routine checks).

**N2. Missing `\value` and missing runnable `\examples` on some exports.**
`ArmijoControl` has no `\value`. Exports with `\value` but no `\examples`:
`conductance_field`, `directed_rates`, `terradish_directed`,
`terradish_directed_algorithm`, `terradish_kron_reduce`,
`terradish_kron_reduce_tiled`. **Fix:** add `@return` to `ArmijoControl` and a
small runnable (or `\donttest`) `@examples` to each.

**N3. `\dontrun{}` used where `\donttest{}` is likely correct.** 12 Rd files
(`covariance_response_power`, `gaussian_scale_summary`,
`gaussian_smoothed_loglinear_conductance`, `mlpe_response_change`,
`NewtonRaphsonControl`, `plot.terradish`, `terradish_assess_settings`,
`terradish_hierarchical`, `terradish_parameters`, `terradish_results`,
`terradish_scale_optim`, `wishart_covariates`). `\dontrun` is only for
genuinely-unrunnable code; slow-but-valid examples belong in `\donttest`.
Reviewers routinely challenge unjustified `\dontrun`.

**N4. Submission artifacts.** No `cran-comments.md` (needed at submission).
`NEWS.md` top entry says `0.0.44 (dev)`; drop the `(dev)` suffix and set a real
release version (0.0.44 is unusual for a first CRAN submission but is allowed).
No `inst/WORDLIST`; add one after `spelling::spell_check_package()` so domain
terms (Laplacian, MLPE, Wishart, conductance, terra) do not NOTE.

**N5. Installed size / vendored amgcl.** Installed size 6.1 MB, of which
`include/` is 1.5 MB: the full vendored `inst/include/amgcl/` header library
(MPI/CUDA/VexCL backends included). The `inst/THIRD_PARTY_NOTICES/amgcl-LICENSE.md`
is present (good). **Fix:** ship only the amgcl headers `src/` actually compiles
against, to trim the tarball and reduce reviewer scrutiny.

**N6. `RandomFields` reference is dead (Suggests hygiene).**
`R/package_access.R:37-45` references `RandomFields`, which is archived on CRAN
and not in Suggests. It is `requireNamespace`-guarded so it will not error, but
it is dead code. **Fix:** remove those helpers or move them to `dev/`.

**Minor.** `1:length()`/`1:nrow()` in a few spots (`mlpe.R:136`,
`radish_optimize.R:1156,1458`); prefer `seq_len`/`seq_along`. DESCRIPTION could
add the McCullagh (2009) generalized-Wishart DOI to the Description field and
ORCIDs to `Authors@R` (both optional, Consensus-verify the DOI first).

---

## Part 2: Correctness findings

The directed engine was checked hard and is **correct**: the forward hitting-time
system, the reverse-mode transpose-solve adjoint (signs, transpose, chain rule,
absorber indexing, R/C++ agreement), the gamma=0 reduction to a reversible
generator, the drift intercept-only reduction to `wishart_covariance`, and the
drift curvature term all verified by hand and against the tests. The real issues
are at the hierarchical model's boundary behavior.

**C1. Hierarchical model breaks silently on the `tau = 0` boundary (HIGH,
verified).** `terradish_algorithm` multiplies both its gradient and its full
`(theta, u)` Hessian by `(1 - subproblem$boundary)` (`R/radish_algorithm.R:900-901`).
`boundary` is TRUE whenever the profiled measurement model lands on `tau = 0`, a
legitimate and common outcome (weak or no resistance signal, IBD-like data). The
hierarchical code never inspects `boundary`, so on that boundary:
- `hres$hessian` is all zeros, so `H_pen` has a zero theta-block and
  `vcov_full <- .safe_invert(...)` (`R/hierarchical_conductance.R:327`) inverts a
  singular matrix, silently falling back to `ginv`. `theta_se` becomes
  meaningless with no warning.
- `edf_field` collapses to 0 (`hierarchical_conductance.R:343-345`,
  `solve(H_pen[uidx,uidx], Huu_lik)` with `Huu_lik = 0`), so `df` undercounts and
  **AIC is wrong**.
- During fitting, the field gradient wrt `u` is zeroed, leaving only the penalty,
  which drives `u -> 0`; a single step onto `tau = 0` **erases the fitted field**.
- `logml` is distorted near the boundary, corrupting tau^2 selection.
This is exactly the regime the Scenario 3 SLiM validation kept hitting
(`tau -> 0`). The existing tests use strong-IBR data and never exercise it.
**Fix:** detect `boundary` in the hierarchical path and either warn and refuse to
report SE/AIC, or fit the field with `nonnegative = FALSE` (allow the
measurement model off the `tau = 0` boundary) for the field-estimation step, and
compute `edf`/`vcov` from the un-zeroed likelihood Hessian. Add a regression test
on weak-signal (near-`tau=0`) data.

**C2. tau^2 REML selection can report a grid-boundary maximum as interior (MED).**
`hierarchical_conductance.R:312` takes `which.max(logml)` over `tau2_grid`
(default `10^seq(-2, 2, 9)`) with no check that the argmax is interior. A monotone
profile pins tau^2 to the grid edge and reports it silently. **Fix:** warn when
the argmax is at the first or last grid point and suggest widening `tau2_grid`.

**C3. Laplace marginal likelihood mixes dimensions (LOW, intentional but worth a
comment).** `hierarchical_conductance.R:127-138` uses the full `(p+m)` `H_pen`
log-determinant while the prior normalizer is `m`-dimensional. The test suite
asserts this exact behavior, so it is a deliberate "joint evidence over (theta,u)"
choice; tau^2 *ranking* is essentially unaffected. Document that `logML` is not a
clean `m`-dimensional evidence so it is not misread across models with different
`p`.

**Minor.** `terradish_directed` with `nu = NULL` on a Wishart model wastes the
whole optimization on sentinel values before erroring cleanly; add a fast up-front
`nu` check.

---

## Part 3: Efficiency and larger landscapes (ranked)

Core pattern: each objective/gradient/Hessian evaluation forms the reduced
Laplacian `Qn` ((N-1)x(N-1), N = #cells), Cholesky-factorizes it, solves
`G = Qn^-1 Zn` for the k focal columns, builds the kxk resistance covariance `E`,
fits the 2-parameter nuisance subproblem, then computes the Hessian with p
additional multi-RHS solves (one per conductance parameter). Ranked by
scalability impact:

1. **AMG crossover is set so high it never triggers (HIGH).**
   `R/radish_algorithm.R:139-193`: `auto_direct_max_vertices = 750000`,
   `auto_amg_min_vertices = 1500000`. So "auto" uses direct sparse Cholesky up to
   1.5M cells, whose factor is O(N log N) nonzeros and O(N^{3/2}) flops, the
   memory wall. AMG-PCG (`src/amg_solver.cpp`) is O(N) memory and is already wired
   into both forward and gradient solves and warm-started across iterations.
   **Fix:** lower the direct->AMG crossover to roughly 5e5-1e6 cells. Also the
   guard at `:176` reverts to direct when `n_rhs > 64`, which defeats AMG for
   `k > 64` focal sites; allow block/per-column AMG for large `k`. **This is the
   primary enabler for 1e6+ cells.**

2. **Supernodal symbolic factor is not cached; re-analyzed per fit (HIGH).**
   The graph stores only a `simplicial_ldl` Cholesky template
   (`R/radish_graph.R:266`), and `update(template, Qn)` (numeric-only refactor,
   good) is used only for that mode. At scale the code selects `supernodal_ll`
   (>= 50000 vertices), for which no template exists, so a full supernodal
   symbolic analysis + AMD ordering (the single most expensive step) is redone
   whenever the reuse chain breaks. **Fix:** precompute and cache the supernodal
   symbolic factor at graph construction (alongside the simplicial one) and add
   `perm = TRUE` (AMD) to both templates. Expected 30-60% off each large-graph
   evaluation's setup.

3. **Newton recomputes the exact Hessian every step; use BFGS at scale (HIGH).**
   `R/radish_optimize.R:133-135`: with `optimizer = "newton"` (default when
   p <= 3), every step solves p extra RHS blocks for the Hessian, on top of
   objective+gradient, plus a fresh factorization per line-search trial. The final
   fit already re-evaluates the exact Hessian once (`:1113`), so per-step Hessians
   are pure convergence-acceleration overhead. **Fix:** make the optimizer choice
   size-aware (prefer BFGS when N is large regardless of p); removes p multi-RHS
   solves per step.

4. **Dense k-pipeline: redundant `tG` copy and O(k^3) work per subproblem
   iteration (MED).** `radish_algorithm.R:746-749` forms `tG <- t(G)` (an explicit
   (N-1)xk transpose copy) and `generalized_wishart.R:178` does a full
   `eigen(SigInvW)` every subproblem Newton iteration. **Fix:** avoid the `tG`
   copy (the C++ helpers accept `G`); factor `Sigma = tau E + e^sigma I` once per
   `(tau,sigma)` and replace the eigendecomposition-based generalized inverse of
   the mean-projected operator with a rank-one (Sherman-Morrison) update, since
   the projector is `I - 11'Sigma^-1/(1'Sigma^-1 1)`.

5. **Parallel Hessian re-spawns a PSOCK cluster and serializes the solver state
   every evaluation (MED).** `radish_algorithm.R:591-609`: `makeCluster` +
   `clusterEvalQ(library(terradish))` + serialize `state` (embedding G, tG, Zn,
   the factor) to each worker, then `stopCluster`, inside each Hessian eval.
   Cluster startup and O(N*k) serialization can cost more than the solve.
   **Fix:** create the cluster once per fit and reuse it, or push the
   per-parameter backpropagation (already C++: `backpropagate_laplacian_to_conductance`)
   into OpenMP and drop the IPC entirely.

6. **Constant structure rebuilt per iteration (MED).** The design matrix `X`, the
   RHS `Zn`, and integer edge-index coercions are rebuilt each evaluation though
   they do not depend on `theta`. **Fix:** cache them on the graph/model; recompute
   only `exp(X theta)`.

7. **Triple dense (N-1)x(p*k) copies at the R/C++ boundary; all double; 8-byte
   indices (MED).** `cholmod_direct.cpp:112-140` and `amg_solver.cpp:210-223` copy
   RHS in and solution out; `.terradish_solver_solve` forces `as.matrix(rhs)`.
   **Fix:** solve in place where CHOLMOD allows; use 32-bit indices for N < 2^31;
   consider single precision on the AMG/PCG path (resistance distances rarely need
   15 digits). Cuts peak solve memory 2-3x.

8. **Laplacian rebuilt from triplets each evaluation (LOW).** `radish.cpp:122-174`
   reassembles from triplets (O(nnz log nnz)) though the sparsity pattern is
   constant; the in-place `Q@x[]` update path exists but is never taken because
   `edge_pairs` is always present. **Fix:** precompute the CSC pattern once and
   recompute only the value vector via an edge->nnz map.

9. **Kron/landmark reduction is not wired into the solver; directed path
   factorizes once per absorber (LOW now, HIGH later).** `kron_reduction.R` is a
   standalone primitive (its Schur-complement adjoint is unimplemented), so it
   does not speed `terradish()` today; the "landmark" path that runs reduces the
   focal set k, not the graph N. `directed_sparse_lu.cpp:107-134` builds a
   separate SparseLU factorization per absorber (k full unsymmetric factorizations
   per evaluation). **Fix (longer term):** implement the Schur-complement adjoint
   so the exact Kron reduction becomes a real backend (order-of-magnitude win for
   large-N/moderate-k); for the directed path, share one base factorization across
   absorbers via low-rank updates or batched multi-RHS.

**The two changes that most directly unblock 1e6-cell landscapes are #1 (lower the
AMG crossover so the O(N)-memory solver is actually used) and #2 (cache the
supernodal symbolic factor).** #3 (BFGS at scale) then removes the per-step
Hessian multiplier. #9a (Schur/Kron adjoint) is the highest-leverage longer-term
algorithmic win.

---

## Recommended sequence

1. **Decide B1 (landgraph):** submit landgraph to CRAN, or vendor its 7 helpers
   and drop the dependency. Nothing submits until this is settled. (User decision.)
2. **Quick compliance sweep (low risk, mechanical):** B2 (set.seed guards), M1
   (LICENSE), M2 (xz recompress), M3 (imports), M4 (.Rbuildignore), M5 (src
   artifacts). Regenerate NAMESPACE and man/.
3. **Docs sweep:** N1 (donttest slow examples), N2 (\value + examples), N3
   (dontrun->donttest), N4 (cran-comments, WORDLIST, NEWS version).
4. **Fix C1 (hierarchical tau=0 boundary)** and add the weak-signal regression
   test; add the C2 grid-boundary warning.
5. **Re-check** on a clean toolchain: `devtools::check_win_devel()` +
   `rhub::rhub_check()` (confirms the PDF/pandoc items are environment-only and
   catches platform-specific issues).
6. **Efficiency:** land #1 and #2 first (biggest large-landscape wins), then #3;
   schedule #9a (Kron/Schur adjoint) as a dedicated effort.

---

## Resolution log (branch `phase-c-cran-readiness`, 2026-07-02)

All requested fixes were applied and verified with `R CMD check --as-cran`
(built from a clean `git archive`; `--no-vignettes` because vignette execution
needs `landgraph` and the precomputed `.rds`). Commits: `11b312a`, `ec7dd28`,
plus the example-wrapping follow-up.

**Fixed and verified:**
- **B2** `set.seed()` guarded (save/restore `.Random.seed` on exit) in
  `terradish_cv`, `simulate_covariance_response`, `simulate.radish`.
- **M1** `LICENSE` corrected. Note: the BSD_3_clause template *requires* the
  `ORGANIZATION` field; removing it created a "License stub records with
  missing fields" NOTE, so it was restored with current holders. The
  DESCRIPTION meta-information check is now OK.
- **M2** `melip.RData` xz-compressed (310 KB -> 274 KB).
- **M3** declared `grDevices`, `methods::as`, `stats::formula`, `stats::optim`;
  the "no visible global function" NOTE is gone.
- **M4** `.Rbuildignore` extended (`.Rhistory`, `.agents`, `.codex`, `figure`,
  `vignettes/*.rds`). **M5** `src/*.dll` added to `src/.gitignore` (artifacts
  were already build-ignored).
- **C1** `terradish_hierarchical` now detects the `tau = 0` boundary and returns
  NA inference with a warning instead of inverting a singular Hessian; the
  degenerate fit is flagged in `print`. **C2** warns on `tau2` grid-edge maxima.
  **C3** documents that `logml` is a joint `(theta, u)` evidence. Fast `nu`
  check added to `terradish_directed`.
- **Docs** `@return` added to `ArmijoControl`; `cran-comments.md` and
  `inst/WORDLIST` created; NEWS updated; the nine examples exceeding ~5s wrapped
  in `\donttest` (the example-timing NOTE is cleared; `--run-donttest` passes).

**Final check status: 1 ERROR, 4 WARNINGs, 3 NOTEs, all non-package or expected:**
- WARNING (incoming feasibility): "New submission" + `landgraph` not yet on CRAN
  (submitted; resolves when it lands) + "no prebuilt vignette index" (artifact of
  `--no-vignettes`).
- 2 vignette WARNINGs: artifacts of the `--no-vignettes` build; a full
  vignette-building check is needed once `landgraph` is installable.
- PDF-manual WARNING + ERROR: local MiKTeX `xkeyval`/graphics fault, not an Rd
  problem. The HTML manual builds fine.
- 3 NOTEs: `pandoc` missing (README/NEWS), and MiKTeX detritus
  (`terradish-manual.tex`, `lastMiKTeXException`). All local-toolchain.
- INFO: installed size 6.0 MB (1.5 MB is the vendored `amgcl` headers; optional
  trim, item N5).
- **Tests pass** under `R CMD check` (`testthat.R` OK).

**Not done (deliberate):**
- Blanket `\dontrun` -> `\donttest`: several `\dontrun` blocks reference
  undefined objects (illustrative snippets, e.g. `NewtonRaphsonControl`) and are
  correctly `\dontrun`; converting them would break the check. Left as-is.
- Runnable examples for six exports lacking them (`conductance_field`,
  `directed_rates`, `terradish_directed(_algorithm)`, `terradish_kron_reduce`,
  `terradish_kron_reduce_tiled`): not a check NOTE; add self-contained
  `\donttest` examples when convenient.
- N5 amgcl trim; the efficiency roadmap (Part 3).

**Before submitting:** run `devtools::check_win_devel()` and `rhub::rhub_check()`
on a clean toolchain (with `pandoc` and a working TeX) to confirm the PDF/pandoc
items are environment-only, run a full vignette-building check once `landgraph`
is installable, reconcile `inst/WORDLIST` with `spelling::spell_check_package()`,
and update `cran-comments.md` with the final results.
