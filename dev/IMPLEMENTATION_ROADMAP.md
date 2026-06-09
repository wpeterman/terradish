# terradish population-genetics integration — implementation roadmap

**Purpose.** This is the durable, returnable reference for the multi-tier effort
to integrate more explicit population genetics into terradish, following the
feasibility assessment in `dev/FEEMS_integration_feasibility.md`. It records the
plan, the math, the validation protocol, the SLiM strategy, the code entry
points, and a living progress log. If work is interrupted, **start here.**

Branch: `dev`. Companion source docs:
`Research/Manuscripts/OSU/SpatialResistance/FEEMS_terradish_critical_review.docx`
and `…/FEEMS_terradish_technical_appendix.pdf`.

---

## 0. Working agreement (the per-tier cycle)

For each tier, in order, do not advance until the current one is **green**:

1. **Build** — implement the feature (R + C++ as needed), `devtools::load_all()`
   clean.
2. **Test** — unit tests (testthat) following existing conventions in
   `tests/testthat/`. Include the reduction-to-base-case test (new model must
   reproduce the existing model at the appropriate parameter values).
3. **Validate** — three-rung ladder (see §3):
   (a) analytic derivatives vs `numDeriv`;
   (b) self-consistency parameter recovery (simulate from the model, refit);
   (c) **SLiM** simulation recovery of a known biological truth.
4. **Document** — roxygen on all exported functions, regenerate `man/` with
   `roxygenise()`, add/extend a vignette, and update this roadmap's progress log
   (§6) + the relevant tier's "status".

Constraints carried from prior sessions unless the user says otherwise:
- Do not fabricate citations or interpretations.
- Do not commit or push unless explicitly asked.

---

## 1. Engine constraints (why the plan is shaped this way)

From reading the actual implementation (see feasibility doc §1). These are the
load-bearing facts:

- **F1 — Symmetric, SPD-only solver.** `src/radish.cpp::assemble_reduced_laplacian`
  builds edge weight `w = conductance(i) + conductance(j)` (the **sum**, not a
  mean), symmetric by construction; factored by `Matrix::Cholesky`/CHOLMOD. All
  backends (simplicial/supernodal LL/LDL, AMG, IC/Jacobi PCG) are SPD. No
  non-symmetric solver exists.
- **F2 — Optimizer is for low-dimensional `theta`.** `terradish_algorithm()`
  forms a **dense** Hessian via one batched solve per parameter
  (`idx <- seq_along(theta)`). Fine for length(theta) ≤ ~20; wrong for 10³–10⁵.
- **F3 — Conductance models are closures** returning `conductance = exp(Xθ)` +
  derivative callbacks. Clean extension seam, subject to F2's dimensionality.
- **F4 — Measurement models carry a thin scalar nuisance vector.**
  `wishart_covariance`: `Σ = τE + exp(σ)I`. `wishart_covariates` already extends
  additively: `Σ = τE + Σ_k λ_k z_k z_kᵀ + exp(σ)I`, fit by low-dim Newton in
  `radish_subproblem`. **This is the seam for richer likelihoods (Tiers 1).**

Two kinds of "asymmetry" (decisive for Tier 3):
- **Reversible** (heterogeneous node size/density): generator similar to a
  symmetric matrix via `D^{1/2}GD^{-1/2}`; SPD machinery reusable.
- **Non-reversible** (directional advection: downstream/downhill/downwind): no
  SPD structure; needs new non-symmetric linear algebra + transpose adjoints +
  an asymmetric distance object. This is the high-demand feature and the deepest
  change.

---

## 2. The measurement-model contract (for Tiers 1–2)

A `terradish_measurement_model` is a function
`g(E, S, phi, nu, gradient, hessian, partial, nonnegative, validate)`:

- **`phi` missing** → return `list(phi = <start>, lower = <vec>, upper = <vec>)`.
- **`phi` present** → return a list with: `objective` (neg loglik), `fitted`,
  `boundary` (logical, e.g. `nonnegative && tau == 0`), `gradient` (d/dphi),
  `hessian` (d²/dphi²), `gradient_E` (d/dE), `partial_E` (∂(d/dE)/∂phi),
  `partial_S` (∂(gradient)/∂S, lower-tri layout), and closures `jacobian_E`,
  `jacobian_S` for reverse-mode AD through E and S.

`radish_subproblem()` profiles out `phi` by `BoxConstrainedNewton` using
`lower`/`upper`, then does one `partial = TRUE` call to get the AD closures.
`terradish_algorithm()` consumes `gradient_E`/`partial_*`/`jacobian_*` to
backprop through the Laplacian to `theta`.

**Generic engine already present:** `.wishart_covariate_finish()` (in
`R/wishart_covariates.R`) is generic over a list `B` of per-parameter basis
matrices with `Σ = Σ_nm B[[nm]]`:
`dPhi[nm] = sum(B[[nm]] * grad_Sigma)`, `ddPhi[i,j] = sum(B[[i]] * dgrad[[j]])`,
plus a single curvature term `ddPhi["sigma","sigma"] += dPhi["sigma"]`. Reuse it.

---

## 3. Validation ladder (every tier)

1. **Analytic vs numerical derivatives.** Use the package's `numDeriv` path:
   each measurement model supports `validate = TRUE`; for new models add a test
   that compares `gradient`/`hessian` to `numDeriv::grad`/`hessian` of
   `objective`, and `partial_S`/`jacobian_*` where applicable. Tolerance ~1e-5.
2. **Self-consistency recovery.** Use `simulate_covariance_response()` /
   `covariance_response_power()` to draw data from the *new* generative model
   with known parameters, refit, confirm recovery and calibrated intervals.
3. **SLiM recovery of biological truth.** Simulate a spatially explicit
   population with the relevant process (Tier 1: spatially varying local density
   / N_e; Tier 2: an unmodeled smooth conductance feature; Tier 3: directional
   gene flow), summarize to the genetic response terradish consumes (covariance
   via `cov_from_genetic_data` / distance via Fst), fit, and confirm the new
   machinery recovers the planted signal and improves on the base model.

**SLiM conventions (fill in as built):** scripts under `dev/slim/`; each script
parameterized and seeded; an R harness (`dev/slim/run_*.R`) that calls SLiM,
reads tree-sequence/VCF output, builds the terradish response, and asserts
recovery. Record SLiM version + command lines in §6.

---

## 4. The tiers

### Tier 1 — Density/drift surface + Wishart-as-default  **[status: build, rungs 1–2, docs GREEN; SLiM scenario 1 PASS, scenario 2 deferred to recapitation]**

**Goal.** (a) Replace the scalar Wishart nugget `exp(σ)I` with a covariate-driven
per-site diagonal `diag(exp(Zγ))` — a drift / effective-size surface that
deconfounds local drift (∝ 1/N_e) from movement. (b) Document
Wishart-on-covariance as the recommended default and tighten `nu` guidance
(largely already shipped — mostly docs).

**Math (covariance model).**
`Σ = τ E + diag(n)`, `n_i = exp((Zγ)_i)`, where `Z` is `n_sites × q` site-level
covariates **including an intercept** and `γ` is length `q`. The intercept term
≡ the old `sigma`, so `Z = [1]` (intercept only) recovers `wishart_covariance`
exactly. Parameter vector `phi = (τ, σ, γ_1, …)` with `σ` the intercept and `γ`
the covariate slopes; the "diagonal parameters" are `phi[-1]`, aligned with the
columns of `Z`.

Per-parameter basis matrices:
- `B[[τ]] = E`
- `B[[d]] = diag(n · Z[,d])` for each diagonal parameter `d` (since
  `∂n_i/∂(γ or σ)_d = n_i Z_{id}`).

Gradient: `dPhi[d] = sum(diag(grad_Sigma) · n · Z[,d])` (only the diagonal of
`grad_Sigma` matters).

Hessian needs the **Σ-curvature** term the generic finish omits for non-`sigma`
diagonal params. Since `∂²Σ/∂phi_a∂phi_b = diag(n · Z[,a] · Z[,b])` for diagonal
params a,b:
`curvature[a,b] = sum(diag(grad_Sigma) · n · Z[,a] · Z[,b]) = (Zᵀ diag(w) Z)[a,b]`
with `w = n · diag(grad_Sigma)`; rows/cols for `τ` are 0. This **reduces** to the
existing `ddPhi["sigma","sigma"] += dPhi["sigma"]` when `Z = [1]`. ✓

The distance (`generalized_wishart`) model uses the same `B`; there a per-site
diagonal is exactly the EEMS/FEEMS within-deme diversity term `½(q_i + q_j)`,
so the drift surface is meaningful for both response types.

Bounds: `τ ≥ 0`; `σ, γ` unconstrained (nugget `= exp(·) > 0` for any real value).

**Interpretation.** `nugget_i ∝ local drift ∝ 1/N_e`. `γ_j > 0` ⇒ covariate j
raises drift / lowers effective size. Document the `N_e ∝ 1/nugget` reading.

**Implementation plan.**
- Refactor `.wishart_covariate_finish()` (+ `.wishart_covariate_covariance`,
  `.wishart_covariate_generalized`) to accept an optional `curvature` function;
  default `NULL` preserves the exact current `sigma` behavior (regression-tested
  by existing `test-wishart-covariates.R`).
- New file `R/wishart_drift_covariates.R`: factory
  `wishart_drift_covariates(x, coords, model, scale, intercept = TRUE)` building
  `Z`; `.wishart_drift_fit()` building `Σ`, `B`, and the curvature closure, then
  delegating to the (refactored) covariance/generalized finishers. Include a
  `subsetter` attribute for CV, mirroring `wishart_covariates`.
- Naming decision (revisit if desired): `wishart_drift_covariates`, params
  `tau`, `sigma` (baseline log-nugget), `gamma_<covariate>`.

**Tests** (`tests/testthat/test-wishart-drift-covariates.R`):
reduction to `wishart_covariance` at `γ = 0` / intercept-only; interface shapes;
`numDeriv` gradient/Hessian/partial checks; recovery via
`simulate_covariance_response` extended with a per-site nugget; end-to-end
`terradish()` fit on melip fixture.

**SLiM validation:** spatially varying local density → spatially varying N_e;
confirm the drift surface recovers the planted N_e gradient and beats the
scalar-nugget model by AIC.

**Open items:** extend `simulate_covariance_response()` to accept a per-site
nugget vector (currently scalar `sigma`); decide whether to also expose drift
covariates through `wishart_covariates` (unified `Σ = τE + Σλ_k K_k + diag(exp(Zγ))`)
— defer to after Tier 1 is green.

### Tier 2 — Hierarchical conductance field `log c = Xθ + u`  **[status: BUILT + tested + validated + docs; SLiM rung-3 via recapitation pending]**

**HARDENED (2026-06-08):** `R/hierarchical_conductance.R` →
`terradish_hierarchical()` (distinct pathway, per user request) + `conductance_field()`
extractor + `print`/`summary`/`coef` S3 methods (all exported; NAMESPACE regenerated).
- API: `terradish_hierarchical(formula, data, conductance_model, measurement_model,
  nu, field_resolution=G, tau2="reml"|numeric, tau2_grid, eps, ...)`.
- Fit: penalized L-BFGS over `(θ,u)` (gradient-only via the adjoint; φ warm-started
  across iterations); coarse G×G piecewise-constant field `Z_u` + GMRF `Q=L_coarse+εI`.
- τ²: Laplace empirical-Bayes marginal likelihood over a grid (clean interior
  optimum on melip at τ²≈0.32: logML −7695→−7634.9→−7681 across 0.01..100).
- **Tests GREEN** (`tests/testthat/test-hierarchical-conductance.R`, 13 checks):
  field-basis/GMRF well-formed; penalized gradient vs numDeriv; reduction to
  terradish at τ²→0 (`max|u|`<0.05, θ matches); recovery of an unmapped feature
  (`cor(field,blob)=0.89`) with θ protected (forestcover −0.63 w/ field vs −0.43
  no-field, true −0.60); `conductance_field()` raster; S3 methods.
- **Docs:** roxygen + man pages (`terradish_hierarchical.Rd`, `conductance_field.Rd`);
  vignette `vignettes/hierarchical-conductance.Rmd` (precompute pattern via
  `vignette-hierarchical.rds`).
- Validation scripts: `dev/proto_tier2.R`, `dev/check_tier2_api.R`.
- **TODO:** SLiM rung-3 for Tier 2 (unmapped-feature recovery in a real
  simulation) — do with recapitation (task #7); v2 refinements: REML over a
  continuous τ² (1-D optim), optional `u⟂X` projection, smooth (non-piecewise)
  field basis, field plotting S3 `plot()` method.

#### original prototype notes:

**Architecture (v1), the key realizations:**
- **Avoid the dense Hessian (F2) by switching the optimizer, not the algebra.**
  The gradient of the penalized objective w.r.t. *all* conductance parameters is
  ONE adjoint solve (`crossprod(df__dtheta_matrix, dl_dC)`), independent of how
  many parameters there are. So put `(θ, u)` in a single design and optimize with
  **L-BFGS** (gradient-only) instead of Newton (Hessian). `terradish_algorithm(...,
  gradient=TRUE, hessian=FALSE, partial=FALSE)` already returns exactly the
  likelihood + gradient needed per L-BFGS step; the package ships `bfgs.R`.
- **Field representation (coarse, piecewise-constant for v1).** Combined design
  `D = [X | Z_u]`, conductance `= exp(D %*% c(θ, u))`. `Z_u` (N_fine × m) is the
  indicator mapping each active cell to one of m coarse super-cells (a G×G coarse
  grid over the raster extent). Build the conductance model directly from `D`
  (same form as `.smooth_loglinear_conductance_from_matrix`). Coarseness is both
  tractability and the identifiability safeguard.
- **GMRF penalty.** `u ~ N(0, τ²(L_coarse + εI)⁻¹)`, `L_coarse` = Laplacian of the
  coarse-grid rook adjacency. Penalized objective
  `J = nll(θ,u) + (1/2τ²) uᵀ(L_coarse+εI) u`; gradient adds `(1/τ²)(L_coarse+εI)u`
  to the `u` block. The εI ridge makes the prior proper and anchors the field
  level (the conductance intercept is non-identifiable).
- **τ² for v1:** small grid + cross-validation (or held-out pairs) to pick τ²;
  validate the *mechanism*. Full REML/Laplace (needs `log|H_uu|` via CHOLMOD) is a
  v2 refinement.
- **Reduces to terradish** as τ²→0 (field shrinks to 0) — a required test.
- **Identifiability safeguards:** coarse field; informative τ² (not maximized
  freely); optional `u ⟂ col(X)` projection.

Prototype plan: (1) build `D`, `Z_u`, `Q_coarse` on the melip surface; (2)
penalized L-BFGS over `(θ,u)` at fixed τ²; (3) numDeriv-check the penalized
gradient; (4) recovery: simulate data with an UNMAPPED conductance feature, show
`u` localizes it and `θ` is protected vs the no-field fit. Then harden into
`R/hierarchical_conductance.R` + `terradish_hierarchical()`.

**PROTOTYPE VALIDATED (`dev/proto_tier2.R`, 2026-06-08):** the architecture works.
- Penalized gradient (adjoint likelihood grad + GMRF penalty grad) vs numDeriv =
  **1.6e-5**. L-BFGS over `(θ,u)` converges; no dense Hessian needed.
- τ²→0 collapses the field (`max|u|`=0.029) → reduces to plain terradish.
- Recovery on melip with an UNMAPPED Gaussian blob in log-conductance:
  `cor(field, true blob) = 0.887`; θ protected — forestcover **−0.631** (true
  −0.60) WITH the field vs **−0.432** (biased ~30% toward 0) for the no-field
  fit; altitude 0.537 (true 0.50). This is the omitted-variable-bias fix in
  action.
- Coarse field = G×G super-cells (here G=6 → m=33 occupied), piecewise-constant
  `Z_u`; `Q=L_coarse+εI`. `loglinear_from_matrix()` builds the combined `[X|Z_u]`
  conductance model; reused `terradish_algorithm(gradient=TRUE, hessian=FALSE)`.
- **TODO (harden):** `R/hierarchical_conductance.R` factory +
  `terradish_hierarchical()`; τ² selection (grid + CV / held-out pairs);
  identifiability (optional `u⟂X`, smoothness); S3 methods (coef/summary/plot of
  the field); tests (numDeriv, reduction, recovery); SLiM validation (best with
  recapitation); docs/vignette.

---

#### (original Tier 2 notes)

**Goal.** Add a coarse, spatially smooth Gaussian residual field `u` to the
log-conductance so the model is mechanistic where covariates explain the data and
flexible where they do not; `u` is a deliverable (map of unexplained structure).

**Math.** `log c = Xθ + u`, `u ~ N(0, τ²(L_grid + εI)⁻¹)` (GMRF = FEEMS
smoothness penalty). Penalized objective
`J = −log p(data | R(θ,u), φ) + (1/2τ²) uᵀ(L_grid+εI) u`; estimate `τ²` by
REML/empirical Bayes with a Laplace approximation around `û`.

**Why not free (F2).** Forward map and gradient reuse the existing
sparse-Cholesky Laplacian solve, but `u` adds one parameter per (coarse) field
cell, which breaks the dense per-parameter Hessian loop. **New machinery
required:** a sparse / Gauss–Newton / Fisher-scoring optimizer for `u` that uses
the sparse GMRF precision and never forms a dense `m×m` Hessian; plus the
log-determinant of the penalized Hessian (CHOLMOD can supply it) for the τ²
outer loop.

**Safeguards (identifiability).** Coarser field than covariate grid; informative
prior on τ²; optional `u ⟂ col(X)` projection.

**Validation:** SLiM with an unmapped barrier/corridor → confirm `θ` is protected
from omitted-variable bias and `u` localizes the unexplained feature.

**Entry points (anticipated):** new conductance model (closure adding `Zu` to the
linear predictor, where `Z` maps coarse field to cells) — but the field cannot
ride the F2 Hessian path; likely a new fitting routine beside
`terradish_algorithm()` / `terradish_optimize` that alternates `θ` (dense) and
`u` (sparse), with τ² in an outer loop.

### Tier 3 — Directional, non-reversible generator (covariate-driven)  **[status: design APPROVED; Phase-0 gate PASSED; building Phases 1–5]**

> **Phase-0 gate PASSED (`dev/proto_tier3_phase0.R`)**: directional flow IS
> identifiable from symmetric commute-time distances on a bounded lattice (sign +
> magnitude; γ̂ sign correct 6/6 under noise). Build justified. Numerical caveat:
> strong |γ| → exponential against-flow hitting times → ill-conditioned solves;
> use moderate-γ regime / log-scale / iterative fallback. Next: Phase 1 (directed
> generator model + edge-gradient helper), Phase 2 (Eigen::SparseLU engine +
> transpose adjoint), Phase 3 (API), Phase 4 (SLiM scen3 + identifiability sweep
> on coalescent data), Phase 5 (docs).


> Full design in **`dev/TIER3_DESIGN.md`**: directed generator
> `G_{a→b}=exp(s_{ab}ᵀθ + d_{ab}ᵀγ_dir)` → symmetric **commute time** (so the
> measurement models need NO change — the genetic response is symmetric) →
> `Eigen::SparseLU` (already a dependency) for the non-symmetric forward +
> transpose-adjoint solves, as a separate `terradish_directed_algorithm()` engine.
> Validated via SLiM Scenario 3 (asymmetric migration) on the recapitation
> pipeline. Central risk is *scientific* (is direction identifiable from symmetric
> distances?), so the plan gates the C++ build behind a cheap Phase-0 R prototype.


**Goal.** Represent directional gene flow (downstream/downhill/downwind,
source–sink) — the field's highest-demand feature.

**Approach.** Covariate-driven directional rates `G_{a→b} = exp(z_{ab}ᵀ γ)` with
`z_{ab}` directional covariates (slope, flow accumulation, wind) so asymmetry is
**forced by data** (few parameters, F2-friendly) rather than free per-edge.

**New machinery (the deep part).**
- Non-symmetric sparse solver (sparse LU via UMFPACK/KLU, or GMRES/BiCGSTAB);
  a new C++ path beside the SPD one.
- Transpose solves for the adjoint (self-adjointness of the resistance-distance
  gradient is lost when `G` is non-symmetric).
- A new distance object: hitting/commute or expected coalescence times from the
  generator, which are **asymmetric** (`T_ij ≠ T_ji`). Decide explicitly how to
  model `T_ij ≠ T_ji` (model the asymmetry vs. symmetrize) — the Gower transform
  and both Wishart models currently assume symmetric `R`/`Σ`.

**Scope.** Treat as a parallel engine / possible successor, not a flag on the SPD
path. Best pursued after Tiers 1–2 prove the lower-risk pop-gen content.

**Validation:** SLiM with asymmetric migration (e.g. downstream-biased) →
confirm directional model recovers the bias and that the symmetric model is
demonstrably biased (reproduce Lundgren & Ralph 2019 / Thomaz et al. 2019 in
miniature).

---

## 5. Code entry-point map

| Concern | File / symbol |
|---|---|
| Laplacian assembly (symmetric, `w=c_a+c_b`) | `src/radish.cpp::assemble_reduced_laplacian` |
| Reverse-mode edge→conductance | `src/radish.cpp::backpropagate_laplacian_to_conductance` |
| Solver dispatch (SPD only) | `R/radish_algorithm.R::.terradish_solver_setup/.solve` |
| Forward + adjoint to `theta` | `R/radish_algorithm.R::terradish_algorithm` |
| Conductance closures (`exp(Xθ)`) | `R/radish_conductance_model.R` |
| Nuisance subproblem (Newton on `phi`) | `R/radish_subproblem.R` |
| Scalar-nugget Wishart (`τE+exp(σ)I`) | `R/wishart_covariance.R` |
| Generic Wishart finish over basis `B` | `R/wishart_covariates.R::.wishart_covariate_finish` |
| IBE kernels (off-diagonal) | `R/wishart_covariates.R` |
| **Drift/density surface (Tier 1)** | `R/wishart_drift_covariates.R` (new) |
| Analytic simulation of covariance data | `R/genetic_distance.R::simulate_covariance_response` |
| Power / recovery harness | `R/covariance_response_power.R` |
| Top-level fit | `R/radish_optimize.R` (`terradish`) |

---

### Task #7 — Recapitation pipeline (tree-seq + pyslim + msprime)  **[status: DONE — pipeline built; scen1 validates drift surface; scen2 yields a scope finding; Tier 3 foundation ready]**

**Scenario 2 (continuous) outcome:** even recapitated + MAF-filtered,
`cor(diag(S),habitat)=+0.33` (wrong sign), robustly across all designs. This is a
**scope finding, not a bug**: the drift surface is a discrete-deme model;
continuous-space sampling confounds the covariance diagonal with local
relatedness + sampling scale, so it does not invert continuous-space density.
Documented in `wishart_drift_covariates` roxygen (`\strong{Scope.}`) and
`dev/slim/VALIDATION_SUMMARY.md`. Continuous-space density estimation belongs to
SBI/mapNN methods (feasibility Directions 2/6).


conda env: `C:\Users\peterman.73\AppData\Local\anaconda3\envs\slim\python.exe`
(msprime/tskit/pyslim/numpy). Pipeline components (all in `dev/slim/`):
- `scen1_ts.slim`, `scen2_ts.slim` — SLiM with `initializeTreeSeq()`, mutation
  rate 0, short forward burn-in, `treeSeqOutput()` (no mutation registry → fast).
- `recap_sample.py` — `pyslim.recapitate(ancestral_Ne, recombination_rate)` →
  `msprime.sim_mutations(BinaryMutationModel)` → sample per subpopulation
  (`--mode population`) or per nearest focal point via recorded individual
  locations (`--mode spatial --focal`) → write per-deme derived-allele counts
  `Y` (demes×loci) + haploid `N` + present deme ids.
- `run_scen1_recap.R`, `run_scen2_recap.R` — orchestrate SLiM→Python→R,
  MAF-filter, `cov_from_biallelic`, fit drift vs scalar.

**KEY FINDINGS:**
- Pipeline is functional end-to-end (SLiM→recap→mutations→covariance→fit).
- **MAF filtering is essential.** Recapitated data has a realistic SFS (many rare
  variants); `cov_from_biallelic`'s `1/sqrt(p(1-p))` standardization upweights
  rare alleles, inflating `diag(S)` to ~1000 and swamping the signal. A standard
  MAF≥0.05 filter restores an O(1)-scaled covariance and the drift signal.
  (Worth surfacing in user docs/QC guidance.)
- **Recapitated Scenario 1 (discrete) validates the drift surface:** with low
  migration (0.005) + long burn-in (4000) + MAF≥0.05, `cor(diag(S),cov)=-0.50`,
  `gamma_cov=-0.095` (correct sign), AIC prefers drift. Confirms the Tier 1 drift
  surface on *realistic* deep-coalescent data (third independent confirmation
  after analytic rungs 1-2 and forward-only scen1 cor=-0.85).
- **IBR vs drift trade-off:** strengthening the drift signal (low migration /
  large size contrast) makes demes nearly independent → little IBR covariance →
  `tau=0`. A joint IBR+drift scenario with `tau>0` needs a separate
  conductance(migration) gradient distinct from the size(drift) gradient. Not
  essential for validating the new features (IBR estimation is terradish's
  established core), but noted as a richer future scenario.

---

## 6. Progress log (append-only)

- **2026-06-08 (Tier 2 + recapitation, latest)** — Tier 2
  `terradish_hierarchical()` hardened/tested(13)/validated/documented. Recapitation
  pipeline built on conda `envs/slim`; drift surface validated on realistic
  discrete-deme data (scen1 recap `cor=-0.50`, AIC-preferred); continuous-space
  scen2 → documented scope finding (drift surface ≠ continuous-density estimator),
  added to `wishart_drift_covariates` roxygen + VALIDATION_SUMMARY. Tasks #1–5, #7
  complete; #6 (Tier 3) is next. Nothing committed.

- **2026-06-08** — Created `dev` branch. Wrote feasibility assessment
  (`dev/FEEMS_integration_feasibility.md`). Reconnaissance of engine confirmed
  F1–F4 and the measurement-model contract. Wrote this roadmap. Began Tier 1:
  derived the drift-surface math (basis matrices + curvature term reducing to the
  existing `sigma` special case), confirmed `.wishart_covariate_finish` is the
  reusable engine. Toolchain verified (R 4.6.0; devtools/numDeriv/testthat/roxygen2).
- **2026-06-08 (cont.)** — Tier 1 BUILD complete: refactored
  `.wishart_covariate_finish`/`_covariance`/`_generalized` to accept an optional
  `curvature` (default `NULL` = exact prior behavior); added
  `R/wishart_drift_covariates.R` (factory + `.wishart_drift_fit` +
  `.wishart_drift_design`). Validation rungs 1–2 GREEN:
  (1) numDeriv — gradient/Hessian/`gradient_E` match to 1e-9 (covariance) and
  1e-9 (generalized); full `terradish_algorithm` θ-gradient vs finite difference
  = 2.9e-7; intercept-only and γ=0 reduce to `wishart_covariance` to machine eps.
  (2) recovery — γ̂ mean 0.688 (true 0.70) over 5 draws; AIC prefers drift 5/5
  vs scalar nugget. No regressions (wishart-covariance 54, wishart-covariates 33,
  wishart-drift-covariates 31 tests pass). `man/wishart_drift_covariates.Rd`
  generated; NAMESPACE export added; DESCRIPTION untouched.
  Validation scripts: `dev/check_tier1.R`, `dev/check_regress.R`, `dev/check_e2e.R`.
  - **NOTE / latent bug found (not mine):** `terradish_algorithm(..., validate=TRUE)`
    does NOT thread `nu` into its recursive numerical-check calls, so the built-in
    numerical validation is broken for ALL nu-requiring (Wishart) measurement
    models (`grad_Sigma` becomes length-0 → "non-conformable"). Worked around by
    finite-differencing the profiled objective manually. Worth a small dedicated
    fix later (thread `nu` through the `validate` block).
  - **Test execution note:** `testthat::test_file()` after
    `devtools::load_all()` SEGFAULTS on Windows (DLL reload) — pre-existing, not
    code-related. Use `Rscript -e 'devtools::test(filter="...")'` instead.
  - **TODO Tier 1:** SLiM rung-3 validation (SLiM not yet located on this
    machine — see below); docs (vignette section + Wishart-default guidance);
    optional: extend `simulate_covariance_response()` to accept a per-site nugget
    vector.
- **2026-06-08 (docs + SLiM setup)** — Tier 1 DOCS drafted (task #4):
  added "Modeling a drift / effective-size surface" section to
  `vignettes/wishart-covariance.Rmd` (live demo; smoke-tested end-to-end on the
  full 37-site melip surface — γ̂ = −0.598 vs true −0.6, τ̂ = 1.003, conductance
  forestcover 1.10 / altitude −0.515, `summary()` shows the drift `phi` table with
  curvature-corrected CIs). D5: fixed a stale microsatellite `nu` recommendation in
  the vignette's "role of `nu`" section to match the conservative locus-count
  guidance. Added `@seealso` cross-refs (`wishart_covariance`, `generalized_wishart`,
  `wishart_covariates` now point to `wishart_drift_covariates`); regenerated man pages.
  - **SLiM RUNG-3 SETUP** — SLiM 5.2 at `C:\msys64\mingw64\bin\slim.exe`. The
    canonical pyslim/msprime/tskit pipeline is NOT available here (no conda; msys
    Python 3.14 lacks the packages and is too new for wheels), so validation uses
    **SLiM-direct allele-frequency output** + finite sampling in R (no Python).
    SLiM 5.x API notes: `individual.genomes` → `haplosomes`; `outputVCF()` is not a
    Haplosome method → switched to `sim.mutationFrequencies()` CSV output.
    Files: `dev/slim/scen1_stepping_stone.slim` (WF 2D stepping-stone, deme size
    varies along a gradient → varying Ne, uniform migration) and
    `dev/slim/run_scen1.R` (runs SLiM, binomial-samples, `cov_from_biallelic`,
    fits drift vs scalar; expects γ_cov < 0 and AIC preferring drift).
    User wants BOTH scenarios: (1) stepping-stone varying K [running], and
    (2) continuous-space nonWF habitat-driven density [TODO].
  - **TODO Tier 1:** finish SLiM scenario 1 analysis (running in background) +
    build/run scenario 2 (continuous space); then mark Tier 1 fully green.
- **2026-06-08 (scenario 2 + Python)** — Built scenario 2
  (`dev/slim/scen2_continuous_space.slim` + `run_scen2.R`): continuous-space nonWF
  following the user's canonical template (Gaussian competition via
  `localPopulationDensity`, mate choice, natal Gaussian dispersal via
  `pointDeviated`), habitat gradient `HAB0 + HABSLOPE*x` drives local density →
  local Ne. Focal sampling by Voronoi assignment to an FGRID×FGRID grid, moved
  into focal subpops via `takeMigrants`, per-focal `mutationFrequencies` output;
  per-focal sample sizes handled in R (`cov_from_biallelic` N as per-population
  vector). SLiM 5.x note: `%`/`/` yield float → must `asInteger()` array indices.
  Probed OK; full run launched in background **in parallel with scenario 1**.
  - **PYTHON RECOMMENDATION (open):** forward-only neutral burn-in is the runtime
    bottleneck (scenario 1 is slow). Recommended the user install a conda-forge
    env (`python=3.12` + `msprime tskit pyslim numpy scikit-allel pandas`) to
    enable the canonical **tree-seq + recapitation** pipeline (~100× faster, true
    equilibrium ancestry, matches their template, needed for Tier 3 asymmetric
    validation). If installed, switch both harnesses to tree-seq recording →
    `pyslim.recapitate` → `msprime.sim_mutations` → sample → covariance.
- **2026-06-08 (SLiM rung-3 results)** — Full results + interpretation in
  `dev/slim/VALIDATION_SUMMARY.md`. **Scenario 1 (discrete stepping-stone): PASS** —
  `cor(diag(S),cov)=-0.852`, `gamma_cov=-0.147` (correct negative sign),
  AIC prefers drift (27908<27961). Caveats (forward-only artifact, not estimator):
  `tau=0` (off-diagonal IBR not captured; young spatially-clustered mutations
  don't match a smooth resistance distance) and muted slope magnitude.
  **Scenario 2 (continuous space): inconclusive** — even with interior focal
  points + larger neighbourhoods, `diag(S)`~hundreds, `cor(diag(S),habitat)≈0.07`,
  signal swamped by the forward-only young-mutation / high-Fst artifact. A clean
  continuous-space (and joint IBR+drift) test needs **recapitation** (Python).
  Net: drift surface validated (rungs 1–2 exact; scenario 1 biological PASS);
  scenario 2 deferred to the recapitation pipeline.
  - Had to make scenario 1 leaner (5×5 demes, L=4e5, 2500 gens) — the original
    6×6/4000-gen WF run ballooned the mutation registry and never finished.
  - **DECISION POINT for user:** proceed to Tier 2 now (revisit scenario 2 +
    sharpen via recapitation later, which Tier 3 needs anyway) vs install Python
    first and clean up scenario 2 before advancing.
