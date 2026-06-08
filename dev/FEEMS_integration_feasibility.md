# Feasibility assessment: explicit population-genetic integration in terradish

**Context.** This note responds to the two companion documents in
`Research/Manuscripts/OSU/SpatialResistance/`
(`FEEMS_terradish_critical_review.docx` and `FEEMS_terradish_technical_appendix.pdf`),
which propose seven directions for evolving terradish toward FEEMS-style and
more explicitly population-genetic models. It evaluates each proposed direction
*against the actual terradish engine* (not against an idealized version), with
particular attention to (1) accounting for asymmetric gene flow, flagged as the
field's highest-demand feature, and (2) deeper population-genetic integration
generally.

The verdict in one sentence: the reviews are conceptually excellent and the
shared-mathematics framing is correct, but the appendix systematically
*understates the engineering distance* between "the gradient flows through the
same solve" and "terradish can fit this," and the single most-demanded feature
(directional gene flow) maps onto the single most invasive change to the
numerical core.

---

## 1. What the engine actually is (the binding constraints)

Everything below follows from four facts about the current implementation. These
are the load-bearing walls; any proposal has to be read against them.

**F1 — The Laplacian is symmetric and the entire solver stack is SPD-only.**
`src/radish.cpp::assemble_reduced_laplacian` builds the reduced Laplacian with a
single edge weight per pair, `w = conductance(i) + conductance(j)` — symmetric in
`(i,j)` by construction. The R layer wraps this with `forceSymmetric` and factors
it with `Matrix::Cholesky` / CHOLMOD; the alternative backends (simplicial and
supernodal LL/LDL, smoothed-aggregation AMG, IC/Jacobi PCG) are *all* SPD
methods. There is no sparse LU, no GMRES/BiCGSTAB, nowhere in the stack that can
factor a non-symmetric operator.

**F2 — The optimizer is built for a handful of conductance parameters.**
`terradish_algorithm()` computes the gradient with one batched solve, but the
Hessian is formed by looping `idx <- seq_along(theta)` and solving the reduced
system with `n_demes × length(theta)` right-hand sides, then assembling a *dense*
`length(theta) × length(theta)` Hessian that the Newton/quasi-Newton drivers
invert. This is excellent for `theta` of length 2–20. It is structurally wrong
for `theta` of length 10³–10⁵.

**F3 — Conductance models are closures, and that is the clean extension seam.**
A `terradish_conductance_model` is a function of `theta` returning
`conductance = exp(Xθ)` plus first/second derivative callbacks
(`df__dtheta`, `d2f__dtheta_dtheta`, …). Anything expressible as
"per-cell conductance = smooth positive function of parameters, with
derivatives" slots in here without touching the solver — *subject to F2's
dimensionality limit.*

**F4 — The measurement models carry a thin, scalar nuisance vector.**
`wishart_covariance` models `Σ = τE + exp(σ)I`: one IBR scale `τ`, one scalar
log-nugget `σ`. `wishart_covariates` already extends this additively with
outer-product IBE kernels `Σ = τE + Σ_k λ_k z_k z_kᵀ + exp(σ)I`. The
nuisance subproblem (`radish_subproblem`) fits these by Newton in low dimension.
This is the seam for richer *likelihoods* (drift surfaces, per-node variance).

---

## 2. The two kinds of "asymmetry" (the most important distinction for this project)

"Accounting for asymmetry in movement" is treated in the literature as one
feature. For implementation it is two completely different problems, and
conflating them is the main risk in scoping this work.

### 2a. Reversible asymmetry (tractable on today's solver)

A Markov chain can have direction-dependent edge *rates* (`m_ij ≠ m_ji`) and
still satisfy detailed balance with respect to a non-uniform stationary
distribution `π` — for example when migration is higher out of large populations
than out of small ones. Such a generator `G` is *similar to a symmetric matrix*
via `D^{1/2} G D^{-1/2}` with `D = diag(π)`. The symmetrized operator is SPD and
can be factored with the **existing CHOLMOD/AMG stack**.

What this costs in code: a new assembly (the symmetrized generator, not
`c_a + c_b`) and a back-transform of the solution; the *solver* and its adjoint
are reusable. This is essentially the population-genetic content of Direction 2
(a density / effective-size surface) and is feasible in the near term.

What this does **not** buy: it cannot represent genuinely directional gene flow.
Detailed balance means there is still an undirected potential underneath; you get
heterogeneous *magnitudes*, not true downstream/downhill/downwind bias.

### 2b. Non-reversible asymmetry (the high-demand feature; successor-method scale)

Directional advection — downstream in a riverscape, downhill, downwind,
source→sink — breaks detailed balance. The generator is genuinely non-symmetric,
has no SPD structure, and **cannot be touched by any current solver (F1).**
Implementing it requires, at minimum:

1. **A non-symmetric sparse linear-algebra backend** — sparse LU (UMFPACK/KLU)
   or non-symmetric Krylov (GMRES/BiCGSTAB) with preconditioning. This is a new
   numerical core sitting beside the SPD one, not a patch to it.
2. **Transpose solves for the adjoint.** The clean gradient identity
   `∂R_ij/∂w_e = −[b_eᵀ L⁺ (e_i − e_j)]² ≤ 0` is a *self-adjointness* result that
   holds only because `L` is symmetric. With a non-symmetric generator the
   reverse-mode pass needs solves against both `N` and `Nᵀ` (the fundamental
   matrix and its transpose), so the elegant single-factorization economy is
   partially lost.
3. **A new distance object.** The target is no longer a symmetric resistance
   distance. Hitting/commute or expected coalescence times from a non-reversible
   generator give an *asymmetric* matrix, `T_ij ≠ T_ji`. The Gower double-centering
   `Σ = −½ J R J` and the Wishart measurement models assume a symmetric `R`/`Σ`.
   You must either model the asymmetric `T` directly (a real modelling question,
   not just numerics) or symmetrize `(T + Tᵀ)/2` and forfeit much of the point.

This is the SAMC (Fletcher et al. 2022) / FRAME (Shen et al. 2025) territory the
appendix points to (Section H). The appendix is right that "the same adjoint
approach applies" *in principle*; it understates that the approach now runs on
linear algebra terradish does not contain.

**Strategic consequence.** The field's loudest request lands on the most
invasive change. Two corollaries shape the realistic path:

- *Make asymmetry covariate-driven, not free.* Free directional edge rates double
  the edge parameters in an inverse problem that already barely identifies the
  symmetric ones (Section 4). The identifiable version is a *directional covariate*
  with one or a few coefficients — `G_{a→b} = exp(z_{ab}ᵀ γ)`, where `z_{ab}` is
  slope, flow accumulation, or wind from `a` to `b`, so `G_{a→b} ≠ G_{b→a}` is
  forced by data, not estimated per edge. This keeps `theta` small (F2-friendly)
  even though the operator is non-symmetric.
- *Treat it as a separate engine / successor method,* not a measurement-model
  option bolted onto the SPD path.

---

## 3. Direction-by-direction feasibility

Ordered by the reviews' own numbering. "Effort" is relative to the current
codebase; "new machinery" flags whether the work is reachable from F1–F4 or
needs something genuinely new.

### Direction 1 — Hierarchical conductance surface `log c = Xθ + u` (covariate-FEEMS hybrid)

**Verdict: the right idea; feasible; but NOT "free," contrary to the appendix.**

This is the strongest recommendation in both documents and I agree it is the
highest-value *tractable* target. It nests terradish (`u = 0`) and FEEMS (drop
`Xθ`), and the residual field `u` is a genuine deliverable (a map of what the
covariates failed to explain) that directly attacks the omitted-covariate bias
behind the ResistanceGA cost-recovery failures (Daniel et al. 2024; Graves et al.
2013).

Where the appendix is right: the *forward* conductance map `exp(Xθ + u)` and the
*gradient* of the penalized likelihood do flow through the same sparse-Cholesky
Laplacian solve. The GMRF prior `u ~ N(0, τ²(L_grid + εI)⁻¹)` is exactly the
smoothness penalty, and `u = 0`/`τ²→0` recovers the current model.

Where the appendix is too optimistic ("No new numerical machinery is required
beyond what terradish already has"):

- **F2 breaks.** `u` adds one parameter *per field cell*. The Hessian loop solves
  once per parameter and assembles a dense Hessian; with a per-cell field that is
  thousands of solves per Newton step and a dense Hessian of size `(p + m)²`.
  Intractable as written. The fix is real new code: a sparse / Gauss–Newton /
  Fisher-scoring optimizer for `u` that exploits the GMRF precision
  (`L_grid + εI` is sparse) and never forms a dense `m × m` Hessian.
- **Variance components are new.** Estimating `τ²` by REML / empirical Bayes with
  a Laplace approximation requires the log-determinant of the penalized Hessian in
  `u` — a large sparse determinant that nothing in terradish currently computes.
  (CHOLMOD can supply it from the factor; the plumbing does not exist yet.)
- **Identifiability.** `θ` and `u` compete for the same signal. The reviews'
  safeguards are correct and necessary: an informative prior on `τ²`, a *coarser*
  field than the covariate grid, and optionally constraining `u ⟂ col(X)`.

Net: a focused but real project — a new latent-field optimizer plus
variance-component estimation — reusing the existing forward/adjoint solve. Worth
doing; budget it as such.

### Direction 2 — Density / effective-size surface (deconfound movement from abundance)

**Verdict: the most tractable genuinely population-genetic extension. Recommend first.**

Conductance currently conflates "how readily an organism moves through a cell"
with "how many organisms a cell supports." A node-level effective-size surface
`log N_i = Zγ` separates them. In the covariance likelihood this is a
*node-specific drift / nugget* term: replace the scalar `exp(σ)I` in
`wishart_covariance` with `diag(exp(Zγ))` (or `τE + diag(exp(Zγ))`).

Why it is easy relative to everything else:

- `Σ` stays SPD; `E` is unchanged; the expensive Laplacian solve is untouched.
- It adds parameters only to the low-dimensional nuisance subproblem (F4), where
  `γ` has a handful of entries — no collision with F2.
- It is mechanistic and interpretable in exactly terradish's idiom (covariates on
  a surface), unlike FEEMS's free node variances.

Management payoff is direct and is the review's strongest practical argument: a
low-similarity region caused by a density trough calls for different action than
one caused by a movement barrier. This is real population genetics (local drift ∝
1/N_e) for modest cost, and it is the natural bridge to the *reversible* asymmetry
of Section 2a.

### Direction 3 — Directional, non-reversible generators (the asymmetry request)

**Verdict: highest demand, deepest surgery. Successor-engine, covariate-driven.**

See Section 2b in full. This is where the loudest user demand and the largest
implementation cost coincide. It is achievable, and the gradient theory
generalizes, but it needs a non-symmetric solver, transpose-based adjoints, and a
reconsidered (asymmetric) distance/measurement model. Scope it as a parallel
engine with covariate-driven directional rates (few parameters), not as a flag on
the SPD path. It is the right *medium-term research program*, not a near-term
release feature.

### Direction 4 — Target a coalescence time, not a commute time

**Verdict: partly already absorbed; the novel part collapses into Direction 3.**

Under symmetric migration at equilibrium, expected coalescence time is affine in
resistance distance — which the measurement models' `α + β·d` already absorb. So
the "lighter route" buys little beyond what terradish does implicitly. The
"heavier route" — computing expected coalescence times from a (possibly
non-reversible) generator — *is* the asymmetric-distance problem of Direction 3
and inherits its costs. Treat Direction 4 as a motivation for Direction 3 rather
than as a separable task.

### Direction 5 — Fit the data, not a distance matrix

**Verdict: largely already shipped. Low-cost to elevate; real work only in extensions.**

terradish already has `generalized_wishart` (Wishart on the distance matrix) and
`wishart_covariance` + `wishart_covariates` (Wishart on the covariance). The
recommendation to "make the generative likelihood the default rather than an
option" is mostly *documentation and defaults*, not new math — and is a good,
cheap win. The genuinely new work is extending the generative target to
identity-by-descent or site-frequency-spectrum summaries, which is substantial
and overlaps with Directions 6–7.

### Directions 6 & 7 — Simulation-based / genealogy-based inference; time-stratified (IBD) surfaces

**Verdict: out of scope for terradish-the-package; this is a successor tool.**

These abandon the analytic link function (SBI on spatial coalescent simulations;
covariate-parameterized conductance fit to long-vs-short IBD blocks for
contemporary connectivity). They are the most scientifically ambitious and
squarely address the timescale-mismatch problem the review rightly calls the
field's most serious practical issue. But they share almost no machinery with the
current package. terradish's durable contribution in that world is the
*interpretable covariate parameterization of conductance and density* as the
thing being inferred; that idea should be carried into a new project rather than
retrofitted here.

---

## 4. The constraint that governs all of it: the ill-posed inverse

Both reviews are blunt that genetic data constrain the conductance surface
*weakly* (Graves et al. 2013; Daniel et al. 2024: optimized costs depart from
truth even when prediction is good). This is the single most important thing to
hold in view when adding flexibility, and it reorders the priorities above:

- Every new degree of freedom must be *paid for* — by a strong prior or by new
  information. Direction 1 works precisely because the GMRF prior pays for the
  field. Direction 3's *free* directional rates are not paid for and would
  amplify non-identifiability; its *covariate-driven* form is, which is why that
  is the only realistic version.
- "Added flexibility without added information will not help, and may hurt" (the
  review's own warning) is the acceptance test. A density surface (Direction 2)
  adds information structure (drift ∝ 1/N_e on the diagonal, a distinct signal
  from the off-diagonal IBR). Asymmetric edge weights, in general, do not add
  information — they add parameters — unless directional covariates supply it.

This argues for sequencing by *information added per parameter*, not by
conceptual ambition.

---

## 5. Recommended sequencing

| Tier | Work | Reuses F1 solver? | New machinery | Effort | Pop-gen value |
|------|------|-------------------|---------------|--------|---------------|
| 1 | **D5 defaults**: make Wishart-on-covariance the documented default; tighten `nu` guidance | yes | none | low | medium (already present) |
| 1 | **D2 density/drift surface**: `diag(exp(Zγ))` nugget in `wishart_covariance` | yes | small (nuisance gradients) | low–med | high |
| 2 | **D1 hierarchical field** at coarse resolution | forward/adjoint yes | **new**: sparse latent-field optimizer + REML/Laplace for τ² | med–high | high |
| 3 | **D3 directional generator** (covariate-driven rates) | **no** | **new**: non-symmetric solver + transpose adjoint + asymmetric distance/measurement model | high | high (the demand) |
| — | D4 | folds into D3 | — | — | — |
| 4 | **D6/D7** SBI, genealogy, IBD time-stratification | no | successor tool | very high | very high (timescale) |

**Concrete first steps (Tier 1, buildable now):**

1. *Density/drift surface.* Add a measurement-model variant where the nugget is
   `diag(exp(Zγ))` instead of `exp(σ)I`. Entry points: `R/wishart_covariance.R`
   (the `Σ = τE + exp(σ)I` construction and its `phi` gradients) and
   `radish_subproblem`. `E` and the Laplacian solve are untouched; this is
   contained to F4. Deliver alongside a simulation that recovers a known `γ`.
2. *Wishart-as-default.* Documentation/vignette work plus a recommendation in
   `terradish()`'s guidance; no engine change.

**Before any Tier-2/3 work:** a small identifiability study on the melip data and
a simulated landscape — how much does a coarse field, or a single directional
covariate, actually change recovered `θ` and predictive accuracy? Given Section 4,
this should gate the larger investments.

---

## 6. Corrections / cautions on the source documents (for the manuscript)

These are accuracy notes, not disagreements with the thrust:

- **Edge-weight function.** The appendix (Section F) says `w_e = h(c_a, c_b)`
  "for example their mean or their conductance combined in series." The actual
  implementation uses the **sum**, `w_e = c_a + c_b`. Minor, but worth stating
  correctly since the gradient algebra depends on it.
- **"No new numerical machinery" (Direction 1 / Section G).** As above, the
  forward and gradient computations reuse the solve, but the *optimizer* and
  *variance-component estimation* are genuinely new. The claim is true of the
  linear algebra per evaluation and false of the fitting procedure; the
  manuscript should make that distinction explicit.
- **"The same adjoint approach applies" (non-reversible generator / Section H).**
  Correct in principle, but it now requires transpose solves against a
  non-symmetric factorization — i.e., the adjoint is reusable as a *concept*, not
  as terradish's current *code*.
- **Asymmetric distance + Wishart.** The Gower transform and both Wishart models
  assume a symmetric `R`/`Σ`. Any non-reversible target needs an explicit
  decision about modelling `T_ij ≠ T_ji`; the appendix does not address this and
  it is the crux of making Direction 3 statistically coherent.
- No fabricated citations or interpretations are introduced here; all references
  are those already in the two source documents.

---

*Prepared on branch `dev`. Nothing committed. Companion source documents:*
*`Research/Manuscripts/OSU/SpatialResistance/FEEMS_terradish_critical_review.docx`,*
*`…/FEEMS_terradish_technical_appendix.pdf`.*
