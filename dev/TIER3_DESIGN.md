# Tier 3 design: directional, non-reversible conductance

**Status: DESIGN FOR REVIEW — no code yet.** Companion to
`dev/IMPLEMENTATION_ROADMAP.md` (§4 Tier 3) and `dev/FEEMS_integration_feasibility.md`
(§2 "two kinds of asymmetry", §3 Direction 3). This document answers the four
questions raised before committing to the build: (1) where the non-symmetric
solver goes, (2) the covariate-driven directional-rate parameterization, (3) the
asymmetric-distance measurement-model question, and (4) the SLiM validation
design. It adds the math, the adjoint, a phased plan, decision points, and an
honest identifiability assessment.

Tier 3 is the package's highest-demand feature (directional gene flow:
downstream / downhill / downwind / source–sink) and its most invasive change: it
breaks the symmetric-positive-definite (SPD) assumption the entire current solver
stack rests on (engine constraint **F1**).

---

## 0. The key realization that makes Tier 3 tractable and clean

The genetic data are **symmetric**: observed genetic distance / covariance
`S_ij = S_ji`. Under a structured coalescent with an **asymmetric** backward
migration matrix, the *expected pairwise coalescence time* is still symmetric
(`E[T_coal(i,j)] = E[T_coal(j,i)]` — coalescence of two lineages is symmetric
under relabeling), even though the lineages move under asymmetric rates.

Two consequences shape the whole design:

1. **The asymmetry lives entirely in the forward map** `θ, γ → R(θ, γ)`, not in
   the likelihood. We map the asymmetric generator to a **symmetric** model
   distance (a directed commute time; §3). Therefore **the measurement models
   need no change** (question 3 resolved — see §4).
2. **Directionality may be only weakly identifiable from symmetric distances**
   (Lundgren & Ralph 2019 show asymmetry leaves a distinctive but subtle
   signature). The validation must probe this honestly (§6, §7).

---

## 1. The mathematical object: directed generator → commute time

### 1.1 What replaces the symmetric Laplacian

Current engine: symmetric Laplacian `L = D − W`, edge weight `w_e = c_a + c_b`,
resistance distance `R_ij = (e_i−e_j)ᵀ L⁺ (e_i−e_j)`, covariance
`E = −½ J R J`. (See `src/radish.cpp::assemble_reduced_laplacian`,
`graph_rhs_crossprod`; `R/radish_algorithm.R::terradish_algorithm`.)

Tier 3: a **continuous-time Markov generator** `G(γ)` on the graph nodes with
directional off-diagonal rates and zero row sums:

> `G_{ab} = rate(a → b) = exp(z_{ab}ᵀ γ)` for neighbours `a,b`;  `G_{aa} = −Σ_b G_{ab}`.

`G` is generally **non-symmetric** (`G_{ab} ≠ G_{ba}`). Circuit theory / the
symmetric Laplacian is the special case where `G` is reversible (§2 gives the
parameterization that makes this nesting explicit).

### 1.2 The distance: directed hitting → symmetric commute time

For a target deme `j` made **absorbing**, the vector of expected hitting times
`h_{·→j}` from all transient nodes solves a single sparse linear system

> `Q_j h_{·→j} = −1`,

where `Q_j` is `G` with row/column `j` removed (the sub-generator on transient
states). This is the spatial-absorbing-Markov-chain view (Fletcher et al. 2022);
the fundamental matrix `N = (I − P_sub)⁻¹` of the embedded chain is the
discrete-time analogue. One solve per focal absorber gives hitting times to that
focal from every node, hence the full `n × n` directed hitting-time matrix
`H_ij = h(i → j)` in **n sparse solves** — the same order as the current
resistance-distance computation, just non-symmetric.

Because the data are symmetric, we use the **commute time** (symmetric by
construction):

> `R^dir_ij = h(i → j) + h(j → i)`.

- **Symmetric** ⇒ feeds the existing measurement models / Gower transform
  unchanged.
- **Reduces to resistance distance**: when `G` is reversible, commute time
  `= 2m · R_ij` (the appendix's `C_ij = 2m R_ij`), so Tier 3 with `γ_dir = 0`
  reproduces the current model up to the known constant absorbed by the
  measurement model's scale (`β`/`τ`).
- **Captures directionality**: the individual hitting times `h(i→j) ≠ h(j→i)`
  encode the asymmetric flow even though their sum is symmetric.

### 1.3 What we are NOT doing in Tier 3a (and why)

The *exact* structured-coalescent coalescence time on a fine raster requires a
CTMC over **pairs** of node-locations — `O(d²)` states (`d` = raster cells, e.g.
22 000² ≈ 5×10⁸). Intractable on terradish's grids. That full coalescent target
is **Tier 3b** (research / successor; deme-scale only, à la Lundgren & Ralph),
documented as future work, not part of this build.

So Tier 3a = **directed commute time from a covariate-parameterized non-reversible
generator** (single-lineage, `d` states, sparse non-symmetric solves). This is
the implementable "accounts for asymmetry in movement" feature and the standard
SAMC formulation; commute time approximates coalescence time exactly under
reversibility and heuristically (directionally-aware) otherwise.

---

## 2. Covariate-driven directional-rate parameterization (question 2)

`z_{ab}` are **edge-and-direction-specific** covariates (`z_{ab} ≠ z_{ba}`).
Decompose the log-rate into a symmetric and an antisymmetric part:

> `log G_{a→b} = s_{ab}ᵀ θ  +  d_{ab}ᵀ γ_dir`,  with  `s_{ab} = s_{ba}` (symmetric)  and  `d_{ab} = −d_{ba}` (antisymmetric).

- **Symmetric part** `s_{ab}ᵀθ` — the usual conductance term built from node
  covariates (e.g. `s_{ab} = ½(x_a + x_b)`), recovering the current model.
- **Antisymmetric part** `d_{ab}ᵀγ_dir` — directional covariates, e.g.:
  - elevation drop `d_{ab} = h_a − h_b` (downhill easier),
  - flow accumulation / stream direction (downstream > upstream),
  - prevailing wind projected onto the edge.

Properties:
- **Parsimonious & identifiable**: a *handful* of `γ_dir` coefficients drive all
  asymmetry (not free per-edge rates, which would be hopelessly ill-posed given
  weak genetic identifiability). This is the feasibility doc's central
  recommendation — *covariate-driven, not free*.
- **Nests the symmetric model exactly**: `γ_dir = 0 ⇒ G_{a→b} = G_{b→a}`
  (reversible) ⇒ current terradish.
- **Sign-interpretable**: e.g. `γ_dir > 0` on elevation drop ⇒ net downhill bias.

**New abstraction.** The conductance-model closure (`R/radish_conductance_model.R`,
cell→conductance) is generalized to a **generator model**: given `(θ, γ_dir)`,
return the directed edge rates `G_{a→b}`, `G_{b→a}` and their derivatives. This
parallels the existing factory but operates at the directed-edge level.

**New graph helper.** A function to build per-directed-edge antisymmetric
covariates from a raster, e.g. `edge_gradient(raster, graph)` returning
`d_{ab} = value_a − value_b` for each directed edge (elevation, flow-potential,
wind potential). Plugs into `conductance_surface()` output (`graph$edge_pairs`,
`graph$vertex_coordinates`).

---

## 3. Non-symmetric solver placement in the dispatch (question 1)

### 3.1 Backend: `Eigen::SparseLU` (already a dependency)

`src/` already links **Eigen via RcppEigen** (38 refs). `Eigen::SparseLU`:
- factors a general (non-symmetric) sparse matrix once,
- provides **both** `solve` (forward, against `Q`) and `transpose().solve`
  (adjoint, against `Qᵀ`) from the **same** factorization — exactly the two
  solves Tier 3 needs (§5),
- no new package dependency (CHOLMOD/AMGCL stay for the SPD path; Armadillo lacks
  SuperLU in this build, so we do not use `arma::spsolve`).

For very large graphs, an iterative fallback (`Eigen::BiCGSTAB` /
`Eigen::GMRES` with `IncompleteLUT` preconditioning) can be added later; direct
`SparseLU` covers terradish's stated domain (<1M cells) for v1.

### 3.2 Placement: a **parallel engine**, not a flag on the SPD path

Recommended: a new `terradish_directed_algorithm()` paralleling
`terradish_algorithm()`, plus C++ helpers in `src/`, rather than overloading the
existing solver dispatch (`.terradish_solver_setup` / `.terradish_solver_solve`).
Rationale:
- The forward model differs throughout (assemble directed sub-generator, absorbing
  solves for hitting times, symmetrize to commute time) — not just the linear
  solver.
- Isolating the non-symmetric algebra from the validated SPD path **eliminates
  regression risk** to Tiers 0–2.
- Matches the feasibility doc's "separate engine / successor, not a flag"
  guidance.

What is **reused unchanged**: the optimizer drivers (Newton/BFGS,
`R/radish_optimize.R`), the measurement models (§4), the nuisance subproblem,
and the conductance-surface graph object.

New C++ (parallel to the named symmetric functions):
- `assemble_directed_subgenerator(rate_ab, rate_ba, edge_pairs, absorbing)` →
  sparse `Q` (cf. `assemble_reduced_laplacian`).
- hitting-time RHS (`−1`) and the absorbing-state bookkeeping (cf. `graph_rhs`).
- `backpropagate_generator_to_rates(...)` adjoint accumulation (cf.
  `backpropagate_laplacian_to_conductance`), using `Qᵀ` solves.

(Exact reduced/absorbing index bookkeeping to be verified against the current
`reduced_index` / `demes` conventions during the build.)

---

## 4. The asymmetric-distance measurement-model question (question 3) — RESOLVED

Because we map the asymmetric generator to a **symmetric** commute-time distance
`R^dir` (§1.2), and the genetic response `S` is symmetric, **no measurement model
changes are required**. `leastsquares`, `mlpe`, `generalized_wishart`, and
`wishart_covariance` all consume a symmetric distance/covariance; the Gower
transform `Σ = −½ J R^dir J` applies as-is. Tier 1's `wishart_drift_covariates`
and Tier 2's hierarchical field compose with it too.

Why not a genuinely asymmetric distance? There is no symmetric→asymmetric
information in the data to fit it to (S is symmetric), so an asymmetric model
distance is unidentifiable by construction. Symmetrizing (commute time) is the
correct and only well-posed choice. The directional information is recovered
through `γ_dir` shaping the symmetric commute time, **not** through an asymmetric
response. (True asymmetric *coalescent* modelling is Tier 3b.)

---

## 5. The adjoint / gradient (transpose solves)

Hitting times `h` satisfy `Q h = −1`. For likelihood `ℓ` with sensitivity
`dℓ/dh` (from the measurement model via the commute-time chain rule):

> `dℓ/dγ_k = −[Q⁻ᵀ (dℓ/dh)]ᵀ (∂Q/∂γ_k) h`.

So: one **transpose solve** `Qᵀ a = dℓ/dh` per focal contrast, then accumulate
`dℓ/dγ_k = −aᵀ (∂Q/∂γ_k) h` over edges (sparse: `∂Q/∂γ_k` is nonzero only on
edges, and `∂G_{a→b}/∂γ = G_{a→b} · z_{ab}`). This is reverse-mode AD, identical
in structure to the current `backpropagate_laplacian_to_conductance` but against
`Qᵀ` instead of the (self-adjoint) `L`. `Eigen::SparseLU` supplies the `Qᵀ` solve
from the forward factorization, so the elegant single-factorization economy is
preserved (the self-adjointness that the symmetric case enjoyed is replaced by
the transpose solve, at no extra factorization cost).

Hessian: for Newton, either finite-difference the gradient (γ is low-dimensional)
or extend the per-parameter solve loop (as in the symmetric engine). Given `γ`
is small (a few covariates), BFGS on the gradient is a safe v1 (as in Tier 2).

**Reversible fast path (optional):** when `γ_dir = 0` the generator is reversible
and the existing SPD machinery applies; the engine can detect this and fall back
to the Cholesky path.

---

## 6. Identifiability — the central scientific risk (be honest about it)

Directionality must be inferred from a **symmetric** genetic summary. Asymmetric
migration does leave a distinctive signature in coalescence times (Lundgren &
Ralph 2019; Thomaz et al. 2019), so `γ_dir` is **not** trivially unidentifiable —
but it may be **weak**, especially:
- when the symmetric conductance term can mimic part of the pattern,
- at migration-drift equilibrium (the asymmetry signature can be subtle),
- with few focal sites.

This is the dominant risk and must gate the build (Phase 1 below). Possible
outcomes, all worth reporting honestly:
- `γ_dir` is identifiable with realistic data → Tier 3a delivers the feature.
- `γ_dir` is only weakly identifiable from symmetric distances alone → document
  the limitation and the auxiliary data that would help (directional covariates
  with strong priors; IBD-asymmetry or allele-frequency-cline signals — which
  point back to Directions 4/6/7 of the feasibility review).

Either way the result is a genuine contribution: a tractable directional model
*and* a clear statement of when symmetric genetic data can and cannot resolve
direction.

---

## 7. SLiM validation design (question 4)

Built on the now-complete recapitation pipeline (`dev/slim/`, task #7).

**Scenario 3 — asymmetric (downstream-biased) migration.**
- `scen3_ts.slim`: WF stepping-stone (or continuous) on a lattice with an
  **elevation/flow gradient**; migration asymmetric along it
  (`m(a→b) = m0·exp(+β_dir·drop)`, `m(b→a) = m0·exp(−β_dir·drop)`), set via
  `setMigrationRates` with unequal rates. Tree-seq output.
- Recapitate + mutations + sample (existing `recap_sample.py`); MAF-filter; build
  the covariance.
- **Fits to compare:**
  1. *Symmetric* terradish (current) — expected to misfit / mislocate the barrier
     under asymmetry (reproducing Lundgren & Ralph "circuit theory misleads",
     in miniature).
  2. *Directed* Tier 3 with the elevation **directional covariate** — expected to
     recover `γ_dir` with the correct sign and improve fit (AIC / cross-validated
     prediction) over the symmetric model.
- **Identifiability probe (the key test):** sweep `β_dir` from 0 (symmetric) to
  strong; measure at what asymmetry strength `γ_dir` becomes reliably recovered
  (sign + magnitude) and the directed model is AIC-preferred. Report the
  detectability threshold honestly.
- **Ground truth check:** the coalescent pipeline can also report the realized
  asymmetry of lineage movement to confirm the planted signal.

**Validation ladder (as for Tiers 1–2):** (1) numDeriv on the `γ` gradient
(transpose-solve adjoint) on a small lattice; (2) self-consistency recovery
(simulate from the directed model, refit); (3) SLiM Scenario 3 as above.

---

## 8. Phased implementation plan

Mirrors the de-risking approach that worked for Tiers 1–2 (prototype → harden).

- **Phase 0 — R prototype + identifiability gate (do first, low cost).** Pure-R
  directed generator on a small lattice; hitting times via `Matrix` sparse LU;
  commute time; fit with `generalized_wishart`; numDeriv-check the `γ` gradient;
  run the §7 identifiability probe on *simulated-from-model* data. **Decision
  gate:** proceed only if `γ_dir` is recoverable at plausible asymmetry strengths.
- **Phase 1 — directed generator model + graph helpers.** `directed generator`
  factory (`log G = s·θ + d·γ_dir`) and `edge_gradient()` raster helper.
- **Phase 2 — C++ engine.** `Eigen::SparseLU` assemble/solve/transpose-adjoint
  (`assemble_directed_subgenerator`, `backpropagate_generator_to_rates`);
  `terradish_directed_algorithm()`; reversible fast-path fallback.
- **Phase 3 — top-level API + optimizer wiring.** `terradish_directed()` (BFGS on
  the `γ` gradient; reuse measurement models + nuisance subproblem); S3 methods;
  directional-effect summary.
- **Phase 4 — validation.** numDeriv + recovery + SLiM Scenario 3 (+ the
  identifiability sweep).
- **Phase 5 — docs.** Vignette (directional gene flow; when it is and isn't
  identifiable), man pages, roadmap update.

---

## 9. Decision points for review (before building)

1. **Distance object.** Directed **commute time** (recommended: symmetric,
   tractable, reduces to resistance distance) vs. pursuing the true structured
   coalescent (Tier 3b; intractable on fine rasters — deme-scale only).
2. **Scope.** Build **Tier 3a** (directed commute time) now; document Tier 3b
   (coalescent) as future. Agree?
3. **Solver.** `Eigen::SparseLU` (recommended; already a dependency; gives the
   transpose solve) vs. an iterative backend now. Recommend SparseLU for v1.
4. **Architecture.** Separate `terradish_directed_algorithm()` engine
   (recommended; isolates non-symmetric algebra) vs. integrating into
   `terradish_algorithm()`.
5. **Gate on the Phase-0 identifiability prototype** before the C++ build?
   Strongly recommended — it is cheap and de-risks the central scientific
   assumption.

---

## 10. Effort & risk summary

| Component | Effort | Risk |
|---|---|---|
| Phase 0 R prototype + identifiability gate | low | **high-value: tests the core assumption** |
| Directed generator model + edge-gradient helper | low–med | low |
| C++ `Eigen::SparseLU` engine + transpose adjoint | **high** | med (new linear-algebra path, but isolated) |
| API + optimizer wiring | med | low (reuses Tier 1–2 patterns) |
| SLiM Scenario 3 + identifiability sweep | med | med (detectability is the open question) |
| Docs | low | low |

**Bottom line.** Tier 3a is feasible and well-scoped: the symmetric-output
realization (§0) keeps the measurement models untouched (question 3 resolved),
the covariate-driven `γ_dir` keeps it identifiable-in-principle and parsimonious
(question 2), `Eigen::SparseLU` supplies the non-symmetric forward+transpose
solves with no new dependency (question 1), and the recapitation pipeline gives a
ready asymmetric-migration validation (question 4). The dominant risk is
*scientific, not engineering* — whether direction is recoverable from symmetric
genetic distances — so the plan gates the invasive C++ build behind a cheap R
identifiability prototype (Phase 0).
