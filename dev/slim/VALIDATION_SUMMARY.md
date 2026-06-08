# Tier 1 SLiM validation summary (drift / effective-size surface)

SLiM 5.2 (`C:\msys64\mingw64\bin\slim.exe`). No Python/tree-sequence stack
available, so these runs are **forward-only** with direct allele-frequency
output and finite sampling in R (`cov_from_biallelic`). See the recapitation
note at the bottom.

## Three-rung validation status

- **Rung 1 (analytic, numDeriv):** PASS — gradient/Hessian/`gradient_E` to 1e-9
  (covariance and generalized models); full `terradish_algorithm` θ-gradient vs
  finite difference = 2.9e-7; intercept-only / γ=0 reduce to `wishart_covariance`
  to machine epsilon. (`dev/check_tier1.R`, `dev/check_e2e.R`.)
- **Rung 2 (self-consistency recovery):** PASS — simulating from the model and
  refitting recovers γ (mean 0.688 vs true 0.70 over 5 draws; vignette demo
  γ̂ = −0.598 vs true −0.6) and AIC prefers the drift model 5/5.
- **Rung 3 (SLiM, biological signal):** see below.

## Scenario 1 — discrete stepping-stone, varying deme size  → **PASS**

5×5 WF lattice, deme size `250·exp(0.7·cov)` (≈7× Nₑ contrast), uniform
neighbour migration; `dev/slim/scen1_stepping_stone.slim` + `run_scen1.R`.

| quantity | value | expectation | verdict |
|---|---|---|---|
| `cor(diag(S), cov)` | **−0.852** | negative (more Nₑ → less drift) | ✓ strong |
| `cor(diag(S), deme size)` | −0.784 | negative | ✓ |
| `gamma_cov` (drift slope) | **−0.147** | negative | ✓ correct sign |
| AIC drift vs scalar | 27908 < 27961 | drift preferred | ✓ |

The drift surface correctly detects the per-deme effective-size gradient and is
selected over a scalar nugget. Two honest caveats, both from the forward-only
simulation (not the estimator):
- `tau = 0`: the off-diagonal IBR structure is not captured — forward-only
  variation is dominated by young, spatially-clustered mutations whose
  covariance does not match a smooth resistance distance. So this run validates
  the **diagonal (drift)** part; it is not a joint IBR+drift test.
- The recovered slope magnitude is muted (observed `diag(S)` contrast ≈1.8× vs
  the true ≈7× Nₑ contrast) because the young-mutation artifact inflates the
  baseline diagonal and compresses the relative gradient.

## Scenario 2 — continuous space, habitat-driven density  → **inconclusive (needs recapitation)**

nonWF xy, `localPopulationDensity` competition, natal Gaussian dispersal,
habitat gradient drives local density → local Nₑ; interior focal grid;
`dev/slim/scen2_continuous_space.slim` + `run_scen2.R`.

Even after fixing edge effects (interior focal points) and enlarging
neighbourhoods (σ_d=0.12, K=2800): `diag(S)` ~ 180–350, `cor(diag(S), habitat)`
≈ 0.07 (no clean signal), `tau = 0`. The density→drift signal is **swamped** by
the forward-only young-mutation / high-Fst artifact (each new variant arises at
a point and spreads slowly via σ_d, producing pathological structure with no
deep coalescent history). A small-sample ascertainment confound
(`cor(diag(S), nsamp) ≈ 0.41`) is also present.

This is a **simulation-design limitation**, not an estimator failure: the
estimator is proven by rungs 1–2 and scenario 1. A clean continuous-space test
requires deep ancestry, i.e. **recapitation**.

## Recapitation note (the fix, when Python is available)

The forward-only burn-in is both slow and biased toward young variants. With a
conda-forge env (`python=3.12` + `msprime tskit pyslim numpy scikit-allel`),
switch both scenarios to: tree-seq recording (mutation rate 0) →
`pyslim.recapitate(ancestral_Ne, recombination_rate)` →
`msprime.sim_mutations(rate, keep=True)` → sample at focal coords → covariance.
This gives equilibrium ancestry, ~100× faster runs, a clean low/moderate-Fst
covariance (`diag(S)` ~ O(1)), `tau > 0` (so the joint IBR+drift model is
testable), and a sharp density→drift gradient in scenario 2. It is also the
prerequisite for Tier 3 (coalescent-based asymmetric-migration validation).

## Recapitation pipeline (added once conda Python was available)

Pipeline: SLiM tree-seq (`scen1_ts.slim`, `scen2_ts.slim`) → `recap_sample.py`
(`pyslim.recapitate` → `msprime.sim_mutations` → per-subpop or spatial-Voronoi
sampling) → R (`run_scen1_recap.R`, `run_scen2_recap.R`). Built and functional.

- **MAF filtering is essential** with recapitated data: the realistic SFS has
  many rare variants, and `cov_from_biallelic`'s `1/sqrt(p(1-p))` standardization
  upweights them (`diag(S)`~1000 unfiltered). MAF≥0.05 restores an O(1) covariance.
- **Scenario 1 (discrete), recapitated:** with low migration (0.005) + long
  burn-in (4000) + MAF≥0.05, `cor(diag(S),cov)=-0.50`, `gamma_cov=-0.095`
  (correct sign), AIC prefers drift. A third independent confirmation of the
  drift surface, now on realistic deep-coalescent data.
- **Scenario 2 (continuous), recapitated:** `cor(diag(S),habitat)=+0.33` — wrong
  sign, robustly across all designs (forward/edge-fixed/large-neighbourhood/
  recapitated+MAF). **Interpretation (a real finding, not a bug):** sampling a
  fixed N nearest individuals from denser (high-habitat) regions draws a
  spatially concentrated, more-related sample → larger covariance diagonal, which
  opposes and overwhelms the neighbourhood-size→drift effect (sample size also
  rises 40→60 with habitat). In continuous space the covariance diagonal tracks
  local relatedness / sampling scale, not a clean 1/Ne.
- **`tau=0` (IBR vs drift trade-off):** strengthening drift (low migration /
  large size contrast) makes demes nearly independent, so there is little IBR
  covariance to detect. A joint IBR+drift scenario with `tau>0` needs a separate
  conductance(migration) gradient distinct from the size(drift) gradient — a
  richer future scenario; not required to validate the new feature.

## Bottom line

The drift / effective-size surface is **validated for deme-structured data**:
exact on analytic and self-consistency rungs, and confirmed on the discrete-deme
biological scenario under BOTH forward-only (`cor=-0.85`) and realistic
recapitated (`cor=-0.50`) simulation. **Scope finding:** it is a discrete-deme
model and does not invert continuous-space density variation, where the genetic
covariance diagonal is confounded by local relatedness and sampling scale — the
continuous-space density problem belongs to simulation-based / mapNN-style
methods (feasibility review, Directions 2/6). This scope should be stated in the
user docs.
