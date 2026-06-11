# Scenario 4: DRAGON rung-3 validation (asymmetric migration + density gradient)

Extends the scen3 recapitation pipeline to plant **both** a directional gene-flow
bias and a **density / Ne gradient correlated with the same elevation axis**, the
worst-case confound the DRAGON robustness sweep identified. It tests, on real
coalescent data, that the directional coefficient is recovered only when the drift
surface is modeled jointly.

## Files

- `scen4_ts.slim` - SLiM 5 WF stepping stone. Asymmetric (downhill) migration
  `MIG*exp(BETA_DIR*(elev_k-elev_i))` plus deme sizes `N_i = N0*exp(KAPPA_N*elev_i)`.
  `BETA_DIR=0` -> symmetric; `KAPPA_N=0` -> uniform sizes (i.e. scen3).
- `recap_sample.py` - unchanged; recapitate + mutate + sample per deme.
- `dragon_engine.py` - self-contained DRAGON engine (covariate directed
  structured-coalescent forward map, analytic adjoint gradient, Wishart fit).
  Depends only on numpy/scipy.
- `run_scen4_recap.py` - orchestrator: SLiM -> recap -> three DRAGON fits
  (symmetric, direction-only, direction+drift) -> verdict.

## Run (from the terradish repo root, conda `slim` Python)

```
& "C:/Users/peterman.73/AppData/Local/anaconda3/envs/slim/python.exe" dev/slim/run_scen4_recap.py
```

Environment overrides: `SLIM`, `SLIM_FORCE=TRUE`, `DIM`, `N0`, `MIG`, `BETA_DIR`,
`KAPPA_N`, `BURNIN`, `NSAMP`, `SEED`. Defaults: DIM=6 (36 demes), BETA_DIR=0.8,
KAPPA_N=0.6.

## Expected result (the confound test)

1. **Direction detected:** `LRT direction vs sym > 3.84`.
2. **Drift detected:** `LRT drift (joint vs dir) > 3.84` (because scen4 plants a
   size gradient; on scen3, which is uniform, this LRT is ~2.8, n.s.).
3. **Direction-only is biased:** `|gamma_dir|` from the direction-only fit is
   inflated relative to the joint fit, because the omitted density gradient runs
   along the directional covariate. The joint fit is the trustworthy one.

`gamma_dir` and `kappa` are in standardized-elevation units (the fit's covariate),
so compare signs and the direction-only-vs-joint shift, not the raw magnitude to
`BETA_DIR`/`KAPPA_N`.

## Validation status

The fit harness was smoke-tested end-to-end on the existing scen3 recap output:
it reproduces the de-risk exactly (symmetric NLL -4595.1, direction -4640.4,
`LRT=90.7`, `gamma_dir=-0.65`) and correctly finds drift **not** needed on
uniform-size scen3 (`LRT_joint_vs_dir=2.8`, `kappa=-0.07`). Running `scen4_ts.slim`
supplies the density gradient that should make the drift term significant and
expose the direction-only bias.
