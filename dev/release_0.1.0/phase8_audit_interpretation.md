# Final regression audit interpretation

The adapted scripts retain the original audit designs wherever those designs
apply to the core interface. Their source paths and SHA-256 hashes are in
`phase8_audits/manifest.json`. Outputs are saved separately from the original
audit. Successful script completion is not itself a statistical claim.

## Acceptance evidence

- Centered covariance, uncentered covariance, and the corresponding distance
  response give matching contrast-Wishart estimates in the five clustered
  simulations. The five-start optimizer experiment reaches the same optimum;
  Newton and BFGS spline fits agree. The near-zero MLPE correlation experiment
  recovers the least-squares limit.
  Its outer optimizer converges (projected gradient 2.08e-8), while the inner
  nuisance optimizer reports a stall near rho = 1.69e-9. That diagnostic is
  retained; matching the least-squares limit is not a claim that every inner
  optimization has a clean exit. The reported natural-scale rho is zero to
  displayed precision, with unavailable boundary uncertainty.
- Gaussian impulse checks recover the correct origin on all four even/odd
  grid combinations, with errors from 6e-15 to 2.13e-13. The external smoothing
  units experiment is retained as a units check, not evidence that a different
  kernel is numerically identical to the native Gaussian kernel.
- Spline prediction retains identical values on unchanged cells (maximum
  log-ratio zero, correlation one). Layer scaling independently spans zero to
  one for both layers.
- The synthetic crop example gives mean resistance increases of 22.9% and
  15.2% with a two-cell buffer and 1.5% and 0.8% with ten cells. These are
  case-specific examples. For melip, the mean increase is 1.4%, maximum 7.2%,
  and correlation 0.9986. Three duplicate-cell pairs have zero resistance and
  are excluded from relative changes; their presence is reported explicitly.
- Both Wishart formulations agree with independent contrast-Wishart log
  density differences within 1e-8. The MLPE objective agrees with the dense
  Gaussian calculation within 1e-8. The default kernel is the double-centered
  unscaled absolute difference, and reconstructs that pairwise difference.
- Non-nested measurement-model comparisons are rejected. The one-component
  boundary LRT uses the chi-bar-square approximation for positive statistics;
  a statistic of exactly zero returns p = 1. Wald intervals near a boundary
  remain an approximation and can extend below zero.
- Gaussian CV scores match manual scoring with internal-scale parameters
  within 1e-6. Passing external-scale parameters can change scores. A changed
  response is rejected when resuming a checkpoint.

## Measurement-extension limitations

The four-seed `ext/check5` run uses the new absolute-difference default.
Covariance and distance Wishart versions agree. The added environmental term
has lambda values 0.01290, 0.03640, 0.00848, and 0.00746; its AIC differences
are -0.636, -15.910, 0.830, and 1.083. Thus reformulating the kernel does not
make AIC a calibrated mechanism test. The term can absorb misspecification or
sampling variation. This result is retained rather than replaced by a claim
that the false association has disappeared.

The eight-replicate noise-covariate CV audit uses identical folds for the two
models in each replicate. Fixed nuisance scoring favors the extension in 3/8
Wishart and 2/8 MLPE replicates. Reprofiled scoring favors it in 8/8 for both.
With a 1e-6 absolute reporting tolerance, a raw positive
difference of 7.85e-13 in Wishart replicate seven is a numerical tie. The
unrounded scores and original greater-than-zero tally remain in the raw log.
This gives the same classification as the earlier audit's scaled machine-
precision tie rule; no other comparison is close to zero.
The reprofiled comparison is reconstructed from separate calls solely to
diagnose the old failure mode; the core API rejects comparing different
measurement models in one reprofiled call.

Fixed scoring has zero failed Wishart folds and one failed MLPE fold across
the model fits. Reprofiled scoring has one failed Wishart fold and nine failed
MLPE folds. Differences use only common successful folds. Failed baseline
scores are separately recorded and do not invalidate an otherwise valid
trained-model density. Eight replicates cannot establish a false-positive
rate, and these results do not support causal interpretation of an IBE term.

## Script corrections and provenance

The first CV attempt failed because its checkpoint-error assertion expected a
different phrase and because it tried the now-prohibited multi-model
reprofiled call. `phase8_cv_v2` is the corrected final receipt. The initial
extension diagnostics mistakenly assumed that scale=TRUE was still the
default; `phase8_corrections` uses the explicit current default and asserts
the kernel identity. The original logs remain available. No package code was
changed to make these audit adaptations pass.
