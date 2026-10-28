# Wishart warm-start nonconvergence fixture

Captured October 28, 2026 from the IBE/IBR vignette's positive statistical
benchmark, at the zero conductance-coefficient start. The response is a
Gaussian/Wishart covariance example, not additional simulated genetics.

The fixture stores a 24-site graph covariance, response covariance, paired
hydroperiod differences, nuisance start and bounds, and the existing strict
nuisance control. No graph solve or SLiM run is needed to reproduce the failure.
The initial nuisance run returned nonconvergence; the same model with its
default nuisance start converged. The regression verifies recovery without
relaxing tolerances, suppressing public warnings, or treating a failed retry
as converged. Original and recovered parameters are retained as provenance.

Captured by dev/vignette_consistency_20261028/diagnose_ibe_capture.R, using the
combined development package and the prototype of the targeted retry fix.

