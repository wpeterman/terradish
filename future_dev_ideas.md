future_dev_ideas
================

Updated: 2026-04-23

Saved response: anisotropy and directional movement in terradish
----------------------------------------------------------------

The short answer is yes for anisotropy, maybe for true directionality, but those are not the same problem.

Right now terradish is built around an undirected, symmetric graph. The graph is created from undirected cell adjacency in R/radish_graph.R, the reduced operator is forced symmetric in R/radish_algorithm.R, and the C++ core uses edge weights of the form w_ij = c_i + c_j in src/radish.cpp. On top of that, the main distance models collapse the latent covariance to symmetric pairwise distances using only the lower triangle in R/leastsquares.R and R/mlpe.R, and the Wishart-style models explicitly symmetrize inputs in R/generalized_wishart.R and R/wishart_covariance.R. So the current architecture is much closer to anisotropy than to true asymmetry.

Ranked ideas, balancing scientific payoff and implementation feasibility:

1. Symmetric edge-based anisotropy
Potential: high. Feasibility: high.
This is the best near-term extension. Instead of one conductance per vertex, we would model conductance per edge using orientation, slope magnitude, river class, current alignment, or edge-specific covariates. That could capture movement is easier along this axis than across it while keeping w_ij = w_ji, so the Laplacian stays symmetric and most of the solver stack survives.

2. Anisotropic diffusion or tensor conductance on rasters
Potential: high. Feasibility: moderate.
This is a more formal version of edge anisotropy: let movement depend on local direction, not just habitat value. For example, east-west and north-south movement could have different weights, or movement aligned with a river corridor could be favored. Still reversible and still compatible with SPD solvers, but it needs a more thoughtful edge-design-matrix layer.

3. Direction-aware but symmetrized metrics
Potential: moderate. Feasibility: moderate.
We could compute forward and backward movement costs separately, then feed a symmetrized summary into existing measurement models, like average, min, max, or geometric mean of d(i->j) and d(j->i). This would let us inject directional ecology without immediately rewriting the measurement layer, but it throws away some of the directional information.

4. Ordered-pair measurement models for directional data
Potential: high when data support it. Feasibility: moderate.
If the response were directional, such as assignment, parentage, mark-recapture, migrant counts, or telemetry-derived movement, terradish could support ordered pairs. The selected-pair work already gives part of the scaffolding. But this does not solve directional inference from symmetric matrices like FST; it mainly becomes powerful when the data themselves are directional.

5. Directed random-walk or CTMC backend
Potential: very high. Feasibility: low to moderate.
This is the real route to asymmetric movement. Instead of an undirected Laplacian, we would estimate a directed generator with q_ij != q_ji, for example downstream bias, current-driven transport, or one-way barriers. That is the cleanest scientific answer, but it is a genuine backend rewrite: no more Cholesky or CG assumptions, no longer an SPD matrix, and the latent summary would no longer be the same symmetric resistance distance.

6. BGR-style directed network covariance
Potential: very high for rivers. Feasibility: low.
White et al. 2020 is the clearest conceptual bridge here: keep the network framing, but move toward a bidirectional geneflow model with network covariance and a generalized Wishart distance likelihood. terradish already has graph-based covariance ideas and generalized Wishart machinery, so there is some architectural overlap. But because the current core is symmetric, this is still a substantial methodological extension, not a light feature add.

7. Advection-diffusion PDE model on the raster
Potential: high. Feasibility: low.
This would be the most general spatial formulation for currents, wind, or flow fields: symmetric diffusion plus directional advection. It is elegant, but much farther from the present package than the options above.

Main caution:
If the observed response is a symmetric pairwise genetic distance matrix, then true one-way movement is only weakly identifiable unless the latent process itself becomes directional. Just adding flow direction as another ordinary raster covariate under the current symmetric solver will not really estimate asymmetric movement; it will mostly estimate spatial heterogeneity correlated with direction.

Practical priority order:

1. Edge-based symmetric anisotropy first.
2. Then a directed backend for true asymmetry.
3. Then a BGR-style river-network implementation if we want to push directly into the White et al. space.


Design memo: anisotropy and directional movement in terradish
-------------------------------------------------------------

Purpose
-------

This memo sketches a development path for adding anisotropy and, eventually,
true asymmetric movement to terradish while preserving as much of the current
optimization and solver machinery as possible.


1. Current architectural constraints
------------------------------------

Current terradish assumptions that matter:

- Graph construction is undirected.
- The reduced graph operator is symmetric positive definite after grounding.
- Edge weights are assembled from vertex conductance as w_ij = c_i + c_j.
- The default latent object is a symmetric covariance or resistance surface.
- The main measurement models assume symmetric pairwise responses.

Practical implication:

- Symmetric anisotropy can fit within the current architecture.
- True directionality requires a new latent process, not just a new formula.


2. Recommended roadmap
----------------------

Phase A. Symmetric anisotropy on the existing solver stack
Priority: highest

Goal:
Allow conductance to vary with edge orientation and edge-specific covariates
while preserving a symmetric Laplacian.

What changes:

- Introduce an edge-level design matrix.
- Support edge covariates such as:
  - orientation class: horizontal, vertical, diagonal
  - continuous bearing features: dx, dy, abs(dx), abs(dy), cos(theta), sin(theta)
  - slope magnitude
  - alignment with a supplied flow or corridor vector field
  - edge averages, edge differences, or edge minima/maxima of raster values
- Define edge conductance directly:
  log w_ij = z_ij^T beta
  with w_ij = w_ji guaranteed by construction.

Why this is attractive:

- Keeps the operator SPD.
- Reuses direct, AMG, and PCG solvers.
- Reuses current reverse-mode derivative structure with modest extension.
- Works for rivers, coastlines, ridges, and movement corridors.

Expected user benefit:

- Can test whether movement is easier along versus across a feature.
- Can quantify directional preference without needing a fully directed model.


Phase B. Semi-directional summaries built on directed costs
Priority: medium

Goal:
Explore directional ecology before full directed inference by computing
directional movement costs and then reducing them to a symmetric pairwise
quantity for current measurement models.

Examples of symmetric summaries:

- mean(d_ij_forward, d_ji_backward)
- geometric mean
- min or max
- weighted combination tuned to biology

Why this is useful:

- Faster scientific prototyping.
- Lets us compare what information is lost by symmetrization.

Main downside:

- Not a true directional model.
- Identifiability remains limited for symmetric genetic distances.


Phase C. True directed random-walk backend
Priority: high scientifically, lower immediately

Goal:
Replace the undirected resistance kernel with a directed generator for movement.

Concept:

- Build a continuous-time Markov chain or directed graph generator Q where
  q_ij != q_ji.
- Parameterize directional bias with edge-level predictors.
- Derive directional movement summaries from the generator, for example:
  - mean first passage times
  - commute-time-like directed analogs
  - stationary covariance or related SAR-like structures

What breaks relative to current terradish:

- No symmetric Laplacian assumption.
- No Cholesky-based direct solve.
- No SPD AMG/PCG assumptions in the current implementation.
- Existing measurement models based on symmetric distances need adaptation or
  replacement.

What likely replaces them:

- Sparse LU or iterative nonsymmetric solvers.
- Ordered-pair likelihoods when the response is directional.
- New covariance construction for Wishart-like or SAR-like models.


3. Minimum viable anisotropic terradish
---------------------------------------

This is the most realistic first implementation.

Core idea:
Keep `conductance_surface()` for node and graph creation, but add an optional
edge-covariate layer that is derived once and then reused during optimization.

Candidate API ideas:

- `conductance_surface(..., edge_features = c("orientation", "slope"))`
- `edge_covariates(surface, features = ...)`
- `anisotropic_loglinear_conductance(formula, x_edge)`
- `edge_loglinear_conductance(formula, x_edge)`

Potential formula style:

- `~ forestcover_mean + altitude_mean + orient_diag`
- `~ slope_abs + flow_align`
- `~ forestcover_mean * flow_align`

Internal representation:

- Store edge-level data frame `data$edge_x`
- Store `edge_pairs`
- Conductance model returns one value per edge rather than one value per vertex

Key implementation choice:

- Do not try to shoehorn anisotropy into the existing vertex-only conductance
  parameterization if it makes derivatives messy.
- It is cleaner to support a second conductance-model family for edges.


4. Concrete implementation options for Phase A
----------------------------------------------

Option A1. Orientation-only anisotropy

Definition:

- Separate coefficients for horizontal, vertical, and diagonal edges.
- Can also interact these with raster summaries.

Example:

- `log w_ij = beta0 + beta_h * I(horizontal) + beta_v * I(vertical) + beta_d * I(diagonal)`

Pros:

- Very easy to implement.
- Good proof of concept.
- Tests whether grid-direction anisotropy matters.

Cons:

- Coarse and grid-dependent.


Option A2. Edge summaries from raster covariates

Definition:

- Build edge covariates from raster values at the two incident cells.

Examples:

- mean: (x_i + x_j) / 2
- absolute difference: |x_i - x_j|
- gradient-aligned change: signed change along edge direction
- minimum, maximum, product

Pros:

- Natural extension of current raster workflow.
- Useful for slope, elevation, moisture, flow accumulation.

Cons:

- Need to decide which summaries are biologically meaningful.


Option A3. Alignment with a vector field

Definition:

- User supplies a two-layer raster or edge-level vector field, such as local
  current direction, prevailing wind, or downstream direction.
- For each edge, compute alignment between edge direction and the field.

Example covariates:

- dot product with unit flow vector
- absolute alignment
- signed alignment if later extended to directional models

Pros:

- Strong fit for river and marine systems.
- Sets up a natural bridge to directed movement later.

Cons:

- Requires careful unit and projection handling.


5. How derivatives could still work in Phase A
----------------------------------------------

Current reverse-mode differentiation already propagates:

- measurement model gradient with respect to E
- E with respect to the reduced Laplacian
- Laplacian with respect to conductance
- conductance with respect to theta

For symmetric edge conductance, the extension is conceptually straightforward:

- Replace vertex conductance vector c with edge conductance vector w
- Build dL/dw directly rather than through vertex aggregation
- Use edge design matrix Z so dw/dtheta is simple under log-link models

This should be substantially easier than a fully directed backend.


6. True directionality: what a new backend would need
-----------------------------------------------------

Directed model form:

- `log q_ij = z_ij^T beta + a_ij^T delta`
- where z_ij are symmetric habitat terms and a_ij are directional terms

Directional terms could include:

- downstream indicator
- signed slope
- alignment with current
- barrier passability in one direction only

Latent summaries to consider:

- directional hitting times
- pairwise migration rates
- stationary covariance from a directed SAR-like model
- source-sink asymmetry summaries

Likely package architecture:

- Keep current terradish as the symmetric engine
- Add a parallel experimental family rather than silently overloading existing
  functions

Candidate naming:

- `terradish_directional()`
- `directed_conductance_surface()`
- `directional_measurement_model()`

This separation would reduce scientific ambiguity and avoid hidden changes to
established behavior.


7. Measurement-model implications
---------------------------------

If the observed response is symmetric:

- least squares and MLPE can still be used for symmetric summaries only
- generalized Wishart and Wishart covariance would need careful theoretical
  work if the latent covariance is no longer symmetric in the same way

If the observed response is directional:

- ordered-pair models become much more appropriate
- recent selected-pair infrastructure could be extended to ordered pairs
- a Poisson, multinomial, or binomial movement model may be more natural than a
  distance-regression model

Bottom line:

- symmetric genetic distances support anisotropy well
- directional responses are much better for true asymmetry


8. River-network-specific path
------------------------------

This is where the White et al. 2020 BGR concept is most relevant.

Potential river-focused development path:

1. Add a network-style graph constructor or allow a user-supplied edge list.
2. Support edge covariates measured on stream reaches rather than raster cells.
3. Implement symmetric anisotropic stream movement first:
   - reach slope
   - stream order
   - confluence effects
   - thermal gradient
4. Explore a directed generator with explicit downstream bias.
5. Revisit generalized Wishart or SAR-style covariance once the directed latent
   process is defined.

This is likely cleaner than trying to force rivers fully into a raster-only
representation.


9. Validation plan
------------------

Phase A validation:

- isotropic limit reproduces current terradish results
- numerical derivative checks for new edge-conductance models
- simulation recovery for known anisotropy coefficients
- graph SPD checks under all supported formulas
- benchmark runtime against current isotropic model

Phase C validation:

- small directed toy graphs with known hitting times
- agreement with analytical results on simple chains
- sensitivity to one-way barriers and directional bias
- simulation study for identifiability under symmetric versus directional data


10. Recommended immediate next steps
------------------------------------

Recommended order of work:

1. Prototype edge-level anisotropy on the current raster graph.
2. Start with orientation-only and slope-based edge covariates.
3. Build one new conductance-model factory for edge conductance.
4. Reuse current measurement models unchanged at first.
5. Add simulation tools to verify coefficient recovery.
6. Only after that, scope a directed backend.

If a first concrete development ticket is needed, it should be:

"Add an experimental edge-based anisotropic conductance family that keeps the
graph symmetric and works with existing leastsquares, mlpe, and Wishart-style
measurement models."


11. Open questions worth revisiting later
-----------------------------------------

- Should anisotropy live inside `conductance_surface()` or in a separate
  `edge_covariates()` helper?
- Should edge conductance models coexist with vertex conductance models or
  replace them internally?
- How should users supply vector fields for alignment-based anisotropy?
- Should river-network support be raster-first, network-first, or both?
- Which directional summaries are defensible for symmetric genetic responses?
- Can generalized Wishart be extended cleanly to a directed latent process, or
  should a different likelihood family be used?


Saved response: FFT landscape metrics and Gaussian scale framework
------------------------------------------------------------------

Reviewed latest multiScaleR on GitHub:

- Repository: wpeterman/multiScaleR
- Version inspected: 0.6.18
- Commit: 4ba958d, 2026-05-03
- Commit message: Reconcile landscape metric merge fixes

The relevant multiScaleR additions are:

- FFT convolution for raster smoothing.
- Explicit covariate specifications with `msr_vars()`, `kernel_var()`, and
  `landscape_var()`.
- FFT projection of categorical landscape metrics, including composition,
  edge, and adjacency metrics.

The short answer for terradish is:

Yes, this is feasible, but the best integration depends on whether the goal is
fixed-radius landscape metrics or true Gaussian scale optimization.


1. Near-term feasible path: FFT metric rasters as preprocessing
---------------------------------------------------------------

The cleanest first integration is to let multiScaleR compute fixed-radius
landscape metric rasters, then pass those numeric rasters into standard
terradish conductance models.

Workflow:

1. Define metric covariates with `multiScaleR::msr_vars()` and
   `multiScaleR::landscape_var()`.
2. Project them across the raster with `kernel_scale.raster(..., fft = TRUE)`.
3. Use the projected metric rasters in `conductance_surface()`.
4. Fit with `loglinear_conductance()` or `linear_conductance()`.

Why this works well:

- By the time terradish sees the data, each metric is a numeric raster layer.
- No change is needed to terradish's inner gradient, Hessian, or solver code.
- Existing formula handling, model comparison, and measurement models still
  apply.
- This approach supports classic landscape metrics whose definitions are based
  on a hard circular radius.

Main limitation:

- The radius is fixed outside the terradish optimizer unless handled by an outer
  optimization such as `terradish_scale_optim()`.


2. Caution: hard-radius metrics are not smooth in sigma
------------------------------------------------------

Classic landscape metrics computed inside a circle change discontinuously when
the radius changes enough for cells or edge pairs to enter or leave the buffer.

This is a poor match for the current `gaussian_smoothed_loglinear_conductance()`
framework because that framework assumes:

- smoothed covariates vary smoothly with sigma,
- first derivatives with respect to sigma exist and are stable,
- second derivatives can be propagated through the conductance model.

Therefore, directly treating a hard-radius landscape metric radius as an inner
Newton or BFGS parameter is mathematically awkward. It may still be useful in
an outer grid or coordinate search, but it should not be presented as the same
thing as Gaussian smoothing.


3. Better long-term path: Gaussian-weighted landscape metrics
------------------------------------------------------------

For joint scale optimization inside terradish, a better approach is to define
Gaussian-weighted analogs of selected landscape metrics.

This would keep the spirit of landscape metrics while matching the existing
Gaussian kernel framework:

- Use Gaussian weights rather than a hard circular inclusion radius.
- Use FFT convolution to evaluate local weighted counts or proportions across
  the whole raster.
- Derive first and second sigma derivatives from the Gaussian kernel, as the
  current Gaussian conductance code already does for continuous covariates.

This is conceptually aligned with the existing implementation in
`R/gaussian_scale_conductance.R`, which already:

- FFT-convolves raster values and masks,
- normalizes by local valid-cell kernel weight,
- computes first and second derivatives of the Gaussian kernel with respect to
  sigma,
- propagates those derivatives through log-linear conductance formulas.


4. Metric feasibility ranking
-----------------------------

High feasibility:

- Gaussian-weighted class proportions.
- Gaussian-weighted binary habitat proportion.
- Gaussian-weighted Shannon diversity style metrics.
- Gaussian-weighted Simpson diversity style metrics.

These can be built from FFT-convolved class-indicator rasters. The resulting
class probabilities are smooth functions of sigma, so derivatives are feasible.

Moderate feasibility:

- Edge density.
- Total edge.
- PLADJ.
- Contagion.

These can be represented using edge-pair indicator rasters and FFT convolution,
as multiScaleR now does for projection. The challenge is that derivative
bookkeeping grows with the number of classes and ordered or unordered class
pairs.

Low feasibility for inner Gaussian optimization:

- Patch richness.
- Relative patch richness.
- Landscape shape index.
- Connected-component or patch-identity metrics.
- Metrics based on class presence/absence thresholds.
- Metrics involving minimum perimeter or discrete topology.

These contain discontinuities or discrete operations that do not naturally
support smooth derivatives with respect to sigma.


5. Suggested terradish development design
-----------------------------------------

Recommended staged design:

Phase A. Fixed-radius metric raster support

- Document and test the workflow using multiScaleR-generated metric rasters.
- Possibly add a small terradish helper or vignette showing:
  - derive FFT metric rasters,
  - scale them with `scale_covariates()`,
  - fit them in terradish,
  - compare against ordinary raster covariates.
- Keep this outside `gaussian_smoothed_loglinear_conductance()`.

Phase B. Outer optimization for hard-radius metrics

- Use `terradish_scale_optim()` with a custom `scale_fun` that returns one
  landscape metric raster for a candidate radius.
- Treat the radius search as a nonsmooth outer problem.
- Prefer grid search or robust coordinate search over derivative-based inner
  optimization.

Phase C. Native Gaussian-weighted metric family

- Add an experimental Gaussian landscape metric smoother.
- Start with class proportions and diversity metrics.
- Reuse the existing FFT and derivative machinery from
  `gaussian_smoothed_loglinear_conductance()`.
- Only support metrics whose values and derivatives are well-defined.


6. Candidate API ideas
----------------------

For fixed-radius preprocessing:

- `landscape_metric_rasters(covariates, specs, radius, engine = "multiScaleR")`
- or simply a vignette using multiScaleR directly.

For native Gaussian metrics:

- `gaussian_landscape_vars()`
- `gaussian_landscape_metric_conductance()`
- `gaussian_smoothed_loglinear_conductance(..., metric_vars = NULL)`

Possible variable specification:

- `gaussian_metric_var("landcover", metric = "shdi")`
- `gaussian_metric_var("forest", metric = "proportion", class = 1)`

Important design rule:

- Keep hard-radius landscape metrics and Gaussian-weighted landscape metrics
  distinct in names, documentation, and interpretation.


7. Practical recommendation
---------------------------

First implementation target:

Use multiScaleR's FFT landscape metric projection as a preprocessing workflow
for terradish.

First native Gaussian target:

Gaussian-weighted class proportions and diversity metrics.

Avoid initially:

- Optimizing classic hard-radius landscape metric radii inside
  `gaussian_smoothed_loglinear_conductance()`.
- Supporting every landscapemetrics-style metric at once.
- Including patch/topology metrics in derivative-based sigma optimization.

Bottom line:

FFT metric projection is immediately useful for terradish as a preprocessing
tool. Native joint sigma optimization is also promising, but it should be built
around Gaussian-weighted metric analogs rather than hard-radius landscape
metrics.
