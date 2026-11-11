# Deprecated radish names

terradish grew out of the radish package. The original function and
class names still work so that older scripts and saved fits keep
running, but they are deprecated. New code should use the `terradish`
names in the table below.

## Details

|  |  |
|----|----|
| **Legacy name** | **Use instead** |
| [`radish()`](https://wpeterman.github.io/terradish/reference/legacy_radish_fit_wrapper.md) | [`terradish()`](https://wpeterman.github.io/terradish/reference/terradish.md) |
| [`radish_grid()`](https://wpeterman.github.io/terradish/reference/legacy_radish_grid_wrapper.md) | [`terradish_grid()`](https://wpeterman.github.io/terradish/reference/terradish_grid.md) |
| [`radish_distance()`](https://wpeterman.github.io/terradish/reference/legacy_radish_distance_wrapper.md) | [`terradish_distance()`](https://wpeterman.github.io/terradish/reference/terradish_distance.md) |
| [`radish_algorithm()`](https://wpeterman.github.io/terradish/reference/legacy_radish_algorithm_wrapper.md) | [`terradish_algorithm()`](https://wpeterman.github.io/terradish/reference/terradish_algorithm.md) |
| [`terra_radish_scale_optim()`](https://wpeterman.github.io/terradish/reference/legacy_terra_radish_scale_optim.md) | [`terradish_scale_optim()`](https://wpeterman.github.io/terradish/reference/terradish_scale_optim.md) |
| methods for class `"radish"` ([`print()`](https://rdrr.io/r/base/print.html), [`summary()`](https://rspatial.github.io/terra/reference/summary.html), [`coef()`](https://rdrr.io/r/stats/coef.html), [`fitted()`](https://rdrr.io/r/stats/fitted.values.html), [`residuals()`](https://rdrr.io/r/stats/residuals.html), [`simulate()`](https://rdrr.io/r/stats/simulate.html), [`anova()`](https://rdrr.io/r/stats/anova.html), [`logLik()`](https://rdrr.io/r/stats/logLik.html), [`AIC()`](https://rdrr.io/r/stats/AIC.html), [`plot()`](https://rspatial.github.io/terra/reference/plot.html)) | the same methods for class `"terradish"`; see [`terradish_methods`](https://wpeterman.github.io/terradish/reference/terradish_methods.md) and [`plot.terradish`](https://wpeterman.github.io/terradish/reference/plot.terradish.md) |
| [`conductance()`](https://wpeterman.github.io/terradish/reference/conductance.md) for class `"radish_graph"` | [`conductance()`](https://wpeterman.github.io/terradish/reference/conductance.md) for class `"terradish_graph"` |
| class `"radish_conductance_model_factory"` | [`terradish_conductance_model_factory`](https://wpeterman.github.io/terradish/reference/terradish_conductance_model_factory.md) |
| class `"radish_measurement_model"` | [`terradish_measurement_model`](https://wpeterman.github.io/terradish/reference/terradish_measurement_model.md) |

The wrappers
[`radish()`](https://wpeterman.github.io/terradish/reference/legacy_radish_fit_wrapper.md),
[`radish_grid()`](https://wpeterman.github.io/terradish/reference/legacy_radish_grid_wrapper.md),
[`radish_distance()`](https://wpeterman.github.io/terradish/reference/legacy_radish_distance_wrapper.md),
and
[`radish_algorithm()`](https://wpeterman.github.io/terradish/reference/legacy_radish_algorithm_wrapper.md)
pass their arguments unchanged to the replacement and signal a
deprecation warning.
[`terra_radish_scale_optim()`](https://wpeterman.github.io/terradish/reference/legacy_terra_radish_scale_optim.md)
is an alias with the same arguments as
[`terradish_scale_optim`](https://wpeterman.github.io/terradish/reference/terradish_scale_optim.md).
Objects created by terradish carry both the new and the legacy class
names, so either set of methods dispatches on them.

Arguments, defaults, and returned values are documented on the
replacement pages. Full compatibility with every historical radish entry
point is not promised.

## See also

[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md),
[`terradish-package`](https://wpeterman.github.io/terradish/reference/terradish-package.md)
