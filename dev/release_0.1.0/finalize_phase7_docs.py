from pathlib import Path
import re

p = Path('dev/check/phase7-package')
words = '''absdiff comparator estimand exponentiate gower layerwise natively
predictively PSD recalibration reoptimizing Reprofiling restandardize unclamped
Unclamped Unsampled'''.split()
f = p / 'inst/WORDLIST'
s = f.read_text(encoding='utf-8').rstrip() + '\n'
existing = set(s.splitlines())
f.write_text(s + ''.join(w + '\n' for w in words if w not in existing), encoding='utf-8')
f = p / 'README.md'
s = f.read_text(encoding='utf-8').replace("landgraph's coherent `gower` diagonal", 'the coherent `gower` diagonal from landgraph')
s = s.replace('projected_sites <- terra::project(sites, "EPSG:5070")', '''lonlat <- terra::crds(terra::project(sites, "EPSG:4326"))
local_crs <- sprintf("+proj=aeqd +lat_0=%f +lon_0=%f +datum=WGS84 +units=m",
                     mean(lonlat[, 2]), mean(lonlat[, 1]))
projected_sites <- terra::project(sites, local_crs)''')
s = s.replace('geographic = dist(terra::crds(projected_sites))', 'geographic_100km = dist(terra::crds(projected_sites)) / 100000')
s = s.replace('The IBE:IBR ratio expresses', 'Geographic differences in this example are in units of 100 km. The IBE:IBR ratio expresses')
f.write_text(s, encoding='utf-8')
# Keep the optional simulation driver but disconnect it from vignette results.
f = p / 'vignettes/precompute.R'
s = f.read_text(encoding='utf-8').replace('# Pre-computes slow simulation results for package vignettes.', '''# Optional exploratory power-study driver. The current vignettes do not load
# these artifacts; they contain self-contained examples and explicit recipes.
# Results are conditional on the generating model and supplied effective nu.
# Archived outputs from older versions are not validation of the current code.''')
s = s.replace('# Output files (committed to source):', '# Optional local outputs (not included in the built package):')
s = s.replace('maxit = 10', 'maxit = 100')
s = s.replace('    nu                = nu,', '    nu                = nu,\n    nu_fit            = nu,')
s = s.replace('  nu         = 100,', '  nu         = 100,\n  nu_fit     = 100,')
f.write_text(s, encoding='utf-8')
f = p / 'inst/benchmarks/compare-radish-terradish.R'
s = f.read_text(encoding='utf-8').replace('maxit = 20L', 'maxit = 100L')
s = s.replace('# avoid namespace and S3-method conflicts.', '''# avoid namespace and S3-method conflicts. Parallel derivatives are experimental;
# compare convergence and likelihoods before interpreting elapsed-time ratios.''')
s = s.replace('        df = fit$df,', '''        df = fit$df,
        convergence = if (is.list(fit$convergence)) fit$convergence$code else NA_real_,''')
s = s.replace('        df = x$df[[1]],', '''        df = x$df[[1]],
        convergence = x$convergence[[1]],''')
f.write_text(s, encoding='utf-8')
# Extract only executable README demonstrations, excluding installation calls.
blocks = re.findall(r'```r\n(.*?)```', (p / 'README.md').read_text(encoding='utf-8'), re.S)
code = '\n\n'.join(b for b in blocks if 'install_github' not in b)
Path('dev/release_0.1.0/phase7_readme_code.R').write_text(code + '\n', encoding='utf-8')
