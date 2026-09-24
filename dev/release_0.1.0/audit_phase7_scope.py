"""Scan the deliverable source, excluding compiler and rendered artifacts."""
from pathlib import Path
import re

root = Path('dev/check/phase7-package')
pattern = re.compile(r'terradish_directed|directed_rates|edge_gradient|edge_flow|terradish_hierarchical|conductance_field|wishart_drift_covariates|(^|[^g_])linear_conductance|pair_subset|terradish_multiscale|radish_multiscale|terradish_kron_reduce|kron_reduce|block_cg|pcg_jacobi|"pcg"|terradish_cv[^_]|terradish_cv_replicates|radish_cv|cv_model_selection|terradish_results|terradish_parameters|radish_parameters|hierarchical-conductance|directional-conductance|out-of-core|exact_refine = FALSE')
files = [root / 'README.md']
for directory in ('R', 'src', 'tests', 'vignettes', 'inst'):
    files.extend(p for p in (root / directory).rglob('*') if p.is_file()
                 and p.suffix in ('.R', '.r', '.Rmd', '.cpp', '.h', '.hpp', '.md', '.txt'))
hits = []
exceptions = []
for p in files:
    for number, line in enumerate(p.read_text(encoding='utf-8', errors='replace').splitlines(), 1):
        if pattern.search(line):
            entry = f'{p.relative_to(root)}:{number}:{line}'
            # The plan's substring expression also matches retained CV names.
            cleaned = re.sub(r'\.?terradish_cv_[A-Za-z_]+', '', line)
            if not pattern.search(cleaned):
                exceptions.append('retained CV identifier: ' + entry)
            elif p.name in ('radish.cpp', 'RcppExports.cpp', 'RcppExports.R') and 'block_cg' in line:
                exceptions.append('S9 shared-file compiled helper, not NAMESPACE-exported: ' + entry)
            elif p.name == 'test-core-scope.R' and 'exact_refine = FALSE' in line:
                exceptions.append('test asserting rejection of removed option: ' + entry)
            else:
                hits.append(entry)
print(f'{len(exceptions)} literal matches classified as retained identifiers, S9 helpers, or rejection tests.')
print('\n'.join(exceptions))
print('\n'.join(hits) if hits else f'Clean: {len(files)} source/document files scanned.')
assert 'export(block_cg' not in (root / 'NAMESPACE').read_text()
raise SystemExit(bool(hits))
