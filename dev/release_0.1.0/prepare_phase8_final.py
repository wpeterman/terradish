from pathlib import Path
import shutil

source = Path('dev/check/phase8-package')
target = Path('dev/check/phase8-final-package')
assert not target.exists()
shutil.copytree(source, target, ignore=shutil.ignore_patterns(
    '*.o', '*.dll', '*.so', '*.a', '*.gcda', '*.gcno', '*.html', 'Rplots.pdf', '*_cache', '*_files'))
p = target / 'vignettes/spline-conductance.Rmd'
s = p.read_text(encoding='utf-8')
s = s.replace('## Observed vs. fitted genetic distances', '## Observed genetic distance versus fitted resistance')
s = s.replace('***Observed vs. fitted pairwise Fst for all candidate models.*** Each point is a pair of sampling sites. The diagonal (dashed) represents a perfect fit.',
    '***Observed pairwise Fst versus fitted resistance for all candidate models.*** Each point is a pair of sampling sites. The blue line is a descriptive least-squares regression.')
s = s.replace('All four panels use the same data and MLPE correlation structure. Systematic deviations from the diagonal (e.g. all large-distance pairs systematically under-predicted) indicate lack of fit that could arise from the conductance formula, measurement model, response construction, or other misspecification.',
    'All four panels use the same data and MLPE correlation structure. Their horizontal axes show resistance, whose scale can differ across fitted conductance models. The blue regression lines summarize the plotted relationship; they do not account for dependence among pairs. Curvature or other systematic departures warrant checking the conductance formula, measurement model, and response construction. Use likelihood diagnostics and cross-validation for model comparison.')
assert 'diagonal (dashed)' not in s
p.write_text(s, encoding='utf-8')
p = target / 'NEWS.md'
p.write_text(p.read_text(encoding='utf-8').replace(
    '* Recorded final core regression audits and release validation evidence.',
    '* Recorded final core regression audits and release validation evidence.\n* Corrected the spline fit-plot caption to describe its resistance axis and\n  descriptive regression line.'), encoding='utf-8')
# The full numerical tests and coverage apply unchanged to this candidate.
for directory in ('R', 'src', 'tests'):
    for file in source.joinpath(directory).rglob('*'):
        if file.is_file() and file.suffix in ('.R', '.cpp', '.h'):
            assert file.read_bytes() == target.joinpath(file.relative_to(source)).read_bytes()
print('Final candidate retains identical R, C++, and test sources.')
