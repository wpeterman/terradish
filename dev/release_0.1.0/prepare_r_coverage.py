from pathlib import Path
import shutil

source = Path('dev/check/phase8-final-package')
target = Path('dev/check/phase8-r-coverage')
assert not target.exists()
shutil.copytree(source, target, ignore=shutil.ignore_patterns(
    '*.o', '*.dll', '*.so', '*.a', '*.gcda', '*.gcno', '*.html', 'Rplots.pdf', '*_cache', '*_files'))
p = Path('dev/release_0.1.0/report_phase8_coverage.R')
p.with_name('report_phase8_r_coverage.R').write_text(
    p.read_text(encoding='utf-8').replace('phase8_coverage', 'phase8_r_coverage')
      .replace('phase8_uncovered', 'phase8_r_uncovered'), encoding='utf-8')
