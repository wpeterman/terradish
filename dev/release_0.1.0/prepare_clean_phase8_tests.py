from pathlib import Path
import shutil

source = Path('dev/check/phase8-package')
target = Path('dev/check/phase8-source-tests')
assert not target.exists()
shutil.copytree(source, target, ignore=shutil.ignore_patterns(
    '*.o', '*.dll', '*.so', '*.a', '*.gcda', '*.gcno', '*.html', 'Rplots.pdf', '*_cache', '*_files'))
p = Path('dev/release_0.1.0/validate_phase8_tests.R')
code = p.read_text(encoding='utf-8').replace('phase8-package', 'phase8-source-tests')
code = code.replace('phase8_tests.rds', 'phase8_tests_clean.rds')
p.with_name('validate_phase8_tests_clean.R').write_text(code, encoding='utf-8')
