from pathlib import Path
import shutil

source = Path('dev/check/phase8-package')
target = Path('dev/check/phase8-coverage')
assert not target.exists()
shutil.copytree(source, target, ignore=shutil.ignore_patterns(
    '*.o', '*.dll', '*.so', '*.a', '*.gcda', '*.gcno', '*.html', 'Rplots.pdf', '*_cache', '*_files'))
