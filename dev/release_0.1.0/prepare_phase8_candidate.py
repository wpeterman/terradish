from pathlib import Path
import shutil

source = Path('dev/check/phase7-package')
target = Path('dev/check/phase8-package')
assert not target.exists()
shutil.copytree(source, target, ignore=shutil.ignore_patterns(
    '*.o', '*.dll', '*.so', '*.a', 'Rplots.pdf', '*.html', '*_cache', '*_files'))
p = target / 'DESCRIPTION'
p.write_text(p.read_text(encoding='utf-8').replace('Version: 0.0.56', 'Version: 0.0.57'), encoding='utf-8')
p = target / 'NEWS.md'
p.write_text('''terradish 0.0.57
----------------
* Added fast analytic-derivative and AMG agreement tests that run on CRAN.
* Recorded final core regression audits and release validation evidence.

''' + p.read_text(encoding='utf-8'), encoding='utf-8')
