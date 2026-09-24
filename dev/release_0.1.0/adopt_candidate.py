"""Copy checked package files, excluding build products and unrelated state."""
from pathlib import Path
import shutil
import sys

phase = sys.argv[1]
assert phase in ('5', '6', '7', '8')
root = Path.cwd()
candidate = root / (sys.argv[2] if len(sys.argv) > 2 else f'dev/check/phase{phase}-package')
assert candidate.resolve().is_relative_to((root / 'dev/check').resolve())
files = [candidate / name for name in ('DESCRIPTION', 'NEWS.md', 'NAMESPACE', 'README.md')]
for folder in ('R', 'man', 'tests'):
    files.extend(p for p in (candidate / folder).rglob('*')
                 if p.is_file() and (folder != 'tests' or p.suffix == '.R'))
if phase in ('6', '7', '8'):
    files.extend([candidate / '.Rbuildignore', candidate / 'inst/WORDLIST'])
    files.extend((candidate / 'man-roxygen').glob('*.R'))
if phase in ('7', '8'):
    files.extend(candidate / name for name in ('CITATION.cff', '.zenodo.json', 'cran-comments.md'))
    files.extend((candidate / 'vignettes').glob('*.Rmd'))
    files.append(candidate / 'vignettes/precompute.R')
    files.extend((candidate / 'inst/examples').glob('*.R'))
    files.append(candidate / 'inst/benchmarks/compare-radish-terradish.R')
for src in files:
    dst = root / src.relative_to(candidate)
    assert dst.resolve().is_relative_to(root.resolve())
    if not dst.exists() or src.read_bytes() != dst.read_bytes():
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        print(dst.relative_to(root))
