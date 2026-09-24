"""Copy checked package files, excluding build products and unrelated state."""
from pathlib import Path
import shutil
import sys

phase = sys.argv[1]
assert phase in ('5', '6')
root = Path.cwd()
candidate = root / 'dev/check' / f'phase{phase}-package'
files = [candidate / name for name in ('DESCRIPTION', 'NEWS.md', 'NAMESPACE', 'README.md')]
for folder in ('R', 'man', 'tests'):
    files.extend(p for p in (candidate / folder).rglob('*')
                 if p.is_file() and (folder != 'tests' or p.suffix == '.R'))
if phase == '6':
    files.extend([candidate / '.Rbuildignore', candidate / 'inst/WORDLIST'])
    files.extend((candidate / 'man-roxygen').glob('*.R'))
for src in files:
    dst = root / src.relative_to(candidate)
    assert dst.resolve().is_relative_to(root.resolve())
    if not dst.exists() or src.read_bytes() != dst.read_bytes():
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        print(dst.relative_to(root))
