from pathlib import Path
import difflib

old = Path('dev/check/phase5-package/R')
new = Path('dev/check/phase6-package/R')
paragraphs = []
for path in sorted(new.glob('*.R')):
    previous = (old / path.name).read_text(encoding='utf-8').splitlines()
    current = path.read_text(encoding='utf-8').splitlines()
    for line in difflib.unified_diff(previous, current):
        if line.startswith("+#'") and not line.startswith("+#' @"):
            paragraphs.append(line[3:].strip())
for path in sorted(Path('dev/check/phase6-package/man-roxygen').glob('*.R')):
    if 'fit' not in path.name:
        paragraphs.extend(line[2:].strip() for line in path.read_text(encoding='utf-8').splitlines()
                          if line.startswith("#'") and not line.startswith("#' @"))
Path('dev/release_0.1.0/phase6_prose.txt').write_text('\n'.join(paragraphs), encoding='utf-8')
