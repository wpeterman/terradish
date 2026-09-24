from pathlib import Path
import re
from html.parser import HTMLParser

pkg = Path('dev/check/phase7-package')
parts = []
for p in [pkg / 'README.md', *(pkg / 'vignettes').glob('*.Rmd')]:
    text = re.sub(r'```.*?```', '', p.read_text(encoding='utf-8'), flags=re.S)
    parts.append(f'FILE {p.name}\n{text}')
Path('dev/release_0.1.0/phase7_prose.txt').write_text('\n\n'.join(parts), encoding='utf-8')

class PlainText(HTMLParser):
    def __init__(self):
        super().__init__()
        self.parts = []
        self.skip = 0
    def handle_starttag(self, tag, attrs):
        if tag in ('style', 'script'): self.skip += 1
    def handle_endtag(self, tag):
        if tag in ('style', 'script'): self.skip -= 1
    def handle_data(self, data):
        if not self.skip: self.parts.append(data)

for folder in ('phase7-render-first', 'phase7-render-remaining-v2'):
    for p in (Path('dev/check') / folder).glob('*.html'):
        parser = PlainText()
        parser.feed(p.read_text(encoding='utf-8'))
        output = '\n'.join(parser.parts)
        target = Path('dev/check/phase7-render-text')
        target.mkdir(exist_ok=True)
        (target / (p.stem + '.txt')).write_text(output, encoding='utf-8')
        for line in output.splitlines():
            if 'Convergence:' in line:
                print(p.stem, line)
