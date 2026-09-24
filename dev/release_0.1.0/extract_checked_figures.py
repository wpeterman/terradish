from pathlib import Path
from html.parser import HTMLParser
import base64
import json
import sys

source = Path(sys.argv[1]) if len(sys.argv) > 1 else Path('dev/check/phase8/terradish.Rcheck/terradish/doc')
output = Path(sys.argv[2]) if len(sys.argv) > 2 else Path('dev/check/phase8-figure-review')
output.mkdir(exist_ok=True)
index = []

class Images(HTMLParser):
    def handle_starttag(self, tag, attrs):
        if tag != 'img':
            return
        attrs = dict(attrs)
        data = attrs.get('src', '')
        if data.startswith('data:image/png;base64,'):
            name = f'{self.document}-{self.number:02d}.png'
            (output / name).write_bytes(base64.b64decode(data.split(',', 1)[1]))
            index.append(dict(document=self.document, figure=self.number, file=str(output / name), alt=attrs.get('alt', '')))
            self.number += 1

for p in sorted(source.glob('*.html')):
    parser = Images()
    parser.document = p.stem
    parser.number = 1
    parser.feed(p.read_text(encoding='utf-8'))
(output / 'index.json').write_text(json.dumps(index, indent=2), encoding='utf-8')
for row in index:
    print(row['file'], row['alt'][:120])
