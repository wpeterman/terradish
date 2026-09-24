from pathlib import Path
import json

root = Path('dev/check/phase8-package')
p = root / 'CITATION.cff'
p.write_text(p.read_text(encoding='utf-8').replace('version: 0.0.56', 'version: 0.0.57'), encoding='utf-8')
p = root / '.zenodo.json'
metadata = json.loads(p.read_text(encoding='utf-8'))
metadata['version'] = '0.0.57'
p.write_text(json.dumps(metadata, indent=2) + '\n', encoding='utf-8')
p = root / 'cran-comments.md'
p.write_text(p.read_text(encoding='utf-8').replace('0.0.56', '0.0.57').replace(
    'Final .56 source-archive and as-cran checks are pending.',
    'Final .57 source-archive and as-cran checks are pending.'), encoding='utf-8')
