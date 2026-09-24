from pathlib import Path
import hashlib
import json
import sys

root = Path.cwd()
final = len(sys.argv) > 1 and sys.argv[1] == 'final'
candidate = root / ('dev/check/phase8-final-package' if final else 'dev/check/phase8-package')
archive = root / ('dev/check/phase8-final/terradish_0.0.57.tar.gz' if final else 'dev/check/phase8/terradish_0.0.57.tar.gz')
files = [candidate / name for name in ('DESCRIPTION', 'NAMESPACE', 'NEWS.md', 'README.md')]
for directory, suffixes in [('R', {'.R'}), ('src', {'.cpp', '.h'}),
                            ('man', {'.Rd'}), ('tests', {'.R'}), ('vignettes', {'.Rmd'})]:
    files.extend(p for p in (candidate / directory).rglob('*') if p.is_file() and p.suffix in suffixes)
record = dict(archive=str(archive.relative_to(root)),
              archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),
              candidate_files={str(p.relative_to(candidate)).replace('\\', '/'):
                  hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)})
target = root / ('dev/release_0.1.0/phase8_final_fingerprint.json' if final else 'dev/release_0.1.0/phase8_fingerprint.json')
target.write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')
print(record['archive_sha256'])
