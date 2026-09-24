from pathlib import Path
import hashlib
import json

folder = Path('dev/release_0.1.0/phase8_audits')
p = folder / 'manifest.json'
rows = json.loads(p.read_text(encoding='utf-8'))
source_root = Path('C:/Users/peterman.73/OneDrive - The Ohio State University/Research/R_packages/terradish/review_work_20260924/package_audit')
known = {Path(row['adapted']).as_posix() for row in rows}
for rel, original in [('cv/e3_converted.R', 'cv/e3.R'), ('cv/e2_modes.R', 'cv/e2.R')]:
    target = folder / rel
    if target.as_posix() not in known:
        source = source_root / original
        rows.append(dict(source=str(source), source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                         adapted=str(target)))
for row in rows:
    target = Path(row['adapted'])
    row['adapted_sha256'] = hashlib.sha256(target.read_bytes()).hexdigest()
p.write_text(json.dumps(rows, indent=2) + '\n', encoding='utf-8')
print(f'Recorded {len(rows)} original/adapted source pairs.')
