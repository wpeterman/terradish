"""Create an isolated source snapshot and retire only unreferenced guide assets."""
from pathlib import Path
import hashlib
import json
import re
import shutil
import subprocess

repo = Path(__file__).resolve().parent.parent
validation = repo / 'dev' / 'pond_vignette_20261013' / 'validation'
validation.mkdir(parents=True, exist_ok=True)
guides = ['getting-started', 'gaussian-scale-optimization', 'categorical-conductance',
          'model-comparison', 'spline-conductance', 'ibe-ibr-workflow']
image_pattern = re.compile(r'precomputed/[A-Za-z0-9_.-]+[.]png')
current_refs = set()
for path in (repo / 'vignettes').glob('*.Rmd'):
    current_refs.update(image_pattern.findall(path.read_text(encoding='utf-8')))
retired = []
for guide in guides:
    old = subprocess.run(['git', 'show', f'HEAD:vignettes/{guide}.Rmd'],
                         cwd=repo, text=True, capture_output=True)
    if old.returncode:
        continue
    for reference in sorted(set(image_pattern.findall(old.stdout)) - current_refs):
        path = (repo / 'vignettes' / reference).resolve()
        if not path.is_relative_to((repo / 'vignettes' / 'precomputed').resolve()):
            raise RuntimeError(f'Asset outside intended directory: {path}')
        if path.exists():
            retired.append({'file': str(path.relative_to(repo)),
                            'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
            path.unlink()
(validation / 'retired_assets.json').write_text(json.dumps(retired, indent=2), encoding='utf-8')

source = validation / 'source' / 'terradish'
if source.exists():
    raise RuntimeError('Snapshot already exists; use its recorded check, or a fresh validation directory.')
source.mkdir(parents=True)
files = subprocess.check_output(['git', 'ls-files', '-z'], cwd=repo).decode().split('\0')
for folder in ['R', 'man', 'data', 'inst', 'tests', 'vignettes']:
    files.extend(str(p.relative_to(repo)) for p in (repo / folder).rglob('*') if p.is_file())
manifest = []
for name in sorted(set(files)):
    if not name or name.replace('\\', '/').startswith(('dev/', 'src/')):
        # Compiled objects are never source inputs; add tracked C++ below.
        continue
    path = repo / name
    if not path.is_file():
        continue
    if path.suffix in {'.dll', '.o', '.so', '.a', '.obj'}:
        continue
    destination = source / name
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(path, destination)
    manifest.append({'file': name.replace('\\', '/'),
                     'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
for name in subprocess.check_output(['git', 'ls-files', '-z', 'src'], cwd=repo).decode().split('\0'):
    if name and (repo / name).is_file():
        destination = source / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(repo / name, destination)
        manifest.append({'file': name, 'sha256': hashlib.sha256((repo / name).read_bytes()).hexdigest()})
(validation / 'source_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
print(json.dumps({'source': str(source), 'files': len(manifest), 'retired_assets': len(retired)}))
