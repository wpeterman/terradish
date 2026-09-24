from pathlib import Path
import tarfile
import re
import sys

final = len(sys.argv) > 1 and sys.argv[1] == 'final'
archive = Path('dev/check/phase8-final/terradish_0.0.57.tar.gz' if final else
               'dev/check/phase8-doc-outputs/terradish_0.0.57.tar.gz')
out = Path('dev/check/phase8-final-caption-docs' if final else 'dev/check/phase8-final-docs')
out.mkdir(exist_ok=True)
with tarfile.open(archive) as package:
    for member in package.getmembers():
        if member.name.startswith('terradish/inst/doc/') and member.name.endswith('.html'):
            data = package.extractfile(member).read()
            target = out / Path(member.name).name
            target.write_bytes(data)
            text = data.decode('utf-8')
            images = len(re.findall(r'<img\b', text))
            print(target.name, 'images:', images)
            if target.stem in ('spline-conductance', 'model-comparison'):
                assert images > 0, f'Missing computed figures in {target.name}'
                assert 'data:image/png' in text
            if final and target.stem == 'spline-conductance':
                assert 'The diagonal (dashed)' not in text
                assert 'descriptive least-squares regression' in text
print('Ordinary-build spline/model-comparison figures are present.')
