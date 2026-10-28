"""Add the same preprocessing qualification to source and executed guides.

This changes prose only. Retain the exact executed-source hashes and record
the publication hashes separately; code and numerical outputs are unchanged.
"""
from pathlib import Path
import hashlib
import json
import re

repo = Path(__file__).resolve().parent.parent
guides = ['getting-started', 'gaussian-scale-optimization', 'model-comparison',
          'spline-conductance', 'ibe-ibr-workflow']
qualification = (
    'For the SLiM SNP response, marker ascertainment and pooled allele-frequency\n'
    'weights use all retained animals. Held-out site scores therefore evaluate\n'
    'transfer conditional on that complete response construction. They do not\n'
    'validate end-to-end SNP discovery or training-only genetic standardization.\n'
    'That target requires estimating marker filters and frequency weights within\n'
    'training data, then applying the same definition to test animals before\n'
    'scoring. The provided CV interface operates on one fixed response matrix.\n\n')

def chunks(text):
    return '\n'.join(re.findall(r'^```\{r[^\n]*\}\n(.*?)^```', text, re.M | re.S))

receipts = []
for guide in guides:
    source = repo / 'inst' / 'vignette-source' / f'{guide}.Rmd'
    target = repo / 'vignettes' / f'{guide}.Rmd'
    before = source.read_text(encoding='utf-8')
    code_before = chunks(before)
    changed = before
    if qualification not in changed:
        location = changed.index('## Quick reference')
        changed = changed[:location] + qualification + changed[location:]
    published = target.read_text(encoding='utf-8')
    if qualification not in published:
        location = published.index('## Quick reference')
        published = published[:location] + qualification + published[location:]
    assert chunks(changed) == code_before
    if guide == 'getting-started':
        old = ('**Isolation by resistance (IBR)** represents the landscape as a graph of\n'
               'connected raster cells. Cells have relative conductance, with larger values\n'
               'making connections through them stronger. The graph determines resistance\n'
               'distances and a corresponding covariance among sampling locations. A\n'
               'measurement model links that graph quantity to the observed genetic\n'
               'summary. Conductance is a fitted model quantity; interpreting it as\n'
               'movement, habitat suitability, or a causal effect requires additional\n'
               'evidence.')
        new = ('**Isolation by resistance (IBR)** links genetic patterns to a graph of\n'
               'connected raster cells. Larger cell conductance strengthens neighboring\n'
               'connections. The graph yields resistance distances and covariance among\n'
               'sampling locations, and a measurement model links those quantities to the\n'
               'observed genetic summary. Conductance is a fitted model quantity. Evidence\n'
               'about movement, habitat suitability, or causal effects requires information\n'
               'beyond the fitted surface.')
        if old in changed:
            assert old in published
            changed = changed.replace(old, new)
            published = published.replace(old, new)
    source.write_text(changed, encoding='utf-8', newline='\n')
    target.write_text(published, encoding='utf-8', newline='\n')
    receipts.append({'guide': guide,
       'executed_source_md5': hashlib.md5(before.encode()).hexdigest(),
       'publication_source_md5': hashlib.md5(changed.encode()).hexdigest(),
       'executed_chunk_sha256': hashlib.sha256(code_before.encode()).hexdigest(),
       'published_chunk_sha256': hashlib.sha256(chunks(changed).encode()).hexdigest(),
       'changes': 'Identical prose qualification in executable and precomputed source; no code/output change.'})
(repo / 'dev' / 'pond_vignette_20261013' / 'prose_publication_receipts.json').write_text(
    json.dumps(receipts, indent=2), encoding='utf-8')
print(json.dumps({'guides': len(receipts), 'all_chunk_hashes_equal':
    all(x['executed_chunk_sha256'] == x['published_chunk_sha256'] for x in receipts)}))
