"""Create the portable reading index and audit its executed mathematical sources."""
from pathlib import Path
import base64
import hashlib
import html
import json
import re
import shutil

repo = Path(__file__).resolve().parent.parent
review = Path('C:/Users/peterman.73/Downloads/terradish_pond_vignettes_review_20261013')
receipts = json.loads((review / 'EXECUTION_RECEIPTS.json').read_text(encoding='utf-8'))
titles = {
    'getting-started': 'Getting Started with terradish',
    'gaussian-scale-optimization': 'Gaussian Scale Optimization',
    'model-comparison': 'Model Comparison and Selection in terradish',
    'spline-conductance': 'Spline Conductance Models',
    'ibe-ibr-workflow': 'Isolation by Environment and Resistance',
    'categorical-conductance': 'Categorical Landscape Covariates',
}
gaussian = receipts['gaussian-scale-optimization']['profile']
validation = []
for guide in titles:
    source = (repo / 'inst/vignette-source' / f'{guide}.Rmd').read_text(encoding='utf-8')
    displayed = (repo / 'vignettes' / f'{guide}.Rmd').read_text(encoding='utf-8')
    assert '#> Error:' not in displayed
    assert '\\eqn{' not in source and '\\deqn{' not in source
    assert '\u2014' not in source
    prose = re.sub(r'^```.*?^```[^\n]*', '', source, flags=re.M | re.S)
    assert len(re.findall(r'(?<!\\)\$', prose)) % 2 == 0
    math = re.findall(r'\$\$(.*?)\$\$|(?<!\$)\$([^$\n]+?)\$(?!\$)', prose, re.S)
    for alternatives in math:
        expression = ''.join(alternatives)
        level = 0
        for char in expression:
            if char == '{': level += 1
            if char == '}': level -= 1
            assert level >= 0
        assert level == 0
    document = (review / f'{guide}.html').read_text(encoding='utf-8')
    assert '<html' in document and '<title>' in document
    assert 'precomputed/' not in re.sub(r'<code.*?</code>', '', document, flags=re.S)
    assert 'data:image/png;base64,' in document or guide == 'ibe-ibr-workflow'
    validation.append({'guide': guide, 'math_expressions': len(math),
        'source_sha256': hashlib.sha256(source.encode()).hexdigest(),
        'html_sha256': hashlib.sha256((review / f'{guide}.html').read_bytes()).hexdigest()})

cv_rows = []
for guide, x in receipts.items():
    for name in ['cv', 'cv_cov']:
        if name in x:
            for row in x[name]['results']:
                assert not row['error'] and row['nuisance_iterations'] == 0
                cv_rows.append({'guide': guide, **row})
assert len(cv_rows) == 42
(review / 'DOCUMENT_QA.json').write_text(json.dumps(validation, indent=2), encoding='utf-8')
summary = {
    'guides': list(titles), 'cv_scores': len(cv_rows),
    'snp_loci': 2500, 'microsatellite_loci': 24, 'retained_animals': 480,
    'ponds': 24, 'diploids_per_pond': 20, 'raw_captures': 780,
    'gaussian_scale_m': gaussian['estimate'], 'profile_interval_m': gaussian['interval'],
    'source_test_assertions': 1114, 'source_test_failures': 0,
    'source_test_warnings': 196, 'source_session_skips': 1,
    'archive_check': 'in progress',
    'render_qa': 'HTML structural and mathematical-source checks plus local scientific-figure inspection; remote browser layout not verified',
    'branch': 'codex/amphibian-vignette', 'committed': False,
}
(review / 'SUMMARY.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
for name in ['vignette_cv_scores.csv', 'vignette_cv_summary.csv',
             'vignette_fit_validation.csv', 'prose_publication_receipts.json']:
    shutil.copy2(repo / 'dev/pond_vignette_20261013' / name, review / name)
shutil.copy2(repo / 'dev/amphibian_example/INDIVIDUAL_REDESIGN_20261013.md', review / 'INDIVIDUAL_REDESIGN.md')
shutil.copy2(repo / 'inst/extdata/ohio-ponds/PROVENANCE.md', review / 'PROVENANCE.md')

md = f'''# Pond vignette development review

October 13, 2026. Six guides now provide a shared educational sequence on the
real Ohio landscape. The main dataset has 24 hypothetical forest ponds,
20 retained diploids per pond, 2,500 SNPs, and 24 matched microsatellite loci.
All 42 displayed model-fold predictive scores execute successfully.

## Read the guides in this order

| Guide | Teaching purpose |
|---|---|
| [Getting started](getting-started.html) | Raw genotypes, sample and pond order, covariance geometry, measurement model, coefficients, and matched microsatellites |
| [Gaussian context](gaussian-scale-optimization.html) | Spatial scale, profile uncertainty, kernel distances, and fixed-versus-estimated prediction |
| [Model comparison](model-comparison.html) | Multivariate candidates, multiple starts, information sensitivity, nested tests, and validation |
| [Splines](spline-conductance.html) | Response shape, stored bases, supported ranges, and the cost of extra flexibility |
| [IBE and IBR](ibe-ibr-workflow.html) | MLPE paired differences, PSD covariance kernels, and biological exchange as distinct constructions |
| [Categorical covariates](categorical-conductance.html) | Factor rasters, reference coding, conductance ratios, and a declared class-effect benchmark |

Each guide includes annotated code, executed output, interpretation, parameter
units, a quick-reference workflow, and links to the related guides. Heavy
analyses are precomputed for ordinary package builds; the executable sources
remain bundled. The real grid and all ponds are retained in the analysis.

## What the executed examples show

- The Gaussian estimate is {gaussian['estimate']:.1f} m against a 500 m
  generating context. Its conditional information-40 profile interval is
  {gaussian['interval']['lower']:.1f} to {gaussian['interval']['upper']:.1f} m.
- Fixed and estimated forest context predict better than the uniform graph.
  The known fixed context is an oracle comparator and slightly outperforms
  the estimated scale on these held-out ponds.
- Forest alone leads the information criteria. Adding slope improves the
  summed predictive score by only about 0.75; its generating direct effect
  is zero. Adding elevation hurts one spatial block.
- The spline's AIC advantage is less than one unit, but its held-out score
  is about 30 units below the linear model. It does not justify a nonlinear
  ecological story in this linear-generating history.
- The IBE positive control recovers the declared 0.10 component with both
  constructions, while retaining their different likelihoods and constraints.
  This is an independent-Gaussian-feature benchmark, not a SNP simulation.
- A binary unsmoothed class model on the actual SNP response drives its
  contrast toward an extreme value and stalls. That failed attempt is
  preserved and not interpreted. The categorical coding lesson uses a
  separate known-class-effect covariance benchmark on the same landscape.

The SNP validation target is conditional on globally ascertained markers and
pooled allele-frequency weights. It does not validate end-to-end discovery or
training-only genetic standardization. Wishart information is an assumption,
not the marker count or the 40 allele copies sampled at each pond.

## Terrestrial redesign

[The redesign protocol](INDIVIDUAL_REDESIGN.md) separates recruitment capacity
from successful dispersal. It proposes a smaller real forest domain, explicit
age and sex, and separate records of movement proposal, survival, arrival,
reproduction, and surviving offspring. The next step is an ecology-only pilot
before new ancestry, marker generation, or genetic fitting. Rates remain
hypothetical. The existing capacity-only histories and their mixed validation
results are preserved.

The design is informed by the verified [PLOS ONE study](https://doi.org/10.1371/journal.pone.0062184)
and [Molecular Ecology study](https://doi.org/10.1111/mec.12747), while keeping
the Ohio setting and species extension separate from fitted field calibration.

## Validation and package scope

The full source suite passed 1,114 assertions, with no failures, 196 warning
receipts matching the repair baseline, and one source-session parallel skip
reserved for installed-package validation. The new data, public optimizer
forwarding, and checkpoint cases are included. An isolated source-archive
check is running; its final result will be added here.

The guides use the suite's existing format and educational conventions.
Prose received the required editing pass and vocabulary audit. Cadence flags
remain in some tables and repeated parameter explanations, where consistent
wording keeps units and model roles unambiguous. No em dashes or British
spelling were detected.

HTML and mathematical-source checks pass, and the generated scientific figures
were visually inspected locally. The browser preview is unavailable on this
host, so rendered phone layout and MathJax display remain unverified. HTML
files embed their figures; mathematical rendering uses the document's MathJax
loader. Native phone summaries use readable notation and scientific images.

Development remains on `codex/amphibian-vignette`. No commit, version bump,
merge, push, or release replacement was made. [Data provenance](PROVENANCE.md)
records the landscape, ancestry, panels, and limitations. The archive check is
for this development branch, not a replacement submission artifact.
'''
(review / 'REVIEW.md').write_text(md, encoding='utf-8')
links = ''.join(f'<li><a href="{g}.html">{html.escape(t)}</a></li>' for g, t in titles.items())
sections = [f'<h1>Pond vignette development review</h1><p>October 13, 2026. Six executed educational guides with matched genetic data on a real Ohio landscape.</p><ol>{links}</ol>',
            '<h2>Shared data and scale</h2><p>24 forest ponds; 20 retained diploids per pond; 2,500 SNPs and 24 matched microsatellite loci. The 500 m generating context is recovered near 480 m. The conditional profile spans approximately 311 to 697 m.</p>',
            '<h2>Model choice and interpretation</h2><p>The spline slightly improves AIC but predicts worse than the linear model. A slope extension gains only about 0.75 predictive units, despite a zero direct generating effect. The guides keep that uncertainty visible.</p>',
            '<h2>Terrestrial redesign</h2><p><a href="INDIVIDUAL_REDESIGN.md">Read the proposed protocol</a>: recruitment capacity and successful dispersal are separate processes. An ecology-only pilot precedes new genomic runs.</p>',
            '<p><a href="REVIEW.md">Complete review and validation notes</a> | <a href="PROVENANCE.md">Provenance</a></p>']
for filename, caption in [('gaussian-scale-optimization-profile-plot-1.png', 'Conditional profile support for forest context'),
                          ('model-comparison-cv-plot-1.png', 'Predictive scores of multivariate candidates'),
                          ('spline-conductance-marginal-1.png', 'Spline association and model interval')]:
    data = (repo / 'vignettes/precomputed' / filename).read_bytes()
    sections.append(f'<figure><img src="data:image/png;base64,{base64.b64encode(data).decode()}" alt="{caption}"><figcaption>{caption}</figcaption></figure>')
document = '<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>Pond vignette development review</title><style>body{max-width:850px;margin:24px auto;padding:0 18px;font:17px/1.55 system-ui;color:#172329}h1,h2{line-height:1.2}a{color:#006790}img{max-width:100%;height:auto}figure{margin:25px 0}figcaption{font-size:15px;color:#40545e}</style><body>' + ''.join(sections) + '</body></html>'
(review / 'REVIEW.html').write_text(document, encoding='utf-8')
print(json.dumps({'guides': len(validation), 'math_expressions': sum(x['math_expressions'] for x in validation),
                  'cv_scores': len(cv_rows), 'review': str(review / 'REVIEW.html')}))
