"""
Scan-pull workflow. Understands the folder layout Knee-Segmentation-Pipeline
writes for every pull on the drive:

    <pull>/
      sorted/mrn_code_mapping.csv      code -> MRN, knee side (laterality), ...
      sorted/<code>/                   one folder per scan
      cleaned_segs/<code>.nrrd         manually cleaned segmentations (some pulls:
                                       cleaned_segmentations/)
      dixon_segmentations/segmentations/<code>.nrrd   Dixon model output
      mapss_segmentations/segmentations/<code>.nrrd   MAPSS model output
      kneemeasure.toml                 optional settings for this pull (below)
      kneemeasure_review.json          human landmark reviews (written by `kneemeasure review`)
      kneemeasure/                     <- results are written here

Segmentation source per scan ('best', the default): the cleaned segmentation
if there is one, otherwise the Dixon model output. Which one was used is
recorded per case (segmentation_source column). 'cleaned', 'dixon' or
'mapss' use only that source.

kneemeasure.toml (all keys optional):

    name = "Controls"                  # display name (default: folder name)
    source = "best"                    # best | cleaned | dixon | mapss
    [exclude]                          # scans to leave out, with the reason
    6YGMS = "tibia mostly outside the scan"
"""
import os
import tomllib

from . import review as review_store
from .runner import run_cohort, compare_runs, write_csv

MAPPING = os.path.join('sorted', 'mrn_code_mapping.csv')
CLEANED_DIRS = ('cleaned_segs', 'cleaned_segmentations')
MODEL_DIRS = {'dixon': os.path.join('dixon_segmentations', 'segmentations'),
              'mapss': os.path.join('mapss_segmentations', 'segmentations')}
SOURCES = ('best', 'cleaned', 'dixon', 'mapss')
RESULTS_DIR = 'kneemeasure'


def _nrrds(d):
    if not os.path.isdir(d):
        return {}
    return {f[:-5]: os.path.join(d, f) for f in sorted(os.listdir(d))
            if f.lower().endswith('.nrrd') and not f.startswith('.')}


def load_config(pull_dir):
    p = os.path.join(pull_dir, 'kneemeasure.toml')
    if not os.path.exists(p):
        return {}
    with open(p, 'rb') as f:
        cfg = tomllib.load(f)
    unknown = set(cfg) - {'name', 'source', 'exclude'}
    if unknown:
        raise ValueError(f'{p}: unknown setting(s) {sorted(unknown)}')
    return cfg


def inspect_pull(pull_dir, source=None):
    """Everything needed to run a pull, without running it. Returns a dict
    with 'eligible' False and a 'reason' if the pull can't be measured."""
    pull_dir = os.path.abspath(pull_dir)
    cfg = load_config(pull_dir)
    source = source or cfg.get('source', 'best')
    if source not in SOURCES:
        raise ValueError(f'source must be one of {SOURCES}, got {source!r}')
    spec = {'path': pull_dir, 'name': cfg.get('name', os.path.basename(pull_dir)), 'source': source,
            'exclude': dict(cfg.get('exclude', {})), 'laterality': os.path.join(pull_dir, MAPPING),
            'out': os.path.join(pull_dir, RESULTS_DIR), 'eligible': False, 'reason': None,
            'review_file': os.path.join(pull_dir, review_store.REVIEW_FILE)}

    if not os.path.exists(spec['laterality']):
        spec['reason'] = f'no {MAPPING}'
        return spec
    with open(spec['laterality'], encoding='utf-8-sig') as f:
        header = f.readline().lower()
    if 'laterality' not in header:
        spec['reason'] = f'{MAPPING} has no laterality column'
        return spec

    cleaned = {}
    for d in CLEANED_DIRS:
        cleaned.update(_nrrds(os.path.join(pull_dir, d)))
    models = {k: _nrrds(os.path.join(pull_dir, v)) for k, v in MODEL_DIRS.items()}
    files, used = {}, {}
    if source in ('best', 'cleaned'):
        for code, path in cleaned.items():
            files[code], used[code] = path, 'cleaned'
    if source == 'best':
        for code, path in models['dixon'].items():
            if code not in files:
                files[code], used[code] = path, 'dixon model'
    if source in ('dixon', 'mapss'):
        for code, path in models[source].items():
            files[code], used[code] = path, f'{source} model'
    if not files:
        spec['reason'] = f'no {source} segmentations found'
        return spec

    spec.update(eligible=True, files=[files[c] for c in sorted(files)],
                case_info={c: {'segmentation_source': used[c]} for c in files},
                counts={k: sum(1 for v in used.values() if v == k) for k in sorted(set(used.values()))})
    return spec


def run_pull(pull_dir, out=None, source=None, jobs=4, rerun=False, quiet=False, **kw):
    spec = inspect_pull(pull_dir, source)
    if not spec['eligible']:
        raise ValueError(f"{spec['name']}: can't be measured -- {spec['reason']}")
    out = out or spec['out']
    reviews = review_store.load(spec['review_file'])
    entered = {c: ((v.get('laterality') or {}).get('side'), 'entered in app') for c, v in reviews['cases'].items()}
    rows, flags = run_cohort(spec['files'], out, laterality=spec['laterality'], carry=['mrn'],
                             exclude=spec['exclude'], case_info=spec['case_info'], jobs=jobs,
                             cohort_name=spec['name'], rerun=rerun, quiet=quiet,
                             reviews=reviews, sides=entered, **kw)
    spec['out'] = out
    return spec, rows, flags


def run_all(root, report=None, jobs=4, rerun=False, only=None, quiet=False, extras=False):
    """Run every measurable pull under root and print a status table.
    report: a folder to also write status.csv/.md and a cross-pull
    comparison to (nothing is written outside the pulls otherwise).
    Returns the list of pull specs."""
    root = os.path.abspath(root)
    out = os.path.abspath(report) if report else None
    specs, status = [], []
    for name in sorted(os.listdir(root)):
        d = os.path.join(root, name)
        if not os.path.isdir(d) or name.startswith(('.', '_')) or os.path.abspath(d) == out:
            continue
        if not os.path.isdir(os.path.join(d, 'sorted')):
            continue
        if only and name not in only:
            continue
        spec = inspect_pull(d)
        row = {'pull': spec['name'], 'folder': name, 'results': spec['out'] if spec['eligible'] else ''}
        if not spec['eligible']:
            row['status'] = f"skipped: {spec['reason']}"
            print(f"{name}: skipped -- {spec['reason']}")
            status.append(row)
            continue
        print(f"{name}: {len(spec['files'])} segmentation(s) "
              f"({', '.join(f'{v} {k}' for k, v in spec['counts'].items())})")
        try:
            spec, rows, flags = run_pull(d, jobs=jobs, rerun=rerun, quiet=quiet, extras=extras)
        except Exception as e:
            row['status'] = f'FAILED: {e}'
            status.append(row)
            continue
        measured = [r for r in rows if not r.get('error')]
        cats = {}
        for f in flags:
            cats[f['category']] = cats.get(f['category'], 0) + 1
        row.update(status='measured', cases=len(rows), errors=len(rows) - len(measured),
                   sources=', '.join(f'{v} {k}' for k, v in spec['counts'].items()),
                   excluded=cats.get('excluded', 0), not_segmented=cats.get('missing segmentation', 0),
                   failed_measurements=cats.get('failed', 0), segmentation_warnings=cats.get('segmentation', 0),
                   tt_out_of_window=sum(1 for r in rows if r.get('TT_possibly_out_of_window') == 1),
                   outliers=cats.get('outlier', 0))
        status.append(row)
        specs.append(spec)

    cols = ['pull', 'status', 'cases', 'sources', 'excluded', 'not_segmented', 'errors', 'failed_measurements',
            'segmentation_warnings', 'tt_out_of_window', 'outliers', 'results', 'folder']
    print()
    for r in status:
        extra = '' if r.get('status') != 'measured' else (
            f" -- {r['cases']} scans ({r['sources']}), {r['failed_measurements']} failed measurement(s), "
            f"{r['tt_out_of_window']} with the tuberosity probably out of view")
        print(f"{r['pull']}: {r['status']}{extra}")
    if out:
        os.makedirs(out, exist_ok=True)
        write_csv(os.path.join(out, 'status.csv'), status, cols)
        _write_status_md(os.path.join(out, 'status.md'), status, cols)
        if len(specs) >= 2:
            compare_runs({s['name']: s['out'] for s in specs}, os.path.join(out, 'comparison'))
    return specs


def _write_status_md(path, status, cols):
    show = [c for c in cols if c != 'folder']
    lines = ['# KAMINO status', '', '| ' + ' | '.join(show) + ' |', '|' + '---|' * len(show)]
    for r in status:
        lines.append('| ' + ' | '.join(str(r.get(c, '')) for c in show) + ' |')
    lines += ['', "Each pull's results are in its own `kneemeasure/` folder. "
              '`comparison/` compares the measured pulls.']
    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')
