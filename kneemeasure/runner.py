"""
Cohort runner: segmentations in, a results folder out.

    <out>/
      measurements.csv            one row per case: values, review status, tibia length, notes
      qc/<measurement>_grid.png   paginated QC grids (+ segmentation overview)
      .kneemeasure/               hidden: rerun cache, run_info.json, log.txt

With extras=True (`--extras`) it also writes:

      statistics.csv              descriptive stats (all cases, and by knee side)
      qc_flags.csv                failures, segmentation warnings, laterality problems, outliers
      segmentation_check.csv      per-case segmentation sanity checks
      landmarks.csv               every landmark: reoriented index + LPS mm
      measurement_definitions.csv column, unit, definition
      changes.csv + history/      what changed since the previous run
      summary.md, plots/distributions.png
      qc/cases/<case>.png         one-page QC sheet per flagged case (or all)
      slicer/<case>.mrk.json      landmarks as 3D Slicer markups
"""
import csv
import datetime
import glob
import shutil
import json
import math
import os
import platform
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import nrrd
import numpy as np

from . import __version__, kra, tt_tg, tibial_pcl, stats, plots, qc, slicer
from . import review as review_store
from .cache import ResultCache, cache_key, code_fingerprint, file_sha1
from .laterality import LateralityTable, resolve
from .pipeline import MEASURES, MEASURE_BY_KEY, measure_case, value_columns
from .validation import COLUMNS as VALIDATION_COLUMNS
from .volume import find_segmentations, case_name_from_path, header_laterality

VOLUME_COLUMNS = [c for c in VALIDATION_COLUMNS if c.startswith('vol_')]
TT_DEPENDENT = ('TT_TG_mm', 'sTT_TG_mm', 'TTL', 'TTL_t_mm', 'TTL_T_mm')
META_DIR = '.kneemeasure'
EXTRA_OUTPUTS = ('statistics.csv', 'qc_flags.csv', 'segmentation_check.csv', 'landmarks.csv',
                 'measurement_definitions.csv', 'changes.csv', 'summary.md', 'run_info.json', 'log.txt',
                 'history', 'plots', 'slicer', os.path.join('qc', 'cases'), '.cache')


def _tidy(out_dir, extras):
    """Move the hidden cache/log from older layouts into .kneemeasure/, and
    (without extras) remove extra outputs an earlier run left behind -- all
    of them regenerable with --extras."""
    meta = os.path.join(out_dir, META_DIR)
    os.makedirs(meta, exist_ok=True)
    old_cache, new_cache = os.path.join(out_dir, '.cache'), os.path.join(meta, 'cache')
    if os.path.isdir(old_cache) and not os.path.exists(new_cache):
        os.replace(old_cache, new_cache)
    if extras:
        return
    for name in EXTRA_OUTPUTS:
        p = os.path.join(out_dir, name)
        if os.path.isdir(p):
            shutil.rmtree(p)
        elif os.path.exists(p):
            os.remove(p)


class Log:
    def __init__(self, path, quiet=False):
        self.f = open(path, 'a')
        self.quiet = quiet

    def __call__(self, msg):
        line = f"[{datetime.datetime.now():%H:%M:%S}] {msg}"
        self.f.write(line + '\n')
        self.f.flush()
        if not self.quiet:
            print(line, flush=True)

    def close(self):
        self.f.close()


def _fmt(v, nd=4):
    if v is None:
        return ''
    if isinstance(v, (float, np.floating)):
        if not math.isfinite(v):
            return ''
        return round(float(v), nd)
    return v


def write_csv(path, rows, columns):
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, columns, extrasaction='ignore')
        w.writeheader()
        for r in rows:
            w.writerow({k: _fmt(r.get(k)) for k in columns})


def _git_commit():
    if getattr(sys, 'frozen', False):          # the packaged app: no git checkout (and no console windows)
        return None
    try:
        here = os.path.dirname(os.path.abspath(__file__))
        return subprocess.run(['git', '-C', here, 'rev-parse', '--short', 'HEAD'], capture_output=True,
                              text=True, timeout=5).stdout.strip() or None
    except Exception:
        return None


def plan_cases(inputs, table=None, priority=None, default_side=None, sides=None):
    """[(path, name, side, source, note)] -- laterality resolved per file,
    before anything is measured. sides: {name: (side, source)} decided
    elsewhere (e.g. entered in the review app); these take precedence."""
    plan = []
    sides = sides or {}
    for path in find_segmentations(inputs):
        name = case_name_from_path(path)
        if sides.get(name) and sides[name][0]:
            plan.append((path, name, sides[name][0], sides[name][1], None))
            continue
        try:
            header_side = header_laterality(nrrd.read_header(path))
        except Exception:
            header_side = None
        side, source, note = resolve(name, header_side, table, priority)
        if side is None and default_side:
            side, source, note = default_side, 'command line', note
        plan.append((path, name, side, source, note))
    return plan


def run_cohort(inputs, out_dir, laterality=None, id_columns=None, side_column=None, priority=None,
               carry=(), measure_keys=None, jobs=1, qc_grids=True, case_sheets='flagged',
               slicer_export=True, outlier_z=3.5, cohort_name=None, quiet=False, exclude=(),
               case_info=None, cache=True, rerun=False, reviews=None, sides=None, progress=None, extras=False):
    """exclude: case names to skip, or {name: reason}. case_info: {name:
    {column: value}} extra per-case columns (e.g. segmentation_source).
    cache: reuse per-case results whose segmentation file, side, options
    and code are unchanged (rerun=True ignores the cache). reviews: the
    review-file dict (review.py) whose landmark decisions are applied.
    sides: {name: (side, source)} that override every other laterality
    source. progress: optional callback(done, total, case_name). extras:
    also write the statistics / flags / landmarks / Slicer / summary files.
    Pulls the tidy measurements.csv + QC grids layout otherwise."""
    t_start = time.time()
    exclude_reasons = dict(exclude) if isinstance(exclude, dict) else {e: '' for e in (exclude or ())}
    exclude = set(exclude_reasons)
    case_info = case_info or {}
    info_cols = sorted({k for v in case_info.values() for k in v})
    os.makedirs(out_dir, exist_ok=True)
    _tidy(out_dir, extras)
    meta = os.path.join(out_dir, META_DIR)
    qc_dir = os.path.join(out_dir, 'qc')
    log = Log(os.path.join(meta, 'log.txt'), quiet=quiet)
    if not extras:
        case_sheets, slicer_export = 'none', False
    measures = [MEASURE_BY_KEY[k] for k in measure_keys] if measure_keys is not None else MEASURES
    cols = value_columns(measures)
    units = {c: u for m in measures for c, u, _ in m.columns}

    table = LateralityTable(laterality, id_columns=id_columns, side_column=side_column) if laterality else None
    if table:
        log(f"laterality table {laterality}: {len(table.rows)} rows, id columns {table.id_columns}, "
            f"side column {table.side_column!r}")
    plan = plan_cases(inputs, table, priority, sides=sides)
    log(f"{len(plan)} segmentation file(s) found in {', '.join(inputs)}")
    flags = []
    if exclude:
        skipped = [p[1] for p in plan if p[1] in exclude]
        plan = [p for p in plan if p[1] not in exclude]
        log(f"excluded {len(skipped)} case(s): {', '.join(skipped)}")
        for n in skipped:
            flags.append({'case': n, 'category': 'excluded', 'measurement': '',
                          'detail': exclude_reasons.get(n) or 'excluded by configuration'})

    by_row = {}
    for path, name, side, source, note in plan:
        if note:
            flags.append({'case': name, 'category': 'laterality', 'measurement': '', 'detail': note})
        if table is not None:
            row = table.lookup(name)[1]
            if row is not None:
                by_row.setdefault(id(row), []).append(name)
    for names in by_row.values():
        if len(names) > 1:
            for n in names:
                flags.append({'case': n, 'category': 'laterality', 'measurement': '',
                              'detail': f"several files match the same table row: {', '.join(names)}"})
    if table is not None:
        matched = set(by_row)
        excluded_rows = {id(table.lookup(n)[1]) for n in exclude}
        for row in table.rows:
            if id(row) not in matched and id(row) not in excluded_rows:
                ident = next((row.get(c) for c in table.id_columns if row.get(c)), '?')
                flags.append({'case': ident, 'category': 'missing segmentation', 'measurement': '',
                              'detail': 'in the laterality table but no segmentation file matched'})

    # ---- measure (reusing cached results for unchanged files)
    results = [None] * len(plan)
    shas = [file_sha1(p[0]) for p in plan]
    reviews = reviews or {}
    ovr = [review_store.overrides_for(reviews, p[1]) or None for p in plan]
    keys = [cache_key(h, p[2], measure_keys, qc_grids, slicer_export, o) for h, p, o in zip(shas, plan, ovr)]
    rcache = ResultCache(os.path.join(meta, 'cache')) if cache else None
    work = []
    for i, p in enumerate(plan):
        hit = None if (rerun or rcache is None) else rcache.get(shas[i], keys[i])
        if hit is not None:
            results[i] = dict(hit, path=p[0], name=p[1])
        else:
            work.append((i, p[0], p[2], ovr[i]))
    n_reused = len(plan) - len(work)
    if n_reused:
        log(f"{n_reused} case(s) unchanged since the last run -- reusing their results")
    done = 0

    def _report(i, r):
        nonlocal done
        done += 1
        if rcache is not None and not (r.get('error') or '').startswith('worker crashed'):
            rcache.put(shas[i], keys[i], r)
        ok = sum(1 for m in measures if r['notes'].get(m.key) is None and r['error'] is None)
        elapsed = time.time() - t_start
        eta = elapsed / done * (len(work) - done)
        if progress:
            progress(done, len(work), plan[i][1])
        log(f"[{done}/{len(work)}] {plan[i][1]} ({plan[i][2] or '?'}): {ok}/{len(measures)} measurements OK, "
            f"{len(r.get('warnings', []))} segmentation warning(s), {r.get('seconds', 0):.1f}s  (ETA {eta/60:.1f} min)")

    if jobs > 1 and len(work) > 1:
        with ProcessPoolExecutor(max_workers=jobs) as ex:
            futs = {ex.submit(measure_case, path, side, measure_keys, qc_grids, slicer_export, o): i
                    for i, path, side, o in work}
            for fut in as_completed(futs):
                i = futs[fut]
                try:
                    results[i] = fut.result()
                except Exception as e:
                    results[i] = {'name': plan[i][1], 'path': plan[i][0], 'knee_side': plan[i][2],
                                  'values': {c: None for c in cols}, 'notes': {}, 'validation': {},
                                  'warnings': [], 'payloads': {}, 'landmarks': [], 'error': f'worker crashed: {e}'}
                _report(i, results[i])
    else:
        for i, path, side, o in work:
            results[i] = measure_case(path, side, measure_keys, qc_grids, slicer_export, o)
            _report(i, results[i])
    if rcache is not None:
        rcache.prune(shas)

    # ---- tables
    carry = [c for c in (carry or []) if table is not None and c in table.fieldnames]
    m_rows, v_rows, lm_rows = [], [], []
    for (path, name, side, source, _), r, sha in zip(plan, results, shas):
        trow = table.lookup(name)[1] if table is not None else None
        row = {'case': r.get('name') or name, 'knee_side': r.get('knee_side'), 'laterality_source': source,
               **case_info.get(name, {})}
        for c in carry:
            row[c] = (trow or {}).get(c)
        row.update(r['values'])
        for m in measures:
            row[f'note_{m.key}'] = r['notes'].get(m.key)
        row['segmentation_warnings'] = len(r.get('warnings', []))
        row['tibia_length_mm'] = (r.get('validation') or {}).get('tibia_length_mm')
        notes = [f"{m.key}: {r['notes'][m.key]}" for m in measures if r['notes'].get(m.key)]
        row['notes'] = '; '.join(([r['error']] if r.get('error') else []) + notes
                                 + r.get('warnings', []) + r.get('cautions', []))
        row['error'] = r.get('error')
        row['file'] = path
        row['file_sha1'] = sha
        m_rows.append(row)

        v_rows.append({'case': row['case'], 'knee_side': row['knee_side'], **r.get('validation', {}),
                       'warnings': '; '.join(r.get('warnings', [])), 'n_warnings': len(r.get('warnings', []))})
        for measure, label, idx, lps in r.get('landmarks', []):
            lm_rows.append({'case': row['case'], 'knee_side': row['knee_side'], 'measurement': measure,
                            'landmark': label, 'x_index': idx[0], 'y_index': idx[1], 'z_index': idx[2],
                            'L_mm': lps[0], 'P_mm': lps[1], 'S_mm': lps[2]})

        if r.get('error'):
            flags.append({'case': row['case'], 'category': 'error', 'measurement': '', 'detail': r['error']})
        for w in r.get('warnings', []):
            flags.append({'case': row['case'], 'category': 'segmentation', 'measurement': '', 'detail': w})
        for w in r.get('cautions', []):
            flags.append({'case': row['case'], 'category': 'caution', 'measurement': 'tt_tg, ttl', 'detail': w})
        for m in measures:
            note = r['notes'].get(m.key)
            if note and not note.startswith('skipped'):
                failed = all(r['values'].get(c) is None for c, _, _ in m.columns)
                category = 'reviewed' if 'by reviewer' in note else ('failed' if failed else 'note')
                flags.append({'case': row['case'], 'category': category,
                              'measurement': m.key, 'detail': note})

    for f in stats.flag_outliers([r for r in m_rows if not r['error']], cols, threshold=outlier_z):
        flags.append({'case': f['case'], 'category': 'outlier', 'measurement': f['measurement'],
                      'detail': f"value {f['value']:.2f} vs cohort median {f['cohort_median']:.2f} "
                                f"(robust z {f['robust_z']:+.1f})", 'robust_z': f['robust_z']})
    for f in stats.flag_outliers(v_rows, VOLUME_COLUMNS, threshold=outlier_z):
        flags.append({'case': f['case'], 'category': 'segmentation outlier', 'measurement': f['measurement'],
                      'detail': f"{f['value']:.1f} cm3 vs cohort median {f['cohort_median']:.1f} "
                                f"(robust z {f['robust_z']:+.1f})", 'robust_z': f['robust_z']})

    if extras:
        m_cols = (['case', 'knee_side', 'laterality_source'] + info_cols + carry + cols
                  + ['TT_possibly_out_of_window', 'TT_review', 'TPCL_review', 'Sulcus_review', 'tibia_length_mm']
                  + [f'note_{m.key}' for m in measures] + ['segmentation_warnings', 'error', 'file', 'file_sha1'])
    else:
        m_cols = (['case', 'knee_side', 'laterality_source'] + info_cols + carry + cols
                  + ['TT_possibly_out_of_window', 'TT_review', 'TPCL_review', 'Sulcus_review', 'tibia_length_mm', 'notes', 'file'])
    changes = _archive_and_diff(out_dir, m_rows, cols) if extras else []
    write_csv(os.path.join(out_dir, 'measurements.csv'), m_rows, m_cols)
    measured = [r for r in m_rows if not r['error']]
    named = [dict(r, name=r.get('name') or p[1]) for p, r in zip(plan, results)]
    if qc_grids:
        _write_qc(named, flags, measures, qc_dir, case_sheets, log, jobs)
    info = _run_info(t_start, cohort_name, inputs, laterality, priority, exclude_reasons, plan, n_reused,
                     measures, jobs, outlier_z)
    with open(os.path.join(meta, 'run_info.json'), 'w') as f:
        json.dump(info, f, indent=2)
    if not extras:
        log(f"done: {len(plan)} cases in {(time.time() - t_start)/60:.1f} min -> {out_dir}")
        log.close()
        return m_rows, flags

    write_csv(os.path.join(out_dir, 'changes.csv'), changes, ['case', 'change', 'column', 'previous', 'current'])
    write_csv(os.path.join(out_dir, 'segmentation_check.csv'), v_rows, ['case', 'knee_side'] + VALIDATION_COLUMNS)
    write_csv(os.path.join(out_dir, 'landmarks.csv'), lm_rows,
              ['case', 'knee_side', 'measurement', 'landmark', 'x_index', 'y_index', 'z_index',
               'L_mm', 'P_mm', 'S_mm'])
    write_csv(os.path.join(out_dir, 'qc_flags.csv'), flags,
              ['case', 'category', 'measurement', 'detail', 'robust_z'])
    write_csv(os.path.join(out_dir, 'measurement_definitions.csv'),
              [{'column': c, 'unit': u, 'measurement': m.title, 'definition': d}
               for m in measures for c, u, d in m.columns], ['column', 'unit', 'measurement', 'definition'])
    desc = stats.describe(measured, cols, units)
    tt_cols = [c for c in cols if c in TT_DEPENDENT]
    in_window = [r for r in measured if not r.get('TT_possibly_out_of_window')]
    if tt_cols and len(in_window) < len(measured):
        for d in stats.describe(in_window, tt_cols, units):
            d['subgroup'] += ' (TT in window)'
            desc.append(d)
    write_csv(os.path.join(out_dir, 'statistics.csv'), desc,
              ['measurement', 'unit', 'subgroup', 'n_cases', 'n', 'n_missing', 'mean', 'sd', 'median',
               'q1', 'q3', 'iqr', 'min', 'max'])
    log('wrote measurements.csv, statistics.csv, qc_flags.csv, segmentation_check.csv, landmarks.csv')

    os.makedirs(os.path.join(out_dir, 'plots'), exist_ok=True)
    outlier_flags = [f for f in flags if f['category'] == 'outlier']
    plots.distribution_figure(measured, cols, units, outlier_flags,
                              os.path.join(out_dir, 'plots', 'distributions.png'),
                              title=f"{cohort_name or 'Cohort'}: measurement distributions "
                                    f"(red diamond = flagged outlier, |robust z| > {outlier_z})")

    if slicer_export:
        sdir = os.path.join(out_dir, 'slicer')
        os.makedirs(sdir, exist_ok=True)
        for old in glob.glob(os.path.join(sdir, '*.mrk.json')):
            os.remove(old)
        for r in named:
            if r.get('landmarks'):
                slicer.write_markups(os.path.join(sdir, f"{r['name']}.mrk.json"), r['name'], r['landmarks'])
        log('wrote Slicer markups to slicer/')

    with open(os.path.join(out_dir, 'run_info.json'), 'w') as f:
        json.dump(info, f, indent=2)

    _write_summary(os.path.join(out_dir, 'summary.md'), cohort_name, plan, m_rows, measures, desc, flags, info,
                   changes, info_cols)
    log(f"done: {len(plan)} cases in {(time.time() - t_start)/60:.1f} min -> {out_dir}")
    log.close()
    return m_rows, flags


def _write_qc(named, flags, measures, qc_dir, case_sheets, log, jobs=1):
    """QC grids, one per measurement, drawn in parallel; pages whose cases
    didn't change are kept (see qc.build_grid)."""
    os.makedirs(qc_dir, exist_ok=True)
    for old in glob.glob(os.path.join(qc_dir, 'cases', '*.png')):
        os.remove(old)
    keys = ['overview'] + [m.key for m in measures]
    tasks = [([qc.grid_inputs(r, key) for r in named], key, qc_dir) for key in keys]
    if jobs > 1:
        with ProcessPoolExecutor(max_workers=min(jobs, len(tasks))) as ex:
            written = list(ex.map(qc.build_grid, *zip(*tasks)))
    else:
        written = [qc.build_grid(*t) for t in tasks]
    keep = {os.path.basename(p) for paths in written for p in paths}
    for old in glob.glob(os.path.join(qc_dir, '*_grid*.png')):          # pages no longer produced
        if os.path.basename(old) not in keep:
            os.remove(old)
            if os.path.exists(os.path.join(qc_dir, f'.{os.path.basename(old)}.key')):
                os.remove(os.path.join(qc_dir, f'.{os.path.basename(old)}.key'))
    for paths in written:
        log(f"QC grid: {', '.join(os.path.basename(p) for p in paths)}")
    if case_sheets != 'none':
        flagged_cases = {f['case'] for f in flags if f['category'] not in ('missing segmentation',)}
        sheet_dir = os.path.join(qc_dir, 'cases')
        os.makedirs(sheet_dir, exist_ok=True)
        n_sheets = 0
        for r in named:
            if case_sheets == 'all' or r['name'] in flagged_cases:
                qc.build_case_sheet(r, os.path.join(sheet_dir, f"{r['name']}.png"),
                                    ['overview'] + [m.key for m in measures])
                n_sheets += 1
        log(f"{n_sheets} per-case QC sheet(s) in qc/cases/")


def _run_info(t_start, cohort_name, inputs, laterality, priority, exclude_reasons, plan, n_reused, measures,
              jobs, outlier_z):
    return {
        'kneemeasure_version': __version__, 'git_commit': _git_commit(),
        'run_started': datetime.datetime.fromtimestamp(t_start).isoformat(timespec='seconds'),
        'run_seconds': round(time.time() - t_start, 1), 'cohort_name': cohort_name,
        'inputs': inputs, 'laterality_table': laterality, 'laterality_priority': priority,
        'excluded_cases': exclude_reasons,
        'n_cases': len(plan), 'n_reused_from_cache': n_reused, 'code_fingerprint': code_fingerprint(),
        'measures': [m.key for m in measures], 'jobs': jobs,
        'parameters': {'KRA_SIGN': kra.KRA_SIGN, 'STT_TG_SIGN': tt_tg.STT_TG_SIGN,
                       'tibial_pcl_plateau_area_fraction': tibial_pcl.find_tibial_pcl.__defaults__[0],
                       'reference_slice_min_depth_mm': 2.0, 'outlier_robust_z': outlier_z},
        'python': sys.version.split()[0], 'platform': platform.platform(),
        'command': ' '.join(sys.argv),
    }


def _archive_and_diff(out_dir, m_rows, cols, tol=0.01):
    """Copy the previous measurements.csv into history/ and list what
    changed: added / removed cases, changed segmentation files, and values
    that moved by more than tol (or appeared / disappeared)."""
    prev_path = os.path.join(out_dir, 'measurements.csv')
    if not os.path.exists(prev_path):
        return []
    stamp = 'previous'
    try:
        with open(os.path.join(out_dir, 'run_info.json')) as f:
            stamp = json.load(f).get('run_started', stamp).replace(':', '')
    except Exception:
        pass
    os.makedirs(os.path.join(out_dir, 'history'), exist_ok=True)
    shutil.copy2(prev_path, os.path.join(out_dir, 'history', f'measurements_{stamp}.csv'))
    with open(prev_path, newline='') as f:
        prev = {r['case']: r for r in csv.DictReader(f)}
    cur = {r['case']: r for r in m_rows}
    out = []
    for case in sorted(set(cur) - set(prev)):
        out.append({'case': case, 'change': 'added'})
    for case in sorted(set(prev) - set(cur)):
        out.append({'case': case, 'change': 'removed'})
    for case in sorted(set(cur) & set(prev)):
        p, c = prev[case], cur[case]
        if p.get('file_sha1') and p['file_sha1'] != c.get('file_sha1'):
            out.append({'case': case, 'change': 'segmentation file changed'})
        if (p.get('knee_side') or None) != (c.get('knee_side') or None):
            out.append({'case': case, 'change': 'knee side changed', 'column': 'knee_side',
                        'previous': p.get('knee_side'), 'current': c.get('knee_side')})
        for col in cols:
            a = p.get(col) or None
            b = c.get(col)
            if a is None and b is None:
                continue
            if a is None or b is None or abs(float(a) - float(b)) > tol:
                out.append({'case': case, 'change': 'value changed', 'column': col, 'previous': a,
                            'current': b})
    return out


def _describe_inputs(inputs):
    """Folders (with file counts) instead of a long list of files."""
    groups = {}
    for p in inputs:
        d = p if os.path.isdir(p) else os.path.dirname(p)
        groups[d] = groups.get(d, 0) + (0 if os.path.isdir(p) else 1)
    return '; '.join(f'{d} ({n} files)' if n else d for d, n in groups.items())


def _write_summary(path, cohort_name, plan, m_rows, measures, desc, flags, info, changes=(), info_cols=()):
    L = [f"# {cohort_name or 'Knee measurement run'}", '',
         f"- Run: {info['run_started']} ({info['run_seconds'] / 60:.1f} min), kneemeasure {info['kneemeasure_version']}"
         + (f", repo commit {info['git_commit']}" if info['git_commit'] else ''),
         f"- Inputs: {_describe_inputs(info['inputs'])}",
         f"- Laterality: {info['laterality_table'] or 'NRRD headers only'}",
         f"- Cases: {len(plan)}"
         + (f" ({info['n_reused_from_cache']} unchanged since the last run, reused)" if info.get('n_reused_from_cache') else ''),
         ]
    for col in info_cols:
        counts = {}
        for r in m_rows:
            counts[r.get(col)] = counts.get(r.get(col), 0) + 1
        L.append(f"- {col}: " + ', '.join(f"{k}: {v}" for k, v in sorted(counts.items(), key=str)))
    if changes:
        kinds = {}
        for c in changes:
            kinds[c['change']] = kinds.get(c['change'], 0) + 1
        L.append('- Changes since the previous run (changes.csv): '
                 + ', '.join(f"{v} {k}" for k, v in kinds.items()))
    elif os.path.exists(os.path.join(os.path.dirname(path), 'history')):
        L.append('- No changes since the previous run.')
    L += ['', '## Success per measurement', '',
         '| measurement | measured | failed / skipped |', '|---|---|---|']
    for m in measures:
        first = m.columns[0][0]
        ok = sum(1 for r in m_rows if r.get(first) is not None)
        L.append(f"| {m.title} | {ok} | {len(m_rows) - ok} |")
    L += ['', '## Descriptive statistics (all cases)', '',
          '| column | unit | n | mean | SD | median | IQR | min | max |', '|---|---|---|---|---|---|---|---|---|']

    def f(v, nd=2):
        return '' if v is None else f'{v:.{nd}f}'
    for d in desc:
        if d['subgroup'] == 'all':
            L.append(f"| {d['measurement']} | {d['unit']} | {d['n']} | {f(d['mean'])} | {f(d['sd'])} | "
                     f"{f(d['median'])} | {f(d['iqr'])} | {f(d['min'])} | {f(d['max'])} |")
    cats = {}
    for fl in flags:
        cats.setdefault(fl['category'], []).append(fl)
    L += ['', '## QC flags', '']
    if not flags:
        L.append('None.')
    for cat in ('error', 'excluded', 'laterality', 'missing segmentation', 'segmentation', 'segmentation outlier',
                'failed', 'outlier', 'caution', 'reviewed', 'note'):
        if cat not in cats:
            continue
        L += [f"### {cat} ({len(cats[cat])})", '']
        for fl in cats[cat]:
            meas = f" [{fl['measurement']}]" if fl.get('measurement') else ''
            L.append(f"- **{fl['case']}**{meas}: {fl['detail']}")
        L.append('')
    L += ['## Caveats', '',
          '- Automated values have not yet been validated against manual measurements in this cohort.',
          '- KRA sign convention (+ = tibial external rotation) has not been checked against manual KRA.',
          '- Review qc/ grids (and qc/cases/ for flagged cases) before using the numbers.']
    with open(path, 'w') as fh:
        fh.write('\n'.join(L) + '\n')


def compare_runs(groups, out_dir, measure_keys=None, tt_in_window_only=True):
    """groups: {name: measurements.csv path or run folder}. With
    tt_in_window_only, TT-TG / sTT-TG / TTL values of cases flagged as a
    likely out-of-window tuberosity are left out of the comparison."""
    os.makedirs(out_dir, exist_ok=True)
    measures = [MEASURE_BY_KEY[k] for k in measure_keys] if measure_keys is not None else MEASURES
    cols = value_columns(measures)
    units = {c: u for m in measures for c, u, _ in m.columns}
    data = {}
    for name, p in groups.items():
        path = os.path.join(p, 'measurements.csv') if os.path.isdir(p) else p
        with open(path, newline='') as fh:
            data[name] = [r for r in csv.DictReader(fh) if not r.get('error')]
        if tt_in_window_only:
            for r in data[name]:
                if r.get('TT_possibly_out_of_window') == '1':
                    for c in TT_DEPENDENT:
                        r[c] = ''
    table, pairwise = stats.compare_groups(data, cols, units)
    names = list(data)
    tcols = ['measurement', 'unit']
    for g in names:
        tcols += [f'{g}_n', f'{g}_mean', f'{g}_sd', f'{g}_median', f'{g}_iqr']
    if len(names) == 2:
        tcols += ['mean_difference', 'welch_t_p', 'welch_t_p_holm', 'mann_whitney_p', 'mann_whitney_p_holm',
                  'cohens_d', 'hedges_g']
    else:
        tcols += ['anova_p', 'anova_p_holm', 'kruskal_p', 'kruskal_p_holm']
    write_csv(os.path.join(out_dir, 'comparison.csv'), table, tcols)
    if pairwise:
        write_csv(os.path.join(out_dir, 'comparison_pairwise.csv'), pairwise,
                  ['measurement', 'group_a', 'group_b', 'mean_a', 'mean_b', 'welch_t_p', 'welch_t_p_holm_pairs',
                   'mann_whitney_p', 'mann_whitney_p_holm_pairs', 'cohens_d', 'hedges_g'])
    plots.comparison_figure(data, cols, units, os.path.join(out_dir, 'comparison.png'))
    with open(os.path.join(out_dir, 'comparison_info.json'), 'w') as fh:
        json.dump({'groups': groups, 'n': {g: len(v) for g, v in data.items()},
                   'tt_dependent_columns_exclude_out_of_window': tt_in_window_only,
                   'created': datetime.datetime.now().isoformat(timespec='seconds')}, fh, indent=2)
    return table, pairwise


def validate_run(run, manual, id_column, pairs, out_dir, side_column=None, exclude=()):
    """Compare a run's measurements.csv against a manual sheet.
    pairs: [(auto_column, manual_column), ...]. TT-dependent measurements
    are also reported with likely out-of-window tuberosity cases excluded."""
    from . import agreement
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(run, 'measurements.csv') if os.path.isdir(run) else run
    with open(path, newline='') as fh:
        auto_rows = [r for r in csv.DictReader(fh) if not r.get('error') and r['case'] not in set(exclude)]
    index = agreement.load_manual(manual, id_column, side_column)
    sheet_cols = list(next(iter(index.values()))[0]) if index else []
    by_stripped = {c.strip().lower(): c for c in sheet_cols}
    resolved = []
    for auto_col, man_col in pairs:
        actual = man_col if man_col in sheet_cols else by_stripped.get(man_col.strip().lower())
        if actual is None:
            raise ValueError(f'manual sheet has no column {man_col!r}')
        resolved.append((auto_col, actual))
    pairs = resolved

    matched, per_case = [], []
    for r in auto_rows:
        row, how = agreement.match(index, r['case'], r.get('knee_side'), side_column)
        per_case.append({'case': r['case'], 'knee_side': r.get('knee_side'), 'matched': how})
        if row is not None:
            matched.append((r, row))

    results, case_rows, plot_pairs = [], [], {}
    for auto_col, man_col in pairs:
        subsets = [('all', matched)]
        if auto_col in TT_DEPENDENT:
            subsets.append(('TT in window', [(a, m) for a, m in matched if a.get('TT_possibly_out_of_window') != '1']))
        for subset, rows in subsets:
            xy = []
            for a, m in rows:
                x, y = agreement._num(a.get(auto_col)), agreement._num(m.get(man_col))
                if x is not None and y is not None:
                    xy.append((x, y, a['case']))
                    if subset == 'all':
                        case_rows.append({'case': a['case'], 'knee_side': a.get('knee_side'),
                                          'automated_column': auto_col, 'manual_column': man_col,
                                          'automated': x, 'manual': y, 'difference': x - y,
                                          'TT_possibly_out_of_window': a.get('TT_possibly_out_of_window')})
            res = agreement.agreement(xy)
            results.append({'automated_column': auto_col, 'manual_column': man_col, 'subset': subset, **res})
            label = auto_col if sum(a == auto_col for a, _ in pairs) == 1 else f'{auto_col} vs {man_col[:24]}'
            plot_pairs[f'{label} ({subset})' if len(subsets) > 1 else label] = xy

    cols = ['automated_column', 'manual_column', 'subset', 'n', 'auto_mean', 'manual_mean', 'bias', 'sd_diff',
            'loa_low', 'loa_high', 'mae', 'pearson_r', 'pearson_p', 'spearman_rho', 'icc_2_1', 'icc_3_1',
            'icc_2_1_if_sign_flipped', 'sign_check']
    write_csv(os.path.join(out_dir, 'agreement.csv'), results, cols)
    write_csv(os.path.join(out_dir, 'agreement_cases.csv'), case_rows,
              ['case', 'knee_side', 'automated_column', 'manual_column', 'automated', 'manual', 'difference',
               'TT_possibly_out_of_window'])
    write_csv(os.path.join(out_dir, 'agreement_matching.csv'), per_case, ['case', 'knee_side', 'matched'])
    plots.agreement_figure(plot_pairs, os.path.join(out_dir, 'agreement.png'))
    return results, per_case
