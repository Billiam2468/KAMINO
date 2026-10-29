"""
A review/measurement project: a segmentations folder, optionally a scans
folder, and where results and human input are kept.

    <root>/kneemeasure_review.json   human input: knee sides entered in the app and landmark
                                     reviews (never deleted by a rebuild)

The app's results (measurements, QC grids, the rerun cache) are kept in the
user's cache folder, one folder per <root> (see results_dir), not next to the
scans: they're regenerable, and the app offers them as downloads. The
command line's `pull` / `run` still write <pull>/kneemeasure/ (pulls.py).

<root> is the scan pull folder when the segmentations live inside one (an
ancestor folder containing sorted/), otherwise the segmentations folder's
parent.

Knee side per scan, first found wins:
  1. entered in the app (stored in kneemeasure_review.json)
  2. the segmentation's NRRD header (laterality:=L)
  3. the scan's NIfTI description (laterality=L)
  4. a mapping sheet next to the folders (sorted/mrn_code_mapping.csv), if any
Scans with none of these are asked for in the app.
"""
import csv
import hashlib
import json
import os
import sys

import nrrd

from . import review as review_store
from .laterality import LateralityTable
from .scans import find_scan, laterality_from_nifti
from .volume import header_laterality, case_name_from_path, parse_side

RECENT_FILE = os.path.join(os.path.expanduser('~'), '.config', 'kneemeasure', 'recent.json')


def cache_root():
    """The user's cache folder for kneemeasure (regenerable files only)."""
    home = os.path.expanduser('~')
    if sys.platform == 'win32':
        base = os.environ.get('LOCALAPPDATA') or os.path.join(home, 'AppData', 'Local')
    elif sys.platform == 'darwin':
        base = os.path.join(home, 'Library', 'Caches')
    else:
        base = os.environ.get('XDG_CACHE_HOME') or os.path.join(home, '.cache')
    return os.path.join(base, 'kneemeasure')


def results_dir(root):
    """Where the app keeps one project's results: <cache>/projects/<name>-<hash of the path>."""
    root = os.path.abspath(root)
    tag = hashlib.sha1(root.encode('utf-8')).hexdigest()[:10]
    return os.path.join(cache_root(), 'projects', f'{os.path.basename(root) or "root"}-{tag}')
SHEET_NAMES = ('mrn_code_mapping.csv',)


def pull_root(path):
    """The scan-pull folder containing path (has sorted/), or None."""
    p = os.path.abspath(path)
    for _ in range(4):
        if os.path.isdir(os.path.join(p, 'sorted')):
            return p
        parent = os.path.dirname(p)
        if parent == p:
            break
        p = parent
    return None


def suggest_scans_dir(segs_dir):
    root = pull_root(segs_dir)
    if root and os.path.isdir(os.path.join(root, 'sorted')):
        return os.path.join(root, 'sorted')
    return None


def find_sheet(*dirs):
    for d in dirs:
        if not d:
            continue
        for cand in (d, os.path.join(d, 'sorted')):
            for name in SHEET_NAMES:
                p = os.path.join(cand, name)
                if os.path.exists(p):
                    return p
    return None


class Project:
    def __init__(self, segs_dir=None, scans_dir=None, files=None, root=None, name=None, sheet='auto'):
        """Either segs_dir (every .nrrd in it) or files ({code: path}, e.g.
        a pull's mixed cleaned/model selection) must be given."""
        if files is None:
            if not segs_dir or not os.path.isdir(segs_dir):
                raise ValueError(f'not a folder: {segs_dir}')
            files = {case_name_from_path(f): os.path.join(segs_dir, f) for f in sorted(os.listdir(segs_dir))
                     if f.lower().endswith('.nrrd') and not f.startswith('.')}
            if not files:
                raise ValueError(f'no .nrrd segmentations in {segs_dir}')
        self.files = dict(sorted(files.items()))
        first_dir = segs_dir or os.path.dirname(next(iter(self.files.values())))
        self.segs_dir = segs_dir or first_dir
        self.root = root or pull_root(first_dir) or os.path.dirname(os.path.abspath(first_dir))
        self.name = name or os.path.basename(self.root)
        self.scans_dir = scans_dir or None
        self.out = results_dir(self.root)
        self.review_file = os.path.join(self.root, review_store.REVIEW_FILE)
        self.sheet = find_sheet(self.root, self.scans_dir) if sheet == 'auto' else sheet
        self.scans = {c: find_scan(self.scans_dir, c) for c in self.files}
        self.refresh_sides()

    def refresh_sides(self):
        table = None
        if self.sheet:
            try:
                table = LateralityTable(self.sheet)
            except Exception:
                table = None
        entered = {c: (v.get('laterality') or {}).get('side')
                   for c, v in review_store.load(self.review_file)['cases'].items()}
        self.sides = {}
        for code, path in self.files.items():
            side, source = entered.get(code), 'entered in app'
            if not side:
                try:
                    side, source = header_laterality(nrrd.read_header(path)), 'segmentation header'
                except Exception:
                    side = None
            if not side and self.scans.get(code):
                side, source = laterality_from_nifti(self.scans[code]), 'scan header'
            if not side and table is not None:
                side, source = table.lookup(code)[0], 'mapping sheet'
            self.sides[code] = (side, source if side else None)

    def set_side(self, code, side):
        side = parse_side(side)
        if code not in self.files or side is None:
            raise ValueError('unknown scan or side')
        data = review_store.load(self.review_file)
        review_store.set_laterality(data, code, side)
        review_store.save(self.review_file, data)
        self.refresh_sides()

    def table(self, reviews=None):
        reviews = reviews if reviews is not None else review_store.load(self.review_file)['cases']
        rows = []
        for code in self.files:
            side, source = self.sides[code]
            t = (reviews.get(code) or {}).get('tuberosity')
            rows.append({'case': code, 'side': side, 'side_source': source, 'has_scan': bool(self.scans[code]),
                         'review': t['status'] if t else None})
        return rows

    def to_measure(self):
        """The scans a run would measure: all of them before the first run,
        then those added or whose segmentation changed since. (Review
        decisions are applied by rebuilding the results, not counted here.)"""
        m = os.path.join(self.out, 'measurements.csv')
        if not os.path.exists(m):
            return list(self.files)
        with open(m, newline='') as f:
            done = {r.get('case') for r in csv.DictReader(f)}
        t = os.path.getmtime(m)
        return [c for c, p in self.files.items() if c not in done or os.path.getmtime(p) > t]

    def summary(self):
        return {'name': self.name, 'root': self.root, 'segs_dir': self.segs_dir, 'scans_dir': self.scans_dir,
                'out': self.out, 'sheet': self.sheet, 'n_cases': len(self.files),
                'n_scans': sum(1 for v in self.scans.values() if v),
                'n_missing_side': sum(1 for s, _ in self.sides.values() if not s),
                'measured': os.path.exists(os.path.join(self.out, 'measurements.csv')),
                'n_to_measure': len(self.to_measure())}


# ---- folder browsing (local machine only)

def list_dir(path):
    path = os.path.abspath(os.path.expanduser(path or '~'))
    if not os.path.isdir(path):
        raise ValueError(f'not a folder: {path}')
    dirs = []
    try:
        names = sorted(os.listdir(path), key=str.lower)
    except PermissionError:
        names = []
    for n in names:
        p = os.path.join(path, n)
        if n.startswith('.') or not os.path.isdir(p):
            continue
        dirs.append({'name': n, 'path': p})
    files = [n for n in names if os.path.isfile(os.path.join(path, n))]
    parent = os.path.dirname(path)
    return {'path': path, 'parent': parent if parent != path else None, 'dirs': dirs,
            'n_nrrd': sum(1 for f in files if f.lower().endswith('.nrrd')),
            'n_nii': sum(1 for f in files if f.lower().endswith(('.nii', '.nii.gz'))),
            'n_scan_folders': sum(1 for d in dirs[:400] if any(os.path.exists(os.path.join(d['path'], s))
                                                              for s in ('water.nii.gz', 'water.nii')))}


def shortcuts():
    out = [{'name': 'Home', 'path': os.path.expanduser('~')}]
    for base in ('/run/media', '/media', '/mnt', '/Volumes'):
        if os.path.isdir(base):
            for d in sorted(os.listdir(base)):
                p = os.path.join(base, d)
                if os.path.isdir(p):
                    subs = [os.path.join(p, s) for s in sorted(os.listdir(p))] if base in ('/run/media', '/media') else [p]
                    for s in subs:
                        if os.path.isdir(s):
                            out.append({'name': os.path.basename(s), 'path': s})
    return out


def load_recent(limit=None):
    """Recently opened folders, newest first; with limit, only up to that
    many whose segmentations folder still exists (for the setup page)."""
    try:
        with open(RECENT_FILE) as f:
            recent = json.load(f)
    except Exception:
        return []
    if limit is None:
        return recent
    return [r for r in recent if os.path.isdir(r.get('segs_dir') or '')][:limit]


def remember(segs_dir, scans_dir):
    recent = [r for r in load_recent() if r.get('segs_dir') != segs_dir]
    recent.insert(0, {'segs_dir': segs_dir, 'scans_dir': scans_dir})
    try:
        os.makedirs(os.path.dirname(RECENT_FILE), exist_ok=True)
        with open(RECENT_FILE, 'w') as f:
            json.dump(recent[:10], f, indent=1)
    except Exception:
        pass
