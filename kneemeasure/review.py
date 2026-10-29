"""
Human input: review decisions for automatically placed landmarks, and knee
sides entered in the review app.

Stored as one JSON file per pull (<pull>/kneemeasure_review.json), kept
OUTSIDE the regenerable results folder so deleting or rebuilding results
never loses review work. Every measurement run applies it, so the command
line and the review app always give the same numbers.

Format:

    {"version": 1,
     "cases": {
       "N31XL": {
         "laterality": {"side": "L", "time": "..."},   # only if entered in the app
         "tpcl": {"status": "replaced", "method": "slice", "z": 31, "theta": 4.2,
                  "medial_pt": [..], "lateral_pt": [..], "auto": {"z": 29, "theta": 1.1}, ...},
         "sulcus": {"status": "replaced", "method": "points", "z": 88, "sulcus": [x, y],
                    "peaks": [[x, y], [x, y]], "auto": {"z": 90, "sulcus": [x, y]}, ...},
         "tuberosity": {
           "status": "replaced",            # accepted | replaced | not_in_view | misaligned
           "point": [212.5, 98.25, 14],     # reoriented voxel index (x, y, z; x, y can be fractional when
                                            #   clicked on the finer axial series); only for accepted / replaced
           "source": "axial",               # only when clicked on the axial series (z may then be negative:
                                            #   below the segmentation, which TT-TG / TTL allow -- they use x, y)
           "point_lps_mm": [..., ..., ...], # same point in scanner LPS mm, for reference
           "auto_point": [205, 101, 6],     # what the algorithm proposed when reviewed
           "reviewer": "BL", "time": "2026-10-28T22:40:00", "note": "",
           "file_sha1": "..."               # segmentation the review was made on
         }}}}

accepted and replaced both pin the reviewed point, so later code changes
can't silently move a point a human already checked. misaligned: the
tuberosity is visible on the axial series, but the segmentation doesn't sit
on the bones there (the patient moved between the two series), so a point
clicked on the axial can't be combined with the segmentation's landmarks;
TT-TG / sTT-TG / TTL are left blank, to be measured manually.
"""
import datetime
import json
import os
import shutil
import tempfile

TUBEROSITY_STATUSES = ('accepted', 'replaced', 'not_in_view', 'misaligned')
REVIEW_FILE = 'kneemeasure_review.json'


def load(path):
    if not path or not os.path.exists(path):
        return {'version': 1, 'cases': {}}
    with open(path) as f:
        data = json.load(f)
    data.setdefault('cases', {})
    return data


def save(path, data):
    """Atomic write, so a crash mid-save can't corrupt review work."""
    d = os.path.dirname(os.path.abspath(path))
    fd, tmp = tempfile.mkstemp(dir=d, prefix='.review_', suffix='.json')
    with os.fdopen(fd, 'w') as f:
        json.dump(data, f, indent=1, sort_keys=True)
    os.replace(tmp, path)


def set_tuberosity(data, case, status, point=None, auto_point=None, reviewer='', note='',
                   file_sha1=None, point_lps_mm=None, source=None):
    if status not in TUBEROSITY_STATUSES:
        raise ValueError(f'status must be one of {TUBEROSITY_STATUSES}')
    if status in ('accepted', 'replaced') and point is None:
        raise ValueError(f'{status} needs a point')
    entry = {'status': status, 'reviewer': reviewer, 'note': note,
             'time': datetime.datetime.now().isoformat(timespec='seconds')}
    has_point = status in ('accepted', 'replaced')
    if point is not None and has_point:
        entry['point'] = [round(float(point[0]), 2), round(float(point[1]), 2), int(point[2])]  # x, y may be sub-voxel
    if point_lps_mm is not None and has_point:
        entry['point_lps_mm'] = [round(float(v), 2) for v in point_lps_mm]
    if source and has_point:
        entry['source'] = source
    if auto_point is not None:
        entry['auto_point'] = [int(v) for v in auto_point]
    if file_sha1:
        entry['file_sha1'] = file_sha1
    data.setdefault('cases', {}).setdefault(case, {})['tuberosity'] = entry
    return entry


def set_laterality(data, case, side):
    """Knee side entered by a person in the app (overrides headers/sheets)."""
    data.setdefault('cases', {}).setdefault(case, {})['laterality'] = {
        'side': side, 'time': datetime.datetime.now().isoformat(timespec='seconds')}


def overrides_for(data, case):
    """The landmark overrides for one case, as passed to measure_case (the
    laterality entry is not a landmark and is left out)."""
    c = dict(data.get('cases', {}).get(case, {}) or {})
    c.pop('laterality', None)
    return c


TPCL_STATUSES = ('accepted', 'replaced')


def set_tpcl(data, case, status, line, auto=None, method='auto', reviewer='', note='', file_sha1=None):
    """line: {'z', 'medial_pt', 'lateral_pt', 'theta'} -- pinned so later code
    changes can't move a line a person already checked. method: auto (the
    automatic line accepted), slice (recomputed on a chosen slice) or points
    (two clicked points)."""
    if status not in TPCL_STATUSES:
        raise ValueError(f'status must be one of {TPCL_STATUSES}')
    entry = {'status': status, 'method': method, 'z': int(line['z']), 'theta': round(float(line['theta']), 4),
             'medial_pt': [round(float(v), 2) for v in line['medial_pt']],
             'lateral_pt': [round(float(v), 2) for v in line['lateral_pt']],
             'reviewer': reviewer, 'note': note, 'time': datetime.datetime.now().isoformat(timespec='seconds')}
    if auto:
        entry['auto'] = {'z': int(auto['z']), 'theta': round(float(auto['theta']), 4)}
    if file_sha1:
        entry['file_sha1'] = file_sha1
    data.setdefault('cases', {}).setdefault(case, {})['tpcl'] = entry
    return entry


SULCUS_STATUSES = ('accepted', 'replaced')


def set_sulcus(data, case, status, ref, auto=None, method='auto', reviewer='', note='', file_sha1=None):
    """ref: {'z', 'sulcus', 'peak1', 'peak2'} -- the trochlear reference
    slice's sulcus and two facet peaks, pinned like the TPCL. method: auto,
    slice (re-detected on a chosen slice) or points (three clicked points)."""
    if status not in SULCUS_STATUSES:
        raise ValueError(f'status must be one of {SULCUS_STATUSES}')
    entry = {'status': status, 'method': method, 'z': int(ref['z']),
             'sulcus': [round(float(v), 2) for v in ref['sulcus']],
             'peaks': [[round(float(v), 2) for v in ref[k]] for k in ('peak1', 'peak2')],
             'reviewer': reviewer, 'note': note, 'time': datetime.datetime.now().isoformat(timespec='seconds')}
    if auto:
        entry['auto'] = {'z': int(auto['z']), 'sulcus': [int(v) for v in auto['sulcus']]}
    if file_sha1:
        entry['file_sha1'] = file_sha1
    data.setdefault('cases', {}).setdefault(case, {})['sulcus'] = entry
    return entry


def clear_review(data, case, key):
    """Remove one landmark's decision ('tuberosity', 'tpcl' or 'sulcus')."""
    c = data.get('cases', {}).get(case)
    if c:
        c.pop(key, None)
        if not c:
            data['cases'].pop(case)


LANDMARK_KEYS = ('tuberosity', 'tpcl', 'sulcus')


def start_over(path):
    """Remove every landmark decision in a review file (knee sides entered
    in the app are kept), after copying the file to
    kneemeasure_review.backup-<time>.json next to it. Returns (number of
    decisions removed, backup path or None when there was nothing to remove)."""
    data = load(path)
    n = sum(1 for c in data['cases'].values() for k in LANDMARK_KEYS if k in c)
    if not n:
        return 0, None
    stamp = datetime.datetime.now().strftime('%Y-%m-%d_%H%M%S')
    backup = os.path.join(os.path.dirname(os.path.abspath(path)), f'kneemeasure_review.backup-{stamp}.json')
    shutil.copy2(path, backup)
    for case in list(data['cases']):
        for k in LANDMARK_KEYS:
            clear_review(data, case, k)
    save(path, data)
    return n, backup


def clear_tuberosity(data, case):
    c = data.get('cases', {}).get(case)
    if c:
        c.pop('tuberosity', None)
        if not c:
            data['cases'].pop(case)
