"""
Per-case orchestration: load one segmentation, check it, compute every
measurement from shared landmarks, and capture the small QC payloads and
landmark coordinates the cohort runner needs -- so the full volume can be
dropped as soon as the case is done.
"""
import time
import traceback

from . import kra, tt_tg, ptltr, ttl, bump_height, cdi, trochlear_morphology, validation
from .landmarks import Landmarks
from .volume import load_case, header_laterality, index_to_lps


class Measure:
    def __init__(self, key, title, module, needs_side=True):
        self.key, self.title, self.module, self.needs_side = key, title, module, needs_side

    @property
    def columns(self):
        return self.module.COLUMNS


# Order here = column order in measurements.csv and QC grid order.
MEASURES = [
    Measure('kra', 'Knee rotation angle (KRA)', kra),
    Measure('tt_tg', 'TT-TG / sTT-TG', tt_tg),
    Measure('ptltr', 'PTLTR', ptltr),
    Measure('ttl', 'Tibial tubercle lateralization (TTL)', ttl),
    Measure('trochlear', 'Trochlear morphology (LTI, sulcus angle, facets)', trochlear_morphology),
    Measure('bump_height', 'Trochlear bump height', bump_height, needs_side=False),
    Measure('cdi', 'Caton-Deschamps index (CDI)', cdi, needs_side=False),
]
MEASURE_BY_KEY = {m.key: m for m in MEASURES}


def value_columns(measures=MEASURES):
    return [c for m in measures for c, _, _ in m.columns]


def measure_case(path, knee_side=None, measure_keys=None, qc=True, landmarks=True, overrides=None):
    """Run the whole suite on one segmentation file.

    knee_side: 'L'/'R', or None to use the file's 'laterality' header
    field. overrides: reviewed landmarks for this case (review.py). Side-dependent measurements are skipped (with a note) if no
    side is known. Never raises: a failure anywhere is recorded in the
    result's notes/'error'.

    Returns a plain dict (picklable, no volumes):
        name, path, knee_side, values {column: value}, notes {measure: note},
        validation {column: value}, warnings [segmentation problems],
        cautions [measurement-level concerns], payloads {measure: ..},
        landmarks [(measure, label, (x, y, z), (L, P, S))], error, seconds
    """
    t0 = time.time()
    measures = [MEASURE_BY_KEY[k] for k in measure_keys] if measure_keys is not None else MEASURES
    out = {'name': None, 'path': path, 'knee_side': knee_side, 'values': {}, 'notes': {},
           'validation': {}, 'warnings': [], 'cautions': [], 'payloads': {}, 'landmarks': [], 'error': None}
    try:
        case = load_case(path)
    except Exception as e:
        out['error'] = f'could not load segmentation: {e}'
        out['seconds'] = time.time() - t0
        return out
    out['name'] = case['name']
    if knee_side is None:
        knee_side = header_laterality(case['header'])
        out['knee_side'] = knee_side

    try:
        out['validation'], out['warnings'] = validation.check_case(case)
        if qc:
            out['payloads']['overview'] = validation.overview_payload(case)
    except Exception as e:
        out['warnings'].append(f'validation failed: {e}')

    lm = Landmarks(case, knee_side, overrides)
    for m in measures:
        for col, _, _ in m.columns:
            out['values'][col] = None
        if m.needs_side and knee_side not in ('L', 'R'):
            out['notes'][m.key] = 'skipped: no laterality'
            continue
        try:
            r = m.module.compute(case, lm)
            out['values'].update(r['values'])
            out['notes'][m.key] = r.get('note')
            if qc:
                out['payloads'][m.key] = m.module.qc_payload(case, lm, r)
            if landmarks:
                for label, idx in m.module.landmark_points(case, lm, r):
                    out['landmarks'].append((m.key, label, tuple(float(v) for v in idx),
                                             index_to_lps(case['header'], *idx)))
        except Exception as e:
            out['notes'][m.key] = f'ERROR {type(e).__name__}: {e}'
            out.setdefault('tracebacks', {})[m.key] = traceback.format_exc()

    _tuberosity_window_check(out, (overrides or {}).get('tuberosity'))
    tpcl_review = (overrides or {}).get('tpcl')
    out['values']['TPCL_review'] = tpcl_review['status'] if tpcl_review else 'auto'
    sulcus_review = (overrides or {}).get('sulcus')
    out['values']['Sulcus_review'] = sulcus_review['status'] if sulcus_review else 'auto'
    out['seconds'] = time.time() - t0
    return out


# Validated on the ORS cohort against the user's own "tuberosity out of window"
# tracking table: every in-window case had >= 38 mm of tibia segmented; a
# tibia shorter than this caught 12/15 out-of-window cases with 0/30 false flags.
# (Testing whether the detected point sits on its search limit did far worse:
# 5/15 caught, 6 false.)
TT_WINDOW_MIN_TIBIA_MM = 37.0


def _tuberosity_window_check(out, review=None):
    """Flag cases where the tibial tuberosity is probably outside the scan's
    field of view: the detector still returns a "most anterior" point, but
    it isn't the tuberosity, so TT-TG / sTT-TG / TTL are unreliable. A human
    review of the tuberosity overrides this guess."""
    out['values']['TT_review'] = 'auto'
    if review:
        out['values']['TT_review'] = review['status'].replace('_', ' ') + (
            ' on axial' if review.get('source') == 'axial' else '')
        out['values']['TT_possibly_out_of_window'] = int(review['status'] in ('not_in_view', 'misaligned'))
        return
    length = out['validation'].get('tibia_length_mm')
    out['values']['TT_possibly_out_of_window'] = None
    if not length:
        return
    out_of_window = length < TT_WINDOW_MIN_TIBIA_MM
    out['values']['TT_possibly_out_of_window'] = int(out_of_window)
    if out_of_window:
        out['cautions'].append(
            f'only {length:.0f} mm of tibia in the scan (< {TT_WINDOW_MIN_TIBIA_MM:.0f} mm) -- the tibial '
            f'tuberosity is probably outside the field of view, so TT-TG / sTT-TG / TTL are unreliable')
