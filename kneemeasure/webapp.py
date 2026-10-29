"""
The KAMINO app:  python -m kneemeasure app

A small web server on this computer (standard library only -- nothing extra
to install) that opens in the browser. Nothing leaves the machine.

1. Setup: pick a segmentations folder and (optionally) a scans folder with
   the built-in folder browser. Each scan's knee side is found automatically
   (see project.py) or entered with L/R buttons.
2. Measure: runs every measurement, with progress. Results go to
   <root>/kneemeasure/ (see project.py for <root>).
3. QC, one step per landmark -- the tibial tuberosity, the tibial posterior
   condylar line and the trochlear sulcus (reference slice). Each starts with
   a prescreen grid where you pick the scans to look at closely (the rest get
   the automatic landmark accepted), then a review of those over the MRI when
   a scan is available. For the tuberosity, the axial series next to the
   scan (sorted/<code>/axial.nii.gz, usually the PD axial) can be shown
   instead of the Dixon: it usually reaches further down the tibia (see
   scans.AxialDisplay). Decisions are saved immediately to
   <root>/kneemeasure_review.json and applied by every later run.

`python -m kneemeasure review <pull>` opens a scan pull straight into review.
"""
import base64
import csv
import functools
import getpass
import json
import os
import queue
import threading
import time
import traceback
import webbrowser
from collections import OrderedDict
from http.server import ThreadingHTTPServer, BaseHTTPRequestHandler
from urllib.parse import urlparse, unquote, parse_qs

import numpy as np

from . import __version__, GITHUB_REPO
from . import review as review_store
from . import tt_tg, ttl, kra, trochlear_morphology
from .tibial_pcl import tpcl_at_slice, tpcl_from_points
from .trochlea import sulcus_on_slice
from .cache import file_sha1
from .landmarks import Landmarks
from .pipeline import TT_WINDOW_MIN_TIBIA_MM
from .project import Project, list_dir, shortcuts, load_recent, remember, suggest_scans_dir, find_sheet
from .runner import run_cohort
from .scans import DisplayScan, AxialDisplay, find_axial
from .volume import load_case, index_to_lps

STATIC = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'webapp')


class App:
    def __init__(self, reviewer=''):
        self.reviewer = reviewer or _login_name()
        self.last_request = time.time()      # for quitting once no browser tab is open (serve's idle_exit_minutes)
        self.project = None
        self.lock = threading.Lock()
        self.cases = OrderedDict()          # code -> loaded case (small LRU)
        self.rows, self.tibia_mm = {}, {}
        self.job = {'running': False, 'done': 0, 'total': 0, 'current': '', 'message': '', 'error': None}
        self.case_info_extra = None          # extra per-case result columns (e.g. a pull's segmentation_source)
        self.loading = {}                    # code -> Event while a case is being loaded
        self.thumbs = {s: {} for s in STEP_KEYS}  # step -> code -> prescreen crops (kept even if the case is evicted)
        self.max_cases = _cache_size()
        self.prefetch_q = queue.Queue()
        for _ in range(LOADERS):             # loading is mostly numpy/zlib work, which runs in parallel
            threading.Thread(target=self._prefetch_worker, daemon=True).start()

    # ---- project
    def open_project(self, segs_dir=None, scans_dir=None, project=None):
        p = project or Project(segs_dir=segs_dir, scans_dir=scans_dir or None)
        with self.lock:
            self.project = p
            self.cases.clear()
            for t in self.thumbs.values():
                t.clear()
        if segs_dir:
            remember(segs_dir, scans_dir)
        self._load_results()
        return self.project_state()

    def project_state(self):
        p = self.project
        if p is None:
            return {'project': None}
        reviews = review_store.load(p.review_file)['cases']
        return {'project': p.summary(), 'cases': p.table(reviews), 'job': self.job}

    def _require(self):
        if self.project is None:
            raise ValueError('open a segmentations folder first')
        return self.project

    def _load_results(self):
        p = self.project
        self.rows, self.tibia_mm = {}, {}
        m = os.path.join(p.out, 'measurements.csv')
        if os.path.exists(m):
            with open(m, newline='') as f:
                self.rows = {r['case']: r for r in csv.DictReader(f)}
        self.tibia_mm = {c: _num(r.get('tibia_length_mm')) for c, r in self.rows.items()}

    def set_side(self, code, side):
        p = self._require()
        p.set_side(code, side)
        with self.lock:
            self.cases.pop(code, None)
        return self.project_state()

    # ---- measuring
    def start_measure(self):
        p = self._require()
        missing = [c for c, (s, _) in p.sides.items() if not s]
        if missing:
            raise ValueError(f"set the knee side for {len(missing)} scan(s) first: {', '.join(missing[:8])}")
        if self.job['running']:
            return self.job

        def progress(done, total, name):
            self.job.update(done=done, total=total, current=name)

        def work():
            try:
                run_cohort(list(p.files.values()), p.out, sides=p.sides, cohort_name=p.name, jobs=MEASURE_JOBS,
                           quiet=True, reviews=review_store.load(p.review_file), progress=progress,
                           case_info=self.case_info_extra)
                self._load_results()
                self.job.update(running=False, message=f'Measured {len(p.files)} scans', error=None)
            except Exception as ex:
                traceback.print_exc()
                self.job.update(running=False, message='', error=str(ex))

        self.job = {'running': True, 'done': 0, 'total': len(p.files), 'current': '',
                    'message': 'Measuring...', 'error': None}
        threading.Thread(target=work, daemon=True).start()
        return self.job

    # ---- review
    def review_state(self, step='tt'):
        p = self._require()
        key = STEP_KEYS[step]
        rv = review_store.load(p.review_file)['cases']
        cases = []
        for code in p.files:
            t = (rv.get(code) or {}).get(key)
            row = self.rows.get(code, {})
            if step == 'tt':
                flag = self._flag(code, t)
                reason = 'maybe out of view' if flag else None
            elif step == 'sulcus':
                a = _num(row.get('sulcus_angle_deg'))
                lo, hi = SULCUS_USUAL_DEG
                reason = None if t or not self.rows else (
                    'no automatic sulcus' if a is None else (f'unusual angle ({a:.0f}°)' if not lo <= a <= hi else None))
                flag = bool(reason)
            else:
                k = _num(row.get('KRA_deg'))
                no_line = bool(self.rows) and k is None and 'tibial PCL not found' in (row.get('notes') or '')
                reason = None if t else ('no automatic line' if no_line else
                                         (f'unusual KRA ({k:+.0f}°)' if k is not None and abs(k) > KRA_UNUSUAL_DEG else None))
                flag = bool(reason)
            cases.append({'case': code, 'side': p.sides[code][0],
                          'status': t['status'] if t else 'unreviewed',
                          'flag_out_of_view': flag, 'flag': flag, 'flag_reason': reason,
                          'tt_tg': _num(row.get('TT_TG_mm')), 'kra': _num(row.get('KRA_deg')),
                          'ttl': _num(row.get('TTL')), 'sulcus_angle': _num(row.get('sulcus_angle_deg')),
                          'tibia_mm': self.tibia_mm.get(code)})
        return {'pull': p.name, 'reviewer': self.reviewer, 'cases': cases, 'job': self.job, 'step': step,
                'measured': bool(self.rows), 'loaded': len(self.cases)}

    def _flag(self, code, reviewed):
        t = self.tibia_mm.get(code)
        return bool(t is not None and t < TT_WINDOW_MIN_TIBIA_MM and not reviewed)

    def _case(self, code):
        """The loaded case (cached). If the prefetcher is already loading it,
        wait for that instead of loading it twice."""
        p = self._require()
        if code not in p.files:
            raise KeyError(code)
        while True:
            with self.lock:
                if code in self.cases:
                    self.cases.move_to_end(code)
                    return self.cases[code]
                ev = self.loading.get(code)
                if ev is None:
                    ev = self.loading[code] = threading.Event()
                    break
            ev.wait()
        try:
            return self._load(code)
        finally:
            with self.lock:
                self.loading.pop(code, None)
            ev.set()

    def _load(self, code):
        p = self.project
        side = p.sides[code][0]
        if not side:
            raise ValueError(f'{code}: set the knee side on the setup page first')
        case = load_case(p.files[code], patella=False)
        lm = Landmarks(case, side)
        auto = lm.tuberosity_auto
        lm.pbcl, lm.tpcl_auto, lm.reference_slice_auto  # the landmarks the live values need, computed once
        troch = _trochlea_store(case)
        for k in ('bone', 'cartilage', 'patella'):  # not needed any more -- keep memory per scan small
            case[k] = None
        zs = np.where(case['tibia'].any(axis=(0, 1)))[0]
        ys = np.where(case['labels'].any(axis=(0, 2)))[0]
        xs = np.where(case['labels'].any(axis=(1, 2)))[0]
        scan, scan_note = None, None
        if p.scans.get(code):
            try:
                scan = DisplayScan(p.scans[code], case['header'], case['labels'].shape)
                scan_note = scan.note
            except Exception as ex:
                scan_note = f'scan not shown: {ex}'
        axial, axial_note = None, None
        if find_axial(p.scans.get(code)):
            try:
                axial = AxialDisplay(find_axial(p.scans[code]), case['header'], case['labels'].shape)
            except Exception as ex:
                axial_note = f'axial scan not shown: {ex}'
        entry = {'case': case, 'side': side, 'auto': auto, 'lm': lm, 'sha1': file_sha1(p.files[code]),
                 'tibia_z': (int(zs.min()), int(zs.max())) if len(zs) else None,
                 'y_extent': (int(ys.min()), int(ys.max())) if len(ys) else None,
                 'x_extent': (int(xs.min()), int(xs.max())) if len(xs) else None,
                 'troch': troch, 'scan': scan, 'scan_note': scan_note, 'axial': axial, 'axial_note': axial_note}
        thumbs = {'tt': _tuberosity_thumb(entry), 'tpcl': _tpcl_thumb(entry), 'sulcus': _sulcus_thumb(entry)}
        with self.lock:
            self.cases[code] = entry
            for step, t in thumbs.items():
                self.thumbs[step][code] = t
            while len(self.cases) > self.max_cases:
                self.cases.popitem(last=False)
        return entry

    def prefetch(self, codes, everything=False):
        """Load these cases in the background, in this order, so opening
        them is instant. A new request replaces the previous queue (the
        reviewer moved on). Normally no more cases than fit in memory are
        queued; everything=True (prescreen) loads all of them once so every
        thumbnail exists, letting the oldest loaded cases be dropped again."""
        codes = [c for c in (codes or []) if self.project and c in self.project.files
                 and self.project.sides[c][0]]
        if everything:
            codes = [c for c in codes if c not in self.thumbs['tt']]
        else:
            codes = codes[:self.max_cases - 1]
        while True:
            try:
                self.prefetch_q.get_nowait()
            except queue.Empty:
                break
        for c in codes:
            self.prefetch_q.put((self.project, c))
        return {'queued': len(codes)}

    def _prefetch_worker(self):
        while True:
            project, code = self.prefetch_q.get()
            if project is not self.project:
                continue
            try:
                self._case(code)
            except Exception:
                pass

    def case_info(self, code, step='tt'):
        e = self._case(code)
        c, auto, sp = e['case'], e['auto'], e['case']['spacing']
        reviews = review_store.load(self.project.review_file)['cases'].get(code) or {}
        t = reviews.get('tuberosity')
        tp = e['lm'].tpcl_auto
        ref = e['lm'].reference_slice_auto
        tibia_mm = (e['tibia_z'][1] - e['tibia_z'][0] + 1) * sp[2] if e['tibia_z'] else None
        now = self.recompute(code)
        return {
            'case': code, 'side': e['side'], 'shape': list(c['labels'].shape), 'spacing': list(sp),
            'tibia_z': e['tibia_z'], 'y_extent': e['y_extent'], 'x_extent': e['x_extent'],
            'auto_point': list(auto['point']) if auto.get('point') else None,
            'auto_note': auto.get('note'),
            'search_window_z': [auto.get('distal_boundary_z'), auto.get('plateau_boundary_z')],
            'review': reviews.get(STEP_KEYS[step]), 'tibia_length_mm': tibia_mm,
            'tpcl_auto': _line(tp) if tp.get('theta') is not None else None,
            'tpcl_auto_note': tp.get('note'),
            'sulcus_auto': _sulcus(ref, e) if ref.get('z') is not None else None,
            'sulcus_auto_note': ref.get('note'),
            'trochlea_z': [e['troch'][0], e['troch'][0] + e['troch'][1].shape[2] - 1] if e['troch'] else None,
            'pbcl': _line(e['lm'].pbcl) if e['lm'].pbcl.get('theta') is not None else None,
            'flag_out_of_view': bool(tibia_mm is not None and tibia_mm < TT_WINDOW_MIN_TIBIA_MM and not t),
            'sulcus_review': now['sulcus'] if reviews.get('sulcus') else None,
            'values': now['values'],
            'has_scan': e['scan'] is not None, 'scan_note': e['scan_note'],
            'axial_z_below': e['axial'].z_below if e['axial'] is not None else None, 'axial_note': e['axial_note'],
        }

    def slice_data(self, code, plane, index, src='dixon'):
        """One slice of the segmentation with the scan under it. src='axial':
        the axial series instead of the Dixon, with the axial view going
        down to z = -z_below and the sagittal one padded below the
        segmentation (see AxialDisplay)."""
        if plane not in ('axial', 'sagittal'):
            raise ValueError(plane)
        e = self._case(code)
        labels = e['case']['labels']
        scan = e['axial'] if src == 'axial' else e['scan']
        below = e['axial'].z_below if src == 'axial' and scan is not None else 0
        axis = 2 if plane == 'axial' else 0
        index = int(np.clip(index, -below if plane == 'axial' else 0, labels.shape[axis] - 1))
        if plane == 'sagittal':
            sl = np.pad(labels[index, :, :], ((0, 0), (below, 0)))
        else:
            sl = labels[:, :, index] if index >= 0 else np.zeros(labels.shape[:2], dtype=labels.dtype)
        out = {'shape': list(sl.shape), 'data': _bytes(sl), 'scan': None, 'scan_shape': None}
        if scan is not None:
            img = scan.slice(plane, index)       # the axial: finer than the segmentation (AxialDisplay.factor)
            if src != 'axial':
                img = _upsample(img)
            out['scan'], out['scan_shape'] = _bytes(img), list(img.shape)
        return out

    # ---- prescreen
    def thumbs_since(self, have, step='tt'):
        p = self._require()
        have = set(have or [])
        with self.lock:
            thumbs = self.thumbs[step]
            new = {c: t for c, t in thumbs.items() if c not in have and c in p.files}
            ready = set(thumbs)
        pending = [c for c in p.files if c not in ready and p.sides[c][0]]
        return {'thumbs': new, 'pending': pending, 'state': self.review_state(step)}

    def accept_auto(self, codes, step='tt', note='accepted at prescreen'):
        """Accept the automatic landmark for many scans at once (scans
        already decided, or without an automatic landmark, are left)."""
        p = self._require()
        key = STEP_KEYS[step]
        done, skipped = [], []
        with self.lock:
            data = review_store.load(p.review_file)
            for code in codes or []:
                t = self.thumbs[step].get(code)
                if (code not in p.files or not t or not t['auto']
                        or (data['cases'].get(code) or {}).get(key)):
                    skipped.append(code)
                    continue
                if step == 'tt':
                    review_store.set_tuberosity(data, code, 'accepted', point=t['auto'], reviewer=self.reviewer,
                                                note=note, auto_point=t['auto'], file_sha1=t['sha1'],
                                                point_lps_mm=t['auto_lps'])
                elif step == 'sulcus':
                    review_store.set_sulcus(data, code, 'accepted', t['auto'], auto=t['auto'], method='auto',
                                            reviewer=self.reviewer, note=note, file_sha1=t['sha1'])
                else:
                    review_store.set_tpcl(data, code, 'accepted', t['auto'], auto=t['auto'], method='auto',
                                          reviewer=self.reviewer, note=note, file_sha1=t['sha1'])
                done.append(code)
            review_store.save(p.review_file, data)
        return {'accepted': done, 'skipped': skipped}

    def tpcl_preview(self, code, z):
        """The line the 'use this slice' button would give, for drawing."""
        e = self._case(code)
        line = tpcl_at_slice(e['case'], e['side'], z)
        return {'line': _line(line) if line else None}

    def decide_tpcl(self, code, status, z=None, points=None, note=''):
        e = self._case(code)
        auto = e['lm'].tpcl_auto
        if status == 'accepted':
            if auto.get('theta') is None:
                raise ValueError('there is no automatic line to accept')
            line, method = auto, 'auto'
        elif points:
            if len(points) != 2:
                raise ValueError('click exactly two points')
            line, method = tpcl_from_points(e['case'], e['side'], points[0], points[1], z), 'points'
        else:
            line, method = tpcl_at_slice(e['case'], e['side'], z), 'slice'
            if line is None:
                raise ValueError("this slice doesn't show two posterior condyles -- click two points instead")
        with self.lock:
            data = review_store.load(self.project.review_file)
            review_store.set_tpcl(data, code, 'accepted' if status == 'accepted' else 'replaced', line,
                                  auto=auto if auto.get('theta') is not None else None, method=method,
                                  reviewer=self.reviewer, note=note, file_sha1=e['sha1'])
            review_store.save(self.project.review_file, data)
        return self.recompute(code)

    def sulcus_preview(self, code, z):
        """The sulcus and facet peaks the 'use this slice' button would give."""
        e = self._case(code)
        r = _sulcus_at(e, z)
        return {'sulcus': _sulcus(r, e) if r else None}

    def decide_sulcus(self, code, status, z=None, points=None, note=''):
        e = self._case(code)
        auto = e['lm'].reference_slice_auto
        if status == 'accepted':
            if auto.get('z') is None:
                raise ValueError('there is no automatic sulcus to accept')
            ref, method = auto, 'auto'
        elif points:
            if len(points) != 3:
                raise ValueError('click three points: the sulcus and the two facet peaks')
            # any order: the middle one from side to side is the sulcus, the outer two the facet peaks
            # (which is medial / lateral follows the knee side, as for the automatic peaks)
            a, s, b = sorted(([float(v) for v in pt[:2]] for pt in points), key=lambda p: p[0])
            if not a[0] < s[0] < b[0]:
                raise ValueError('the three points need to be side by side, not one above another')
            ref, method = {'z': int(z), 'sulcus': s, 'peak1': a, 'peak2': b}, 'points'
        else:
            ref, method = _sulcus_at(e, z), 'slice'
            if ref is None:
                raise ValueError("no sulcus found on this slice -- click the three points instead")
        with self.lock:
            data = review_store.load(self.project.review_file)
            review_store.set_sulcus(data, code, 'accepted' if status == 'accepted' else 'replaced', ref,
                                    auto=auto if auto.get('z') is not None else None, method=method,
                                    reviewer=self.reviewer, note=note, file_sha1=e['sha1'])
            review_store.save(self.project.review_file, data)
        return self.recompute(code)

    def decide(self, code, status, point=None, note='', source=None):
        """source='axial': the point was clicked on the axial series, so it
        may lie below the segmentation (negative z)."""
        e = self._case(code)
        if status == 'accepted':
            if not e['auto'].get('point'):
                raise ValueError('there is no automatic point to accept')
            point, source = e['auto']['point'], None
        has_point = point is not None and status in ('accepted', 'replaced')
        if has_point and not (-(e['axial'].z_below if source == 'axial' and e['axial'] else 0)
                              <= point[2] < e['case']['labels'].shape[2]):
            raise ValueError('the point is outside the scan')
        lps = index_to_lps(e['case']['header'], *point) if has_point else None
        with self.lock:
            data = review_store.load(self.project.review_file)
            review_store.set_tuberosity(
                data, code, status, point=point if has_point else None, reviewer=self.reviewer, note=note,
                auto_point=e['auto'].get('point'), file_sha1=e['sha1'], point_lps_mm=lps,
                source=source if has_point else None)
            review_store.save(self.project.review_file, data)
        return self.recompute(code)

    def clear(self, code, step='tt'):
        with self.lock:
            data = review_store.load(self.project.review_file)
            review_store.clear_review(data, code, STEP_KEYS[step])
            review_store.save(self.project.review_file, data)
        return self.recompute(code)

    def recompute(self, code):
        """TT-TG / sTT-TG / TTL / KRA / trochlear values with the case's current
        reviews, from the already-loaded case and landmarks (same functions a
        full run uses)."""
        e = self._case(code)
        overrides = review_store.overrides_for(review_store.load(self.project.review_file), code)
        lm = Landmarks(e['case'], e['side'], overrides)
        for k in ('pbcl', 'tpcl_auto', 'reference_slice_auto', 'tuberosity_auto'):
            lm.__dict__[k] = getattr(e['lm'], k)
        values, notes = {}, {}
        for key, module in (('tt_tg', tt_tg), ('ttl', ttl), ('kra', kra), ('trochlea', trochlear_morphology)):
            r = module.compute(e['case'], lm)
            values.update(r['values'])
            notes[key] = r.get('note')
        tp, ref = lm.tpcl, lm.reference_slice
        return {'values': {k: values.get(k) for k in ('TT_TG_mm', 'sTT_TG_mm', 'TTL', 'KRA_deg', 'sulcus_angle_deg',
                                                      'LTI_deg', 'trochlear_depth_mm')},
                'notes': notes, 'review': (overrides or {}).get('tuberosity'),
                'tpcl_review': (overrides or {}).get('tpcl'), 'sulcus_review': (overrides or {}).get('sulcus'),
                'tpcl': _line(tp) if tp.get('theta') is not None else None,
                'sulcus': _sulcus(ref, e) if ref.get('z') is not None else None}


STEP_KEYS = {'tt': 'tuberosity', 'tpcl': 'tpcl', 'sulcus': 'sulcus'}
LOADERS = max(2, min(6, (os.cpu_count() or 4) // 2))   # background threads loading scans for review
MEASURE_JOBS = max(2, min(8, (os.cpu_count() or 4) // 2))   # processes measuring in parallel
KRA_UNUSUAL_DEG = 15.0      # prescreen pre-selects scans whose KRA is beyond +-this
SULCUS_USUAL_DEG = (140.0, 170.0)   # ... and whose sulcus angle is outside this (all pulls so far: 144-164)


def _line(r):
    return {'z': int(r['z']), 'theta': float(r['theta']),
            'medial_pt': [float(v) for v in r['medial_pt']], 'lateral_pt': [float(v) for v in r['lateral_pt']]}


def _sulcus(r, entry):
    """A reference-slice result as drawn by the app: sulcus, facet peaks
    labelled medial / lateral (LEFT knee: lateral = larger x) and the angle."""
    sp = entry['case']['spacing']
    a, b = sorted(([float(v) for v in r['peak1'][:2]], [float(v) for v in r['peak2'][:2]]), key=lambda p: p[0])
    lateral, medial = (b, a) if entry['side'] == 'L' else (a, b)
    s = [float(v) for v in r['sulcus'][:2]]
    return {'z': int(r['z']), 'sulcus': s, 'peak1': [float(v) for v in r['peak1'][:2]],
            'peak2': [float(v) for v in r['peak2'][:2]], 'medial_pt': medial, 'lateral_pt': lateral,
            'angle': trochlear_morphology.compute_sulcus_angle_deg(s, medial, lateral, sp[:2])}


def _trochlea_store(case):
    """The femur bone (bit 1) and cartilage (bit 2) over the cartilage's
    slices, so a slice a reviewer picks can be re-detected after the full
    masks are dropped. (z0, array) or None."""
    zs = np.where(case['cartilage'].any(axis=(0, 1)))[0]
    if not len(zs):
        return None
    z0, z1 = int(zs.min()), int(zs.max()) + 1
    return z0, (case['bone'][:, :, z0:z1].astype(np.uint8) | (case['cartilage'][:, :, z0:z1].astype(np.uint8) << 1))


def _sulcus_at(entry, z):
    """Facet peaks + sulcus re-detected on slice z, or None."""
    if entry['troch'] is None or z is None:
        return None
    z0, arr = entry['troch']
    z = int(z)
    if not z0 <= z < z0 + arr.shape[2]:
        return None
    sl = arr[:, :, z - z0]
    return sulcus_on_slice((sl & 1).astype(bool), (sl & 2).astype(bool), entry['case']['spacing'][0], z)


def _sulcus_thumb(entry, margin_mm=8.0):
    """The axial slice the automatic sulcus was found on, cropped to the
    front half of the femur (the trochlea), with the sulcus and the two
    facet peaks."""
    case, scan, lm = entry['case'], entry['scan'], entry['lm']
    labels, (sx, sy, sz) = case['labels'], case['spacing']
    ref = lm.reference_slice_auto
    auto = _sulcus(ref, entry) if ref.get('z') is not None else None
    troch = entry['troch']
    z = auto['z'] if auto else (troch[0] + troch[1].shape[2] * 2 // 3 if troch else labels.shape[2] // 2)
    fem = labels[:, :, z] == 1
    xs, ys = np.where(fem.any(axis=1))[0], np.where(fem.any(axis=0))[0]
    if len(xs) == 0:
        xs, ys = np.array([0, labels.shape[0] - 1]), np.array([0, labels.shape[1] - 1])
    mx, my = int(margin_mm / sx), int(margin_mm / sy)
    x0, x1 = max(0, xs.min() - mx), min(labels.shape[0], xs.max() + mx + 1)
    y0, y1 = max(0, ys.min() - my), min(labels.shape[1], (ys.min() + ys.max()) // 2 + 1)
    lab = labels[x0:x1, y0:y1, z]
    img = scan.slice('axial', z)[x0:x1, y0:y1] if scan is not None else None
    crop = {'shape': list(lab.shape), 'data': _b64(lab), 'scan': _b64(img) if img is not None else None,
            'sulcus': None}
    if auto:
        crop['sulcus'] = [[p[0] - x0, p[1] - y0] for p in (auto['medial_pt'], auto['sulcus'], auto['lateral_pt'])]
    return {'side': entry['side'], 'spacing': [sx, sy, sz], 'auto': auto, 'sha1': entry['sha1'],
            'axial': crop, 'z': z}


def _tpcl_thumb(entry, margin_mm=10.0):
    """The axial tibia slice the automatic TPCL was taken from, cropped to
    the tibia, with the line's two condyle points (for drawing the line)."""
    case, scan, lm = entry['case'], entry['scan'], entry['lm']
    labels, (sx, sy, sz) = case['labels'], case['spacing']
    tp = lm.tpcl_auto
    line = _line(tp) if tp.get('theta') is not None else None
    z = line['z'] if line else (int(np.mean(entry['tibia_z'])) if entry['tibia_z'] else labels.shape[2] // 2)
    tib = case['tibia'][:, :, z]
    xs, ys = np.where(tib.any(axis=1))[0], np.where(tib.any(axis=0))[0]
    if len(xs) == 0:
        xs, ys = np.array([0, labels.shape[0] - 1]), np.array([0, labels.shape[1] - 1])
    mx, my = int(margin_mm / sx), int(margin_mm / sy)
    x0, x1 = max(0, xs.min() - mx), min(labels.shape[0], xs.max() + mx + 1)
    y0, y1 = max(0, ys.min() - my), min(labels.shape[1], ys.max() + my + 1)
    lab = labels[x0:x1, y0:y1, z]
    img = scan.slice('axial', z)[x0:x1, y0:y1] if scan is not None else None
    crop = {'shape': list(lab.shape), 'data': _b64(lab), 'scan': _b64(img) if img is not None else None,
            'line': None}
    if line:
        crop['line'] = [[line['medial_pt'][0] - x0, line['medial_pt'][1] - y0],
                        [line['lateral_pt'][0] - x0, line['lateral_pt'][1] - y0]]
    return {'side': entry['side'], 'spacing': [sx, sy, sz], 'auto': line, 'sha1': entry['sha1'],
            'axial': crop, 'z': z}


def _tuberosity_thumb(entry, half_mm=35.0, up_mm=55.0):
    """Small crops around the automatic tuberosity point for the prescreen
    grid: the axial slice through it and the sagittal slice through it (that
    one reaching down to the bottom of the scan when it's close, so a tibia
    cut off by the field of view is visible)."""
    case, auto, scan = entry['case'], entry['auto'], entry['scan']
    labels, (sx, sy, sz) = case['labels'], case['spacing']
    nx, ny, nz = labels.shape
    pt = auto.get('point')
    if pt is None:                                   # no automatic point: centre on the tibia
        tib = np.argwhere(case['tibia'])
        pt = tuple(int(v) for v in (tib.mean(axis=0) if len(tib) else (nx // 2, ny // 2, nz // 2)))
    x, y, z = (int(v) for v in pt)
    hx, hy, hz = int(half_mm / sx), int(half_mm / sy), int(half_mm / sz)
    x0, x1 = max(0, x - hx), min(nx, x + hx + 1)
    y0, y1 = max(0, y - hy), min(ny, y + hy + 1)
    z0, z1 = max(0, z - hz), min(nz, z + int(up_mm / sz) + 1)      # more room above: the tubercle's profile

    def crop(plane):
        if plane == 'axial':
            lab = labels[x0:x1, y0:y1, z]
            img = scan.slice('axial', z)[x0:x1, y0:y1] if scan is not None else None
            mark = (x - x0, y - y0)
        else:
            lab = labels[x, y0:y1, z0:z1]
            img = scan.slice('sagittal', x)[y0:y1, z0:z1] if scan is not None else None
            mark = (y - y0, z - z0)
        return {'shape': list(lab.shape), 'data': _b64(lab), 'scan': _b64(img) if img is not None else None,
                'marker': list(mark) if auto.get('point') else None}

    tibia_mm = (entry['tibia_z'][1] - entry['tibia_z'][0] + 1) * sz if entry['tibia_z'] else None
    return {'side': entry['side'], 'spacing': [sx, sy, sz], 'auto_point': list(auto['point']) if auto.get('point') else None,
            'auto': list(auto['point']) if auto.get('point') else None,
            'auto_lps': index_to_lps(case['header'], *auto['point']) if auto.get('point') else None,
            'sha1': entry['sha1'], 'tibia_length_mm': tibia_mm,
            'scan_bottom_visible': z0 == 0 and entry['tibia_z'] is not None and entry['tibia_z'][0] == 0,
            'axial': crop('axial'), 'sagittal': crop('sagittal'),
            'values': None}


def _cache_size(per_case_mb=60, share=0.25, floor=6, ceiling=300):
    """How many loaded scans to keep: a quarter of the free memory, at about
    60 MB per scan (measured: segmentation, tibia mask, display-ready MRI)."""
    try:
        with open('/proc/meminfo') as f:
            avail_kb = next(int(l.split()[1]) for l in f if l.startswith('MemAvailable'))
        return max(floor, min(ceiling, int(avail_kb / 1024 * share / per_case_mb)))
    except Exception:
        return 40


def _upsample(img, factor=3):
    """A Dixon slice made `factor` times finer by cubic (Catmull-Rom)
    interpolation, for display: the browser then only stretches it a little,
    which looks sharper than its own smoothing of the 0.9 mm voxels (blurry)
    or none (blocky). Four taps along each axis, a few ms per slice; fine
    pixel s is centred at (s + 0.5) / factor - 0.5, where the page expects."""
    out = img.astype(np.float32)
    for axis in (0, 1):
        idx, w = _cubic_taps(out.shape[axis], factor)
        out = sum(np.take(out, idx[k], axis=axis) * (w[k][:, None] if axis == 0 else w[k][None, :]) for k in range(4))
    return np.clip(np.rint(out), 0, 255).astype(np.uint8)


@functools.lru_cache(maxsize=16)
def _cubic_taps(n, factor):
    """For n samples made factor times finer: the 4 source indices (edges
    clamped) and Catmull-Rom weights of every output sample, each (4, n * factor)."""
    t = (np.arange(n * factor) + 0.5) / factor - 0.5
    i0 = np.floor(t).astype(int)
    idx, w = [], []
    for k in (-1, 0, 1, 2):
        d = np.abs(t - (i0 + k))
        w.append(np.where(d < 1, 1.5 * d ** 3 - 2.5 * d ** 2 + 1,
                          np.where(d < 2, -0.5 * d ** 3 + 2.5 * d ** 2 - 4 * d + 2, 0)).astype(np.float32))
        idx.append(np.clip(i0 + k, 0, n - 1))
    return np.array(idx), np.array(w)


def _bytes(a):
    return np.ascontiguousarray(a, dtype=np.uint8).tobytes()


def _b64(a):
    return base64.b64encode(_bytes(a)).decode()


def _login_name():
    try:
        return getpass.getuser()
    except Exception:
        return ''


def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def make_handler(app):
    class Handler(BaseHTTPRequestHandler):
        def log_message(self, *args):
            pass

        def _send(self, code, body, ctype='application/json', headers=None):
            data = body if isinstance(body, bytes) else json.dumps(body).encode()
            self.send_response(code)
            self.send_header('Content-Type', ctype)
            for k, v in (headers or {}).items():
                self.send_header(k, v)
            self.send_header('Content-Length', str(len(data)))
            self.send_header('Cache-Control', 'no-store')
            self.end_headers()
            self.wfile.write(data)

        def _route(self, method):
            app.last_request = time.time()
            url = urlparse(self.path)
            q = {k: v[0] for k, v in parse_qs(url.query).items()}
            parts = [unquote(s) for s in url.path.strip('/').split('/') if s]
            try:
                if method == 'GET' and not parts:
                    with open(os.path.join(STATIC, 'index.html'), 'rb') as f:
                        return self._send(200, f.read(), 'text/html; charset=utf-8')
                if parts[:1] != ['api']:
                    return self._send(404, {'error': 'not found'})
                p = parts[1:]
                body = {}
                if method == 'POST':
                    n = int(self.headers.get('Content-Length') or 0)
                    body = json.loads(self.rfile.read(n) or b'{}')
                if method == 'GET':
                    if p == ['fs']:
                        return self._send(200, list_dir(q.get('path')))
                    if p == ['version']:          # also the open page's once-a-minute "still here"
                        return self._send(200, {'app': 'kneemeasure', 'version': __version__, 'repo': GITHUB_REPO})
                    if p == ['places']:
                        return self._send(200, {'shortcuts': shortcuts(), 'recent': load_recent(limit=5)})
                    if p == ['suggest']:
                        segs = q.get('segs_dir', '')
                        scans = suggest_scans_dir(segs)
                        return self._send(200, {'scans_dir': scans, 'sheet': find_sheet(os.path.dirname(segs), scans)})
                    if p == ['project']:
                        return self._send(200, app.project_state())
                    if p == ['job']:
                        return self._send(200, app.job)
                    if p == ['loaded']:
                        return self._send(200, {'loaded': len(app.cases),
                                                'total': len(app.project.files) if app.project else 0})
                    if p == ['state']:
                        return self._send(200, app.review_state(q.get('step', 'tt')))
                    if len(p) == 2 and p[0] == 'case':
                        return self._send(200, app.case_info(p[1], q.get('step', 'tt')))
                    if len(p) == 4 and p[0] == 'case' and p[2] in ('axial', 'sagittal'):
                        # raw bytes (segmentation, then the scan), shapes in the headers: no base64 to decode
                        d = app.slice_data(p[1], p[2], int(p[3]), q.get('src', 'dixon'))
                        shape = lambda s: ','.join(str(v) for v in s) if s else ''
                        return self._send(200, d['data'] + (d['scan'] or b''), 'application/octet-stream',
                                          {'X-Shape': shape(d['shape']), 'X-Scan-Shape': shape(d['scan_shape'])})
                if method == 'POST':
                    if p == ['project']:
                        return self._send(200, app.open_project(body.get('segs_dir'), body.get('scans_dir')))
                    if p == ['prefetch']:
                        return self._send(200, app.prefetch(body.get('cases'), bool(body.get('everything'))))
                    if p == ['thumbs']:
                        return self._send(200, app.thumbs_since(body.get('have'), body.get('step', 'tt')))
                    if p == ['prescreen', 'accept']:
                        return self._send(200, app.accept_auto(body.get('cases'), body.get('step', 'tt')))
                    if len(p) == 3 and p[0] == 'case' and p[2] == 'tpcl_preview':
                        return self._send(200, app.tpcl_preview(p[1], body.get('z')))
                    if len(p) == 3 and p[0] == 'case' and p[2] == 'sulcus_preview':
                        return self._send(200, app.sulcus_preview(p[1], body.get('z')))
                    if p == ['side']:
                        return self._send(200, app.set_side(body.get('case'), body.get('side')))
                    if p in (['measure'], ['rebuild']):
                        return self._send(200, app.start_measure())
                    if len(p) == 3 and p[0] == 'case' and p[2] == 'decide':
                        if body.get('step') == 'sulcus':
                            return self._send(200, app.decide_sulcus(p[1], body.get('status'), body.get('z'),
                                                                     body.get('points'), body.get('note', '')))
                        if body.get('step') == 'tpcl':
                            return self._send(200, app.decide_tpcl(p[1], body.get('status'), body.get('z'),
                                                                   body.get('points'), body.get('note', '')))
                        return self._send(200, app.decide(p[1], body.get('status'), body.get('point'),
                                                          body.get('note', ''), body.get('source')))
                    if len(p) == 3 and p[0] == 'case' and p[2] == 'clear':
                        return self._send(200, app.clear(p[1], body.get('step', 'tt')))
                return self._send(404, {'error': 'not found'})
            except KeyError as e:
                return self._send(404, {'error': f'unknown scan {e}'})
            except Exception as e:
                traceback.print_exc()
                return self._send(400, {'error': str(e)})

        def do_GET(self):
            self._route('GET')

        def do_POST(self):
            self._route('POST')

    return Handler


def serve(project=None, segs_dir=None, scans_dir=None, port=8765, reviewer='', open_browser=True,
          case_info=None, any_port_if_busy=False, idle_exit_minutes=None):
    """Start the app. project: an already-built Project (opens in review);
    segs_dir/scans_dir: open that folder on the setup page; neither: empty
    setup page. any_port_if_busy: use a free port if `port` is taken.
    idle_exit_minutes: stop once no page has been open for that long (the
    page checks in every minute) and nothing is being measured -- for the
    packaged app, which has no window to close."""
    app = App(reviewer)
    app.case_info_extra = case_info
    start = '#setup'
    if project is not None:
        app.open_project(project=project)
        start = '#review'
    elif segs_dir:
        app.open_project(segs_dir, scans_dir)
    try:
        server = ThreadingHTTPServer(('127.0.0.1', port), make_handler(app))
    except OSError:
        if not any_port_if_busy:
            raise
        server = ThreadingHTTPServer(('127.0.0.1', 0), make_handler(app))
    url = f'http://127.0.0.1:{server.server_address[1]}/{start}'
    print(f'KAMINO {__version__} is running at {url}')
    print('Everything stays on this computer. Press Ctrl+C here to stop it.')
    if open_browser:
        threading.Timer(0.5, lambda: webbrowser.open(url)).start()
    if idle_exit_minutes:
        def watch():
            while True:
                time.sleep(20)
                if not app.job['running'] and time.time() - app.last_request > idle_exit_minutes * 60:
                    print(f'No KAMINO page open for {idle_exit_minutes} minutes -- stopping.')
                    server.shutdown()
                    return
        threading.Thread(target=watch, daemon=True).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
