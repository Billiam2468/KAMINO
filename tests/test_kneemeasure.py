"""
Run from the repository folder:

    python -m unittest discover -s tests -v

Unit tests need no data. The regression tests re-measure a few real cases
(not included: patient data) and compare against values that were
cross-checked against the original research scripts (510/510 Controls and
272/272 cartilage-damage values identical); they run only when
KNEEMEASURE_DATA points at the folder holding the scan pulls:

    KNEEMEASURE_DATA="/path/to/Scan Pulls" python -m unittest discover -s tests
"""
import json
import os
import tempfile
import unittest

import numpy as np

from kneemeasure import agreement, kra, stats
from kneemeasure.laterality import LateralityTable, resolve
from kneemeasure.volume import normalize_slicer_layers, parse_side

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.environ.get('KNEEMEASURE_DATA', '')    # the Scan Pulls folder, for the regression tests (skipped without it)


class TestUnits(unittest.TestCase):
    def test_parse_side(self):
        for v, want in (('L', 'L'), ('Left', 'L'), (' right ', 'R'), ('Rt', 'R'), ('x', None), ('', None)):
            self.assertEqual(parse_side(v), want)

    def test_kra_angle_wrap(self):
        self.assertAlmostEqual(kra._signed_angle_diff_deg(89.0, -89.0), -2.0)
        self.assertAlmostEqual(kra._signed_angle_diff_deg(-89.0, 89.0), 2.0)
        self.assertAlmostEqual(kra._signed_angle_diff_deg(10.0, 4.0), 6.0)

    def test_slicer_layers(self):
        # layer 0: femur bone (1), patella bone (2 -> Segment_3), tibia (3 -> Segment_5)
        # layer 1: femur cartilage (value 1 -> Segment_2), overlapping one bone voxel
        data = np.zeros((2, 4, 4, 4), dtype=np.uint8)
        data[0, 0, :, :] = 1
        data[0, 1, :, :] = 2
        data[0, 3, :, :] = 3
        data[1, 0, 0, 0] = 1
        data[1, 2, :, :] = 1
        header = {'space': 'left-posterior-superior', 'dimension': 4,
                  'space directions': [[np.nan] * 3, [1, 0, 0], [0, 1, 0], [0, 0, 1]],
                  'space origin': np.zeros(3), 'kinds': ['list', 'domain', 'domain', 'domain']}
        for i, (name, layer, value) in enumerate((('Segment_1', 0, 1), ('Segment_2', 1, 1),
                                                   ('Segment_3', 0, 2), ('Segment_5', 0, 3))):
            header.update({f'Segment{i}_ID': name, f'Segment{i}_Name': name,
                           f'Segment{i}_Layer': str(layer), f'Segment{i}_LabelValue': str(value)})
        out, h3, notes = normalize_slicer_layers(data, header)
        self.assertEqual(out.shape, (4, 4, 4))
        self.assertEqual(out[0, 1, 1], 1)      # femur bone
        self.assertEqual(out[0, 0, 0], 2)      # cartilage wins the overlap
        self.assertEqual(out[1, 0, 0], 3)      # patella bone
        self.assertEqual(out[2, 0, 0], 2)      # femur cartilage from layer 1
        self.assertEqual(out[3, 0, 0], 5)      # tibia
        self.assertEqual(h3['dimension'], 3)
        self.assertTrue(notes)

    def test_plain_3d_passthrough(self):
        data = np.ones((3, 3, 3), dtype=np.uint8)
        out, h, notes = normalize_slicer_layers(data, {'space': 'left-posterior-superior'})
        self.assertIs(out, data)
        self.assertEqual(notes, [])

    def test_laterality_table(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, 'map.csv')
            with open(p, 'w') as f:
                f.write('code,mrn,laterality,folder_name\nABC,1,L,ABC_2024\nDEF,2,Right,DEF_2023\n')
            t = LateralityTable(p)
            self.assertEqual(t.lookup('ABC')[0], 'L')
            self.assertEqual(t.lookup('DEF_2023')[0], 'R')       # via folder_name
            self.assertEqual(t.lookup('ABC_1')[0], 'L')          # trailing suffix dropped
            self.assertIsNone(t.lookup('XYZ')[0])
            self.assertEqual(resolve('ABC', 'R', t)[1], 'conflict')
            self.assertEqual(resolve('ABC', 'R', t, priority='header')[0], 'R')
            self.assertEqual(resolve('ABC', None, t)[0], 'L')

    def test_holm(self):
        adj = stats.holm([0.01, 0.04, 0.03, None])
        self.assertAlmostEqual(adj[0], 0.03)
        self.assertAlmostEqual(adj[2], 0.06)
        self.assertAlmostEqual(adj[1], 0.06)
        self.assertIsNone(adj[3])

    def test_icc_perfect_and_offset(self):
        x = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        icc2, icc3 = agreement.icc(x, x)
        self.assertAlmostEqual(icc2, 1.0)
        icc2, icc3 = agreement.icc(x, x + 2.0)
        self.assertAlmostEqual(icc3, 1.0)       # consistency ignores a constant offset
        self.assertLess(icc2, 1.0)              # absolute agreement doesn't

    def test_outlier_flags(self):
        rows = [{'case': str(i), 'v': v} for i, v in enumerate([10, 11, 9, 10, 12, 10, 11, 9, 40])]
        flags = stats.flag_outliers(rows, ['v'])
        self.assertEqual([f['case'] for f in flags], ['8'])


class TestWorkflow(unittest.TestCase):
    def _pull(self, d, cleaned=(), dixon=(), mapping='code,mrn,laterality\nAAAAA,1,L\nBBBBB,2,R\nCCCCC,3,L\n'):
        os.makedirs(os.path.join(d, 'sorted'))
        with open(os.path.join(d, 'sorted', 'mrn_code_mapping.csv'), 'w') as f:
            f.write(mapping)
        for sub, codes in (('cleaned_segs', cleaned), (os.path.join('dixon_segmentations', 'segmentations'), dixon)):
            if codes:
                os.makedirs(os.path.join(d, sub))
            for c in codes:
                open(os.path.join(d, sub, f'{c}.nrrd'), 'w').close()

    def test_pull_best_source(self):
        from kneemeasure.pulls import inspect_pull
        with tempfile.TemporaryDirectory() as d:
            self._pull(d, cleaned=['AAAAA'], dixon=['AAAAA', 'BBBBB'])
            with open(os.path.join(d, 'kneemeasure.toml'), 'w') as f:
                f.write('name = "Test pull"\n[exclude]\nCCCCC = "bad scan"\n')
            spec = inspect_pull(d)
            self.assertTrue(spec['eligible'])
            self.assertEqual(spec['name'], 'Test pull')
            self.assertEqual(spec['case_info']['AAAAA']['segmentation_source'], 'cleaned')
            self.assertEqual(spec['case_info']['BBBBB']['segmentation_source'], 'dixon model')
            self.assertEqual(spec['exclude'], {'CCCCC': 'bad scan'})
            self.assertEqual(len(inspect_pull(d, source='cleaned')['files']), 1)

    def test_pull_without_laterality(self):
        from kneemeasure.pulls import inspect_pull
        with tempfile.TemporaryDirectory() as d:
            self._pull(d, dixon=['AAAAA'], mapping='code,mrn\nAAAAA,1\n')
            spec = inspect_pull(d)
            self.assertFalse(spec['eligible'])
            self.assertIn('laterality', spec['reason'])

    def test_cache_roundtrip_and_key(self):
        from kneemeasure.cache import ResultCache, cache_key
        with tempfile.TemporaryDirectory() as d:
            c = ResultCache(d)
            k = cache_key('abc', 'L', None, True, True)
            c.put('abc', k, {'values': {'KRA_deg': 1.0}})
            self.assertEqual(c.get('abc', k)['values']['KRA_deg'], 1.0)
            self.assertIsNone(c.get('abc', cache_key('abc', 'R', None, True, True)))   # side changed
            c.prune(['other'])
            self.assertIsNone(c.get('abc', k))

    def test_change_report(self):
        from kneemeasure.runner import _archive_and_diff, write_csv
        with tempfile.TemporaryDirectory() as d:
            write_csv(os.path.join(d, 'measurements.csv'),
                      [{'case': 'A', 'knee_side': 'L', 'KRA_deg': 1.0, 'file_sha1': 'x'},
                       {'case': 'B', 'knee_side': 'R', 'KRA_deg': 2.0, 'file_sha1': 'y'}],
                      ['case', 'knee_side', 'KRA_deg', 'file_sha1'])
            changes = _archive_and_diff(d, [{'case': 'A', 'knee_side': 'L', 'KRA_deg': 1.5, 'file_sha1': 'z'},
                                            {'case': 'C', 'knee_side': 'L', 'KRA_deg': 0.0, 'file_sha1': 'w'}],
                                        ['KRA_deg'])
            kinds = sorted((c['case'], c['change']) for c in changes)
            self.assertEqual(kinds, [('A', 'segmentation file changed'), ('A', 'value changed'),
                                     ('B', 'removed'), ('C', 'added')])
            self.assertTrue(os.listdir(os.path.join(d, 'history')))


class TestReview(unittest.TestCase):
    def test_review_store_roundtrip(self):
        from kneemeasure import review
        with tempfile.TemporaryDirectory() as d:
            f = os.path.join(d, review.REVIEW_FILE)
            data = review.load(f)
            review.set_tuberosity(data, 'A', 'replaced', point=(10.4, 20.6, 3), auto_point=(9, 20, 3), reviewer='x')
            review.set_tuberosity(data, 'B', 'not_in_view', reviewer='x')
            review.save(f, data)
            back = review.load(f)
            self.assertEqual(back['cases']['A']['tuberosity']['point'], [10.4, 20.6, 3])        # x, y kept sub-voxel
            self.assertNotIn('point', back['cases']['B']['tuberosity'])
            with self.assertRaises(ValueError):
                review.set_tuberosity(back, 'C', 'replaced')            # needs a point
            review.clear_tuberosity(back, 'A')
            self.assertNotIn('A', back['cases'])

    def test_tuberosity_on_axial_and_misaligned(self):
        from kneemeasure import review
        from kneemeasure.landmarks import Landmarks
        data = review.load(None)
        review.set_tuberosity(data, 'A', 'replaced', point=(40, 12, -9), source='axial')   # below the segmentation
        review.set_tuberosity(data, 'B', 'misaligned', point=(40, 12, -9), source='axial')
        self.assertEqual(data['cases']['A']['tuberosity']['source'], 'axial')
        self.assertNotIn('point', data['cases']['B']['tuberosity'])
        self.assertNotIn('source', data['cases']['B']['tuberosity'])
        self.assertEqual(Landmarks({}, 'L', review.overrides_for(data, 'A')).tuberosity['point'], (40, 12, -9))
        self.assertIsNone(Landmarks({}, 'L', review.overrides_for(data, 'B')).tuberosity['point'])

    def test_landmark_override(self):
        from kneemeasure.landmarks import Landmarks
        lm = Landmarks({}, 'L', {'tuberosity': {'status': 'not_in_view'}})
        self.assertIsNone(lm.tuberosity['point'])
        lm = Landmarks({}, 'L', {'tuberosity': {'status': 'replaced', 'point': [1, 2, 3]}})
        self.assertEqual(lm.tuberosity['point'], (1, 2, 3))
        self.assertEqual(lm.tuberosity['review'], 'replaced')

    def test_tpcl_override_and_store(self):
        from kneemeasure import review
        from kneemeasure.landmarks import Landmarks
        line = {'z': 7, 'theta': 3.5, 'medial_pt': [10.0, 20.0], 'lateral_pt': [40.0, 22.0]}
        data = review.load(None)
        review.set_tpcl(data, 'A', 'replaced', line, auto={'z': 5, 'theta': 1.0}, method='slice', reviewer='x')
        o = review.overrides_for(data, 'A')
        lm = Landmarks({}, 'R', o)
        self.assertEqual(lm.tpcl['z'], 7)
        self.assertEqual(lm.tpcl['review'], 'replaced')
        self.assertAlmostEqual(lm.tpcl['theta'], 3.5)
        review.clear_review(data, 'A', 'tpcl')
        self.assertNotIn('A', data['cases'])

    def test_sulcus_override_and_store(self):
        from kneemeasure import review
        from kneemeasure.landmarks import Landmarks
        ref = {'z': 9, 'sulcus': (30, 12), 'peak1': (10, 5), 'peak2': (50, 6)}
        data = review.load(None)
        review.set_sulcus(data, 'A', 'replaced', ref, auto={'z': 11, 'sulcus': (31, 12)}, method='points')
        lm = Landmarks({}, 'L', review.overrides_for(data, 'A'))
        self.assertEqual(lm.reference_slice['point'], (30.0, 12.0, 9))
        self.assertEqual(lm.reference_slice['review'], 'replaced')
        self.assertIsNone(lm.reference_slice['bone_slice'])                # no masks loaded (the review app)
        review.clear_review(data, 'A', 'sulcus')
        self.assertNotIn('A', data['cases'])

    def test_tpcl_from_points_sides(self):
        from kneemeasure.tibial_pcl import tpcl_from_points
        case = {'spacing': (1.0, 1.0, 1.0)}
        r = tpcl_from_points(case, 'R', (10, 50), (60, 55), 4)
        self.assertEqual(r['medial_pt'], (60.0, 55.0))          # RIGHT knee: medial = larger x
        l = tpcl_from_points(case, 'L', (10, 50), (60, 55), 4)
        self.assertEqual(l['medial_pt'], (10.0, 50.0))
        self.assertAlmostEqual(r['theta'], l['theta'])

    def test_cache_key_changes_with_review(self):
        from kneemeasure.cache import cache_key
        a = cache_key('h', 'L', None, True, True)
        b = cache_key('h', 'L', None, True, True, {'tuberosity': {'status': 'not_in_view'}})
        self.assertNotEqual(a, b)


def _write_nifti(path, data, srow, descrip=b''):
    """Minimal NIfTI-1 writer for tests (int16, sform)."""
    import gzip
    import struct
    hdr = bytearray(352)
    struct.pack_into('<i', hdr, 0, 348)
    struct.pack_into('<8h', hdr, 40, 3, *data.shape, 1, 1, 1, 1)
    struct.pack_into('<h', hdr, 70, 4)
    struct.pack_into('<h', hdr, 72, 16)
    struct.pack_into('<8f', hdr, 76, 1, 1, 1, 1, 0, 0, 0, 0)
    struct.pack_into('<f', hdr, 108, 352.0)
    struct.pack_into('<f', hdr, 112, 1.0)
    hdr[148:148 + len(descrip)] = descrip
    struct.pack_into('<h', hdr, 254, 1)
    struct.pack_into('<12f', hdr, 280, *np.asarray(srow, dtype=float).ravel())
    hdr[344:348] = b'n+1\0'
    with gzip.open(path, 'wb') as f:
        f.write(bytes(hdr) + data.astype('<i2').tobytes(order='F'))


class TestScansAndProject(unittest.TestCase):
    def _seg(self, path, shape=(6, 5, 4), laterality=None):
        import nrrd
        h = {'space': 'left-posterior-superior', 'space directions': np.eye(3), 'space origin': np.zeros(3)}
        if laterality:
            h['laterality'] = laterality
        nrrd.write(path, np.zeros(shape, dtype=np.uint8), h)

    def test_nifti_reader_and_grid(self):
        import nrrd
        from kneemeasure.scans import read_nifti, same_grid, laterality_from_nifti, DisplayScan
        with tempfile.TemporaryDirectory() as d:
            data = np.arange(6 * 5 * 4).reshape(6, 5, 4)
            srow = [[-1, 0, 0, 0], [0, -1, 0, 0], [0, 0, 1, 0]]          # RAS of an LPS identity grid
            _write_nifti(os.path.join(d, 's.nii.gz'), data, srow, b'laterality=R')
            arr, h = read_nifti(os.path.join(d, 's.nii.gz'))
            self.assertTrue(np.array_equal(arr, data))
            self.assertEqual(laterality_from_nifti(os.path.join(d, 's.nii.gz')), 'R')
            self._seg(os.path.join(d, 'c.nrrd'))
            seg_h = nrrd.read_header(os.path.join(d, 'c.nrrd'))
            self.assertTrue(same_grid(h, seg_h)[0])
            shifted = [[-1, 0, 0, 5], [0, -1, 0, 0], [0, 0, 1, 0]]
            _write_nifti(os.path.join(d, 't.nii.gz'), data, shifted)
            self.assertFalse(same_grid(read_nifti(os.path.join(d, 't.nii.gz'))[1], seg_h)[0])
            ds = DisplayScan(os.path.join(d, 's.nii.gz'), seg_h, (6, 5, 4))
            self.assertEqual(ds.slice('axial', 2).shape, (6, 5))

    def test_axial_display_reaches_below(self):
        import nrrd
        from kneemeasure.scans import AxialDisplay, find_axial
        with tempfile.TemporaryDirectory() as d:
            self._seg(os.path.join(d, 'c.nrrd'))                         # LPS identity grid, z = 0..3
            seg_h = nrrd.read_header(os.path.join(d, 'c.nrrd'))
            data = np.broadcast_to(np.arange(1, 8) * 10, (6, 5, 7))      # one value per slice
            _write_nifti(os.path.join(d, 'axial.nii.gz'), data, [[-1, 0, 0, 0], [0, -1, 0, 0], [0, 0, 1, -3]])
            open(os.path.join(d, 'water.nii.gz'), 'wb').close()
            path = find_axial(os.path.join(d, 'water.nii.gz'))
            ax = AxialDisplay(path, seg_h, (6, 5, 4))
            self.assertEqual(ax.z_below, 3)                              # the axial starts at z = -3
            self.assertTrue(np.array_equal(ax.slice('axial', -2), ax.raw[:, :, 1]))
            self.assertTrue(np.array_equal(ax.slice('axial', 2), ax.raw[:, :, 5]))
            self.assertEqual(ax.slice('sagittal', 1).shape, (5, 7))
            self.assertIsNone(find_axial(None))

    def test_project_side_precedence(self):
        from kneemeasure.project import Project
        with tempfile.TemporaryDirectory() as d:
            segs = os.path.join(d, 'segs')
            os.makedirs(segs)
            self._seg(os.path.join(segs, 'AAAAA.nrrd'), laterality='L')
            self._seg(os.path.join(segs, 'BBBBB.nrrd'))
            p = Project(segs_dir=segs)
            self.assertEqual(p.sides['AAAAA'], ('L', 'segmentation header'))
            self.assertEqual(p.sides['BBBBB'], (None, None))
            p.set_side('BBBBB', 'right')
            p.set_side('AAAAA', 'R')                                      # entered in app wins
            self.assertEqual(p.sides['BBBBB'], ('R', 'entered in app'))
            self.assertEqual(p.sides['AAAAA'], ('R', 'entered in app'))
            self.assertFalse(p.out.startswith(d))                       # the app's results stay out of the data folder


class TestPatellaWidth(unittest.TestCase):
    """A synthetic patella: an ellipse 45 x 20 mm (bone) with a cartilage
    band on its posterior surface, between its medial and lateral thirds,
    tilted so the lateral end moves posterior; PBCL along x."""
    SP, A, B, C = 0.9, 22.5, 10.0, 15.0             # voxel mm; bone semi-axes; cartilage half-width (mm)

    def _inside(self, x, y, tilt, scale=1.0):
        """Unrotated (patella-frame) coordinates of mm points: bone, cartilage masks."""
        t = np.radians(tilt)
        xp, yp = x * np.cos(t) + y * np.sin(t), -x * np.sin(t) + y * np.cos(t)
        a, b = self.A * scale, self.B * scale
        bone = (xp / a) ** 2 + (yp / b) ** 2 <= 1
        cart = ~bone & ((xp / a) ** 2 + (yp / (b + 3 * scale)) ** 2 <= 1) & (yp > 0) & (np.abs(xp) <= self.C * scale)
        return bone, cart

    class _Landmarks:
        """The two landmarks patella_width uses: a PBCL along x, and the patella axis (computed when asked)."""
        def __init__(self, case, side):
            self.case, self.knee_side = case, side
            self.pbcl = {'theta': 0.0, 'z': 2, 'medial_pt': (40.0, 100.0), 'lateral_pt': (80.0, 100.0)}

        @property
        def patella_axis(self):
            from kneemeasure.patella_width import find_patella_axis
            return find_patella_axis(self.case, self.knee_side)

        @property
        def patella_effective(self):
            from kneemeasure.patella_width import find_effective_width
            return find_effective_width(self.case, self.pbcl, self.knee_side)

    def _case(self, tilt, side):
        n, nz = 120, 7
        x, y = np.meshgrid((np.arange(n) - n / 2) * self.SP, (np.arange(n) - n / 2) * self.SP, indexing='ij')
        if side == 'R':
            x = -x                                       # lateral = smaller x for a right knee
        labels = np.zeros((n, n, nz), dtype=np.uint8)
        for z, scale in ((1, 0.6), (2, 0.8), (3, 1.0), (4, 0.8), (5, 0.6)):
            bone, cart = self._inside(x, y, tilt, scale)
            labels[:, :, z][cart] = 4
            labels[:, :, z][bone] = 3
        case = {'labels': labels, 'spacing': (self.SP, self.SP, self.SP), 'patella': None}
        return case, self._Landmarks(case, side)

    def _expected_effective(self, tilt):
        """The continuous cartilage band's extent along x at this tilt (dense sampling)."""
        g = np.arange(-30, 30, 0.05)
        x, y = np.meshgrid(g, g, indexing='ij')
        _, cart = self._inside(x, y, tilt)
        return np.ptp(x[cart])

    def test_widths_and_tilt(self):
        from kneemeasure import patella_width
        for side in ('L', 'R'):
            got = {}
            for tilt in (0, 10, 25):
                case, lm = self._case(tilt, side)
                r = patella_width.compute(case, lm)
                v = r['values']
                got[tilt] = v
                with self.subTest(side=side, tilt=tilt):
                    self.assertAlmostEqual(v['patellar_tilt_deg'], tilt, delta=1.0)
                    self.assertAlmostEqual(v["PW_max_mm"], 2 * self.A, delta=1.01 * self.SP)  # within a voxel of 45 mm
                    self.assertAlmostEqual(v['PW_effective_mm'], self._expected_effective(tilt), delta=1.5 * self.SP)
                    self.assertEqual(r['details']['max']['z'], 3)                  # the widest slice, for the QC
                    self.assertEqual(r['details']['effective']['z'], 3)
                    self.assertIsNone(r['note'])
            with self.subTest(side=side):
                pw_max = [got[t]['PW_max_mm'] for t in (0, 10, 25)]
                self.assertLess(max(pw_max) - min(pw_max), 1.5 * self.SP)            # PW_max doesn't change with tilt
                self.assertGreater(got[0]['PW_effective_mm'], got[10]['PW_effective_mm'])
                self.assertGreater(got[10]['PW_effective_mm'], got[25]['PW_effective_mm'])

    def test_edge_to_edge(self):
        from kneemeasure.patella_width import effective_width, max_chord
        m = np.zeros((20, 20), dtype=bool)
        m[5:15, 10] = True                               # a bar 10 voxels long: centre to centre 9
        self.assertAlmostEqual(effective_width(m, (1.0, 1.0), np.array([1.0, 0.0]), 'L')[0], 10.0)
        self.assertAlmostEqual(max_chord(m, (1.0, 1.0), 'L')[0], 10.0)             # along the bar, not its diagonal
        u = np.array([np.cos(np.radians(30)), np.sin(np.radians(30))])
        self.assertAlmostEqual(effective_width(m, (1.0, 1.0), u, 'L')[0], 9 * u[0] + u[0] + u[1])   # + |cos| + |sin| voxels

    def test_flags(self):
        from kneemeasure import patella_width
        case, lm = self._case(0, 'L')
        case['labels'][95:110, 95:110, 2:5] = 3                                      # a separate bony piece, ~25%
        note = patella_width.compute(case, lm)['note']
        self.assertIn('bipartite', note)
        self.assertIn('more than one piece on the PW_max slice', note)
        case, lm = self._case(0, 'L')
        case['labels'][0, 60, 3] = 3                                                 # touches the scan's edge
        self.assertIn('edge of the scan', patella_width.compute(case, lm)['note'])
        case, lm = self._case(0, 'L')
        case['labels'][case['labels'] == 4] = 0
        r = patella_width.compute(case, lm)
        self.assertIn('no patellar cartilage', r['note'])
        self.assertIsNone(r['values']['PW_effective_mm'])
        self.assertIsNotNone(r['values']['PW_max_mm'])


class TestFemoralWidth(unittest.TestCase):
    def _run(self, theta, side='L', at_edge=False):
        from types import SimpleNamespace
        from kneemeasure import femoral_width
        bone = np.zeros((100, 80, 3), dtype=bool)
        bone[15:85, 20:60, 1] = True                    # 70 x 40 voxels of 0.9 mm on the PBCL slice
        if at_edge:
            bone[0:85, 20:60, 1] = True
        case = {'bone': bone, 'spacing': (0.9, 0.9, 0.9)}
        lm = SimpleNamespace(knee_side=side, pbcl={'theta': theta, 'z': 1, 'medial_pt': (20, 59), 'lateral_pt': (80, 59)})
        return femoral_width.compute(case, lm)

    def test_width_edge_to_edge(self):
        self.assertAlmostEqual(self._run(0.0)['values']['FW_mm'], 70 * 0.9)
        t = np.radians(20.0)                                # a box's extent along a line at 20 deg to it
        self.assertAlmostEqual(self._run(20.0, 'R')['values']['FW_mm'], (70 * np.cos(t) + 40 * np.sin(t)) * 0.9)
        r = self._run(0.0, 'L')['details']
        self.assertGreater(r['lateral_edge'][0], r['medial_edge'][0])     # lateral = larger x for a left knee

    def test_flags(self):
        self.assertIsNone(self._run(0.0)['note'])
        self.assertIn('edge of the scan', self._run(0.0, at_edge=True)['note'])


class TestBoneTrochlea(unittest.TestCase):
    """A synthetic trochlea: a V groove (half-width A, depth D) on the femur
    bone, its facets falling away at 45 deg beyond the peaks, under
    cartilage of uniform thickness T. PBCL along x; the tuberosity 10 mm
    lateral of the groove. Expected values come from the geometry: the bony
    points are the V's corners; the cartilage peaks sit T above the bony
    peaks and the cartilage sulcus T / sin(half-angle) above the bony one."""
    SP, A, D, T = 0.25, 12.0, 8.0, 4.0

    def _case(self, side):
        from scipy.ndimage import distance_transform_edt
        n, nz = 200, 5
        x = (np.arange(n) - n / 2) * self.SP
        y = (np.arange(n) - 60) * self.SP                # the bony groove at y = 0; anterior = smaller y
        X, Y = np.meshgrid(x, y, indexing='ij')
        ax = np.abs(X)
        h = np.where(ax <= self.A, self.D * ax / self.A, self.D - (ax - self.A))     # height of the bone surface
        bone2 = (ax <= self.A + self.D) & (Y >= -h) & (Y <= 30)
        cart2 = ~bone2 & (distance_transform_edt(~bone2) * self.SP <= self.T) & (ax <= self.A + 3) & (Y < -h)
        bone = np.zeros((n, n, nz), dtype=bool)
        cart = np.zeros_like(bone)
        bone[:, :, 1:4] = bone2[:, :, None]
        cart[:, :, 1:4] = cart2[:, :, None]
        lateral = 10.0 if side == 'L' else -10.0         # lateral = larger x for a left knee
        tub = (n / 2 + lateral / self.SP, 60 + 40 / self.SP, 1)
        return {'bone': bone, 'cartilage': cart, 'spacing': (self.SP, self.SP, 2.0)}, tub

    def _run(self, side, overrides=None):
        from kneemeasure import trochlear_morphology, tt_tg
        from kneemeasure.landmarks import Landmarks
        case, tub = self._case(side)
        lm = Landmarks(case, side, overrides)
        lm.__dict__['pbcl'] = {'theta': 0.0, 'z': 1, 'medial_pt': (40.0, 150.0), 'lateral_pt': (160.0, 150.0)}
        lm.__dict__['tuberosity'] = {'point': tub, 'note': None}
        t, m = tt_tg.compute(case, lm), trochlear_morphology.compute(case, lm)
        return dict(t['values'], **m['values']), [t['note'], m['note']]

    def test_bone_and_cartilage_geometry(self):
        A, D, T = self.A, self.D, self.T
        half = np.arctan(A / D)                          # the V's half-angle at the sulcus
        rise = T / np.sin(half)                          # cartilage sulcus above the bony one
        d_cart = D + T - rise                            # cartilage peaks T above the bony peaks
        for side in ('L', 'R'):
            v, notes = self._run(side)
            with self.subTest(side=side):
                self.assertEqual(notes, [None, None])
                self.assertAlmostEqual(v['sulcus_angle_bone_deg'], np.degrees(2 * half), delta=1.0)
                self.assertAlmostEqual(v['sulcus_angle_deg'], np.degrees(2 * np.arctan(A / d_cart)), delta=1.0)
                self.assertGreater(v['sulcus_angle_deg'] - v['sulcus_angle_bone_deg'], 3.0)  # not equal, even uniform
                self.assertAlmostEqual(v['LTI_bone_deg'], np.degrees(np.arctan(D / A)), delta=1.0)
                self.assertAlmostEqual(v['LTI_deg'], np.degrees(np.arctan(d_cart / A)), delta=1.0)
                self.assertAlmostEqual(v['trochlear_depth_bone_mm'], D, delta=2 * self.SP)
                self.assertAlmostEqual(v['trochlear_depth_mm'], d_cart, delta=2 * self.SP)
                self.assertAlmostEqual(v['lateral_facet_length_bone_mm'], np.hypot(A, D), delta=2 * self.SP)
                self.assertAlmostEqual(v['lateral_facet_length_mm'], np.hypot(A, d_cart), delta=2 * self.SP)
                # the voxelised V's floor is a short flat run, so the sulcus can sit a voxel off centre (~3%)
                self.assertAlmostEqual(v['facet_asymmetry_bone_ratio'], 1.0, delta=0.05)
                self.assertAlmostEqual(v['facet_asymmetry_ratio'], 1.0, delta=0.05)
                # the groove moves straight back, so TT-TG is the same on both surfaces
                self.assertAlmostEqual(v['TT_TG_bone_mm'], 10.0, delta=self.SP)
                self.assertAlmostEqual(v['TT_TG_mm'], v['TT_TG_bone_mm'], delta=self.SP)
                # ... and sTT-TG differs by how far the cartilage groove sits in front of the bony one
                self.assertAlmostEqual(v['sTT_TG_bone_mm'] - v['sTT_TG_mm'], rise, delta=2 * self.SP)

    def test_reviewed_sulcus(self):
        from kneemeasure.trochlea import find_reference_slice
        case, _ = self._case('L')
        ref = find_reference_slice(case)
        sulcus = {'z': ref['z'], 'sulcus': [float(v) for v in ref['sulcus']],
                  'peaks': [[float(v) for v in ref['peak1']], [float(v) for v in ref['peak2']]]}
        auto, _ = self._run('L')
        for status in ('accepted', 'replaced'):
            v, notes = self._run('L', {'sulcus': dict(sulcus, status=status)})
            with self.subTest(status=status):          # the bone points stay automatic on the reviewed slice
                self.assertAlmostEqual(v['sulcus_angle_bone_deg'], auto['sulcus_angle_bone_deg'])
                self.assertAlmostEqual(v['TT_TG_bone_mm'], auto['TT_TG_bone_mm'])
                self.assertEqual(notes, [None, None])
        # clicked cartilage points on another slice: the bone points move to that slice
        clicked = dict(sulcus, status='replaced', z=ref['z'] - 1, sulcus=[100.0, 70.0], peaks=[[52.0, 30.0], [148.0, 30.0]])
        v, _ = self._run('L', {'sulcus': clicked})
        self.assertNotAlmostEqual(v['sulcus_angle_deg'], auto['sulcus_angle_deg'], places=1)   # cartilage = the clicks
        self.assertAlmostEqual(v['sulcus_angle_bone_deg'], auto['sulcus_angle_bone_deg'])      # same shape on z - 1

    def test_column_order(self):
        from kneemeasure.pipeline import MEASURES, value_columns
        cols = value_columns(MEASURES)
        for c in cols:
            if '_bone_' in c:
                name, unit = c.split('_bone_')
                self.assertEqual(cols[cols.index(c) - 1], f'{name}_{unit}')          # right after its cartilage version


class TestTrochlearWidth(unittest.TestCase):
    """A synthetic knee: the V-groove trochlea of TestBoneTrochlea on slices
    1-3, its anterior cartilage reaching 17, 16 and 15 mm lateral of the
    groove (so the most lateral slice is the lowest, and the TT-TG groove
    slice, the top one, reaches 15 mm), and a band of patellar cartilage
    30 mm wide on slices 6-7, its medial edge 17 mm medial of the groove.
    PBCL along x. Edge to edge adds half a voxel at each end."""
    SP, A, D, T = 0.25, 12.0, 8.0, 4.0

    def _case(self, side, shift=0.0, patella_centre=-2.0, femoral_cartilage=True):
        from scipy.ndimage import distance_transform_edt
        n, nz = 240, 9
        x = (np.arange(n) - n / 2) * self.SP
        y = (np.arange(n) - 100) * self.SP
        X, Y = np.meshgrid(x, y, indexing='ij')
        lateral = X if side == 'L' else -X                # lateral = larger x for a left knee
        ax = np.abs(X)
        h = np.where(ax <= self.A, self.D * ax / self.A, self.D - (ax - self.A))
        bone2 = (ax <= self.A + self.D) & (Y >= -h) & (Y <= 30)
        near = ~bone2 & (distance_transform_edt(~bone2) * self.SP <= self.T) & (Y < -h)
        labels = np.zeros((n, n, nz), dtype=np.uint8)
        bone, cart = np.zeros(labels.shape, dtype=bool), np.zeros(labels.shape, dtype=bool)
        for z, reach in ((1, 5.0), (2, 4.0), (3, 3.0)):
            bone[:, :, z] = bone2
            if femoral_cartilage:
                cart[:, :, z] = near & (ax <= self.A + reach)
        c = patella_centre + shift
        patella = (lateral >= c - 15) & (lateral <= c + 15) & (Y >= -24) & (Y <= -20)
        for z in (6, 7):
            labels[:, :, z][patella] = 4
        labels[bone], labels[cart] = 1, 2
        return {'labels': labels, 'bone': bone, 'cartilage': cart, 'spacing': (self.SP, self.SP, 2.0)}

    def _run(self, case, side):
        from kneemeasure import trochlear_width
        from kneemeasure.landmarks import Landmarks
        lm = Landmarks(case, side)
        lm.__dict__['pbcl'] = {'theta': 0.0, 'z': 1, 'medial_pt': (60.0, 150.0), 'lateral_pt': (180.0, 150.0)}
        return trochlear_width.compute(case, lm)

    def test_values_sides_and_shift(self):
        half = self.SP / 2
        for side in ('L', 'R'):
            with self.subTest(side=side):
                r = self._run(self._case(side), side)
                self.assertIsNone(r['note'])
                self.assertAlmostEqual(r['values']['TW_mm'], (17 + half) + (17 + half))     # the lowest slice
                self.assertAlmostEqual(r['values']['APTO_mm'], (15 + half) + (17 + half))   # the TT-TG groove slice
                self.assertEqual(r['details']['TW_mm']['z'], 1)
                self.assertEqual(r['details']['APTO_mm']['z'], 3)
                self.assertAlmostEqual(r['details']['APTO_mm']['offset_mm'], -6.0)         # 3 slices below the patella's
                for d in (1.0, 2.5):                     # the patella moved laterally by d: both shrink by d
                    moved = self._run(self._case(side, shift=d), side)['values']
                    self.assertAlmostEqual(moved['TW_mm'], r['values']['TW_mm'] - d)
                    self.assertAlmostEqual(moved['APTO_mm'], r['values']['APTO_mm'] - d)

    def test_flags(self):
        r = self._run(self._case('L', patella_centre=35.0), 'L')                  # medial border lateral of the trochlea
        self.assertLessEqual(r['values']['TW_mm'], 0)
        self.assertIn("patella's medial border is lateral of the trochlear edge", r['note'])
        r = self._run(self._case('L', femoral_cartilage=False), 'L')
        self.assertIsNone(r['values']['TW_mm'])
        self.assertIn('no femoral cartilage in the trochlear region', r['note'])
        case = self._case('L')
        edge = 120 + int(16 / self.SP)                    # cut the scan at 16 mm lateral: the cartilage meets its edge
        for k in ('labels', 'bone', 'cartilage'):
            case[k] = case[k][:edge + 1]
        self.assertIn('edge of the field of view', self._run(case, 'L')['note'])


@unittest.skipUnless(os.path.isdir(DATA), f'segmentation data not found at {DATA}')
class TestRegression(unittest.TestCase):
    def test_tpcl_at_slice_matches_method(self):
        from kneemeasure.tibial_pcl import tpcl_at_slice
        from kneemeasure.volume import load_case
        with open(os.path.join(HERE, 'expected_values.json')) as f:
            e = next(iter(json.load(f).values()))
        path = os.path.join(DATA, e['file'])
        if not os.path.exists(path):
            self.skipTest('case missing')
        case = load_case(path)
        zs = np.where(case['tibia'].any(axis=(0, 1)))[0]
        line = tpcl_at_slice(case, e['knee_side'], int(zs.mean()))
        self.assertIsNotNone(line)
        self.assertLess(abs(line['theta']), 45)

    def test_sulcus_on_slice_matches_reference(self):
        from kneemeasure.trochlea import find_reference_slice, sulcus_on_slice
        from kneemeasure.volume import load_case
        with open(os.path.join(HERE, 'expected_values.json')) as f:
            e = next(iter(json.load(f).values()))
        path = os.path.join(DATA, e['file'])
        if not os.path.exists(path):
            self.skipTest('case missing')
        case = load_case(path, patella=False)
        ref = find_reference_slice(case)
        z = ref['z']
        again = sulcus_on_slice(case['bone'][:, :, z], case['cartilage'][:, :, z], case['spacing'][0], z)
        self.assertEqual(again['point'], ref['point'])

    def test_expected_values(self):
        from kneemeasure.pipeline import measure_case
        with open(os.path.join(HERE, 'expected_values.json')) as f:
            expected = json.load(f)
        for case, e in expected.items():
            path = os.path.join(DATA, e['file'])
            if not os.path.exists(path):
                continue
            with self.subTest(case=case):
                r = measure_case(path, e['knee_side'], qc=False, landmarks=False)
                for col, want in e['values'].items():
                    got = r['values'].get(col)
                    if want is None:
                        self.assertIsNone(got, col)
                    else:
                        self.assertAlmostEqual(got, want, places=3, msg=col)


if __name__ == '__main__':
    unittest.main()
