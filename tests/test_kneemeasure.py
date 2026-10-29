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
