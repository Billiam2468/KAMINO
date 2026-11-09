"""
Per-case landmark cache. Several measurements share the same landmarks
(the PBCL feeds KRA, TT-TG, PTLTR and LTI; the Dejour reference slice
feeds TT-TG and LTI; the tuberosity feeds TT-TG and TTL). Each is
deterministic, so it is computed once per case here and shared, instead
of being recomputed by every measurement.

overrides: reviewed landmarks (see review.py) that replace the automatic
ones -- the tibial tuberosity, the tibial posterior condylar line and the
trochlear reference slice's sulcus and facet peaks.
"""
from functools import cached_property

from .patella_width import find_effective_width, find_patella_axis
from .pbcl import find_pbcl
from .tibial_pcl import find_tibial_pcl
from .tibial_tuberosity import find_tibial_tuberosity
from .trochlea import find_reference_slice, find_deepest_cartilage_groove, most_anterior_component, sulcus_on_slice


class Landmarks:
    def __init__(self, case, knee_side, overrides=None):
        self.case = case
        self.knee_side = knee_side
        self.overrides = overrides or {}

    @cached_property
    def pbcl(self):
        return find_pbcl(self.case, self.knee_side)

    @cached_property
    def tpcl_auto(self):
        return find_tibial_pcl(self.case, self.knee_side)

    @cached_property
    def tpcl(self):
        review = self.overrides.get('tpcl')
        if not review:
            return dict(self.tpcl_auto, review='auto')
        return {'theta': float(review['theta']), 'theta_true': None, 'z': int(review['z']),
                'medial_pt': tuple(review['medial_pt']), 'lateral_pt': tuple(review['lateral_pt']),
                'notch_depth_mm': None, 'review': review['status'], 'note': None}

    @cached_property
    def reference_slice_auto(self):
        return find_reference_slice(self.case, min_depth_mm=2.0)

    @cached_property
    def reference_slice(self):
        review = self.overrides.get('sulcus')
        if not review:
            return dict(self.reference_slice_auto, review='auto')
        z = int(review['z'])
        sulcus = tuple(float(v) for v in review['sulcus'])
        bone_sl = cart_sl = None
        if self.case.get('bone') is not None:            # for the QC image (the review app drops these)
            bone_sl = self.case['bone'][:, :, z]
            cart_sl = most_anterior_component(self.case['cartilage'][:, :, z], bone_sl)
        return {'z': z, 'depth_mm': None, 'sulcus': sulcus,
                'peak1': tuple(review['peaks'][0]), 'peak2': tuple(review['peaks'][1]),
                'point': (sulcus[0], sulcus[1], z), 'hull_points': None,
                'bone_slice': bone_sl, 'cartilage_component': cart_sl, 'note': None, 'review': review['status']}

    @cached_property
    def reference_slice_bone(self):
        """The bony sulcus and facet peaks: the reviewer's (a 'sulcus_bone'
        decision, on any slice), else reference_slice_bone_auto."""
        review = self.overrides.get('sulcus_bone')
        if not review:
            return dict(self.reference_slice_bone_auto, review='auto')
        z = int(review['z'])
        sulcus = tuple(float(v) for v in review['sulcus'])
        return {'z': z, 'depth_mm': None, 'sulcus': sulcus,
                'peak1': tuple(review['peaks'][0]), 'peak2': tuple(review['peaks'][1]),
                'point': (sulcus[0], sulcus[1], z), 'hull_points': None, 'note': None, 'review': review['status']}

    @cached_property
    def reference_slice_bone_auto(self):
        """The bony sulcus and facet peaks found automatically on the slice
        the cartilage values use (the same algorithm, on the femur bone):
        the reference slice, or the reviewer's slice when the cartilage
        sulcus was reviewed."""
        ref = self.reference_slice
        if ref['z'] is None:
            return {'z': None, 'note': ref['note']}
        if self.case.get('bone') is None:
            return {'z': None, 'note': 'no bone loaded'}
        z = ref['z']
        found = sulcus_on_slice(self.case['bone'][:, :, z], None, self.case['spacing'][0], z, surface='bone')
        return found or {'z': None, 'note': f'no bony groove on the reference slice (z {z})'}

    @cached_property
    def patella_axis(self):
        """The PW_max line on the patella's widest slice; its lateral end is
        PTLTR's patella point too."""
        return find_patella_axis(self.case, self.knee_side)

    @cached_property
    def patella_effective(self):
        """The PW_effective line (widest patellar cartilage along the PBCL);
        its medial edge is trochlear width's and APTO's patellar point too."""
        return find_effective_width(self.case, self.pbcl, self.knee_side)

    @cached_property
    def deepest_groove(self):
        return find_deepest_cartilage_groove(self.case)

    @cached_property
    def tuberosity_auto(self):
        return find_tibial_tuberosity(self.case, self.knee_side)

    @cached_property
    def tuberosity(self):
        review = self.overrides.get('tuberosity')
        if not review:
            return dict(self.tuberosity_auto, review='auto')
        if review['status'] == 'not_in_view':
            return {'point': None, 'note': 'marked not in view of the scan by reviewer', 'review': 'not in view'}
        if review['status'] == 'misaligned':
            return {'point': None, 'review': 'misaligned',
                    'note': "the axial series and the segmentation don't line up -- measure manually"}
        x, y, z = review['point']
        return {'point': (float(x), float(y), int(z)), 'note': None, 'source': review.get('source'),
                'review': 'accepted' if review['status'] == 'accepted' else 'replaced'}
