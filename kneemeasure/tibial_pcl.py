"""
Tibial posterior condylar line (TPCL): the tibial analogue of the femoral
PBCL, used as the rotational reference axis for KRA and tubercle
lateralization.

One slice, no averaging:

1. The plateau: scanning down from the top of the segmented tibia, the
   first slice whose cross-section reaches plateau_area_fraction of the
   tibia's largest one. Above it the slices cut through the tibial spines
   and the articular surfaces -- small and irregular, with false notches.
2. The slice: from the plateau down, the first (most cranial) slice where
   the posterior notch is clearly present -- the posterior contour has a
   concavity at least min_notch_depth_mm deep between two convex condylar
   lobes. Each lobe must stand at least min_lobe_mm posterior of the
   contour's end on its own side and sit at least min_lobe_mm in from that
   end (so a contour that just curves forward at the bone's edge isn't
   counted as a lobe). The slice must be one piece and mostly convex.
3. The line: the bitangent of the posterior contour across the notch --
   the edge of the contour's posterior convex hull that spans it, touching
   both lobes.

Validated against manual KRA on the ORS cohort (45 knees): ICC(2,1) 0.82,
up from 0.79 for the previous method (the 5 deepest-notch slices,
averaged). Plateau fractions 0.85-0.95 all gave 0.81-0.84.
"""
import numpy as np
import cv2

from .pbcl import posterior_bitangent, _medial_lateral_from_peaks, _theta_from_points
from .volume import true_physical_lp, count_components


def _solidity(mask2d):
    """mask area / convex-hull area: ~0.95 for a clean tibial cross-section
    (its natural posterior notch costs a little), much lower (~0.6) for a
    folded, multi-lobed badly-segmented slice."""
    mask_u8 = mask2d.astype(np.uint8) * 255
    contours, _ = cv2.findContours(mask_u8, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    if not contours:
        return 0.0
    c = max(contours, key=cv2.contourArea)
    area = cv2.contourArea(c)
    hull_area = cv2.contourArea(cv2.convexHull(c))
    return area / hull_area if hull_area > 0 else 0.0


def find_tibial_pcl(case, knee_side, plateau_area_fraction=0.9, min_notch_depth_mm=2.0, min_lobe_mm=3.0,
                     min_solidity=0.85):
    """Returns {'theta', 'theta_true', 'z', 'medial_pt', 'lateral_pt',
    'notch_depth_mm', 'plateau_z'}, or {'theta': None, 'note': ...} on
    failure."""
    tibia = case['tibia']
    spacing_x, spacing_y = case['spacing'][0], case['spacing'][1]
    result = {'theta': None, 'theta_true': None, 'z': None, 'medial_pt': None,
              'lateral_pt': None, 'notch_depth_mm': None, 'plateau_z': None}

    areas = tibia.sum(axis=(0, 1))
    z_nonempty = np.where(areas > 0)[0]
    if len(z_nonempty) == 0:
        result['note'] = 'no tibia bone found'
        return result
    z_top = int(z_nonempty.max())
    plateau_z = next(z for z in range(z_top, -1, -1) if areas[z] >= plateau_area_fraction * areas.max())
    result['plateau_z'] = plateau_z

    for z in range(plateau_z, int(z_nonempty.min()) - 1, -1):
        sl = tibia[:, :, z]
        _, comps = count_components(sl, min_component_size=0)
        if len(comps) != 1 or _solidity(sl) < min_solidity:
            continue
        notch = posterior_bitangent(sl, spacing_x, spacing_y, min_lobe_mm)
        if notch is None or notch['depth_mm'] < min_notch_depth_mm:
            continue
        medial_pt, lateral_pt = _medial_lateral_from_peaks([notch['p1'], notch['p2']], knee_side)
        result.update(theta=_theta_from_points(medial_pt, lateral_pt, spacing_x, spacing_y), z=z,
                      medial_pt=medial_pt, lateral_pt=lateral_pt, notch_depth_mm=notch['depth_mm'])
        med_L, med_P = true_physical_lp(case['header'], medial_pt[0], medial_pt[1], z)
        lat_L, lat_P = true_physical_lp(case['header'], lateral_pt[0], lateral_pt[1], z)
        p1, p2 = sorted([(float(med_L), float(med_P)), (float(lat_L), float(lat_P))], key=lambda p: p[0])
        result['theta_true'] = float(np.degrees(np.arctan2(p2[1] - p1[1], p2[0] - p1[0])))
        return result

    result['note'] = 'no slice below the plateau shows a posterior notch between two condylar lobes'
    return result


def tpcl_at_slice(case, knee_side, z, min_lobe_mm=3.0):
    """The TPCL on one chosen axial slice (a reviewer's pick): the same
    bitangent, without the plateau and notch-depth rules. Returns {'z',
    'medial_pt', 'lateral_pt', 'theta'} or None if this slice doesn't show
    a notch between two condylar lobes."""
    spacing_x, spacing_y = case['spacing'][0], case['spacing'][1]
    z = int(z)
    if not 0 <= z < case['tibia'].shape[2]:
        return None
    notch = posterior_bitangent(case['tibia'][:, :, z], spacing_x, spacing_y, min_lobe_mm)
    if notch is None:
        return None
    medial_pt, lateral_pt = _medial_lateral_from_peaks([notch['p1'], notch['p2']], knee_side)
    return {'z': z, 'medial_pt': medial_pt, 'lateral_pt': lateral_pt,
            'theta': _theta_from_points(medial_pt, lateral_pt, spacing_x, spacing_y)}


def tpcl_from_points(case, knee_side, p1, p2, z):
    """A TPCL through two points a reviewer clicked on slice z; which point
    is medial follows the knee side (RIGHT knee: medial = larger x)."""
    a, b = (float(p1[0]), float(p1[1])), (float(p2[0]), float(p2[1]))
    if a[0] == b[0]:
        raise ValueError('the two points must be side by side, not one above the other')
    medial_pt, lateral_pt = _medial_lateral_from_peaks([a, b], knee_side)
    spacing_x, spacing_y = case['spacing'][0], case['spacing'][1]
    return {'z': int(z), 'medial_pt': tuple(medial_pt), 'lateral_pt': tuple(lateral_pt),
            'theta': _theta_from_points(medial_pt, lateral_pt, spacing_x, spacing_y)}
