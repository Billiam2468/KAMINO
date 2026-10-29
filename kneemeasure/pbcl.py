"""
Femoral posterior bicondylar line (PBCL), per Dejour, Giovannetti de
Sanctis et al. 2025 (KSSTA 33(8):2833-2846), Figure 2: the line tangent to the
subchondral bone at the posterior aspect of both femoral condyles.

On each axial femur slice the bone's POSTERIOR boundary is traced per
column, and the bitangent across the intercondylar notch is taken: the
edge of that contour's posterior convex hull spanning the notch, touching
both condyles (the same construction as the tibial PCL). Both contact
points must be real condylar lobes, standing out from the bone's edge on
their side. Every slice is tried; the one whose bitangent lies most
posterior wins, so the slice and the line are chosen together.

Previously the line joined the two condyles' local maxima (peak tops). On
the ORS cohort the bitangent gives the same agreement with manual KRA,
TT-TG, sTT-TG, PTLTR and LTI (median change 0.1 deg), and it is the
standard definition.

The helpers here (_posterior_contour, posterior_bitangent,
_medial_lateral_from_peaks, _theta_from_points) are shared with the
tibial PCL.

Medial/lateral: axis0 (x) increases toward patient LEFT, so for a RIGHT
knee medial = larger x; for a LEFT knee medial = smaller x.
"""
import numpy as np

from .volume import true_physical_lp


def _posterior_contour(mask2d):
    """For each x-column with tissue present, the most POSTERIOR (max-y)
    pixel."""
    cols = np.where(mask2d.any(axis=1))[0]
    if len(cols) == 0:
        return np.array([]), np.array([])
    y = np.array([np.max(np.where(mask2d[c])[0]) for c in cols], dtype=float)
    return cols, y


def _posterior_hull(px, py):
    """Indices of the posterior (max-y) convex hull chain of contour points
    sorted by x."""
    hull = []
    for i in range(len(px)):
        while len(hull) >= 2:
            o, a = hull[-2], hull[-1]
            if (px[a] - px[o]) * (py[i] - py[o]) - (py[a] - py[o]) * (px[i] - px[o]) >= 0:
                hull.pop()
            else:
                break
        hull.append(i)
    return hull


def posterior_bitangent(mask2d, spacing_x, spacing_y, min_lobe_mm=3.0):
    """The deepest concavity of the posterior contour that sits between two
    convex lobes, and the bitangent across it. Each lobe must stand at
    least min_lobe_mm posterior of the contour's end on its own side and
    sit at least min_lobe_mm in from that end. A gap with no tissue at all
    (two separate condyles) counts as an open notch, depth inf. Returns
    {'p1', 'p2' (the two contact points, (x, y) voxels), 'depth_mm'} or
    None."""
    x, y = _posterior_contour(mask2d)
    if len(x) < 5:
        return None
    px, py = x * spacing_x, y * spacing_y
    best = None
    hull = _posterior_hull(px, py)
    for a, b in zip(hull[:-1], hull[1:]):
        gap = x[b] - x[a] > b - a              # columns with no tissue between the contact points
        if b - a < 2 and not gap:
            continue
        if min(py[a] - py[0], py[b] - py[-1], px[a] - px[0], px[-1] - px[b]) < min_lobe_mm:
            continue
        if gap:
            depth = np.inf
        else:
            dx, dy = px[b] - px[a], py[b] - py[a]
            between = slice(a + 1, b)
            depth = float(np.max(np.abs((px[between] - px[a]) * dy - (py[between] - py[a]) * dx)) / np.hypot(dx, dy))
        if best is None or depth > best['depth_mm']:
            best = {'p1': (float(x[a]), float(y[a])), 'p2': (float(x[b]), float(y[b])), 'depth_mm': depth}
    return best


def _medial_lateral_from_peaks(top2, knee_side):
    """Label two condyle peaks by x-position + knee side: RIGHT knee
    medial = larger x; LEFT knee medial = smaller x."""
    sorted_by_x = sorted(top2, key=lambda p: p[0])
    smaller_x_pt, larger_x_pt = sorted_by_x[0][:2], sorted_by_x[1][:2]
    return (larger_x_pt, smaller_x_pt) if knee_side == 'R' else (smaller_x_pt, larger_x_pt)


def _theta_from_points(pt_a, pt_b, spacing_x, spacing_y):
    """Line angle computed smaller-x -> larger-x (not medial->lateral), so
    it doesn't depend on knee side or wrap near +-180deg -- a line has no
    inherent direction. Converts to physical mm per axis first, so the
    angle is right even when spacing_x != spacing_y."""
    p1, p2 = sorted([pt_a, pt_b], key=lambda p: p[0])
    dx_mm = (p2[0] - p1[0]) * spacing_x
    dy_mm = (p2[1] - p1[1]) * spacing_y
    return float(np.degrees(np.arctan2(dy_mm, dx_mm)))


def find_pbcl(case, knee_side, min_notch_depth_mm=2.0, min_lobe_mm=3.0):
    """Returns {'theta', 'theta_true', 'z', 'medial_pt', 'lateral_pt'}, or
    {'theta': None, ...} on failure.

    'theta' is in the reoriented (x, y) frame -- the angle every relative
    measurement uses. 'theta_true' is the PBCL's actual physical tilt
    relative to true Left/Posterior (they differ for files whose header
    has an oblique in-plane rotation)."""
    bone = case['bone']
    spacing_x, spacing_y = case['spacing'][0], case['spacing'][1]
    result = {'theta': None, 'theta_true': None, 'z': None, 'medial_pt': None, 'lateral_pt': None}

    best_y, best_z, best_med, best_lat = -np.inf, None, None, None
    for z in np.where(bone.any(axis=(0, 1)))[0]:
        notch = posterior_bitangent(bone[:, :, z], spacing_x, spacing_y, min_lobe_mm)
        if notch is None or notch['depth_mm'] < min_notch_depth_mm:
            continue
        medial_pt, lateral_pt = _medial_lateral_from_peaks([notch['p1'], notch['p2']], knee_side)
        y_avg = (medial_pt[1] + lateral_pt[1]) / 2.0
        if y_avg > best_y:
            best_y, best_z, best_med, best_lat = y_avg, int(z), medial_pt, lateral_pt

    if best_z is None:
        return result
    result['theta'] = _theta_from_points(best_med, best_lat, spacing_x, spacing_y)
    result['z'], result['medial_pt'], result['lateral_pt'] = best_z, best_med, best_lat

    med_L, med_P = true_physical_lp(case['header'], best_med[0], best_med[1], best_z)
    lat_L, lat_P = true_physical_lp(case['header'], best_lat[0], best_lat[1], best_z)
    p1, p2 = sorted([(float(med_L), float(med_P)), (float(lat_L), float(lat_P))], key=lambda p: p[0])
    result['theta_true'] = float(np.degrees(np.arctan2(p2[1] - p1[1], p2[0] - p1[0])))
    return result


def pcl_axes(theta_deg):
    """(u_along, u_anterior): unit vectors along a posterior condylar line
    and perpendicular to it (anterior-facing), in the reoriented mm frame.
    Anterior = smaller y, so u_anterior = (sin(theta), -cos(theta))."""
    theta = np.radians(theta_deg)
    return (np.array([np.cos(theta), np.sin(theta)]),
            np.array([np.sin(theta), -np.cos(theta)]))
