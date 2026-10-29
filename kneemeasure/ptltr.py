"""
PTLTR: the mediolateral offset between the patella's lateral margin and
the lateral trochlear ridge, projected onto the PBCL (same "project both
landmarks onto the PBCL" convention as TT-TG). Unsigned.

Landmarks:
- Ridge point: the single point on the LATERAL bony trochlea with the
  greatest anterior projection beyond the PBCL, searched from the proximal
  trochlea down to the condyle split but only on slices that have a real
  (>= 2 mm) cartilaginous sulcus. That gate excludes the proximal shaft
  (whose natural anterior curvature otherwise wins) and the distal notch.
  "Lateral" is lateral of that slice's own sulcus.
- Patella point: the most lateral voxel of the patella bone (largest 3D
  component), tie-broken to the middle of a flat lateral facet.

    PTLTR = |(P_patella - P_ridge) . u_along_PBCL|
"""
import numpy as np

from . import qc_draw
from .pbcl import pcl_axes
from .trochlea import most_anterior_component, find_peaks_and_sulcus
from .volume import count_components

COLUMNS = [
    ('PTLTR_mm', 'mm', 'patella lateral margin to lateral trochlear ridge, along the PBCL (unsigned)'),
]


def find_lateral_patella_point(patella_mask, knee_side):
    """Most lateral voxel of the patella bone -- lateral = smaller x for a
    RIGHT knee, larger x for a LEFT knee. Among voxels tied at the extreme
    x, the one closest to their median (y, z)."""
    xs, ys, zs = np.where(patella_mask)
    if len(xs) == 0:
        return {'point': None, 'note': 'no patella bone found in this case'}
    target_x = xs.min() if knee_side == 'R' else xs.max()
    tied = xs == target_x
    ys_t, zs_t = ys[tied], zs[tied]
    y_med, z_med = np.median(ys_t), np.median(zs_t)
    best = int(np.argmin((ys_t - y_med) ** 2 + (zs_t - z_med) ** 2))
    return {'point': (int(target_x), int(ys_t[best]), int(zs_t[best])), 'note': None}


def find_lateral_bony_apex(case, pbcl_result, knee_side, min_depth_mm=2.0):
    """Lateral trochlear bony apex -- see module docstring. On each slice
    the lateral side is everything lateral of that slice's own cartilage
    sulcus (lateral = smaller x for a RIGHT knee). (Splitting at the
    midpoint of the PBCL's two points instead made the search reach onto
    the medial facet when those points shifted by a few voxels.)"""
    bone, cart = case['bone'], case['cartilage']
    spacing_x = case['spacing'][0]
    spacing_xy = np.array(case['spacing'][:2])
    lateral_sign = -1 if knee_side == 'R' else 1

    _, u_anterior = pcl_axes(pbcl_result['theta'])
    pbcl_anchor_mm = np.array(pbcl_result['medial_pt'], dtype=float) * spacing_xy

    z_nonempty = np.where(cart.any(axis=(0, 1)))[0]
    if len(z_nonempty) == 0:
        return {'point': None, 'note': 'no cartilage found in this case'}

    best = None
    for z in range(int(z_nonempty.max()), int(z_nonempty.min()) - 1, -1):
        bone_sl = bone[:, :, z]
        _, comps = count_components(bone_sl)
        if len(comps) >= 2:
            break
        cart_sl = cart[:, :, z]
        component = most_anterior_component(cart_sl, bone_sl)
        if component.sum() == 0:
            continue
        groove = find_peaks_and_sulcus(component, bone_sl, spacing_x, min_depth_mm=min_depth_mm)
        if groove is None:
            continue
        xs, ys = np.where(bone_sl)
        lateral_mask = np.sign(xs - groove['sulcus'][0]) == lateral_sign
        if not lateral_mask.any():
            continue
        xs_lat, ys_lat = xs[lateral_mask], ys[lateral_mask]
        pts_mm = np.stack([xs_lat, ys_lat], axis=1) * spacing_xy
        offsets = (pts_mm - pbcl_anchor_mm) @ u_anterior
        idx = int(np.argmax(offsets))
        if best is None or offsets[idx] > best['offset']:
            best = {'offset': float(offsets[idx]), 'x': int(xs_lat[idx]), 'y': int(ys_lat[idx]), 'z': z}

    if best is None:
        return {'point': None, 'note': 'no slice with a detectable cartilaginous sulcus found in scan range'}
    return {'point': (best['x'], best['y'], best['z']), 'anterior_offset_mm': best['offset'], 'note': None}


def compute(case, lm):
    pbcl = lm.pbcl
    if pbcl['theta'] is None:
        return {'values': {'PTLTR_mm': None}, 'note': 'PBCL not found'}
    apex = find_lateral_bony_apex(case, pbcl, lm.knee_side)
    if apex['point'] is None:
        return {'values': {'PTLTR_mm': None}, 'note': f"lateral bony apex not found -- {apex['note']}"}
    patella = find_lateral_patella_point(case['patella'], lm.knee_side)
    if patella['point'] is None:
        return {'values': {'PTLTR_mm': None}, 'note': patella['note']}

    spacing_xy = np.array(case['spacing'][:2])
    v = np.array(patella['point'][:2]) * spacing_xy - np.array(apex['point'][:2]) * spacing_xy
    u_along, _ = pcl_axes(pbcl['theta'])
    return {'values': {'PTLTR_mm': abs(float(np.dot(v, u_along)))}, 'note': None,
            'details': {'ridge_point': apex['point'], 'patella_point': patella['point']}}


def qc_payload(case, lm, result, margin_mm=8.0, patella_window_mm=25.0):
    """The patella sits a real ~2-3 cm anterior of the trochlea, so its
    displayed footprint is cropped to a window around the measured lateral
    point to keep the panel framed on what's measured."""
    if result['values']['PTLTR_mm'] is None:
        return None
    d, pbcl = result['details'], lm.pbcl
    sx, sy = case['spacing'][0], case['spacing'][1]
    gz, pz = d['ridge_point'][2], d['patella_point'][2]
    px, py = d['patella_point'][0], d['patella_point'][1]
    fb, fc, pat = case['bone'][:, :, gz], case['cartilage'][:, :, gz], case['patella'][:, :, pz]
    shape = fb.shape

    def _extent(mask):
        xs = np.where(mask.any(axis=1))[0]
        ys = np.where(mask.any(axis=0))[0]
        return xs.min(), xs.max(), ys.min(), ys.max()

    x0f, x1f, y0f, y1f = _extent(fb)
    x0p, x1p, y0p, y1p = _extent(pat)
    win_x, win_y = patella_window_mm / sx, patella_window_mm / sy
    x0p, x1p = max(x0p, int(px - win_x)), min(x1p, int(px + win_x))
    y0p, y1p = max(y0p, int(py - win_y)), min(y1p, int(py + win_y))
    box = (max(min(x0f, x0p) - int(margin_mm / sx), 0), min(max(x1f, x1p) + int(margin_mm / sx), shape[0]),
           max(min(y0f, y0p) - int(margin_mm / sy), 0), min(max(y1f, y1p) + int(margin_mm / sy), shape[1]))
    return {'box': box, 'spacing_xy': (sx, sy),
            'femur_bone': qc_draw.crop(fb, box), 'femur_cart': qc_draw.crop(fc, box),
            'patella': qc_draw.crop(pat, box),
            'ridge_point': d['ridge_point'], 'patella_point': d['patella_point'],
            'pbcl_medial': pbcl['medial_pt'], 'pbcl_lateral': pbcl['lateral_pt'], 'pbcl_theta': pbcl['theta'],
            **result['values']}


def render_panel(ax, p, fontsize=9):
    sx, sy = p['spacing_xy']
    box = p['box']
    X, Y = qc_draw.mesh_mm(box, p['spacing_xy'])
    qc_draw.plot_mask(ax, X, Y, p['femur_bone'], 'Greens', 0.6)
    qc_draw.plot_mask(ax, X, Y, p['femur_cart'], 'Blues', 0.6)
    qc_draw.plot_mask(ax, X, Y, p['patella'], 'Purples', 0.5)

    ridge_mm = np.array(p['ridge_point'][:2]) * (sx, sy)
    patella_mm = np.array(p['patella_point'][:2]) * (sx, sy)
    u, u_anterior = pcl_axes(p['pbcl_theta'])
    pbcl_mid = (np.array(p['pbcl_medial'], dtype=float) * (sx, sy)
                + np.array(p['pbcl_lateral'], dtype=float) * (sx, sy)) / 2.0

    x0, x1, y0, y1 = box
    span = max((x1 - x0) * sx, (y1 - y0) * sy)
    a, b = pbcl_mid - u * span, pbcl_mid + u * span
    ax.plot([a[0], b[0]], [a[1], b[1]], 'r-', linewidth=1.5, alpha=0.8)

    feet = []
    for pt in (ridge_mm, patella_mm):
        foot = pbcl_mid + np.dot(pt - pbcl_mid, u) * u
        feet.append(foot)
        ax.plot([pt[0], foot[0]], [pt[1], foot[1]], color='black', linestyle=':', linewidth=1.5)
    ax.scatter(*ridge_mm, color='red', edgecolor='black', s=90, zorder=6)
    ax.scatter(*patella_mm, color='darkviolet', edgecolor='black', s=90, zorder=6, marker='D')
    ax.annotate('', xy=tuple(feet[1]), xytext=tuple(feet[0]),
                arrowprops=dict(arrowstyle='<->', color='purple', linewidth=2, shrinkA=0, shrinkB=0))
    mid = (feet[0] + feet[1]) / 2.0 + 3.0 * u_anterior
    ax.text(mid[0], mid[1], f"PTLTR = {p['PTLTR_mm']:.1f}mm", color='purple', fontsize=fontsize,
            fontweight='bold', ha='center', va='center')
    qc_draw.finish_axial(ax, box, p['spacing_xy'])


def panel_title(p):
    return f"PTLTR={p['PTLTR_mm']:.1f}mm"


def landmark_points(case, lm, result):
    d = result.get('details')
    if not d:
        return []
    return [('lateral trochlear ridge', d['ridge_point']), ('patella lateral margin', d['patella_point'])]
