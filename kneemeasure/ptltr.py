"""
PTLTR: the mediolateral offset between the patella's lateral margin and
the lateral trochlear ridge, along the posterior condylar line, as in
Mistovich et al. 2018 (patellar tendon-lateral trochlear ridge distance)
except that the patella's edge stands in for the tendon's (the
segmentation has no patellar tendon). Unsigned.

Landmarks, following the paper's method:
- Tendon cut: the first axial slice below the patella bone -- where the
  paper marks the tendon ("the first distal cut demonstrating the patellar
  tendon without bone").
- Ridge point: the apex of the lateral trochlear ridge on the tendon cut
  ("if this cannot be identified on the same cut, this may be marked on the
  most visible cut, then translated"): the femur BONE point lateral of the
  slice's own cartilage sulcus lying furthest anterior of the posterior
  condylar line. A slice qualifies when it has a real (>= 2 mm) cartilage
  sulcus and the condyles haven't split yet; when the tendon cut doesn't,
  the nearest one that does is used (a tie goes to the distal slice). With
  the knee straight the patella often sits above the trochlea (16 of the 43
  ORS knees), and the ridge then comes from the top of the trochlea.
- Posterior condylar line: drawn on the ridge's slice, as in the paper's
  Figure 2B (pbcl.pbcl_on_slice); the knee's usual PBCL (pbcl.find_pbcl,
  the most posterior slice) when that slice has no qualifying notch.
- Patella point: the most lateral voxel of the patella bone (largest 3D
  component), tie-broken to the middle of a flat lateral facet.

    PTLTR = |(P_patella - P_ridge) . u_along_PCL|

On the ORS cohort this agrees with the manual PTLTR about as well as the
earlier ridge (the most anterior lateral point over the whole trochlea,
against the usual PBCL): ICC 0.73 vs 0.75 -- see VALIDATION.md.
"""
import numpy as np

from . import qc_draw
from .pbcl import pcl_axes, pbcl_on_slice
from .trochlea import most_anterior_component, find_peaks_and_sulcus
from .volume import count_components

COLUMNS = [
    ('PTLTR_mm', 'mm', 'patella lateral margin to lateral trochlear ridge at the tendon cut, along the posterior condylar line on that slice (unsigned)'),
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


def find_ridge_slice(case, min_depth_mm=2.0):
    """The slice the lateral trochlear ridge is marked on: the tendon cut,
    or the nearest slice with a usable ridge (see the module docstring).
    Returns {'z', 'sulcus_x' (that slice's sulcus), 'tendon_z', 'note'}."""
    pz = np.where(case['patella'].any(axis=(0, 1)))[0]
    if not len(pz):
        return {'z': None, 'note': 'no patella bone found in this case'}
    tendon_z = int(pz.min()) - 1
    usable = _usable_slices(case, min_depth_mm)
    if usable is None:
        return {'z': None, 'note': 'no cartilage found in this case'}
    if not usable:
        return {'z': None, 'note': 'no slice with a detectable cartilaginous sulcus found in scan range'}
    z = min(usable, key=lambda k: (abs(k - tendon_z), k))
    return {'z': z, 'sulcus_x': usable[z], 'tendon_z': tendon_z, 'note': None}


def _usable_slices(case, min_depth_mm):
    """{z: sulcus x} for the slices from the top of the femoral cartilage
    down to the condyle split that have a cartilage sulcus at least
    min_depth_mm deep (None when there is no cartilage). The depth gate
    leaves out the proximal shaft and the distal notch."""
    bone, cart = case['bone'], case['cartilage']
    z_nonempty = np.where(cart.any(axis=(0, 1)))[0]
    if len(z_nonempty) == 0:
        return None
    out = {}
    for z in range(int(z_nonempty.max()), int(z_nonempty.min()) - 1, -1):
        bone_sl = bone[:, :, z]
        _, comps = count_components(bone_sl)
        if len(comps) >= 2:
            break
        component = most_anterior_component(cart[:, :, z], bone_sl)
        if component.sum() == 0:
            continue
        groove = find_peaks_and_sulcus(component, bone_sl, case['spacing'][0], min_depth_mm=min_depth_mm)
        if groove is not None:
            out[z] = groove['sulcus'][0]
    return out


def ridge_apex(case, z, sulcus_x, pcl, knee_side):
    """Apex of the lateral trochlear ridge on slice z: the femur bone point
    lateral of the sulcus (lateral = smaller x for a RIGHT knee) furthest
    anterior of the posterior condylar line pcl. (Splitting at the
    midpoint of the line's two points instead let the search reach onto
    the medial facet when those points shifted by a few voxels.)"""
    spacing_xy = np.array(case['spacing'][:2])
    xs, ys = np.where(case['bone'][:, :, z])
    lateral = np.sign(xs - sulcus_x) == (-1 if knee_side == 'R' else 1)
    if not lateral.any():
        return None
    _, u_anterior = pcl_axes(pcl['theta'])
    offsets = (np.stack([xs[lateral], ys[lateral]], axis=1) * spacing_xy
               - np.array(pcl['medial_pt'], dtype=float) * spacing_xy) @ u_anterior
    i = int(np.argmax(offsets))
    return {'point': (int(xs[lateral][i]), int(ys[lateral][i]), int(z)), 'anterior_offset_mm': float(offsets[i])}


def compute(case, lm):
    empty = {'PTLTR_mm': None}
    if lm.pbcl['theta'] is None:
        return {'values': empty, 'note': 'PBCL not found'}
    sl = find_ridge_slice(case)
    if sl['z'] is None:
        return {'values': empty, 'note': f"lateral trochlear ridge not found -- {sl['note']}"}
    pcl = pbcl_on_slice(case, sl['z'], lm.knee_side)
    if pcl['theta'] is None:
        pcl = lm.pbcl
    apex = ridge_apex(case, sl['z'], sl['sulcus_x'], pcl, lm.knee_side)
    if apex is None:
        return {'values': empty, 'note': 'lateral trochlear ridge not found -- no bone lateral of the sulcus'}
    patella = find_lateral_patella_point(case['patella'], lm.knee_side)
    if patella['point'] is None:
        return {'values': empty, 'note': patella['note']}

    spacing_xy = np.array(case['spacing'][:2])
    v = np.array(patella['point'][:2]) * spacing_xy - np.array(apex['point'][:2]) * spacing_xy
    u_along, _ = pcl_axes(pcl['theta'])
    return {'values': {'PTLTR_mm': abs(float(np.dot(v, u_along)))}, 'note': None,
            'details': {'ridge_point': apex['point'], 'patella_point': patella['point'], 'tendon_z': sl['tendon_z'],
                        'pcl': {k: pcl[k] for k in ('theta', 'z', 'medial_pt', 'lateral_pt')}}}


def qc_payload(case, lm, result, margin_mm=8.0, patella_window_mm=25.0):
    """The patella sits a real ~2-3 cm anterior of the trochlea, so its
    displayed footprint is cropped to a window around the measured lateral
    point to keep the panel framed on what's measured."""
    if result['values']['PTLTR_mm'] is None:
        return None
    d = result['details']
    pbcl = d['pcl']            # the posterior condylar line PTLTR was measured along (on the ridge's slice)
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
