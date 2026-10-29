"""
TT-TG (tibial tuberosity to trochlear groove) distance, Schoettle-style:
project both landmark points onto the PBCL direction and measure the
offset between the projections. sTT-TG is the same pair of points along
the PBCL-PERPENDICULAR (anteroposterior) axis instead.

Landmarks:
- Groove: the sulcus of the Dejour cartilage reference slice (first slice,
  proximal -> distal, whose cartilage sulcus depth clears 2 mm).
  Schoettle's original method uses the bony groove; this pipeline uses the
  cartilage point by deliberate choice (the bone/cartilage, max-depth /
  first-crossing variants were roughly tied against manual TT-TG; the
  literature first-crossing criterion won).
- Tuberosity: tibial_tuberosity.find_tibial_tuberosity.

With PBCL unit vector u_along and anterior-facing perpendicular u_anterior
(reoriented-frame mm), and v = P_tubercle - P_groove:

    TT-TG  = |v . u_along|
    sTT-TG = STT_TG_SIGN * (v . u_anterior)    (negative = TT posterior to TG)
"""
import numpy as np

from . import qc_draw
from .pbcl import pcl_axes

STT_TG_SIGN = 1

COLUMNS = [
    ('TT_TG_mm', 'mm', 'TT-TG distance along the PBCL'),
    ('sTT_TG_mm', 'mm', 'sagittal TT-TG, perpendicular to the PBCL (- = TT posterior to TG)'),
]


def compute(case, lm):
    empty = {'TT_TG_mm': None, 'sTT_TG_mm': None}
    pbcl = lm.pbcl
    if pbcl['theta'] is None:
        return {'values': empty, 'note': 'PBCL not found'}
    groove = lm.reference_slice
    if groove['z'] is None:
        return {'values': empty, 'note': f"trochlear groove not found -- {groove['note']}"}
    tubercle = lm.tuberosity
    if tubercle['point'] is None:
        return {'values': empty, 'note': f"tibial tuberosity not found -- {tubercle['note']}"}

    spacing_xy = np.array(case['spacing'][:2])
    v = np.array(tubercle['point'][:2]) * spacing_xy - np.array(groove['point'][:2]) * spacing_xy
    u_along, u_anterior = pcl_axes(pbcl['theta'])
    return {'values': {'TT_TG_mm': abs(float(np.dot(v, u_along))),
                       'sTT_TG_mm': STT_TG_SIGN * float(np.dot(v, u_anterior))},
            'note': None,
            'aux': {'groove_z': groove['z'], 'tuberosity_z': tubercle['point'][2]}}


def qc_payload(case, lm, result, margin_mm=15.0):
    if result['values']['TT_TG_mm'] is None:
        return None
    g, t, pbcl = lm.reference_slice, lm.tuberosity, lm.pbcl
    sxy = case['spacing'][:2]
    gz, tz = g['point'][2], t['point'][2]
    fb, fc = case['bone'][:, :, gz], case['cartilage'][:, :, gz]
    tb = qc_draw.mask_slice(case['tibia'], tz)
    box = qc_draw.crop_box([fb, tb], sxy, margin_mm, fb.shape)
    return {'box': box, 'spacing_xy': sxy,
            'femur_bone': qc_draw.crop(fb, box), 'femur_cart': qc_draw.crop(fc, box),
            'tibia': qc_draw.crop(tb, box),
            'peak1': g['peak1'], 'peak2': g['peak2'],
            'groove_point': g['point'], 'tubercle_point': t['point'],
            'pbcl_medial': pbcl['medial_pt'], 'pbcl_lateral': pbcl['lateral_pt'], 'pbcl_theta': pbcl['theta'],
            **result['values']}


def render_panel(ax, p, fontsize=6):
    """Groove slice (bone + cartilage) with the tuberosity slice see-through
    on top, the facet chord, both landmarks, the PBCL axis at its own
    position, each point's projection onto it, the TT-TG bracket along the
    PBCL and the sTT-TG bracket perpendicular to it."""
    sx, sy = p['spacing_xy']
    box = p['box']
    X, Y = qc_draw.mesh_mm(box, p['spacing_xy'])
    qc_draw.plot_mask(ax, X, Y, p['femur_bone'], 'Greens', 0.6)
    qc_draw.plot_mask(ax, X, Y, p['femur_cart'], 'Blues', 0.6)
    qc_draw.plot_mask(ax, X, Y, p['tibia'], 'Oranges', 0.3)

    p1 = np.array(p['peak1']) * (sx, sy)
    p2 = np.array(p['peak2']) * (sx, sy)
    ax.plot([p1[0], p2[0]], [p1[1], p2[1]], 'c--', linewidth=2)

    groove_mm = np.array(p['groove_point'][:2]) * (sx, sy)
    tubercle_mm = np.array(p['tubercle_point'][:2]) * (sx, sy)
    u, u_anterior = pcl_axes(p['pbcl_theta'])
    pbcl_mid = (np.array(p['pbcl_medial'], dtype=float) * (sx, sy)
                + np.array(p['pbcl_lateral'], dtype=float) * (sx, sy)) / 2.0

    x0, x1, y0, y1 = box
    span = max((x1 - x0) * sx, (y1 - y0) * sy)
    a, b = pbcl_mid - u * span, pbcl_mid + u * span
    ax.plot([a[0], b[0]], [a[1], b[1]], 'r-', linewidth=1.5, alpha=0.8)

    for pt in (groove_mm, tubercle_mm):
        foot = pbcl_mid + np.dot(pt - pbcl_mid, u) * u
        ax.plot([pt[0], foot[0]], [pt[1], foot[1]], color='black', linestyle=':', linewidth=1.5)
    ax.scatter(*groove_mm, color='red', edgecolor='black', s=90, zorder=6)
    ax.scatter(*tubercle_mm, color='red', edgecolor='black', s=90, zorder=6, marker='D')

    end = groove_mm + np.dot(tubercle_mm - groove_mm, u) * u
    ax.annotate('', xy=tuple(end), xytext=tuple(groove_mm),
                arrowprops=dict(arrowstyle='<->', color='purple', linewidth=2, shrinkA=0, shrinkB=0))
    mid = (groove_mm + end) / 2.0 + 3.0 * u_anterior
    ax.text(mid[0], mid[1], f"TT-TG = {p['TT_TG_mm']:.1f}mm", color='purple', fontsize=fontsize,
            fontweight='bold', ha='center', va='center')

    end_ap = tubercle_mm - np.dot(tubercle_mm - groove_mm, u_anterior) * u_anterior
    ax.annotate('', xy=tuple(end_ap), xytext=tuple(tubercle_mm),
                arrowprops=dict(arrowstyle='<->', color='darkorange', linewidth=2, shrinkA=0, shrinkB=0))
    mid_ap = (tubercle_mm + end_ap) / 2.0
    ax.text(mid_ap[0], mid_ap[1], f"  sTT-TG = {p['sTT_TG_mm']:+.1f}mm", color='darkorange',
            fontsize=fontsize, fontweight='bold', va='center')
    qc_draw.finish_axial(ax, box, p['spacing_xy'])


def panel_title(p):
    return f"TT-TG={p['TT_TG_mm']:.1f}mm, sTT-TG={p['sTT_TG_mm']:+.1f}mm"


def landmark_points(case, lm, result):
    pts = []
    g = lm.reference_slice
    if g['z'] is not None:
        pts.append(('trochlear groove / reference slice sulcus', g['point']))
    t = lm.tuberosity
    if t['point'] is not None:
        pts.append(('tibial tuberosity', t['point']))
    return pts
