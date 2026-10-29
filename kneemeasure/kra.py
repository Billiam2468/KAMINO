"""
Knee rotation angle (KRA): the angle between the femoral posterior
condylar line (PBCL) and the tibial posterior condylar line (TPCL) on
axial slices -- "Knee Rotation Angle (MRI), +ve = tibial ER relative to
femur" on this project's manual measurement sheets.

Both line angles use the same smaller-x -> larger-x convention, so their
signed difference is wrapped into (-90, 90] (a LINE's angle only means
something mod 180).

Sign: positive = tibia externally rotated relative to the femur. The raw
(x, y) image frame has no knee-side awareness -- a right and left knee are
mirror images in it -- so the same physical rotation has opposite raw
signs on the two sides; side_sign corrects for that. KRA_SIGN is the
single global flip point if validation against manual KRA shows the sign
is backwards. NOT YET VALIDATED against manual values.
"""
import numpy as np

from . import qc_draw

KRA_SIGN = 1

COLUMNS = [
    ('KRA_deg', 'deg', 'Knee rotation angle: TPCL vs PBCL (+ = tibial external rotation)'),
]


def _signed_angle_diff_deg(a, b):
    """(a - b), wrapped into (-90, 90]."""
    return (a - b + 90.0) % 180.0 - 90.0


def compute(case, lm):
    pbcl = lm.pbcl
    if pbcl['theta'] is None:
        return {'values': {'KRA_deg': None}, 'note': 'PBCL not found'}
    tpcl = lm.tpcl
    if tpcl['theta'] is None:
        return {'values': {'KRA_deg': None}, 'note': f"tibial PCL not found -- {tpcl.get('note', '')}".rstrip(' -')}
    side_sign = 1 if lm.knee_side == 'L' else -1
    kra = KRA_SIGN * side_sign * _signed_angle_diff_deg(tpcl['theta'], pbcl['theta'])
    return {'values': {'KRA_deg': kra}, 'note': None,
            'aux': {'pbcl_theta': pbcl['theta'], 'tpcl_theta': tpcl['theta'],
                    'pbcl_z': pbcl['z'], 'tpcl_z': tpcl['z']}}


def qc_payload(case, lm, result, margin_mm=15.0):
    if result['values']['KRA_deg'] is None:
        return None
    pbcl, tpcl = lm.pbcl, lm.tpcl
    sxy = case['spacing'][:2]
    fb = case['bone'][:, :, pbcl['z']]
    fc = case['cartilage'][:, :, pbcl['z']]
    tb = case['tibia'][:, :, tpcl['z']]
    box = qc_draw.crop_box([fb, tb], sxy, margin_mm, fb.shape)
    return {'box': box, 'spacing_xy': sxy,
            'femur_bone': qc_draw.crop(fb, box), 'femur_cart': qc_draw.crop(fc, box),
            'tibia': qc_draw.crop(tb, box),
            'pbcl': (pbcl['medial_pt'], pbcl['lateral_pt'], pbcl['theta']),
            'tpcl': (tpcl['medial_pt'], tpcl['lateral_pt'], tpcl['theta']),
            'kra': result['values']['KRA_deg']}


def render_panel(ax, p, fontsize=8):
    """Femur bone/cartilage at the PBCL slice, tibia at the TPCL slice
    (see-through), and both lines through their own points, extended.
    Drawn in the reoriented mm frame -- the frame the angle is computed
    in."""
    sx, sy = p['spacing_xy']
    box = p['box']
    X, Y = qc_draw.mesh_mm(box, p['spacing_xy'])
    qc_draw.plot_mask(ax, X, Y, p['femur_bone'], 'Greens', 0.6)
    qc_draw.plot_mask(ax, X, Y, p['femur_cart'], 'Blues', 0.6)
    qc_draw.plot_mask(ax, X, Y, p['tibia'], 'Oranges', 0.4)

    x0, x1, y0, y1 = box
    span = max((x1 - x0) * sx, (y1 - y0) * sy)

    def _line(medial_pt, lateral_pt, color):
        p1 = np.array(medial_pt, dtype=float) * (sx, sy)
        p2 = np.array(lateral_pt, dtype=float) * (sx, sy)
        mid = (p1 + p2) / 2.0
        d = p2 - p1
        norm = np.linalg.norm(d)
        u = d / norm if norm > 0 else np.array([1.0, 0.0])
        a, b = mid - u * span, mid + u * span
        ax.plot([a[0], b[0]], [a[1], b[1]], color=color, linewidth=2)
        ax.scatter(*p1, color='yellow', edgecolor='black', s=40, zorder=6)
        ax.scatter(*p2, color=color, edgecolor='black', s=40, zorder=6)

    _line(p['pbcl'][0], p['pbcl'][1], 'red')
    _line(p['tpcl'][0], p['tpcl'][1], 'purple')
    qc_draw.finish_axial(ax, box, p['spacing_xy'])


def panel_title(p):
    return f"KRA={p['kra']:+.1f}deg (PBCL={p['pbcl'][2]:+.1f}, TPCL={p['tpcl'][2]:+.1f})"


def landmark_points(case, lm, result):
    pts = []
    if lm.pbcl['theta'] is not None:
        z = lm.pbcl['z']
        pts += [('PBCL medial', (*lm.pbcl['medial_pt'], z)), ('PBCL lateral', (*lm.pbcl['lateral_pt'], z))]
    if lm.tpcl['theta'] is not None:
        z = lm.tpcl['z']
        pts += [('TPCL medial', (*lm.tpcl['medial_pt'], z)), ('TPCL lateral', (*lm.tpcl['lateral_pt'], z))]
    return pts
