"""
Tibial tubercle lateralization (TTL), per Tensho et al. (J Bone Joint Surg
Am 2015;97(17):1441-8, and 2018;100(9):e58, Fig. 2):

    T = the tibia's mediolateral width measured PARALLEL to the tibial
        posterior condylar line (TPCL), at the TPCL's own slice
    t = distance along that same TPCL direction from T's medial edge to the
        tibial tuberosity point (at the tuberosity's own slice --
        translate-not-rotate, as for TT-TG/KRA)

    TTL = t / T

Medial vs lateral for T's two edges is resolved by the TPCL's own
(knee-side-aware) medial_pt/lateral_pt ordering along the same axis.

Tensho report TTL as a percentage (controls 64.6 +- 2.8% on CT); here it is
a fraction. Their tuberosity point is the tubercle tip where the patellar
tendon fully contacts it; the segmentations have no tendon, so the
tuberosity detector's point is used instead.
"""
import numpy as np

from . import qc_draw

COLUMNS = [
    ('TTL', 'ratio', 'tibial tubercle lateralization t/T'),
    ('TTL_t_mm', 'mm', 'TTL t: medial tibial edge to tuberosity along the TPCL'),
    ('TTL_T_mm', 'mm', 'TTL T: tibial width along the TPCL'),
]


def _direction(theta_deg):
    theta = np.radians(theta_deg)
    return np.array([np.cos(theta), np.sin(theta)])


def compute(case, lm):
    empty = {'TTL': None, 'TTL_t_mm': None, 'TTL_T_mm': None}
    tpcl = lm.tpcl
    if tpcl['theta'] is None:
        return {'values': empty, 'note': f"tibial PCL not found -- {tpcl.get('note', '')}".rstrip(' -')}
    tubercle = lm.tuberosity
    if tubercle['point'] is None:
        return {'values': empty, 'note': f"tibial tuberosity not found -- {tubercle['note']}"}

    spacing_xy = np.array(case['spacing'][:2])
    u = _direction(tpcl['theta'])
    xs, ys = np.where(case['tibia'][:, :, tpcl['z']])
    if len(xs) == 0:
        return {'values': empty, 'note': 'empty tibia bone slice at TPCL z'}
    projections = np.stack([xs * spacing_xy[0], ys * spacing_xy[1]], axis=1) @ u

    med_proj_ref = float(np.dot(np.array(tpcl['medial_pt'], dtype=float) * spacing_xy, u))
    lat_proj_ref = float(np.dot(np.array(tpcl['lateral_pt'], dtype=float) * spacing_xy, u))
    medial_is_min = med_proj_ref <= lat_proj_ref
    idx_min, idx_max = int(np.argmin(projections)), int(np.argmax(projections))
    medial_idx = idx_min if medial_is_min else idx_max
    lateral_idx = idx_max if medial_is_min else idx_min

    medial_edge_proj = float(projections[medial_idx])
    width_mm = abs(float(projections[lateral_idx]) - medial_edge_proj)
    tubercle_proj = float(np.dot(np.array(tubercle['point'][:2], dtype=float) * spacing_xy, u))
    t_mm = abs(tubercle_proj - medial_edge_proj)
    ttl = t_mm / width_mm if width_mm > 0 else None

    return {'values': {'TTL': ttl, 'TTL_t_mm': t_mm, 'TTL_T_mm': width_mm}, 'note': None,
            'details': {'medial_edge_point': (int(xs[medial_idx]), int(ys[medial_idx]), tpcl['z']),
                        'lateral_edge_point': (int(xs[lateral_idx]), int(ys[lateral_idx]), tpcl['z'])}}


def qc_payload(case, lm, result, margin_mm=15.0):
    if result['values']['TTL'] is None:
        return None
    tpcl, t = lm.tpcl, lm.tuberosity
    sxy = case['spacing'][:2]
    wide = case['tibia'][:, :, tpcl['z']]
    narrow = qc_draw.mask_slice(case['tibia'], t['point'][2])
    box = qc_draw.crop_box([wide, narrow], sxy, margin_mm, wide.shape)
    return {'box': box, 'spacing_xy': sxy,
            'tibia_tpcl': qc_draw.crop(wide, box), 'tibia_tt': qc_draw.crop(narrow, box),
            'tpcl_medial': tpcl['medial_pt'], 'tpcl_lateral': tpcl['lateral_pt'], 'tpcl_theta': tpcl['theta'],
            'medial_edge': result['details']['medial_edge_point'][:2],
            'lateral_edge': result['details']['lateral_edge_point'][:2],
            'tubercle_point': t['point'], **result['values']}


def render_panel(ax, p, fontsize=9):
    """Tibia at the TPCL slice (green) with the tuberosity slice
    see-through (orange), the TPCL (red), the T width drawn parallel to the
    TPCL (gold) and the t bracket from the medial edge (purple)."""
    sx, sy = p['spacing_xy']
    box = p['box']
    X, Y = qc_draw.mesh_mm(box, p['spacing_xy'])
    qc_draw.plot_mask(ax, X, Y, p['tibia_tpcl'], 'Greens', 0.6)
    qc_draw.plot_mask(ax, X, Y, p['tibia_tt'], 'Oranges', 0.4)

    x0, x1, y0, y1 = box
    span = max((x1 - x0) * sx, (y1 - y0) * sy)
    u = _direction(p['tpcl_theta'])
    mid = (np.array(p['tpcl_medial'], dtype=float) * (sx, sy)
           + np.array(p['tpcl_lateral'], dtype=float) * (sx, sy)) / 2.0
    a, b = mid - u * span, mid + u * span
    ax.plot([a[0], b[0]], [a[1], b[1]], color='red', linewidth=2)

    u_perp = np.array([-u[1], u[0]])
    med_mm = np.array(p['medial_edge'], dtype=float) * (sx, sy)
    lat_mm = np.array(p['lateral_edge'], dtype=float) * (sx, sy)
    avg_perp = (np.dot(med_mm - mid, u_perp) + np.dot(lat_mm - mid, u_perp)) / 2.0
    med_draw = mid + np.dot(med_mm - mid, u) * u + avg_perp * u_perp
    lat_draw = mid + np.dot(lat_mm - mid, u) * u + avg_perp * u_perp
    ax.plot([med_draw[0], lat_draw[0]], [med_draw[1], lat_draw[1]], color='gold', linewidth=2)
    ax.scatter(*med_draw, color='gold', edgecolor='black', s=60, zorder=6, marker='s')
    ax.scatter(*lat_draw, color='gold', edgecolor='black', s=60, zorder=6, marker='s')

    tt_mm = np.array(p['tubercle_point'][:2], dtype=float) * (sx, sy)
    tt_foot = mid + np.dot(tt_mm - mid, u) * u
    med_foot = mid + np.dot(med_mm - mid, u) * u
    ax.plot([tt_mm[0], tt_foot[0]], [tt_mm[1], tt_foot[1]], color='black', linestyle=':', linewidth=1.5)
    ax.scatter(*tt_mm, color='red', edgecolor='black', s=70, zorder=7, marker='D')
    ax.annotate('', xy=tuple(tt_foot), xytext=tuple(med_foot),
                arrowprops=dict(arrowstyle='<->', color='purple', linewidth=2))
    mb = (med_foot + tt_foot) / 2.0
    ax.text(mb[0], mb[1], f"  t = {p['TTL_t_mm']:.1f}mm", color='purple', fontsize=fontsize,
            fontweight='bold', va='center')
    qc_draw.finish_axial(ax, box, p['spacing_xy'])


def panel_title(p):
    return f"TTL={p['TTL']:.2f} (t={p['TTL_t_mm']:.1f}, T={p['TTL_T_mm']:.1f})"


def landmark_points(case, lm, result):
    d = result.get('details')
    if not d:
        return []
    return [('TTL medial tibial edge', d['medial_edge_point']), ('TTL lateral tibial edge', d['lateral_edge_point'])]
