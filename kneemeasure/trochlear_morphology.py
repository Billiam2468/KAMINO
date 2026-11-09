"""
Trochlear morphology at the Dejour cartilage reference slice (first slice,
cranial -> caudal, with a >= 2 mm cartilage sulcus), per Dejour,
Giovannetti de Sanctis et al., KSSTA 2025;33(8):2833-2846:

- LTI (lateral trochlear inclination): angle between the PBCL (translated,
  not rotated, to the reference slice -- only its angle matters) and the
  sulcus -> lateral-facet-peak line. Positive if the lateral facet peak is
  ANTERIOR to the PBCL.
- Sulcus angle: angle at the sulcus between rays to the two facet peaks
  (PBCL-independent).
- Medial/lateral facet depth: each facet peak's offset from the sulcus
  along the PBCL's anterior-facing perpendicular (positive = facet anterior
  to the sulcus, the normal case).
- Trochlear depth: perpendicular-to-PBCL distance from the sulcus to the
  medial->lateral facet chord, interpolated at the sulcus's own along-PBCL
  position (positive = sulcus deeper than the chord).
- Facet lengths: straight-line sulcus -> facet-peak distances; asymmetry =
  medial / lateral length.

Lateral facet: for a LEFT knee the larger-x peak, for a RIGHT knee the
smaller-x peak.

Bone versions (*_bone_* columns): the same eight values from the bony
sulcus and facet peaks (Landmarks.reference_slice_bone) -- the same
reference slice, PBCL and hull algorithm, on the femur bone surface, or a
reviewer's bony points (which may be on another slice).
"""
import numpy as np

from .pbcl import pcl_axes
from .volume import true_physical_lp

BONE_COLOR = 'magenta'                     # the bony points and lines in the QC grids

COLUMNS = [
    ('LTI_deg', 'deg', 'lateral trochlear inclination (+ = lateral facet anterior to PBCL)'),
    ('LTI_bone_deg', 'deg', 'bony LTI: the same, on the femur bone surface'),
    ('sulcus_angle_deg', 'deg', 'cartilage sulcus angle'),
    ('sulcus_angle_bone_deg', 'deg', 'bony sulcus angle'),
    ('medial_facet_depth_mm', 'mm', 'medial facet height above the sulcus, perpendicular to PBCL'),
    ('medial_facet_depth_bone_mm', 'mm', 'bony medial facet height above the bony sulcus, perpendicular to PBCL'),
    ('lateral_facet_depth_mm', 'mm', 'lateral facet height above the sulcus, perpendicular to PBCL'),
    ('lateral_facet_depth_bone_mm', 'mm', 'bony lateral facet height above the bony sulcus, perpendicular to PBCL'),
    ('trochlear_depth_mm', 'mm', 'sulcus depth below the facet chord, perpendicular to PBCL'),
    ('trochlear_depth_bone_mm', 'mm', 'bony sulcus depth below the bony facet chord, perpendicular to PBCL'),
    ('medial_facet_length_mm', 'mm', 'sulcus to medial facet peak'),
    ('medial_facet_length_bone_mm', 'mm', 'bony sulcus to bony medial facet peak'),
    ('lateral_facet_length_mm', 'mm', 'sulcus to lateral facet peak'),
    ('lateral_facet_length_bone_mm', 'mm', 'bony sulcus to bony lateral facet peak'),
    ('facet_asymmetry_ratio', 'ratio', 'medial / lateral facet length'),
    ('facet_asymmetry_bone_ratio', 'ratio', 'bony medial / lateral facet length'),
]


def _to_mm(pt, spacing):
    return np.array([pt[0] * spacing[0], pt[1] * spacing[1]], dtype=float)


def compute_lti(sulcus, lateral_peak, pbcl_theta_deg, spacing):
    u, v = pcl_axes(pbcl_theta_deg)
    vec = _to_mm(lateral_peak, spacing) - _to_mm(sulcus, spacing)
    return float(np.degrees(np.arctan2(float(np.dot(vec, v)), abs(float(np.dot(vec, u))))))


def compute_sulcus_angle_deg(sulcus, medial_pt, lateral_pt, spacing):
    s = _to_mm(sulcus, spacing)
    v1 = _to_mm(medial_pt, spacing) - s
    v2 = _to_mm(lateral_pt, spacing) - s
    cos_angle = np.clip(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)), -1.0, 1.0)
    return float(np.degrees(np.arccos(cos_angle)))


def compute_facet_depths_mm(sulcus, medial_pt, lateral_pt, pbcl_theta_deg, spacing):
    _, v = pcl_axes(pbcl_theta_deg)
    s = _to_mm(sulcus, spacing)
    return (float(np.dot(_to_mm(medial_pt, spacing) - s, v)),
            float(np.dot(_to_mm(lateral_pt, spacing) - s, v)))


def compute_trochlear_depth_mm(sulcus, medial_pt, lateral_pt, pbcl_theta_deg, spacing):
    u, v = pcl_axes(pbcl_theta_deg)
    s = _to_mm(sulcus, spacing)
    m = _to_mm(medial_pt, spacing) - s
    l = _to_mm(lateral_pt, spacing) - s
    m_u, m_v = float(np.dot(m, u)), float(np.dot(m, v))
    l_u, l_v = float(np.dot(l, u)), float(np.dot(l, v))
    if l_u == m_u:
        return (m_v + l_v) / 2.0
    t = (0.0 - m_u) / (l_u - m_u)
    return m_v + t * (l_v - m_v)


def _facet_lengths_mm(sulcus, medial_pt, lateral_pt, spacing):
    def dist(p1, p2):
        d = _to_mm(p1, spacing) - _to_mm(p2, spacing)
        return float(np.hypot(d[0], d[1]))
    return dist(sulcus, medial_pt), dist(sulcus, lateral_pt)


def bone_column(col):
    """The bone version's column name: '_bone' before the unit (LTI_deg -> LTI_bone_deg)."""
    name, unit = col.rsplit('_', 1)
    return f'{name}_bone_{unit}'


def _surface_values(ref, knee_side, theta, spacing):
    """The eight values from one surface's sulcus and facet peaks (cartilage
    column names), and the medial and lateral peaks."""
    sorted_pts = sorted([ref['peak1'], ref['peak2']], key=lambda p: p[0])
    lateral_pt, medial_pt = (sorted_pts[1], sorted_pts[0]) if knee_side == 'L' else (sorted_pts[0], sorted_pts[1])
    s = ref['sulcus']
    med_depth, lat_depth = compute_facet_depths_mm(s, medial_pt, lateral_pt, theta, spacing)
    med_len, lat_len = _facet_lengths_mm(s, medial_pt, lateral_pt, spacing)
    values = {
        'LTI_deg': compute_lti(s, lateral_pt, theta, spacing),
        'sulcus_angle_deg': compute_sulcus_angle_deg(s, medial_pt, lateral_pt, spacing),
        'medial_facet_depth_mm': med_depth,
        'lateral_facet_depth_mm': lat_depth,
        'trochlear_depth_mm': compute_trochlear_depth_mm(s, medial_pt, lateral_pt, theta, spacing),
        'medial_facet_length_mm': med_len,
        'lateral_facet_length_mm': lat_len,
        'facet_asymmetry_ratio': med_len / lat_len if lat_len != 0 else None,
    }
    return values, medial_pt, lateral_pt


def compute(case, lm):
    empty = {c: None for c, _, _ in COLUMNS}
    ref = lm.reference_slice
    if ref['z'] is None:
        return {'values': empty, 'note': f"reference slice: {ref['note']}"}
    pbcl = lm.pbcl
    if pbcl['theta'] is None:
        return {'values': empty, 'note': 'no valid PBCL slice found'}

    spacing = (case['spacing'][0], case['spacing'][1])
    theta = pbcl['theta']
    values, medial_pt, lateral_pt = _surface_values(ref, lm.knee_side, theta, spacing)
    details = {'z': ref['z'], 'medial_pt': medial_pt, 'lateral_pt': lateral_pt, 'sulcus': ref['sulcus']}
    note = None
    bone = lm.reference_slice_bone
    if bone['z'] is None:
        values.update({bone_column(k): None for k in list(values)})
        note = f"bone version: {bone['note']}"
    else:
        b, b_med, b_lat = _surface_values(bone, lm.knee_side, theta, spacing)
        values.update({bone_column(k): v for k, v in b.items()})
        details['bone'] = {'medial_pt': b_med, 'lateral_pt': b_lat, 'sulcus': bone['sulcus'], 'z': bone['z']}
    return {'values': {c: values[c] for c, _, _ in COLUMNS}, 'note': note, 'details': details}


def qc_payload(case, lm, result, margin=15):
    """Rendered in TRUE physical (L, P) mm so the image matches the real
    scan (the measured angles are relative, so the frame doesn't affect
    them)."""
    d = result.get('details')
    if not d or result['values']['LTI_deg'] is None:
        return None
    ref, header, z = lm.reference_slice, case['header'], d['z']
    bone_sl, comp = ref['bone_slice'], ref['cartilage_component']
    xp, yp = np.where(bone_sl.any(axis=1))[0], np.where(bone_sl.any(axis=0))[0]
    x0, x1 = max(xp.min() - margin, 0), xp.max() + margin
    y0, y1 = max(yp.min() - margin, 0), yp.max() + margin
    X, Y = np.meshgrid(np.arange(x0, x1), np.arange(y0, y1), indexing='ij')
    L, P = true_physical_lp(header, X, Y, z)

    def to_true(pt, at=z):
        tl, tp = true_physical_lp(header, pt[0], pt[1], at)
        return float(tl), float(tp)

    b = d.get('bone')
    return {'L': L, 'P': P, 'bone': bone_sl[x0:x1, y0:y1].copy(), 'cart': comp[x0:x1, y0:y1].copy(),
            'medial': to_true(d['medial_pt']), 'lateral': to_true(d['lateral_pt']),
            'sulcus': to_true(d['sulcus']),
            'bone_points': {k: to_true(b[k], b['z']) for k in ('medial_pt', 'lateral_pt', 'sulcus')} if b else None,
            'bone_z': b['z'] if b and b['z'] != z else None,      # a reviewed bony sulcus on another slice
            **result['values']}


def render_panel(ax, p, fontsize=8):
    """Femur bone (green) and trochlear cartilage (purple) on the reference
    slice. Cartilage points: medial peak yellow, lateral peak orange, sulcus
    cyan diamond, facet chord dashed yellow, LTI line orange. Bone points and
    lines the same shapes in magenta."""
    from matplotlib import colormaps
    from matplotlib.colors import ListedColormap
    for mask, cmap, alpha in ((p['bone'], 'Greens', 0.6), (p['cart'], 'Purples', 0.8)):
        # solid, pre-mixed with white: the mesh is in true physical coordinates (slightly rotated), and
        # see-through cells there leave visible seams where they meet. Bone and cartilage don't overlap.
        c = np.array(colormaps[cmap](1 / 1.5))
        c[:3] = alpha * c[:3] + (1 - alpha)
        m = mask.astype(float)
        ax.pcolormesh(p['L'], p['P'], np.ma.masked_where(m == 0, m), cmap=ListedColormap([c]), shading='nearest')
    ax.set_aspect('equal')
    ax.invert_yaxis()
    med, lat, s = p['medial'], p['lateral'], p['sulcus']
    ax.plot([med[0], lat[0]], [med[1], lat[1]], color='yellow', linewidth=1.5, linestyle='--')
    ax.plot([s[0], lat[0]], [s[1], lat[1]], color='orange', linewidth=2)
    ax.scatter(*lat, color='orange', edgecolor='black', s=60, zorder=5)
    ax.scatter(*med, color='yellow', edgecolor='black', s=60, zorder=5)
    ax.scatter(*s, color='cyan', edgecolor='black', s=60, zorder=5, marker='D')
    b = p.get('bone_points')
    if b:                                            # the bony points and lines, magenta
        bm, bl, bs = b['medial_pt'], b['lateral_pt'], b['sulcus']
        ax.plot([bm[0], bl[0]], [bm[1], bl[1]], color=BONE_COLOR, linewidth=1.5, linestyle='--')
        ax.plot([bs[0], bl[0]], [bs[1], bl[1]], color=BONE_COLOR, linewidth=2)
        ax.scatter(*zip(bm, bl), color=BONE_COLOR, edgecolor='black', s=45, zorder=5)
        ax.scatter(*bs, color=BONE_COLOR, edgecolor='black', s=45, zorder=5, marker='D')
    ax.axis('off')


def panel_title(p):
    asym = f"{p['facet_asymmetry_ratio']:.2f}" if p['facet_asymmetry_ratio'] is not None else 'n/a'
    at = f" (z {p['bone_z']})" if p.get('bone_z') is not None else ''
    bone = (f"bone{at}: LTI={p['LTI_bone_deg']:.1f}deg, sulcus={p['sulcus_angle_bone_deg']:.1f}deg, "
            f"depth={p['trochlear_depth_bone_mm']:.1f}mm" if p.get('LTI_bone_deg') is not None else 'bone: n/a')
    return (f"LTI={p['LTI_deg']:.1f}deg, sulcus={p['sulcus_angle_deg']:.1f}deg\n"
            f"med={p['medial_facet_depth_mm']:.1f}mm, lat={p['lateral_facet_depth_mm']:.1f}mm, asym={asym}\n{bone}")


def landmark_points(case, lm, result):
    d = result.get('details')
    if not d:
        return []
    z = d['z']
    pts = [('reference slice medial facet peak', (*d['medial_pt'], z)),
           ('reference slice lateral facet peak', (*d['lateral_pt'], z))]
    b = d.get('bone')
    if b:
        pts += [('reference slice bony medial facet peak', (*b['medial_pt'], b['z'])),
                ('reference slice bony lateral facet peak', (*b['lateral_pt'], b['z'])),
                ('reference slice bony sulcus', (*b['sulcus'], b['z']))]
    return pts
