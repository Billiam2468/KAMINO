"""
Femoral width (Chen et al., OJSM 2024): the distance between the lines
perpendicular to the femoral posterior condylar line (PBCL) through the most
prominent medial and lateral points of the femoral condyles -- the femur's
extent along the PBCL, measured on the BONE (label 1, largest 3D component),
edge to edge (each extreme voxel's outer corner, as for tibial width and the
patellar widths). 2D, in-plane on one axial slice.

Slice: the PBCL's own slice (pbcl.find_pbcl: of the slices whose posterior
outline has two condylar lobes and a notch, the one whose bitangent lies
most posterior), as tibial width uses the tibial line's slice. Chen takes the
most proximal slice showing the posterior condyles with the Roman arch; that
is the top of the notch, where the femur narrows fast into the metaphysis,
so its width there hinges on where the arch is said to end (a 4 vs 8 mm
notch-depth rule moved it a median 4 mm, over 17 mm on 1 knee in 10). The
PBCL slice is inside the same arch zone, at its distal end, with both
condyles and a well-formed notch; on 85 knees it read a median 2.5 mm
narrower than the widest slice (the epicondyles, ~10 mm higher).

Medial/lateral: x increases toward patient LEFT, so lateral = larger x for
a LEFT knee, smaller x for a RIGHT knee.
"""
import numpy as np

from . import qc_draw
from .patella_width import effective_width
from .pbcl import pcl_axes

COLUMNS = [
    ('FW_mm', 'mm', 'femoral width (Chen 2024): femur bone extent along the PBCL, on the PBCL slice'),
]


def compute(case, lm):
    pbcl = lm.pbcl
    if pbcl['theta'] is None:
        return {'values': {'FW_mm': None}, 'note': 'PBCL not found'}
    z = pbcl['z']
    bone = case['bone'][:, :, z]
    u, _ = pcl_axes(pbcl['theta'])
    width, medial, lateral = effective_width(bone, case['spacing'][:2], u, lm.knee_side)   # extent along u, edge to edge
    note = None
    if bone[[0, -1], :].any() or bone[:, [0, -1]].any():
        note = f'the femur reaches the edge of the scan on the PBCL slice (z {z}): the width may be cut short'
    return {'values': {'FW_mm': width}, 'note': note,
            'details': {'z': z, 'medial_edge': tuple(medial), 'lateral_edge': tuple(lateral),
                        'pbcl': {k: pbcl[k] for k in ('theta', 'medial_pt', 'lateral_pt')}}}


def qc_payload(case, lm, result, margin_mm=10.0):
    d = result.get('details')
    if not d:
        return None
    z, sxy = d['z'], case['spacing'][:2]
    bone, cart = case['bone'][:, :, z], case['cartilage'][:, :, z]
    box = qc_draw.crop_box([bone], sxy, margin_mm, bone.shape)
    return {'box': box, 'spacing_xy': sxy, 'z': z,
            'femur_bone': qc_draw.crop(bone, box), 'femur_cart': qc_draw.crop(cart, box),
            'medial_edge': d['medial_edge'], 'lateral_edge': d['lateral_edge'], 'pbcl': d['pbcl'],
            **result['values']}


def render_panel(ax, p, fontsize=8):
    """Femur bone (green) and cartilage (blue) on the PBCL slice, the PBCL
    through its condyle points (red), and the width drawn on the PBCL
    (gold): the projection of the femur's medial and lateral edges onto it,
    along the dashed perpendiculars through them."""
    sx, sy = p['spacing_xy']
    s = np.array([sx, sy])
    box = p['box']
    X, Y = qc_draw.mesh_mm(box, p['spacing_xy'])
    qc_draw.plot_mask(ax, X, Y, p['femur_bone'], 'Greens', 0.6)
    qc_draw.plot_mask(ax, X, Y, p['femur_cart'], 'Blues', 0.6)
    x0, x1, y0, y1 = box
    span = max((x1 - x0) * sx, (y1 - y0) * sy)

    u, u_ant = pcl_axes(p['pbcl']['theta'])
    m1, m2 = np.array(p['pbcl']['medial_pt'], dtype=float) * s, np.array(p['pbcl']['lateral_pt'], dtype=float) * s
    mid = (m1 + m2) / 2.0
    a, b = mid - u * span, mid + u * span
    ax.plot([a[0], b[0]], [a[1], b[1]], color='red', linewidth=1.5)
    ax.scatter(*m1, color='yellow', edgecolor='black', s=30, zorder=6)
    ax.scatter(*m2, color='red', edgecolor='black', s=30, zorder=6)

    med, lat = np.array(p['medial_edge']) * s, np.array(p['lateral_edge']) * s
    feet = [mid + np.dot(q - mid, u) * u for q in (med, lat)]          # the edges projected onto the PBCL
    for q, f in zip((med, lat), feet):
        top = f + (np.dot(q - f, u_ant) + span * 0.5) * u_ant           # up through the edge, well past the femur
        ax.plot([f[0], top[0]], [f[1], top[1]], color='gold', linewidth=1, linestyle='--', zorder=5)
    ax.plot([feet[0][0], feet[1][0]], [feet[0][1], feet[1][1]], color='gold', linewidth=3, zorder=5)
    qc_draw.finish_axial(ax, box, p['spacing_xy'])


def panel_title(p):
    return f"FW={p['FW_mm']:.1f}mm (z {p['z']})"


def landmark_points(case, lm, result):
    d = result.get('details')
    if not d:
        return []
    return [('femoral width medial edge', (*d['medial_edge'], d['z'])),
            ('femoral width lateral edge', (*d['lateral_edge'], d['z']))]
