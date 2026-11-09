"""
Trochlear width (TW) and axial patellar-trochlear overlap (APTO): how far
the patella's medial border sits medial of the lateral edge of the trochlear
cartilage, along the femoral posterior condylar line (PBCL) -- the "jump
distance" the patella must cross to dislocate (Camp et al., AJSM 2016),
first described as the axial trochlear length / engagement (Guilbert et al.,
OTSR 2013). A position measure, not a size measure. 2D, in-plane, like the
other axial measurements.

Both share the patellar point: the medial end of the PW_effective line
(Landmarks.patella_effective: the medial edge of the patellar cartilage,
label 4, on the PW_effective slice). The two points are usually on
different slices; as in Camp's transfer between images, both are projected
onto the PBCL direction (translated, not rotated), so the slice gap doesn't
enter the arithmetic.

The trochlear point is the lateral edge of the ANTERIOR trochlear
cartilage: on each slice, the most anterior femoral-cartilage piece that
sits on bone (trochlea.most_anterior_component, the piece the sulcus is
found on). These slices also carry cartilage on the back of the condyles,
which reaches further lateral (+5.6 mm median TW on 122 knees without this
rule).

- TW_mm (Camp 2016, as written; Hevesi 2023): the most lateral such point
  across the trochlear range -- the top of the femoral cartilage down to the
  last slice before the femur splits into two condyles (the notch roof).
- APTO_mm (our lab's manual "Axial Patellar-Trochlear Overlap"): the lateral
  edge on the TT-TG groove slice (Landmarks.reference_slice, the reviewer's
  slice when the sulcus was reviewed). Against manual APTO on 84 knees:
  bias +1.9 mm, ICC(2,1) 0.82; TW_mm: +4.1 mm, 0.77 -- its lateral point
  sits near the notch roof, where the edge mostly reflects how wide the
  distal femur is.

Edge to edge, as PW_effective: each point is its voxel's outer corner along
the PBCL direction. Side-aware: + = the patella's medial border is medial
of the trochlear edge (the normal case).
"""
import numpy as np

from . import qc_draw
from .labels import PAT_C
from .patella_width import effective_width
from .pbcl import pcl_axes
from .trochlea import most_anterior_component
from .volume import count_components

COLUMNS = [
    ('TW_mm', 'mm', 'trochlear width (Camp 2016): PW_effective medial edge to the most lateral anterior trochlear '
                    'cartilage across the trochlear range, along the PBCL'),
    ('APTO_mm', 'mm', 'axial patellar-trochlear overlap: the same, to the lateral edge of the anterior trochlear '
                      'cartilage on the TT-TG groove slice'),
]
TW_COLOR, APTO_COLOR, PATELLA_COLOR = 'darkviolet', 'deepskyblue', 'gold'


def _lateral(u, knee_side):
    """Unit vector toward lateral along the PBCL (u points toward larger x: lateral for a LEFT knee)."""
    return u * (1.0 if knee_side == 'L' else -1.0)


def trochlear_range(case):
    """(top, bottom) slices of the trochlea: the top of the femoral cartilage
    down to the last slice before the femur splits into two condyles; None
    without femoral cartilage."""
    cart, bone = case['cartilage'], case['bone']
    zs = np.where(cart.any(axis=(0, 1)))[0]
    if not len(zs):
        return None
    top = int(zs.max())
    for z in range(top, int(zs.min()) - 1, -1):
        if len(count_components(bone[:, :, z])[1]) >= 2:
            return (top, z + 1) if z < top else None
    return top, int(zs.min())


def _trochlear_edge(case, z, u, knee_side):
    """The anterior trochlear cartilage's lateral edge on slice z, in index
    units (outer corner), and whether that voxel is on the scan's border; or
    None without anterior cartilage there."""
    piece = most_anterior_component(case['cartilage'][:, :, z], case['bone'][:, :, z])
    if not piece.any():
        return None
    sxy = np.asarray(case['spacing'][:2], dtype=float)
    _, _, lat = effective_width(piece, sxy, u, knee_side)
    xs, ys = np.nonzero(piece)
    i = int(np.argmax((np.stack([xs, ys], axis=1) * sxy) @ _lateral(u, knee_side)))   # the lateral voxel itself
    nx, ny = piece.shape
    return lat, bool(xs[i] <= 0 or ys[i] <= 0 or xs[i] >= nx - 1 or ys[i] >= ny - 1)


def compute(case, lm):
    empty = {c: None for c, _, _ in COLUMNS}
    pbcl = lm.pbcl
    if pbcl['theta'] is None:
        return {'values': empty, 'note': 'PBCL not found'}
    eff = lm.patella_effective if (case['labels'] == PAT_C).any() else None
    if eff is None:
        return {'values': empty, 'note': 'no patellar cartilage (label 4): no PW_effective medial point'}
    sxy = np.asarray(case['spacing'][:2], dtype=float)
    sz = case['spacing'][2]
    u, _ = pcl_axes(pbcl['theta'])
    lat_dir = _lateral(u, lm.knee_side)
    medial = np.asarray(eff['medial_edge'], dtype=float)
    med_mm = float((medial * sxy) @ lat_dir)
    values, notes = dict(empty), []
    details = {'medial': (*eff['medial_edge'], eff['z']), 'pbcl': {k: pbcl[k] for k in ('theta', 'medial_pt', 'lateral_pt')}}

    def point(name, z, found):
        lat, on_border = found
        values[name] = float((lat * sxy) @ lat_dir) - med_mm
        details[name] = {'z': int(z), 'lateral': (float(lat[0]), float(lat[1]), int(z)),
                         'offset_mm': (int(z) - eff['z']) * sz}           # from the PW_effective slice, + = proximal
        if on_border:
            notes.append(f'{name}: the lateral trochlear point is on the edge of the field of view (z {z})')
        if values[name] <= 0:
            notes.append(f"{name} {values[name]:.1f} mm: the patella's medial border is lateral of the trochlear edge")

    rng = trochlear_range(case)
    edges = {}
    if rng is not None:
        for z in range(rng[0], rng[1] - 1, -1):
            found = _trochlear_edge(case, z, u, lm.knee_side)
            if found is not None:
                edges[z] = found
    if not edges:
        notes.append('no femoral cartilage in the trochlear region')
    else:
        z_tw = max(edges, key=lambda z: float((edges[z][0] * sxy) @ lat_dir))
        point('TW_mm', z_tw, edges[z_tw])

    ref = lm.reference_slice
    if ref['z'] is None:
        notes.append(f"APTO_mm: TT-TG groove slice not found -- {ref['note']}")
    else:
        z = int(ref['z'])
        found = edges.get(z) or _trochlear_edge(case, z, u, lm.knee_side)
        if found is None:
            notes.append(f'APTO_mm: no anterior trochlear cartilage on the TT-TG groove slice (z {z})')
        else:
            point('APTO_mm', z, found)
    return {'values': values, 'note': '; '.join(notes) or None, 'details': details}


def qc_payload(case, lm, result, margin_mm=10.0):
    d = result.get('details') or {}
    z = d['APTO_mm']['z'] if 'APTO_mm' in d else (d['TW_mm']['z'] if 'TW_mm' in d else None)
    if z is None:
        return None
    sxy = case['spacing'][:2]
    bone = case['bone'][:, :, z]
    piece = most_anterior_component(case['cartilage'][:, :, z], bone)
    patella = case['labels'][:, :, d['medial'][2]] == PAT_C
    box = qc_draw.crop_box([bone, patella], sxy, margin_mm, bone.shape)
    return {'box': box, 'spacing_xy': sxy, 'z': z, 'patella_z': d['medial'][2],
            'femur_bone': qc_draw.crop(bone, box), 'trochlear_cart': qc_draw.crop(piece, box),
            'other_cart': qc_draw.crop(case['cartilage'][:, :, z] & ~piece, box),
            'patella_cart': qc_draw.crop(patella, box),
            'medial': d['medial'], 'tw': d.get('TW_mm'), 'apto': d.get('APTO_mm'), 'pbcl': d['pbcl'],
            **result['values']}


def render_panel(ax, p, fontsize=8):
    """The APTO slice (the TT-TG groove slice): femur bone (green), the
    anterior trochlear cartilage (blue) and other femoral cartilage (grey),
    with the PW_effective slice's patellar cartilage see-through on top
    (red), as in Camp's transfer between images. Dashed lines perpendicular
    to the PBCL through the patella's medial edge (gold), the APTO lateral
    point (sky blue) and the TW lateral point (violet; hollow when it's on
    another slice); the two widths as bars along the PBCL in front of the
    femur."""
    sx, sy = p['spacing_xy']
    s = np.array([sx, sy])
    box = p['box']
    X, Y = qc_draw.mesh_mm(box, p['spacing_xy'])
    qc_draw.plot_mask(ax, X, Y, p['femur_bone'], 'Greens', 0.5)
    qc_draw.plot_mask(ax, X, Y, p['other_cart'], 'Greys', 0.4)
    qc_draw.plot_mask(ax, X, Y, p['trochlear_cart'], 'Blues', 0.7)
    qc_draw.plot_mask(ax, X, Y, p['patella_cart'], 'Reds', 0.35)
    x0, x1, y0, y1 = box
    span = max((x1 - x0) * sx, (y1 - y0) * sy)
    u, u_ant = pcl_axes(p['pbcl']['theta'])

    med = np.array(p['medial'][:2]) * s
    pts = [(med, PATELLA_COLOR, True)]
    for key, color in (('apto', APTO_COLOR), ('tw', TW_COLOR)):
        if p.get(key):
            pts.append((np.array(p[key]['lateral'][:2]) * s, color, p[key]['z'] == p['z']))
    front = max(np.dot(q, u_ant) for q, _, _ in pts)                    # the most anterior point
    for q, color, here in pts:
        r0, r1 = q - u_ant * span * 0.3, q + u_ant * span * 0.3
        ax.plot([r0[0], r1[0]], [r0[1], r1[1]], color=color, linewidth=1, linestyle='--')
        ax.scatter(*q, s=45, zorder=6, linewidths=1.5, facecolors=color if here else 'none',
                   edgecolors='black' if here else color)
    for i, (key, color) in enumerate((('apto', APTO_COLOR), ('tw', TW_COLOR))):
        if p.get(key):
            level = front + 3.0 + 3.0 * i
            ends = [np.dot(q, u) * u + level * u_ant for q in (med, np.array(p[key]['lateral'][:2]) * s)]
            ax.plot([ends[0][0], ends[1][0]], [ends[0][1], ends[1][1]], color=color, linewidth=2.5)
    qc_draw.finish_axial(ax, box, p['spacing_xy'])


def panel_title(p):
    f = lambda v: 'n/a' if v is None else f'{v:.1f}mm'
    tz = f" z {p['tw']['z']}" if p.get('tw') else ''
    return f"APTO={f(p['APTO_mm'])} (z {p['z']}), TW={f(p['TW_mm'])}{tz}; patella z {p['patella_z']}"


def landmark_points(case, lm, result):
    d = result.get('details') or {}
    pts = [('TW / APTO patella medial edge (PW_effective)', d['medial'])] if 'medial' in d else []
    if 'TW_mm' in d:
        pts.append(('TW lateral trochlear edge', d['TW_mm']['lateral']))
    if 'APTO_mm' in d:
        pts.append(('APTO lateral trochlear edge', d['APTO_mm']['lateral']))
    return pts
