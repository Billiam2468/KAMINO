"""
Patellar width, two published ways, and lateral patellar tilt -- all 2D, in
plane on the segmentation's own axial slices (which match the PD axial the
raters read), like every other axial measurement here.

- PW_effective (Hevesi et al., AJSM 2023, "effective patellar width"): the
  widest distance between the medial and lateral borders of the patellar
  articular CARTILAGE (label 4), measured parallel to the femoral posterior
  condylar line (PBCL, translated -- not rotated -- to each slice). Every
  slice with cartilage is tried; the widest wins.
- PW_max (Chen et al., OJSM 2024; Staubli): the maximum width of the BONY
  patella (label 3) along the patella's own axis -- the line through its
  medial and lateral tips, found from the longest chord of the slice's
  cross-section (the patella tapers to points medially and laterally, so
  the longest chord runs tip to tip; see max_chord). Every slice is tried;
  the widest wins.
- Patellar tilt (Chen 2024): on the PW_max slice, the angle between the
  PW_max line and the PBCL. Positive = lateral tilt: the patella's lateral
  border is more posterior than the PBCL would put it (the lines converge
  laterally). Side-aware like KRA, so left and right knees read alike.

The Hevesi width is a projection onto the femoral axis, so a tilted patella
reads narrower; PW_max isn't, and the tilt says by how much.

Both widths run edge to edge, as measured by hand: from the outer corner of
the outermost voxels, not their centres (as for TTL's tibial width). Along a
direction at angle a to the grid that adds |cos a| + |sin a| voxels (1 to
1.41 for square voxels), not exactly one.

Each structure is its largest 3D component (stray voxels can't widen it);
the notes say when a second bony piece is big enough to be a bipartite
fragment, when a label splits on its chosen slice, when the patella reaches
the edge of the scan, and when PW_effective > PW_max (it shouldn't: the
cartilage is narrower than the bone and a projection only shortens).

Medial/lateral: x increases toward patient LEFT, so lateral = larger x for
a LEFT knee, smaller x for a RIGHT knee.
"""
import numpy as np
from scipy import ndimage
from scipy.spatial import ConvexHull

from . import qc_draw
from .labels import PAT_B, PAT_C
from .pbcl import pcl_axes
from .volume import largest_component_3d

COLUMNS = [
    ('PW_effective_mm', 'mm', 'effective patellar width (Hevesi 2023): widest medial-lateral extent of the '
                              'patellar cartilage, parallel to the PBCL'),
    ('PW_max_mm', 'mm', 'maximum patellar width (Chen 2024 / Staubli): bony patella, tip to tip along its '
                        'own axis'),
    ('patellar_tilt_deg', 'deg', 'lateral patellar tilt (Chen 2024): PW_max line vs PBCL '
                                 '(+ = lateral border posterior)'),
]
# written to measurements.csv, but not measurements: kept out of the statistics, outlier flags and plots
INFO_COLUMNS = [
    ('PW_effective_slice', 'slice', 'axial slice (z) PW_effective was measured on'),
    ('PW_max_slice', 'slice', 'axial slice (z) PW_max and the tilt were measured on'),
    ('PW_effective_offset_mm', 'mm', "PW_effective slice's distance from the PBCL slice (+ = proximal)"),
    ('PW_max_offset_mm', 'mm', "PW_max slice's distance from the PBCL slice (+ = proximal)"),
]

BIPARTITE_FRACTION = 0.05     # a second bony piece this big (of the largest) may be a bipartite fragment
MAX_CHORD_ANGLE_DEG = 45.0    # a tilt beyond this is worth a look in the QC (a severely tilted patella, or a bad axis)
TIP_MM = 2.0                  # the patella's medial/lateral tips: voxels within this of each end of its axis


def _lateral_sign(knee_side):
    return 1.0 if knee_side == 'L' else -1.0


def effective_width(mask2d, spacing_xy, u, knee_side):
    """Hevesi width of one slice: the extent of mask2d along the unit vector
    u (in mm), edge to edge. Returns (width_mm, medial_edge, lateral_edge),
    the edges in index units (each extreme voxel's outer corner along u)."""
    xs, ys = np.nonzero(mask2d)
    s = np.asarray(spacing_xy, dtype=float)
    proj = (np.stack([xs, ys], axis=1) * s) @ u
    lo, hi = int(np.argmin(proj)), int(np.argmax(proj))
    out = 0.5 * np.sign(u)                                 # corner offset toward +u, in voxels
    lo_edge = np.array([xs[lo], ys[lo]], dtype=float) - out
    hi_edge = np.array([xs[hi], ys[hi]], dtype=float) + out
    width = float((hi_edge * s) @ u - (lo_edge * s) @ u)
    # u points toward larger x (the PBCL's theta runs smaller-x -> larger-x), so +u is lateral for a LEFT knee
    lateral_is_hi = (u[0] >= 0) == (knee_side == 'L')
    return (width, lo_edge if lateral_is_hi else hi_edge, hi_edge if lateral_is_hi else lo_edge)


def max_chord(mask2d, spacing_xy, knee_side, tip_mm=TIP_MM, rounds=2):
    """PW_max of one slice along the patella's own axis: (length_mm,
    medial_end, lateral_end), ends in index units.

    The axis starts as the longest chord between voxel corners. Near the
    long axis a chord's length hardly changes with its angle, so on a voxel
    grid that chord can sit a voxel off the true axis and tip it by several
    degrees; the axis is therefore redrawn through the centres of the two
    tips (the voxels within tip_mm of each end along it), as a reader clicks
    the medial and lateral borders. The width is the extent along that
    axis, edge to edge (each end voxel's outer corner)."""
    s = np.asarray(spacing_xy, dtype=float)
    xs, ys = np.nonzero(mask2d)
    centres = np.stack([xs, ys], axis=1).astype(float)
    corners = np.concatenate([centres + (dx, dy) for dx in (-0.5, 0.5) for dy in (-0.5, 0.5)])
    try:
        hull = corners[ConvexHull(corners * s).vertices]
    except Exception:                                      # degenerate: a single row or column of voxels
        hull = corners
    d = np.linalg.norm((hull[:, None, :] - hull[None, :, :]) * s, axis=2)
    i, j = np.unravel_index(int(np.argmax(d)), d.shape)
    e = (hull[j] - hull[i]) * s
    e /= np.linalg.norm(e)
    for _ in range(rounds):
        p = (centres * s) @ e
        tip_a = centres[p <= p.min() + tip_mm].mean(axis=0) * s
        tip_b = centres[p >= p.max() - tip_mm].mean(axis=0) * s
        if np.allclose(tip_a, tip_b):
            break
        e = (tip_b - tip_a) / np.linalg.norm(tip_b - tip_a)
    half = 0.5 * float(np.abs(e) @ s)                      # an end voxel's outer corner along e
    p = (centres * s) @ e
    lo, hi = p.min() - half, p.max() + half
    mid = (tip_a + tip_b) / 2.0
    a, b = (mid + (lo - mid @ e) * e) / s, (mid + (hi - mid @ e) * e) / s
    lateral_is_b = (b[0] > a[0]) == (knee_side == 'L')
    return float(hi - lo), (a if lateral_is_b else b), (b if lateral_is_b else a)


def tilt_deg(medial_end, lateral_end, spacing_xy, pbcl_theta, knee_side):
    """Signed angle (deg) from the PBCL to the medial->lateral PW_max line;
    + = the lateral end is posterior of the PBCL's direction (lateral tilt)."""
    u, _ = pcl_axes(pbcl_theta)
    side = _lateral_sign(knee_side)
    toward_lateral = u * side                              # u points toward larger x: lateral for a LEFT knee
    d = (np.asarray(lateral_end) - np.asarray(medial_end)) * np.asarray(spacing_xy, dtype=float)
    cross = toward_lateral[0] * d[1] - toward_lateral[1] * d[0]
    return side * float(np.degrees(np.arctan2(cross, float(toward_lateral @ d))))


def _n_pieces(mask2d):
    return ndimage.label(mask2d, structure=np.ones((3, 3)))[1]


def compute(case, lm):
    empty = {c: None for c, _, _ in COLUMNS + INFO_COLUMNS}
    labels, (sx, sy, sz) = case['labels'], case['spacing']
    sxy = (sx, sy)
    bone_raw, cart_raw = labels == PAT_B, labels == PAT_C
    notes = []
    if not bone_raw.any():
        notes.append('no patella bone (label 3)')
    if not cart_raw.any():
        notes.append('no patellar cartilage (label 4)')
    if not bone_raw.any() and not cart_raw.any():
        return {'values': empty, 'note': '; '.join(notes)}

    pieces, n = ndimage.label(bone_raw, structure=np.ones((3, 3, 3)))
    if n > 1:
        sizes = np.sort(np.bincount(pieces.ravel())[1:])[::-1]
        if sizes[1] > BIPARTITE_FRACTION * sizes[0]:
            notes.append(f'second patella bone piece is {100 * sizes[1] / sizes[0]:.0f}% of the largest '
                         f'(bipartite patella?); only the largest is measured')
    either = bone_raw | cart_raw
    if (either[[0, -1], :, :].any() or either[:, [0, -1], :].any() or either[:, :, [0, -1]].any()):
        notes.append('the patella reaches the edge of the scan')

    pbcl = lm.pbcl
    values, details = dict(empty), {'pbcl': None}
    if pbcl['theta'] is None:
        notes.append('PBCL not found: no PW_effective or tilt')
    else:
        details['pbcl'] = {k: pbcl[k] for k in ('theta', 'z', 'medial_pt', 'lateral_pt')}

    # PW_effective: widest cartilage extent along the PBCL, over all slices
    if cart_raw.any() and pbcl['theta'] is not None:
        cart = largest_component_3d(cart_raw)
        u, _ = pcl_axes(pbcl['theta'])
        best = None
        for z in np.where(cart.any(axis=(0, 1)))[0]:
            w, med, lat = effective_width(cart[:, :, z], sxy, u, lm.knee_side)
            if best is None or w > best[0]:
                best = (w, med, lat, int(z))
        w, med, lat, z = best
        values.update(PW_effective_mm=w, PW_effective_slice=z, PW_effective_offset_mm=(z - pbcl['z']) * sz)
        details['effective'] = {'z': z, 'medial_edge': tuple(med), 'lateral_edge': tuple(lat)}
        if _n_pieces(cart_raw[:, :, z]) > 1:
            notes.append(f'patellar cartilage is in more than one piece on the PW_effective slice (z {z})')

    # PW_max: longest bony chord over all slices, and the tilt on that slice
    if bone_raw.any():
        bone = case['patella'] if case.get('patella') is not None else largest_component_3d(bone_raw)
        best = None
        for z in np.where(bone.any(axis=(0, 1)))[0]:
            length, med, lat = max_chord(bone[:, :, z], sxy, lm.knee_side)
            if best is None or length > best[0]:
                best = (length, med, lat, int(z))
        length, med, lat, z = best
        values.update(PW_max_mm=length, PW_max_slice=z)
        details['max'] = {'z': z, 'medial_end': tuple(med), 'lateral_end': tuple(lat)}
        if _n_pieces(bone_raw[:, :, z]) > 1:
            notes.append(f'patella bone is in more than one piece on the PW_max slice (z {z})')
        if pbcl['theta'] is not None:
            tilt = tilt_deg(med, lat, sxy, pbcl['theta'], lm.knee_side)
            values.update(patellar_tilt_deg=tilt, PW_max_offset_mm=(z - pbcl['z']) * sz)
            if abs(tilt) > MAX_CHORD_ANGLE_DEG:
                notes.append(f'patellar tilt {tilt:+.0f} deg: check in the QC that the PW_max line runs from the '
                             f'medial to the lateral border')

    if values['PW_effective_mm'] is not None and values['PW_max_mm'] is not None \
            and values['PW_effective_mm'] > values['PW_max_mm']:
        notes.append(f"PW_effective ({values['PW_effective_mm']:.1f} mm) is larger than PW_max "
                     f"({values['PW_max_mm']:.1f} mm)")
    return {'values': values, 'note': '; '.join(notes) or None, 'details': details}


def qc_payload(case, lm, result, margin_mm=10.0):
    d = result.get('details') or {}
    if 'max' not in d:
        return None
    labels = case['labels']
    zm = d['max']['z']
    ze = d['effective']['z'] if 'effective' in d else None
    bone, cart = labels[:, :, zm] == PAT_B, labels[:, :, zm] == PAT_C
    cart_eff = labels[:, :, ze] == PAT_C if ze is not None and ze != zm else None
    sxy = case['spacing'][:2]
    box = qc_draw.crop_box([bone, cart, cart_eff], sxy, margin_mm, bone.shape)
    return {'box': box, 'spacing_xy': sxy, 'bone': qc_draw.crop(bone, box), 'cart': qc_draw.crop(cart, box),
            'cart_eff': qc_draw.crop(cart_eff, box) if cart_eff is not None else None,
            'max': d['max'], 'effective': d.get('effective'), 'pbcl': d.get('pbcl'), **result['values']}


def render_panel(ax, p, fontsize=8):
    """The PW_max slice: bone (purple), cartilage (red), the PW_effective
    slice's cartilage see-through (orange) when it's another slice; PW_max
    (blue), PW_effective parallel to the PBCL with perpendicular reference
    lines at its ends (gold), and the PBCL's direction (red), translated to
    just behind the patella."""
    sx, sy = p['spacing_xy']
    s = np.array([sx, sy])
    box = p['box']
    X, Y = qc_draw.mesh_mm(box, p['spacing_xy'])
    qc_draw.plot_mask(ax, X, Y, p['bone'], 'Purples', 0.5)
    qc_draw.plot_mask(ax, X, Y, p['cart'], 'Reds', 0.6)
    if p.get('cart_eff') is not None:
        qc_draw.plot_mask(ax, X, Y, p['cart_eff'], 'Oranges', 0.35)
    x0, x1, y0, y1 = box
    span = max((x1 - x0) * sx, (y1 - y0) * sy)

    if p.get('pbcl'):
        u, u_ant = pcl_axes(p['pbcl']['theta'])
        mid = np.array([(x0 + x1) / 2 * sx, (y1 - 2) * sy])            # just behind the patella, in the box
        a, b = mid - u * span, mid + u * span
        ax.plot([a[0], b[0]], [a[1], b[1]], color='red', linewidth=1.5)
        if p.get('effective'):
            med, lat = np.array(p['effective']['medial_edge']) * s, np.array(p['effective']['lateral_edge']) * s
            level = (np.dot(med, u_ant) + np.dot(lat, u_ant)) / 2.0
            ends = [np.dot(q, u) * u + level * u_ant for q in (med, lat)]
            ax.plot([ends[0][0], ends[1][0]], [ends[0][1], ends[1][1]], color='gold', linewidth=2.5)
            for e in ends:
                r0, r1 = e - u_ant * span * 0.25, e + u_ant * span * 0.25
                ax.plot([r0[0], r1[0]], [r0[1], r1[1]], color='gold', linewidth=1, linestyle='--')

    med, lat = np.array(p['max']['medial_end']) * s, np.array(p['max']['lateral_end']) * s
    ax.plot([med[0], lat[0]], [med[1], lat[1]], color='dodgerblue', linewidth=2)
    ax.scatter(*med, color='dodgerblue', edgecolor='black', s=40, zorder=6, marker='s')
    ax.scatter(*lat, color='dodgerblue', edgecolor='black', s=40, zorder=6, marker='D')    # lateral end
    qc_draw.finish_axial(ax, box, p['spacing_xy'])


def panel_title(p):
    f = lambda v, fmt: 'n/a' if v is None else format(v, fmt)
    return (f"PW_eff={f(p['PW_effective_mm'], '.1f')} (z {f(p['PW_effective_slice'], 'd')}), "
            f"PW_max={f(p['PW_max_mm'], '.1f')} (z {p['PW_max_slice']}), tilt={f(p['patellar_tilt_deg'], '+.1f')}")


def landmark_points(case, lm, result):
    d = result.get('details') or {}
    pts = []
    if 'effective' in d:
        z = d['effective']['z']
        pts += [('PW_effective medial edge', (*d['effective']['medial_edge'], z)),
                ('PW_effective lateral edge', (*d['effective']['lateral_edge'], z))]
    if 'max' in d:
        z = d['max']['z']
        pts += [('PW_max medial end', (*d['max']['medial_end'], z)),
                ('PW_max lateral end', (*d['max']['lateral_end'], z))]
    return pts
