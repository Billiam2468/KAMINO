"""
Caton-Deschamps index (CDI), a patellar-height ratio read off ONE true
sagittal slice:

    CDI = (tibial anterosuperior corner -> inferior pole of the patellar
           articular cartilage) / (patellar articular cartilage length)

- Slice: the sagittal slice with the LONGEST patellar cartilage, the
  standard "best sagittal slice" criterion (avoids foreshortened oblique
  cuts through the periphery of the curved patella).
- Cartilage length: straight-line distance between its most superior and
  most inferior points on that slice, using both y and z (the strip is
  often tilted). Ties at an extreme z use the median y.
- Tibial anterosuperior corner: the tibia pixel maximizing (z*sz - y*sy),
  simultaneously the most superior and most anterior point -- the convex
  vertex where the near-vertical anterior cortex meets the near-horizontal
  plateau.

Reads the unfiltered label map (no largest-component filter), matching the
original implementation.
"""
import numpy as np

from .labels import FEM_B, FEM_C, PAT_B, PAT_C, TIB_B
from .qc_draw import plot_sagittal_label

COLUMNS = [
    ('CDI', 'ratio', 'Caton-Deschamps index'),
    ('CDI_numerator_mm', 'mm', 'tibial anterosuperior corner to inferior patellar cartilage pole'),
    ('patella_cartilage_length_mm', 'mm', 'patellar articular cartilage length (CDI denominator)'),
]


def find_best_sagittal_slice(labels, sy, sz):
    cart = labels == PAT_C
    lengths_mm = np.zeros(cart.shape[0])
    points = [None] * cart.shape[0]
    for x in range(cart.shape[0]):
        ys, zs = np.where(cart[x])
        if len(zs) == 0:
            continue
        z_lo, z_hi = zs.min(), zs.max()
        y_inf = float(np.median(ys[zs == z_lo]))
        y_sup = float(np.median(ys[zs == z_hi]))
        lengths_mm[x] = float(np.hypot((y_sup - y_inf) * sy, (z_hi - z_lo) * sz))
        points[x] = ((y_inf, float(z_lo)), (y_sup, float(z_hi)))
    best_x = int(np.argmax(lengths_mm))
    if points[best_x] is None:
        return None
    return {'x': best_x, 'length_mm': float(lengths_mm[best_x]),
            'inferior_point': points[best_x][0], 'superior_point': points[best_x][1]}


def find_tibial_anterosuperior_corner(tibia_slice_mask, sy, sz):
    ys, zs = np.where(tibia_slice_mask)
    if len(ys) == 0:
        return None
    idx = int(np.argmax(zs * sz - ys * sy))
    return (int(ys[idx]), int(zs[idx]))


def _dist_mm(p1, p2, sy, sz):
    return float(np.hypot((p2[0] - p1[0]) * sy, (p2[1] - p1[1]) * sz))


def compute(case, lm):
    empty = {'CDI': None, 'CDI_numerator_mm': None, 'patella_cartilage_length_mm': None}
    sy, sz = case['spacing'][1], case['spacing'][2]
    best = find_best_sagittal_slice(case['labels'], sy, sz)
    if best is None:
        return {'values': empty, 'note': 'no patella cartilage'}
    corner = find_tibial_anterosuperior_corner(case['labels'][best['x']] == TIB_B, sy, sz)
    values = dict(empty, patella_cartilage_length_mm=best['length_mm'])
    note = None
    if corner is None:
        note = 'no tibia on the best patellar-cartilage sagittal slice'
    else:
        num = _dist_mm(corner, best['inferior_point'], sy, sz)
        values['CDI_numerator_mm'] = num
        values['CDI'] = num / best['length_mm'] if num else None
    return {'values': values, 'note': note, 'details': {**best, 'corner': corner}}


def qc_payload(case, lm, result):
    d = result.get('details')
    if not d:
        return None
    return {'spacing_yz': (case['spacing'][1], case['spacing'][2]),
            'slice': case['labels'][d['x']].copy(), 'x': d['x'],
            'inferior_point': d['inferior_point'], 'superior_point': d['superior_point'],
            'corner': d['corner'], **result['values']}


def render_panel(ax, p, fontsize=7):
    """Sagittal slice at the longest patellar cartilage, anterior=left:
    cartilage length (black dashed), tibial corner (cyan star) and the CDI
    numerator (blue dash-dot)."""
    sy, sz = p['spacing_yz']
    sl = p['slice']
    for label, color in ((FEM_B, 'mediumseagreen'), (FEM_C, 'steelblue'), (TIB_B, 'darkorange'),
                         (PAT_B, 'darkviolet'), (PAT_C, 'red')):
        plot_sagittal_label(ax, sl == label, color, sy, sz)
    (yi, zi), (ys_, zs_) = p['inferior_point'], p['superior_point']
    a, b = (yi * sy, zi * sz), (ys_ * sy, zs_ * sz)
    ax.plot([a[0], b[0]], [a[1], b[1]], color='black', linewidth=1.5, linestyle='--', zorder=8)
    ax.scatter(*a, color='white', edgecolor='black', s=25, zorder=9)
    ax.scatter(*b, color='white', edgecolor='black', s=25, zorder=9)
    ax.text((a[0] + b[0]) / 2, (a[1] + b[1]) / 2, f"  {p['patella_cartilage_length_mm']:.1f}mm",
            fontsize=fontsize, fontweight='bold', va='center', zorder=9)
    if p['corner'] is not None:
        c = (p['corner'][0] * sy, p['corner'][1] * sz)
        ax.scatter(*c, color='cyan', edgecolor='black', s=50, marker='*', zorder=10)
        ax.plot([c[0], a[0]], [c[1], a[1]], color='dodgerblue', linewidth=1.5, linestyle='-.', zorder=7)
        ax.text((c[0] + a[0]) / 2, (c[1] + a[1]) / 2, f"  {p['CDI_numerator_mm']:.1f}mm",
                color='dodgerblue', fontsize=fontsize, fontweight='bold', va='center', zorder=9)
    ax.set_aspect('equal')
    ax.invert_xaxis()
    ax.axis('off')


def panel_title(p):
    cdi = f"CDI={p['CDI']:.2f}" if p['CDI'] is not None else 'CDI=n/a'
    return f"x={p['x']}, len={p['patella_cartilage_length_mm']:.1f}mm, {cdi}"


def landmark_points(case, lm, result):
    d = result.get('details')
    if not d:
        return []
    x = d['x']
    pts = [('patellar cartilage inferior pole', (x, *d['inferior_point'])),
           ('patellar cartilage superior pole', (x, *d['superior_point']))]
    if d['corner'] is not None:
        pts.append(('tibial anterosuperior corner', (x, *d['corner'])))
    return pts
