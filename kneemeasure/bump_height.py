"""
Trochlear "bump" height (supratrochlear spur), Dejour-style: how far the
trochlear surface protrudes anterior to the straight line of the femoral
shaft's anterior cortex, on one TRUE SAGITTAL slice. Reported for the bone
surface and the cartilage surface against the same line.

1. Sagittal slice: the x-index (mediolateral position) of the DEEPEST
   cartilage sulcus across the trochlear range (max relative depth -- not
   the first-crossing Dejour slice TT-TG uses, since "the deepest part of
   the groove" is a maximum over the range).
2. Anterior cortex line: the femur bone's anterior-most contour point at
   every z on that slice, fit with an exhaustive RANSAC (every pair of
   points tried, scored in true mm). A candidate must span >= min_span_mm
   along z, reach within anchor_top_mm of the slice's most proximal point
   (so it is the shaft, not a straight-looking run further down) and have
   >= min_inliers points; the winner is refit to its inliers by SVD. No
   qualifying line -> the case fails rather than returning a number from a
   poorly supported line.
3. Bump height: the MAXIMUM signed perpendicular distance of any bone
   (resp. cartilage) voxel on the slice from that line, positive =
   anterior. A bump is a protrusion, so each value is floored at 0.
"""
import numpy as np

from .qc_draw import plot_sagittal_label

COLUMNS = [
    ('bone_bump_mm', 'mm', 'bone trochlear bump height above the anterior cortex line'),
    ('cartilage_bump_mm', 'mm', 'cartilage trochlear bump height above the anterior cortex line'),
]


def _anterior_contour_by_z(bone_sl):
    """(zs, ys): the anterior-most (smallest y) bone point at every z."""
    zs_present = np.where(bone_sl.any(axis=0))[0]
    ys = np.array([np.where(bone_sl[:, z])[0].min() for z in zs_present], dtype=float)
    return zs_present.astype(float), ys


def _ransac_line(zs, ys, spacing_z, spacing_y, residual_mm=0.4, min_span_mm=20.0, anchor_top_mm=8.0,
                  min_inliers=15):
    """Exhaustive-pair RANSAC fit of the anterior femoral cortex line; see
    module docstring. Returns {'direction', 'point', 'inlier_mask'} (unit
    direction + anchor point in (z, y) TRUE mm) or None."""
    pts_mm = np.stack([zs * spacing_z, ys * spacing_y], axis=1)
    n = len(pts_mm)
    if n < 2:
        return None
    z_top_mm = pts_mm[:, 0].max()

    best_inliers, best_count = None, -1
    for i in range(n):
        d = pts_mm[i + 1:] - pts_mm[i]
        norms = np.hypot(d[:, 0], d[:, 1])
        valid = norms > 1e-6
        for j_idx in np.where(valid)[0]:
            d_unit = d[j_idx] / norms[j_idx]
            normal = np.array([-d_unit[1], d_unit[0]])
            residuals = np.abs((pts_mm - pts_mm[i]) @ normal)
            inliers = residuals <= residual_mm
            count = int(inliers.sum())
            if count <= best_count:
                continue
            inlier_z_mm = pts_mm[inliers, 0]
            if (inlier_z_mm.max() - inlier_z_mm.min()) < min_span_mm:
                continue
            if (z_top_mm - inlier_z_mm.max()) > anchor_top_mm:
                continue
            best_count, best_inliers = count, inliers

    if best_inliers is None or best_count < max(2, min_inliers):
        return None

    inlier_pts = pts_mm[best_inliers]
    mean = inlier_pts.mean(axis=0)
    _, _, vt = np.linalg.svd(inlier_pts - mean)
    return {'direction': vt[0], 'point': mean, 'inlier_mask': best_inliers}


def find_anterior_cortex_line(case, x, residual_mm=0.4, min_span_mm=20.0, anchor_top_mm=8.0, min_inliers=15):
    sy, sz = case['spacing'][1], case['spacing'][2]
    zs, ys = _anterior_contour_by_z(case['bone'][x])
    fit = _ransac_line(zs, ys, sz, sy, residual_mm=residual_mm, min_span_mm=min_span_mm,
                        anchor_top_mm=anchor_top_mm, min_inliers=min_inliers)
    if fit is None:
        return {'direction': None,
                'note': f'no shaft-consistent line found spanning at least {min_span_mm:.0f}mm '
                        f'and reaching within {anchor_top_mm:.0f}mm of the proximal edge'}
    return {'direction': fit['direction'], 'point': fit['point'], 'inlier_mask': fit['inlier_mask'],
            'zs': zs, 'ys': ys, 'note': None}


def _signed_perpendicular_distances_mm(direction, point, idx_y, idx_z, spacing_y, spacing_z, sign_ref_idx):
    """Signed perpendicular distance (true mm) from the cortex line to each
    (y, z) index point, positive = ANTERIOR. The SVD direction's sign is
    arbitrary, so the result is flipped if needed to make sign_ref_idx (a
    point known to be posterior of the line) come out negative."""
    dz, dy = direction
    pz, py = point
    y_mm, z_mm = idx_y * spacing_y, idx_z * spacing_z
    raw = dz * (y_mm - py) - dy * (z_mm - pz)
    ref_y_mm, ref_z_mm = sign_ref_idx[0] * spacing_y, sign_ref_idx[1] * spacing_z
    ref_raw = dz * (ref_y_mm - py) - dy * (ref_z_mm - pz)
    return raw if ref_raw < 0 else -raw


def compute_bump_height(case, x, cortex):
    sy, sz = case['spacing'][1], case['spacing'][2]
    bone_ys, bone_zs = np.where(case['bone'][x])
    if len(bone_ys) == 0:
        return {'bone_bump_mm': None, 'cartilage_bump_mm': None, 'note': 'empty bone slice'}
    sign_ref_idx = (bone_ys[np.argmax(bone_ys)], bone_zs[np.argmax(bone_ys)])

    bone_dist = _signed_perpendicular_distances_mm(cortex['direction'], cortex['point'],
                                                    bone_ys, bone_zs, sy, sz, sign_ref_idx)
    bone_idx = int(np.argmax(bone_dist))
    result = {'bone_bump_mm': max(0.0, float(bone_dist[bone_idx])),
              'bone_point': (int(bone_ys[bone_idx]), int(bone_zs[bone_idx])),
              'cartilage_bump_mm': None, 'cartilage_point': None, 'note': None}

    cart_ys, cart_zs = np.where(case['cartilage'][x])
    if len(cart_ys) == 0:
        result['note'] = 'no cartilage found in this slice'
        return result
    cart_dist = _signed_perpendicular_distances_mm(cortex['direction'], cortex['point'],
                                                    cart_ys, cart_zs, sy, sz, sign_ref_idx)
    cart_idx = int(np.argmax(cart_dist))
    result['cartilage_bump_mm'] = max(0.0, float(cart_dist[cart_idx]))
    result['cartilage_point'] = (int(cart_ys[cart_idx]), int(cart_zs[cart_idx]))
    return result


def compute(case, lm, residual_mm=0.4, min_span_mm=20.0):
    empty = {'bone_bump_mm': None, 'cartilage_bump_mm': None}
    groove = lm.deepest_groove
    if groove['z'] is None:
        return {'values': empty, 'note': groove['note']}
    x = groove['point'][0]
    cortex = find_anterior_cortex_line(case, x, residual_mm=residual_mm, min_span_mm=min_span_mm)
    if cortex['direction'] is None:
        return {'values': empty, 'note': cortex['note']}
    bump = compute_bump_height(case, x, cortex)
    return {'values': {'bone_bump_mm': bump['bone_bump_mm'], 'cartilage_bump_mm': bump['cartilage_bump_mm']},
            'note': bump['note'],
            'details': {'x': x, 'cortex': cortex, 'bone_point': bump.get('bone_point'),
                        'cartilage_point': bump.get('cartilage_point')}}


def qc_payload(case, lm, result):
    d = result.get('details')
    if not d or result['values']['bone_bump_mm'] is None:
        return None
    x = d['x']
    return {'spacing_yz': (case['spacing'][1], case['spacing'][2]),
            'bone': case['bone'][x].copy(), 'cartilage': case['cartilage'][x].copy(),
            'cortex': {k: d['cortex'][k] for k in ('direction', 'point', 'inlier_mask', 'zs', 'ys')},
            'bone_point': d['bone_point'], 'cartilage_point': d['cartilage_point'],
            **result['values']}


def render_panel(ax, p, fontsize=8):
    """Sagittal slice (femur bone + cartilage), anterior=left, superior=up,
    with the anterior contour (gray = rejected, orange = shaft inliers),
    the fitted cortex line, and each bump point with its perpendicular to
    the line."""
    sy, sz = p['spacing_yz']
    plot_sagittal_label(ax, p['bone'], 'mediumseagreen', sy, sz)
    plot_sagittal_label(ax, p['cartilage'], 'steelblue', sy, sz)

    c = p['cortex']
    zs, ys, inl = c['zs'], c['ys'], c['inlier_mask']
    ax.scatter(ys[~inl] * sy, zs[~inl] * sz, color='gray', s=8, zorder=5)
    ax.scatter(ys[inl] * sy, zs[inl] * sz, color='orange', s=10, zorder=6)
    dz, dy = c['direction']
    pz, py = c['point']
    z_lo, z_hi = zs.min() * sz, zs.max() * sz
    if abs(dz) > 1e-9:
        t_lo, t_hi = (z_lo - pz) / dz, (z_hi - pz) / dz
        ax.plot([py + t_lo * dy, py + t_hi * dy], [z_lo, z_hi], color='red', linewidth=1.8)

    p0, d = np.array([pz, py]), np.array([dz, dy])
    for idx_point, mm, color in ((p['bone_point'], p['bone_bump_mm'], 'crimson'),
                                 (p['cartilage_point'], p['cartilage_bump_mm'], 'darkmagenta')):
        if idx_point is None or mm is None:
            continue
        pt = np.array([idx_point[1] * sz, idx_point[0] * sy])
        foot = p0 + np.dot(pt - p0, d) * d
        ax.plot([pt[1], foot[1]], [pt[0], foot[0]], color=color, linestyle='--', linewidth=1.5)
        ax.scatter(pt[1], pt[0], color=color, edgecolor='black', s=70, zorder=8)
        ax.text(pt[1], pt[0], f'  {mm:.1f}mm', color=color, fontsize=fontsize, fontweight='bold')
    ax.set_aspect('equal')
    ax.invert_xaxis()


def panel_title(p):
    cart = f"{p['cartilage_bump_mm']:.1f}" if p['cartilage_bump_mm'] is not None else 'n/a'
    return f"bone={p['bone_bump_mm']:.1f}mm, cartilage={cart}mm"


def landmark_points(case, lm, result):
    d = result.get('details')
    if not d:
        return []
    x = d['x']
    pts = []
    if d['bone_point'] is not None:
        pts.append(('bone bump apex', (x, *d['bone_point'])))
    if d['cartilage_point'] is not None:
        pts.append(('cartilage bump apex', (x, *d['cartilage_point'])))
    return pts
