"""
Trochlear landmarks from the femur CARTILAGE: facet peaks + sulcus on an
axial slice, via convex-hull convexity defects of the isolated anterior
cartilage component.

Two slice-selection rules are used by different measurements:

- reference_slice (Dejour reference slice, per Dejour, Giovannetti de
  Sanctis et al., KSSTA 2025;33(8):2833-2846): scanning cranial -> caudal, the FIRST slice whose
  cartilage sulcus depth clears 2 mm -- "the first image with the medial
  facet beginning to form". Used by LTI / sulcus angle / facet depth and
  as the trochlear-groove point for TT-TG.
- deepest_cartilage_groove: the slice of MAXIMUM relative sulcus depth
  across the whole trochlear range. Used to choose the sagittal slice for
  bump height ("the deepest part of the groove").

Both stop at the first slice where the femur bone splits into two
components -- once the condyles separate, we've left the trochlea.

No rotation step is needed: convexity-defect detection on the closed
cartilage contour is orientation-independent.
"""
import numpy as np
import cv2
from scipy import ndimage

from .volume import count_components


# ------------------------------------------------- component isolation --

def _touches_bone(component_mask, bone2d, min_fraction=0.5, tol_px=2):
    """True if this cartilage component genuinely sits on the bone
    surface: for most of its columns, bone is present within tol_px
    pixels directly POSTERIOR to the component's deepest pixel.
    Directional and column-wise, not omnidirectional -- an omnidirectional
    adjacency check can't distinguish real trochlear cartilage (only
    ~15-20% of its pixels ever directly touch bone, since cartilage has
    thickness) from a component floating in an open notch cavity."""
    cols = np.where(component_mask.any(axis=1))[0]
    if len(cols) == 0:
        return False
    hits = 0
    for c in cols:
        ys = np.where(component_mask[c])[0]
        y_deep = ys.max()
        if bone2d[c, y_deep + 1:y_deep + 1 + tol_px].any():
            hits += 1
    return (hits / len(cols)) >= min_fraction


def most_anterior_component(mask2d, bone2d=None):
    """Cartilage at a single slice can include multiple disconnected
    bodies. Restricting to components that actually sit on the bone
    surface, then picking whichever is most ANTERIOR (smallest y),
    isolates the trochlear cartilage specifically."""
    lbl, n = ndimage.label(mask2d, structure=np.ones((3, 3)))
    if n == 0:
        return mask2d
    labels = list(range(1, n + 1))
    if bone2d is not None:
        labels = [c for c in labels if _touches_bone(lbl == c, bone2d)]
        if not labels:
            return np.zeros_like(mask2d)
    if len(labels) == 1:
        return lbl == labels[0]
    min_y_per_component = [np.argwhere(lbl == c)[:, 1].min() for c in labels]
    most_anterior = labels[int(np.argmin(min_y_per_component))]
    return lbl == most_anterior


# --------------------------------------------- convexity-defect detection --

def _faces_anterior(pts, component_mask, bone2d, tol_px=1):
    """True where a contour point has OPEN SPACE (neither cartilage nor
    bone) immediately anterior to it -- a genuine articular-surface point
    always has air in front of it; a point on the deep side of the closed
    contour doesn't."""
    occupied = component_mask | bone2d
    result = np.zeros(len(pts), dtype=bool)
    for i, (x, y) in enumerate(pts):
        y_lo = max(y - tol_px, 0)
        result[i] = not occupied[x, y_lo:y].any()
    return result


def find_peaks_and_sulcus(component_mask, bone2d, spacing_x, min_depth_mm=1.5):
    """Two facet peaks + one sulcus point on an isolated cartilage
    component via convex-hull convexity defects.

    Hull and defects are computed on the FULL closed contour -- cv2 treats
    an open point sequence as implicitly closed, so trimming first would
    fabricate a fake closing edge that can out-rank the real sulcus. A
    defect is only eligible if its far point (sulcus) AND both flanking
    hull points (peaks) face anterior. min_depth_mm is physical, not
    pixels, since voxel spacing varies between scans.

    Returns {'hull_points', 'defects', 'peak1', 'peak2', 'sulcus',
    'contour', 'depth_mm'} or None if no eligible defect exists."""
    mask_u8 = (component_mask.T.astype(np.uint8)) * 255
    contours, _ = cv2.findContours(mask_u8, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    if not contours:
        return None
    contour = max(contours, key=cv2.contourArea)
    if len(contour) < 5:
        return None

    hull_idx = cv2.convexHull(contour, returnPoints=False)
    hull_idx = np.sort(hull_idx, axis=0)
    if len(hull_idx) < 3:
        return None
    hull_points = contour[hull_idx[:, 0]]

    defects = cv2.convexityDefects(contour, hull_idx)
    if defects is None:
        return None

    min_depth_px = min_depth_mm / spacing_x
    pts = contour[:, 0, :]
    far_ok = _faces_anterior(pts[defects[:, 2]], component_mask, bone2d)
    start_ok = _faces_anterior(pts[defects[:, 0]], component_mask, bone2d)
    end_ok = _faces_anterior(pts[defects[:, 1]], component_mask, bone2d)
    eligible = np.where(far_ok & start_ok & end_ok & (defects[:, 3] / 256.0 >= min_depth_px))[0]
    if len(eligible) == 0:
        return None

    depths = defects[eligible, 3] / 256.0
    best = eligible[np.argmax(depths)]
    s_idx, e_idx, f_idx, depth_px = defects[best]
    depth_mm = (depth_px / 256.0) * spacing_x

    return {
        'hull_points': hull_points[:, 0, :],
        'defects': defects,
        'peak1': contour[s_idx][0],
        'peak2': contour[e_idx][0],
        'sulcus': contour[f_idx][0],
        'contour': contour,
        'depth_mm': depth_mm,
    }


# ------------------------------------------------------- slice selection --

def find_reference_slice(case, min_depth_mm=2.0):
    """Dejour reference slice: scan cranial (high z) -> caudal (low z),
    return the FIRST slice whose cartilage sulcus depth clears
    min_depth_mm. Stops (without a result) at the first z where the femur
    bone splits into 2 components.

    Returns {'z', 'depth_mm', 'peak1', 'peak2', 'sulcus', 'point',
    'hull_points', 'bone_slice', 'cartilage_component', 'note': None} on
    success, or {'z': None, 'note': str} on failure. 'point' is the sulcus
    as (x, y, z) -- the trochlear groove point used by TT-TG."""
    bone, cart = case['bone'], case['cartilage']
    spacing_x = case['spacing'][0]

    z_nonempty = np.where(cart.any(axis=(0, 1)))[0]
    if len(z_nonempty) == 0:
        return {'z': None, 'note': 'no cartilage found in this case'}

    for z in range(int(z_nonempty.max()), int(z_nonempty.min()) - 1, -1):
        bone_sl = bone[:, :, z]
        _, comps = count_components(bone_sl)
        if len(comps) >= 2:
            return {'z': None, 'note': f'reached split condyles (z={z}) before finding a valid reference slice'}
        result = sulcus_on_slice(bone_sl, cart[:, :, z], spacing_x, z, min_depth_mm=min_depth_mm)
        if result is not None:
            return result

    return {'z': None, 'note': 'no slice found with a real sulcus above the split boundary'}


def sulcus_on_slice(bone_sl, cart_sl, spacing_x, z, min_depth_mm=0.0):
    """Facet peaks + sulcus on one axial slice (the reference-slice result
    format), or None if this slice has no eligible sulcus. Also used for a
    slice a reviewer picks, with no depth floor."""
    component = most_anterior_component(cart_sl, bone_sl)
    if component.sum() == 0:
        return None
    result = find_peaks_and_sulcus(component, bone_sl, spacing_x, min_depth_mm=min_depth_mm)
    if result is None:
        return None
    return {
        'z': int(z), 'depth_mm': result['depth_mm'],
        'peak1': result['peak1'], 'peak2': result['peak2'], 'sulcus': result['sulcus'],
        'point': (int(result['sulcus'][0]), int(result['sulcus'][1]), int(z)),
        'hull_points': result['hull_points'],
        'bone_slice': bone_sl, 'cartilage_component': component,
        'note': None,
    }


def find_deepest_cartilage_groove(case, min_depth_mm=0.0):
    """The slice of MAXIMUM relative cartilage sulcus depth between the
    proximal trochlea and the condyle split (same per-slice isolation as
    find_reference_slice). min_depth_mm defaults to 0 -- since the max is
    taken over the whole range, a positive floor only risks dropping a
    case entirely rather than changing which slice wins.

    Returns {'point': (x, y, z), 'z', 'depth_mm', 'peak1', 'peak2',
    'sulcus', 'note': None} or {'z': None, 'note': str}."""
    bone, cart = case['bone'], case['cartilage']
    spacing_x = case['spacing'][0]
    z_nonempty = np.where(cart.any(axis=(0, 1)))[0]
    if len(z_nonempty) == 0:
        return {'z': None, 'note': 'no cartilage found in this case'}

    best = None
    for z in range(int(z_nonempty.max()), int(z_nonempty.min()) - 1, -1):
        bone_sl = bone[:, :, z]
        _, comps = count_components(bone_sl)
        if len(comps) >= 2:
            break
        cart_sl = cart[:, :, z]
        component = most_anterior_component(cart_sl, bone_sl)
        if component.sum() == 0:
            continue
        result = find_peaks_and_sulcus(component, bone_sl, spacing_x, min_depth_mm=min_depth_mm)
        if result is None:
            continue
        if best is None or result['depth_mm'] > best['depth_mm']:
            best = {'z': z, 'depth_mm': result['depth_mm'], 'peak1': result['peak1'],
                    'peak2': result['peak2'], 'sulcus': result['sulcus']}

    if best is None:
        return {'z': None, 'note': 'no eligible cartilage slice found in scan range'}

    x, y = int(best['sulcus'][0]), int(best['sulcus'][1])
    return {'point': (x, y, best['z']), **best, 'note': None}
