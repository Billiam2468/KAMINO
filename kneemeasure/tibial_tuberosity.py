"""
Tibial tuberosity (TT) point: the most anterior point of the tuberosity,
used by TT-TG/sTT-TG and tubercle lateralization.

Method (validated against user-identified ground-truth columns):

1. Fit ONE fixed reference plane from the clean bottom-of-shaft slices
   only: skip the first couple mm off the segmentation's distal cut edge,
   pool the next `use_mm` of slices into one 2D point cloud and PCA-fit a
   major/minor axis + reference centroid. A single fixed fit -- per-slice
   refitting drifts toward off-axis anatomy (e.g. the medial condyle).
2. Plane thickness (along the minor axis) = the NARROWEST minor-axis width
   across those bottom slices, plus a +2 mm margin.
3. Extend that rigid plane through the entire tibia. Candidates within
   plateau_margin_mm of the plateau or distal_margin_mm of the distal cut
   edge are excluded (both ends are prone to boundary artifacts).
4. "Most anterior" is rarely a sharp winner: among points within
   tie_tolerance_mm of the anteriority minimum, the most LATERAL is
   chosen (the true tuberosity is lateral of the ridge).
5. Within that column the most-anterior points form a z-plateau; its
   MIDDLE is taken.

All of this runs in TRUE physical (L, P) mm: "most anterior" is an
absolute direction, so the naive reoriented frame would be wrong by any
oblique in-plane rotation baked into the header.

Note: if the scan's field of view ends above the tuberosity, the "most
anterior point" found is not the tuberosity -- check the TT-TG/TTL QC.
"""
import numpy as np

from .volume import true_physical_lp


def _bottom_slices_axis(tibia_mask, header, spacing_z, skip_mm=2.0, use_mm=15.0):
    """(ref_centroid, major_dir, minor_dir, thickness_mm) fit from the
    clean bottom-of-shaft slices, in TRUE physical (L, P) mm, or None if
    there isn't enough clean bottom material to fit."""
    z_nonempty = np.where(tibia_mask.any(axis=(0, 1)))[0]
    if len(z_nonempty) == 0:
        return None
    z_lo = int(z_nonempty.min())
    lo_px, hi_px = skip_mm / spacing_z, (skip_mm + use_mm) / spacing_z
    clean_z = [z for z in z_nonempty if lo_px <= (z - z_lo) <= hi_px]

    pooled, centroids, spans = [], [], []
    for z in clean_z:
        pts = np.argwhere(tibia_mask[:, :, z])
        if len(pts) < 10:
            continue
        L, P = true_physical_lp(header, pts[:, 0], pts[:, 1], z)
        pts_mm = np.stack([L, P], axis=1)
        c = pts_mm.mean(axis=0)
        centroids.append(c)
        pooled.append(pts_mm - c)
    if len(pooled) < 3:
        return None

    pooled_all = np.concatenate(pooled, axis=0)
    cov = np.cov(pooled_all.T)
    eigvals, eigvecs = np.linalg.eigh(cov)
    minor_dir, major_dir = eigvecs[:, 0], eigvecs[:, 1]

    for z in clean_z:
        pts = np.argwhere(tibia_mask[:, :, z])
        if len(pts) < 10:
            continue
        L, P = true_physical_lp(header, pts[:, 0], pts[:, 1], z)
        pts_mm = np.stack([L, P], axis=1)
        proj = (pts_mm - pts_mm.mean(axis=0)) @ minor_dir
        spans.append(proj.max() - proj.min())
    thickness_mm = float(np.min(spans))

    return np.mean(centroids, axis=0), major_dir, minor_dir, thickness_mm


def find_tibial_tuberosity(case, knee_side, skip_mm=2.0, use_mm=15.0, width_margin_mm=2.0,
                            plateau_margin_mm=15.0, distal_margin_mm=5.0, tie_tolerance_mm=1.0):
    """Returns {'point': (x, y, z), 'ref_centroid', 'major_dir',
    'minor_dir', 'thickness_mm', 'plateau_boundary_z', 'distal_boundary_z',
    'note': None} on success, or {'point': None, 'note': str} on failure.
    'point' is in the reoriented index frame (matching case['tibia']).
    Lateral = smaller x for a RIGHT knee, larger x for a LEFT knee."""
    tibia = case['tibia']
    header = case['header']
    spacing_z = case['spacing'][2]

    axis = _bottom_slices_axis(tibia, header, spacing_z, skip_mm, use_mm)
    if axis is None:
        return {'point': None, 'note': 'not enough clean bottom-of-shaft slices to fit a reference plane'}
    ref_centroid, major_dir, minor_dir, thickness_mm = axis
    thickness_mm += width_margin_mm

    z_nonempty = np.where(tibia.any(axis=(0, 1)))[0]
    z_lo, z_hi = int(z_nonempty.min()), int(z_nonempty.max())
    plateau_boundary_z = z_hi - (plateau_margin_mm / spacing_z)
    distal_boundary_z = z_lo + (distal_margin_mm / spacing_z)

    coords = np.argwhere(tibia)
    L, P = true_physical_lp(header, coords[:, 0], coords[:, 1], coords[:, 2])
    coords_mm = np.stack([L, P], axis=1)
    offset = (coords_mm - ref_centroid) @ minor_dir
    in_slab = np.abs(offset) <= thickness_mm / 2.0
    below_plateau_margin = coords[:, 2] <= plateau_boundary_z
    above_distal_margin = coords[:, 2] >= distal_boundary_z
    eligible = in_slab & below_plateau_margin & above_distal_margin

    candidates = coords[eligible]
    candidates_P = P[eligible]
    if len(candidates) == 0:
        return {'point': None,
                'note': 'no tibia bone found inside the fixed reference plane between the plateau and distal margins',
                'plateau_boundary_z': plateau_boundary_z, 'distal_boundary_z': distal_boundary_z}

    min_P = candidates_P.min()
    near_tie = candidates_P <= (min_P + tie_tolerance_mm)
    tied_candidates = candidates[near_tie]
    if knee_side == 'R':
        winning_x = tied_candidates[:, 0].min()  # lateral = smaller x
    else:
        winning_x = tied_candidates[:, 0].max()  # lateral = larger x

    at_column = tied_candidates[tied_candidates[:, 0] == winning_x]
    z_mid = np.median(at_column[:, 2])
    best_i = np.argmin(np.abs(at_column[:, 2] - z_mid))
    point = tuple(int(v) for v in at_column[best_i])

    return {
        'point': point, 'ref_centroid': ref_centroid, 'major_dir': major_dir,
        'minor_dir': minor_dir, 'thickness_mm': thickness_mm,
        'plateau_boundary_z': plateau_boundary_z, 'distal_boundary_z': distal_boundary_z, 'note': None,
    }
