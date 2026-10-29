"""
Segmentation sanity checks, run on every case before (and independent of)
the measurements. They catch the problems that have actually turned up in
this project's cohorts: missing or tiny labels, cartilage segmented without
its bone (or floating off it), multi-piece labels, a tibia cut off too short
to measure, a patella in the wrong place, and Slicer layered exports.

Each check yields per-case columns plus human-readable warnings; cohort-
level outlier checks (e.g. a femur far smaller than the rest of the
cohort) are done afterwards in stats.flag_outliers.
"""
import numpy as np
from scipy import ndimage

from .labels import FEM_B, FEM_C, PAT_B, PAT_C, TIB_B, STANDARD_LABELS, LABEL_NAMES
from .qc_draw import plot_sagittal_label

# A label below this volume is treated as missing (a few stray voxels).
MIN_LABEL_VOLUME_CM3 = {FEM_B: 20.0, FEM_C: 1.0, PAT_B: 2.0, PAT_C: 0.3, TIB_B: 10.0}
# Minimum fraction of a cartilage label within ~1.5 mm of its own bone.
MIN_CARTILAGE_ON_BONE = 0.25
# Tibia shorter than this (superior-inferior) leaves too little to measure from.
# (A tibia under ~37 mm usually means the tuberosity is out of the field of
# view -- that is reported as a measurement caution, see pipeline.py.)
MIN_TIBIA_LENGTH_MM = 20.0
# A label whose largest 3D component holds less than this share of it is fragmented.
MIN_LARGEST_COMPONENT_FRACTION = 0.95

COLUMNS = (
    ['spacing_x_mm', 'spacing_y_mm', 'spacing_z_mm', 'n_slices']
    + [f'vol_{n}_cm3' for n in ('femur_bone', 'femur_cart', 'patella_bone', 'patella_cart', 'tibia')]
    + ['femur_cart_on_bone', 'patella_cart_on_bone', 'tibia_length_mm', 'femur_length_mm',
       'patella_anterior_of_femur_mm', 'n_warnings', 'warnings']
)

_VOL_COL = {FEM_B: 'vol_femur_bone_cm3', FEM_C: 'vol_femur_cart_cm3', PAT_B: 'vol_patella_bone_cm3',
            PAT_C: 'vol_patella_cart_cm3', TIB_B: 'vol_tibia_cm3'}


def _on_bone_fraction(cart, bone, spacing, tol_mm=1.5):
    if not cart.any() or not bone.any():
        return None
    iters = max(1, int(round(tol_mm / min(spacing))))
    near = ndimage.binary_dilation(bone, iterations=iters)
    return float((near & cart).sum() / cart.sum())


def check_case(case):
    """Returns (row_dict, warnings_list)."""
    labels, spacing = case['labels'], case['spacing']
    voxel_cm3 = spacing[0] * spacing[1] * spacing[2] / 1000.0
    row = {'spacing_x_mm': spacing[0], 'spacing_y_mm': spacing[1], 'spacing_z_mm': spacing[2],
           'n_slices': int(labels.shape[2])}
    warnings = list(case.get('notes', []))

    present = set(int(v) for v in np.unique(labels)) - {0}
    unexpected = sorted(present - set(STANDARD_LABELS))
    if unexpected:
        warnings.append(f'unexpected label values {unexpected} (expected 1-5)')

    masks = {}
    for lab in STANDARD_LABELS:
        m = labels == lab
        masks[lab] = m
        vol = float(m.sum() * voxel_cm3)
        row[_VOL_COL[lab]] = vol
        if vol < MIN_LABEL_VOLUME_CM3[lab]:
            warnings.append(f'{LABEL_NAMES[lab]} missing or tiny ({vol:.2f} cm3)')
            continue
        lbl, n = ndimage.label(m, structure=np.ones((3, 3, 3)))
        if n > 1:
            sizes = np.bincount(lbl.ravel())[1:]
            frac = sizes.max() / sizes.sum()
            if frac < MIN_LARGEST_COMPONENT_FRACTION:
                warnings.append(f'{LABEL_NAMES[lab]} fragmented ({n} pieces, largest {frac:.0%})')

    row['femur_cart_on_bone'] = _on_bone_fraction(masks[FEM_C], masks[FEM_B], spacing)
    row['patella_cart_on_bone'] = _on_bone_fraction(masks[PAT_C], masks[PAT_B], spacing)
    for col, name in (('femur_cart_on_bone', 'femur'), ('patella_cart_on_bone', 'patella')):
        v = row[col]
        if v is not None and v < MIN_CARTILAGE_ON_BONE:
            warnings.append(f'{name} cartilage mostly not touching {name} bone ({v:.0%}) -- check labels')

    def _z_extent_mm(m):
        zs = np.where(m.any(axis=(0, 1)))[0]
        return float((zs.max() - zs.min() + 1) * spacing[2]) if len(zs) else 0.0

    row['tibia_length_mm'] = _z_extent_mm(masks[TIB_B])
    row['femur_length_mm'] = _z_extent_mm(masks[FEM_B])
    if 0 < row['tibia_length_mm'] < MIN_TIBIA_LENGTH_MM:
        warnings.append(f"tibia only {row['tibia_length_mm']:.0f} mm long")

    row['patella_anterior_of_femur_mm'] = None
    if masks[PAT_B].any() and masks[FEM_B].any():
        dy = (np.argwhere(masks[FEM_B])[:, 1].mean() - np.argwhere(masks[PAT_B])[:, 1].mean()) * spacing[1]
        row['patella_anterior_of_femur_mm'] = float(dy)
        if dy < 10:
            warnings.append(f'patella not anterior of the femur ({dy:.0f} mm) -- check labels')
        pz = np.where(masks[PAT_B].any(axis=(0, 1)))[0]
        fz = np.where(masks[FEM_B].any(axis=(0, 1)))[0]
        if pz.max() < fz.min() or pz.min() > fz.max():
            warnings.append('patella outside the femur z-range -- check labels')

    row['n_warnings'] = len(warnings)
    row['warnings'] = '; '.join(warnings)
    return row, warnings


def overview_payload(case):
    """Lateral (side-view) projection of every label, for the segmentation
    overview grid -- quick to eyeball for missing/misplaced structures."""
    labels = case['labels']
    proj = {lab: (labels == lab).any(axis=0) for lab in STANDARD_LABELS}
    return {'proj': proj, 'spacing_yz': (case['spacing'][1], case['spacing'][2])}


def render_overview(ax, p, fontsize=7):
    sy, sz = p['spacing_yz']
    for lab, color in ((FEM_B, 'mediumseagreen'), (TIB_B, 'darkorange'), (FEM_C, 'steelblue'),
                       (PAT_B, 'darkviolet'), (PAT_C, 'red')):
        plot_sagittal_label(ax, p['proj'][lab], color, sy, sz)
    ax.set_aspect('equal')
    ax.invert_xaxis()
    ax.axis('off')
