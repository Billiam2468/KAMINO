"""
Segmentation I/O and the shared voxel/physical geometry every measurement
builds on.

Frames used throughout this package:

- The REORIENTED index frame (x, y, z): the stored array permuted/flipped
  by reorient_to_lps so axis0 increases toward patient LEFT, axis1 toward
  POSTERIOR, axis2 toward SUPERIOR. Every landmark is reported in this
  frame. Relative angles/distances (LTI, TT-TG, KRA, ...) are computed in
  it directly, scaled by voxel spacing -- a whole-frame rotation cancels
  out of a relative quantity.
- TRUE physical LPS mm: the full 'space directions' transform, including
  any oblique in-plane rotation baked into the header that reorient_to_lps
  drops (it only handles permutations and sign flips). Needed for ABSOLUTE
  directions ("most anterior" for the tibial tuberosity), for rendering
  axial QC images that look like the real scan, and for exporting landmarks
  to 3D Slicer.

Loading also normalizes 3D Slicer "layered" exports (4D arrays where
overlapping segments are stored in extra layers) into a single 3D
standard-label map -- see normalize_slicer_layers.
"""
import os
import re
import numpy as np
import nrrd
from scipy import ndimage

from .labels import FEM_B, FEM_C, PAT_B, PAT_C, TIB_B, STANDARD_LABELS

LATERALITY_HEADER_KEY = 'laterality'


# ---------------------------------------------------------------- geometry --

def _lps_axis_order(header):
    """(order, signs) so permuting array axes by `order` then flipping
    axis i wherever signs[i]<0 reorients the array so axis0 increases
    along +L, axis1 along +P, axis2 along +S -- necessary because
    different source files store their array axes in different orders
    relative to physical space."""
    space = header.get('space', '').lower()
    assert 'left' in space and 'posterior' in space and 'superior' in space, \
        f"Unsupported space convention: {header.get('space')!r}"
    directions_signed = np.array(header['space directions'], dtype=float)
    directions = np.abs(directions_signed)
    order, signs, used = [], [], set()
    for phys_axis in range(3):
        col = directions[:, phys_axis].copy()
        col[list(used)] = -1
        array_axis = int(np.argmax(col))
        used.add(array_axis)
        order.append(array_axis)
        signs.append(1 if directions_signed[array_axis, phys_axis] >= 0 else -1)
    return order, signs


def reorient_to_lps(data, header):
    order, signs = _lps_axis_order(header)
    reoriented = np.transpose(data, axes=order)
    for phys_axis, sign in enumerate(signs):
        if sign < 0:
            reoriented = np.flip(reoriented, axis=phys_axis)
    return np.ascontiguousarray(reoriented)


def get_voxel_spacing_mm(header):
    """(spacing_x, spacing_y, spacing_z) in mm/voxel along the volume's
    axes AFTER reorient_to_lps."""
    order, _ = _lps_axis_order(header)
    directions = np.array(header['space directions'], dtype=float)
    return tuple(float(np.linalg.norm(directions[a])) for a in order)


def _reoriented_to_raw_index(header, x, y, z):
    order, signs = _lps_axis_order(header)
    raw_size = [int(header['sizes'][a]) for a in range(3)]
    reor_vals = [np.asarray(x, dtype=float), np.asarray(y, dtype=float), np.asarray(z, dtype=float)]
    raw_idx = [None, None, None]
    for k in range(3):
        raw_axis = order[k]
        v = reor_vals[k]
        if signs[k] < 0:
            v = raw_size[raw_axis] - 1 - v
        raw_idx[raw_axis] = v
    return raw_idx


def true_physical_lp(header, x, y, z):
    """Convert reoriented-frame (x, y, z) index point(s) to TRUE physical
    (L, P) mm using the full space-directions matrix. x, y may be scalars
    or arrays of the same shape; z may be a scalar (one slice) or an array
    matching x/y. Returns (L, P) with the shape of x/y."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    if z.shape != x.shape:
        z = np.broadcast_to(z, x.shape)
    directions = np.array(header['space directions'], dtype=float)
    origin = np.array(header['space origin'], dtype=float)
    raw_idx = _reoriented_to_raw_index(header, x, y, z)
    L = origin[0] + raw_idx[0] * directions[0, 0] + raw_idx[1] * directions[1, 0] + raw_idx[2] * directions[2, 0]
    P = origin[1] + raw_idx[0] * directions[0, 1] + raw_idx[1] * directions[1, 1] + raw_idx[2] * directions[2, 1]
    return L, P


def index_to_lps(header, x, y, z):
    """Reoriented-frame (x, y, z) index point -> TRUE physical (L, P, S) mm,
    the coordinate system 3D Slicer markups use (coordinateSystem 'LPS')."""
    directions = np.array(header['space directions'], dtype=float)
    origin = np.array(header['space origin'], dtype=float)
    raw_idx = _reoriented_to_raw_index(header, x, y, z)
    return tuple(float(origin[c] + sum(raw_idx[k] * directions[k, c] for k in range(3))) for c in range(3))


def largest_component_3d(mask):
    """Keep only the largest 3D-connected (26-connectivity) component.
    Stray single-voxel specks disconnected from the real anatomy can sit
    slightly beyond the true surface and win a raw extremal-point search;
    on a clean segmentation this is a no-op."""
    box = ndimage.find_objects(mask.astype(np.uint8))       # label only the mask's bounding box: same result, faster
    if not box:
        return mask
    box = box[0]
    lbl, n = ndimage.label(mask[box], structure=np.ones((3, 3, 3)))
    if n <= 1:
        return mask
    sizes = np.bincount(lbl.ravel())[1:]
    out = np.zeros_like(mask, dtype=bool)
    out[box] = lbl == int(np.argmax(sizes)) + 1
    return out


def count_components(sl, min_component_size=200):
    """2D 8-connected components of a slice mask with more than
    min_component_size pixels. Returns (label_image, component_ids)."""
    lbl, n = ndimage.label(sl, structure=np.ones((3, 3)))
    sizes = ndimage.sum(sl, lbl, range(1, n + 1))
    comps = [i + 1 for i, s in enumerate(sizes) if s > min_component_size]
    return lbl, comps


# ------------------------------------------------- 3D Slicer layered files --

_SEGMENT_KEY = re.compile(r'^Segment(\d+)_(\w+)$')

# Keyword rules for mapping a Slicer segment NAME to a standard label, used
# only when segments carry anatomical names (checked in order).
_NAME_RULES = (
    (FEM_C, ('femur', 'femoral'), ('cart',)),
    (PAT_C, ('patella', 'patellar'), ('cart',)),
    (FEM_B, ('femur', 'femoral'), ()),
    (PAT_B, ('patella', 'patellar'), ()),
    (TIB_B, ('tibia', 'tibial'), ()),
)


def _slicer_segments(header):
    """{segment_index: {'ID', 'Name', 'Layer', 'LabelValue', ...}} from a
    3D Slicer segmentation header, or {} if it has none."""
    segs = {}
    for key, value in header.items():
        m = _SEGMENT_KEY.match(key)
        if m:
            segs.setdefault(int(m.group(1)), {})[m.group(2)] = value
    return segs


def _standard_label_for_segment(seg):
    """Standard label for one Slicer segment: by anatomical name if it has
    one, else by its auto-generated 'Segment_k' ID/name (k = standard label
    -- confirmed on every Slicer-exported file in this project: Segment_1..5
    are femur bone, femur cartilage, patella bone, patella cartilage, tibia).
    Returns (label or None, how)."""
    name = str(seg.get('Name', '')).lower()
    for label, must_any, must_also in _NAME_RULES:
        if any(w in name for w in must_any) and all(w in name for w in must_also):
            return label, 'name'
    for field in ('ID', 'Name'):
        m = re.fullmatch(r'segment_(\d+)', str(seg.get(field, '')).strip().lower())
        if m and int(m.group(1)) in STANDARD_LABELS:
            return int(m.group(1)), 'segment order'
    return None, None


def normalize_slicer_layers(data, header):
    """Returns (data3d, header3d, notes).

    3D Slicer exports a segmentation with OVERLAPPING segments as a 4D
    array (layer, x, y, z): each layer is its own label map, and a segment
    is identified by (Segment{i}_Layer, Segment{i}_LabelValue) -- so the
    same label value can mean different structures in different layers
    (e.g. femur bone = value 1 in layer 0, femur cartilage = value 1 in
    layer 1). This rebuilds one 3D standard-label map from that metadata.
    Where two segments overlap (a few hundred voxels at a bone/cartilage
    interface), cartilage wins.

    Plain 3D files pass through unchanged -- unless their Slicer metadata
    names the segments anatomically with label values that differ from the
    standard ones, in which case they're relabeled the same way."""
    notes = []
    segs = _slicer_segments(header)

    if data.ndim == 3 and not segs:
        return data, header, notes
    if data.ndim not in (3, 4):
        raise ValueError(f'unsupported segmentation array shape {data.shape}')
    if data.ndim == 4 and not segs:
        raise ValueError(f'4D segmentation {data.shape} without 3D Slicer segment metadata')

    mapping, how = {}, set()
    for i, seg in segs.items():
        label, why = _standard_label_for_segment(seg)
        if label is None:
            raise ValueError(f"can't map Slicer segment {seg.get('Name')!r} to a standard label")
        layer = int(seg.get('Layer', 0)) if data.ndim == 4 else 0
        mapping[i] = (layer, int(seg['LabelValue']), label)
        how.add(why)

    identity = data.ndim == 3 and all(v == lab for _, v, lab in mapping.values())
    if identity:
        return data, header, notes

    out_shape = data.shape[1:] if data.ndim == 4 else data.shape
    out = np.zeros(out_shape, dtype=np.uint8)
    for _, (layer, value, label) in sorted(mapping.items(), key=lambda kv: kv[1][2] in (FEM_C, PAT_C)):
        layer_data = data[layer] if data.ndim == 4 else data
        m = layer_data == value
        if label in (FEM_C, PAT_C):
            out[m] = label
        else:
            out[m & (out == 0)] = label

    new_header = {k: v for k, v in header.items()
                  if not (k.startswith('Segment') or k.startswith('Segmentation'))}
    if data.ndim == 4:
        new_header['dimension'] = 3
        new_header['space directions'] = np.array(header['space directions'][1:], dtype=float)
        new_header['kinds'] = ['domain'] * 3
        new_header['sizes'] = np.array(out.shape)
        notes.append(f"converted 4D Slicer layered export ({data.shape[0]} layers) to 3D "
                     f"(label mapping by {' + '.join(sorted(how))})")
    else:
        new_header['sizes'] = np.array(out.shape)
        notes.append(f"relabeled Slicer segments to standard labels (mapping by {' + '.join(sorted(how))})")
    return out, new_header, notes


# ------------------------------------------------------------------- I/O --

SEGMENTATION_EXTENSIONS = ('.nrrd', '.seg.nrrd')


def case_name_from_path(path):
    base = os.path.basename(path)
    for ext in ('.seg.nrrd', '.nrrd'):
        if base.lower().endswith(ext):
            return base[:-len(ext)]
    return os.path.splitext(base)[0]


def find_segmentations(paths):
    """Expand files/directories into a sorted list of segmentation files."""
    found = []
    for p in paths:
        if os.path.isdir(p):
            for f in sorted(os.listdir(p)):
                if f.lower().endswith('.nrrd') and not f.startswith('.'):
                    found.append(os.path.join(p, f))
        elif os.path.isfile(p):
            found.append(p)
        else:
            raise FileNotFoundError(p)
    return found


def read_segmentation(path):
    """Read a segmentation NRRD in its stored (raw) array order, normalized
    to a 3D standard-label map. Returns (data, header, notes)."""
    data, header = nrrd.read(path)
    data, header, notes = normalize_slicer_layers(data, header)
    return data, header, notes


def header_laterality(header):
    """'L'/'R' from a 'laterality' NRRD header field, or None."""
    value = header.get(LATERALITY_HEADER_KEY)
    if value is None:
        return None
    return parse_side(value)


def parse_side(value):
    """'L'/'R' from common spellings (L, R, Left, Right, lt, rt), else None."""
    v = str(value).strip().lower()
    if v in ('l', 'left', 'lt'):
        return 'L'
    if v in ('r', 'right', 'rt'):
        return 'R'
    return None


def load_case(path, patella=True):
    """Load one segmentation into the dict every measurement consumes:

        'labels'    reoriented standard-label array (uint8)
        'bone'      femur bone      -- largest 3D component
        'cartilage' femur cartilage -- largest 3D component
        'tibia'     tibia bone      -- largest 3D component
        'patella'   patella bone    -- largest 3D component
        'header', 'spacing', 'path', 'name', 'notes'

    The largest-component filter matches the original pipeline exactly
    (CDI deliberately reads the unfiltered 'labels' array, as before).
    patella=False skips the patella mask (the review app doesn't need it)."""
    data, header, notes = read_segmentation(path)
    labels = reorient_to_lps(data, header).astype(np.uint8, copy=False)
    return {
        'path': path,
        'name': case_name_from_path(path),
        'labels': labels,
        'bone': largest_component_3d(labels == FEM_B),
        'cartilage': largest_component_3d(labels == FEM_C),
        'tibia': largest_component_3d(labels == TIB_B),
        'patella': largest_component_3d(labels == PAT_B) if patella else None,
        'header': header,
        'spacing': get_voxel_spacing_mm(header),
        'notes': notes,
    }


def write_with_laterality(src_path, dst_path, side):
    """Copy a segmentation to dst_path (may equal src_path) with a
    'laterality' header field set to 'L'/'R'. The array and every other
    header field are written back unchanged."""
    data, header = nrrd.read(src_path)
    header = dict(header)
    header[LATERALITY_HEADER_KEY] = side
    nrrd.write(dst_path, data, header)
