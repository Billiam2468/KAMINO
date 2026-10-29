"""
MRI scans (NIfTI .nii / .nii.gz) for display under the segmentation in the
review app. Measurements never use the scan itself.

A small NIfTI-1 reader (numpy + gzip, no extra dependencies). Scans the
pipeline produces share their segmentation's voxel grid (the segmentation is
predicted on the scan's own grid), so the scan's stored array lines up
voxel-for-voxel with the segmentation's and the segmentation's reorientation
applies to both. A scan on a different grid is resampled per slice through
scanner coordinates (see DisplayScan).
"""
import gzip
import os
import re
import struct

import numpy as np

from .volume import reorient_to_lps

_DTYPES = {2: np.uint8, 4: np.int16, 8: np.int32, 16: np.float32, 64: np.float64, 256: np.int8,
           512: np.uint16, 768: np.uint32}
SCAN_NAMES = ('water.nii.gz', 'water.nii')


def _open(path):
    return gzip.open(path, 'rb') if path.endswith('.gz') else open(path, 'rb')


def read_nifti_header(path):
    with _open(path) as f:
        raw = f.read(352)
    if len(raw) < 348:
        raise ValueError(f'{path}: not a NIfTI file')
    endian = '<' if struct.unpack('<i', raw[:4])[0] == 348 else '>'
    if struct.unpack(endian + 'i', raw[:4])[0] != 348:
        raise ValueError(f'{path}: only NIfTI-1 is supported')
    dim = struct.unpack(endian + '8h', raw[40:56])
    h = {
        'endian': endian,
        'shape': tuple(int(d) for d in dim[1:1 + dim[0]]),
        'datatype': struct.unpack(endian + 'h', raw[70:72])[0],
        'vox_offset': int(struct.unpack(endian + 'f', raw[108:112])[0]),
        'scl_slope': struct.unpack(endian + 'f', raw[112:116])[0],
        'scl_inter': struct.unpack(endian + 'f', raw[116:120])[0],
        'descrip': raw[148:228].split(b'\0')[0].decode('latin-1', 'replace'),
        'sform_code': struct.unpack(endian + 'h', raw[254:256])[0],
        'srow': np.array(struct.unpack(endian + '12f', raw[280:328]), dtype=float).reshape(3, 4),
    }
    return h


def read_nifti(path):
    h = read_nifti_header(path)
    if h['datatype'] not in _DTYPES:
        raise ValueError(f"{path}: unsupported NIfTI datatype {h['datatype']}")
    dt = np.dtype(_DTYPES[h['datatype']]).newbyteorder(h['endian'])
    n = int(np.prod(h['shape']))
    with _open(path) as f:
        f.read(h['vox_offset'])
        data = np.frombuffer(f.read(n * dt.itemsize), dtype=dt, count=n)
    data = data.reshape(h['shape'], order='F').astype(np.float32)
    if h['scl_slope'] not in (0.0, 1.0) or h['scl_inter'] != 0.0:
        data = data * (h['scl_slope'] or 1.0) + h['scl_inter']
    return data, h


def laterality_from_nifti(path):
    """'L'/'R' from a 'laterality=L' tag in the NIfTI description field
    (written by the sorting step), else None."""
    try:
        m = re.search(r'laterality\s*[=:]\s*([LR])\b', read_nifti_header(path)['descrip'], re.I)
    except Exception:
        return None
    return m.group(1).upper() if m else None


def same_grid(nifti_header, nrrd_header, tol_mm=0.05, origin_tol_mm=0.5):
    """True if the NIfTI scan and NRRD segmentation share one voxel grid.
    NIfTI stores RAS world coordinates, NRRD here LPS: flip x and y."""
    shape = tuple(nifti_header['shape'][:3])
    if shape != tuple(int(s) for s in nrrd_header['sizes'][:3]):
        return False, f'scan is {shape}, segmentation is {tuple(int(v) for v in nrrd_header["sizes"][:3])}'
    flip = np.array([-1.0, -1.0, 1.0])
    nii_dirs = (nifti_header['srow'][:, :3] * flip[:, None]).T     # row i = array axis i, LPS
    nii_origin = nifti_header['srow'][:, 3] * flip
    seg_dirs = np.array(nrrd_header['space directions'], dtype=float)
    seg_origin = np.array(nrrd_header['space origin'], dtype=float)
    if np.abs(nii_dirs - seg_dirs).max() > tol_mm or np.abs(nii_origin - seg_origin).max() > origin_tol_mm:
        return False, 'scan and segmentation are positioned differently'
    return True, None


def _lps_affine(nifti_header):
    """(3x3 matrix whose columns are the array axes, origin) in LPS mm."""
    flip = np.array([-1.0, -1.0, 1.0])
    return nifti_header['srow'][:, :3] * flip[:, None], nifti_header['srow'][:, 3] * flip


def _window(data):
    """The scan as 0-255 display values (1st to 99.5th percentile of the
    non-zero voxels)."""
    sample = data[::3, ::3, ::3]                     # the window doesn't need every voxel
    vals = sample[sample > 0]
    lo, hi = (np.percentile(vals, [1, 99.5]) if vals.size else (0.0, 1.0))
    return (np.clip((data - lo) / max(hi - lo, 1e-6), 0, 1) * 255).astype(np.uint8)


def _resample(raw, inv, o, seg_header, x, y, z):
    """raw (a scan's display array) at the segmentation's reoriented voxels
    (x, y, z), by nearest neighbour through scanner coordinates; 0 outside
    the scan."""
    from .volume import _reoriented_to_raw_index
    idx = _reoriented_to_raw_index(seg_header, x, y, z)
    dirs = np.array(seg_header['space directions'], dtype=float)
    origin = np.array(seg_header['space origin'], dtype=float)
    world = origin + sum(idx[k][..., None] * dirs[k] for k in range(3))
    ijk = np.rint((world - o) @ inv.T).astype(int)
    out = np.zeros(np.shape(x), dtype=np.uint8)
    inside = np.all((ijk >= 0) & (ijk < np.array(raw.shape)), axis=-1)
    i = ijk[inside]
    out[inside] = raw[i[:, 0], i[:, 1], i[:, 2]]
    return out


def _scan_index(inv, o, seg_header, x, y, z):
    """The scan's (continuous) array index at the segmentation's reoriented
    voxels (x, y, z), shape (..., 3)."""
    from .volume import _reoriented_to_raw_index
    idx = _reoriented_to_raw_index(seg_header, x, y, z)
    dirs = np.array(seg_header['space directions'], dtype=float)
    origin = np.array(seg_header['space origin'], dtype=float)
    world = origin + sum(idx[k][..., None] * dirs[k] for k in range(3))
    return (world - o) @ inv.T


def _sample_inplane_linear(raw, c):
    """raw at continuous array index c (..., 3): bilinear within each of the
    scan's own slices (array axes 0 and 1), nearest between them -- an axial
    series' slices are ~5 mm apart, too far to blend. 0 outside the scan."""
    i, j, k = c[..., 0], c[..., 1], np.rint(c[..., 2]).astype(int)
    n = raw.shape
    inside = (i >= 0) & (i <= n[0] - 1) & (j >= 0) & (j <= n[1] - 1) & (k >= 0) & (k < n[2])
    i, j, k = i[inside], j[inside], k[inside]
    i0, j0 = np.clip(np.floor(i).astype(int), 0, n[0] - 2), np.clip(np.floor(j).astype(int), 0, n[1] - 2)
    u, v = i - i0, j - j0
    out = np.zeros(c.shape[:-1], dtype=np.uint8)
    out[inside] = np.rint(raw[i0, j0, k] * (1 - u) * (1 - v) + raw[i0 + 1, j0, k] * u * (1 - v)
                          + raw[i0, j0 + 1, k] * (1 - u) * v + raw[i0 + 1, j0 + 1, k] * u * v)
    return out


def _lps_to_reoriented_index(seg_header, pts):
    """Scanner LPS mm points (n, 3) -> the segmentation's reoriented
    (continuous) voxel index; the inverse of volume.index_to_lps."""
    from .volume import _lps_axis_order
    order, signs = _lps_axis_order(seg_header)
    dirs = np.array(seg_header['space directions'], dtype=float)
    origin = np.array(seg_header['space origin'], dtype=float)
    raw = (np.asarray(pts, dtype=float) - origin) @ np.linalg.inv(dirs)
    sizes = [int(s) for s in seg_header['sizes'][:3]]
    out = np.empty_like(raw)
    for k in range(3):
        v = raw[:, order[k]]
        out[:, k] = sizes[order[k]] - 1 - v if signs[k] < 0 else v
    return out


class DisplayScan:
    """A scan prepared for display on a segmentation's reoriented grid.

    Same grid as the segmentation (the normal case): reoriented once, slices
    are plain array slices. Different grid (e.g. a segmentation saved at a
    finer resolution): each requested slice is resampled from the scan by
    nearest neighbour through both files' scanner coordinates."""

    def __init__(self, path, nrrd_header, reoriented_shape):
        data, h = read_nifti(path)
        if data.ndim == 4:
            data = data[..., 0]
        self.raw = _window(data)
        self.header, self.seg_header, self.shape = h, nrrd_header, tuple(reoriented_shape)
        ok, why = same_grid(h, nrrd_header)
        self.same_grid = ok
        self.note = None if ok else f'resampled onto the segmentation grid ({why})'
        self.reoriented = reorient_to_lps(self.raw, nrrd_header) if ok else None
        if ok:
            self.raw = None                          # only the reoriented copy is needed
        else:
            m, o = _lps_affine(h)
            self._inv, self._o = np.linalg.inv(m), o

    def slice(self, plane, index):
        if self.reoriented is not None:
            return self.reoriented[:, :, index] if plane == 'axial' else self.reoriented[index, :, :]
        nx, ny, nz = self.shape
        if plane == 'axial':
            a, b = np.meshgrid(np.arange(nx), np.arange(ny), indexing='ij')
            x, y, z = a, b, np.full(a.shape, index)
        else:
            a, b = np.meshgrid(np.arange(ny), np.arange(nz), indexing='ij')
            x, y, z = np.full(a.shape, index), a, b
        return _resample(self.raw, self._inv, self._o, self.seg_header, x, y, z)


AXIAL_NAME = 'axial.nii.gz'


def find_axial(scan_path):
    """The sorting step's axial series next to a scan
    (sorted/<code>/axial.nii.gz), or None."""
    if not scan_path:
        return None
    p = os.path.join(os.path.dirname(scan_path), AXIAL_NAME)
    return p if os.path.exists(p) else None


class AxialDisplay:
    """The axial series from the sorting step (usually the PD axial, which
    reaches further down the tibia than the Dixon) shown on the
    segmentation's reoriented grid, extended downward as far as the axial
    goes: z_below extra slices, numbered -1, -2, ... under the
    segmentation's slice 0. Resampled by nearest neighbour through scanner
    coordinates, so the segmentation drawn over it sits where the Dixon puts
    the bones -- a patient who moved between the two series shows as an
    offset -- and a point clicked on it is already in the segmentation's
    frame.

    The axial is finer in-plane than the segmentation (0.44-0.66 mm pixels
    vs the Dixon's 0.9 mm), so it's sampled `factor` times finer than the
    segmentation's voxels across each slice (and smoothly): slices come back
    factor times larger along x and y, each fine pixel centred at voxel
    (s + 0.5) / factor - 0.5."""

    MAX_FACTOR = 3

    def __init__(self, path, nrrd_header, reoriented_shape):
        from .volume import get_voxel_spacing_mm
        data, h = read_nifti(path)
        if data.ndim == 4:
            data = data[..., 0]
        self.raw = _window(data)
        self.seg_header, self.shape = nrrd_header, tuple(reoriented_shape)
        m, o = _lps_affine(h)
        self._inv, self._o = np.linalg.inv(m), o
        last = np.array(data.shape[:3]) - 1
        corners = np.array([m @ (last * np.array(c)) + o for c in np.ndindex(2, 2, 2)])
        lowest = _lps_to_reoriented_index(nrrd_header, corners)[:, 2].min()
        self.z_below = max(0, int(np.floor(-lowest)))
        seg_mm = min(get_voxel_spacing_mm(nrrd_header)[:2])
        pixel_mm = min(np.linalg.norm(m[:, 0]), np.linalg.norm(m[:, 1]))
        self.factor = int(min(self.MAX_FACTOR, max(1, np.ceil(seg_mm / pixel_mm - 0.01))))
        self._axial_base = None

    def slice(self, plane, index):
        """axial: slice z = index (may be negative), (nx, ny) * factor;
        sagittal: column x = index over z = -z_below ... nz - 1,
        (ny * factor, nz + z_below)."""
        nx, ny, nz = self.shape
        f = self.factor
        fine = lambda n: (np.arange(n * f) + 0.5) / f - 0.5
        if plane == 'axial':
            if self._axial_base is None:     # the scan index is linear in z: work out the fine grid once
                a, b = np.meshgrid(fine(nx), fine(ny), indexing='ij')
                c0 = _scan_index(self._inv, self._o, self.seg_header, a, b, np.zeros(a.shape))
                c1 = _scan_index(self._inv, self._o, self.seg_header, *(np.zeros(1), np.zeros(1), np.ones(1)))
                c00 = _scan_index(self._inv, self._o, self.seg_header, *(np.zeros(1),) * 3)
                self._axial_base = (c0.astype(np.float32), (c1 - c00)[0].astype(np.float32))
            c0, dz = self._axial_base
            return _sample_inplane_linear(self.raw, c0 + np.float32(index) * dz)
        a, b = np.meshgrid(fine(ny), np.arange(-self.z_below, nz, dtype=float), indexing='ij')
        c = _scan_index(self._inv, self._o, self.seg_header, np.full(a.shape, float(index)), a, b)
        return _sample_inplane_linear(self.raw, c)


def find_scan(scans_dir, code):
    """The scan for one case in a scans folder: <code>/water.nii.gz (the
    pipeline's sorted/ layout), else <code>/<only .nii.gz>, else
    <code>.nii.gz / <code>_water.nii.gz directly in the folder."""
    if not scans_dir:
        return None
    sub = os.path.join(scans_dir, code)
    if os.path.isdir(sub):
        for name in SCAN_NAMES:
            if os.path.exists(os.path.join(sub, name)):
                return os.path.join(sub, name)
        niis = sorted(f for f in os.listdir(sub) if f.endswith(('.nii', '.nii.gz')) and 'fat' not in f.lower())
        if len(niis) == 1:
            return os.path.join(sub, niis[0])
    for name in (f'{code}_water.nii.gz', f'{code}.nii.gz', f'{code}_water.nii', f'{code}.nii'):
        if os.path.exists(os.path.join(scans_dir, name)):
            return os.path.join(scans_dir, name)
    return None
