"""Small drawing helpers shared by the measurement QC panels."""
import numpy as np
from matplotlib.colors import ListedColormap


def crop_box(masks, spacing_xy, margin_mm, shape):
    """(x0, x1, y0, y1) bounding all non-empty 2D masks plus margin_mm."""
    sx, sy = spacing_xy
    xs0, xs1, ys0, ys1 = [], [], [], []
    for m in masks:
        if m is None or not m.any():
            continue
        xs = np.where(m.any(axis=1))[0]
        ys = np.where(m.any(axis=0))[0]
        xs0.append(xs.min()); xs1.append(xs.max()); ys0.append(ys.min()); ys1.append(ys.max())
    x0 = max(min(xs0) - int(margin_mm / sx), 0)
    x1 = min(max(xs1) + int(margin_mm / sx), shape[0])
    y0 = max(min(ys0) - int(margin_mm / sy), 0)
    y1 = min(max(ys1) + int(margin_mm / sy), shape[1])
    return x0, x1, y0, y1


def tuberosity_slice(tibia, z, spacing_z):
    """The tibia's axial slice through the tuberosity at z, and how far
    below the segmentation the point is (mm; None when it's inside). A
    tuberosity marked on the axial series can lie below the segmentation:
    the lowest segmented tibia slice then stands in, to show the point
    against (drawn as an outline, see stand_in_outline)."""
    if 0 <= z < tibia.shape[2]:
        return tibia[:, :, z], None
    zs = np.where(tibia.any(axis=(0, 1)))[0]
    if z >= 0 or not len(zs):
        return np.zeros(tibia.shape[:2], dtype=tibia.dtype), None
    lowest = int(zs.min())
    return tibia[:, :, lowest], (lowest - z) * spacing_z


def stand_in_outline(ax, X, Y, mask, below_mm, fontsize):
    """The lowest segmented tibia slice as a dashed outline, with a note
    that the tuberosity was marked below it (see tuberosity_slice)."""
    if mask.any():
        ax.contour(X, Y, mask.astype(float), levels=[0.5], colors='darkorange', linestyles='--', linewidths=1.2)
    ax.text(0.5, 0.02, f"tibia outline: lowest segmented slice\n(tuberosity marked {below_mm:.0f} mm below it, on the axial scan)",
            transform=ax.transAxes, ha='center', va='bottom', fontsize=fontsize, color='darkorange')


def crop(mask, box):
    x0, x1, y0, y1 = box
    return mask[x0:x1, y0:y1]


def mesh_mm(box, spacing_xy):
    x0, x1, y0, y1 = box
    X, Y = np.meshgrid(np.arange(x0, x1), np.arange(y0, y1), indexing='ij')
    return X * spacing_xy[0], Y * spacing_xy[1]


def plot_mask(ax, X, Y, mask, cmap, alpha=0.6):
    m = mask.astype(float)
    ax.pcolormesh(X, Y, np.ma.masked_where(m == 0, m), cmap=cmap, vmin=0, vmax=1.5,
                  shading='nearest', alpha=alpha)


def finish_axial(ax, box, spacing_xy):
    """Equal aspect, anterior (smaller y) at the top, axes hidden."""
    x0, x1, y0, y1 = box
    sx, sy = spacing_xy
    ax.set_aspect('equal')
    ax.set_xlim(x0 * sx, x1 * sx)
    ax.set_ylim(y1 * sy, y0 * sy)
    ax.axis('off')


def plot_sagittal_label(ax, sl, color, sy, sz):
    """One label of a (y, z) sagittal slice mask as a flat color."""
    if not sl.any():
        return
    ys, zs = np.where(sl)
    y0, y1 = ys.min(), ys.max() + 1
    z0, z1 = zs.min(), zs.max() + 1
    c = sl[y0:y1, z0:z1].astype(float)
    Y, Z = np.meshgrid(np.arange(y0, y1) * sy, np.arange(z0, z1) * sz, indexing='ij')
    ax.pcolormesh(Y, Z, np.ma.masked_where(c == 0, c), cmap=ListedColormap(['none', color]),
                  vmin=0, vmax=1, shading='nearest')
