"""
Static statistics figures (PNG):

- distribution_figure: one small-multiple panel per measurement showing
  every case's value (jittered strip over a thin box), with QC-flagged
  outliers marked and labeled by case.
- comparison_figure: the same per-measurement panels with one box+strip per
  cohort, for between-group comparison.

Colors are a validated categorical palette applied in fixed slot order
(groups are also named on the x axis, so identity is never color-only); the
flagged-outlier color is a reserved status color, always paired with a
case label.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SURFACE = '#fcfcfb'
TEXT = '#0b0b0b'
TEXT_2 = '#52514e'
GRID = '#e4e3df'
SERIES = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948']
CRITICAL = '#d03b3b'
BOX = '#8a8984'


def _style(ax, title, unit):
    ax.set_facecolor(SURFACE)
    ax.set_title(title, fontsize=9, color=TEXT, loc='left')
    if unit:
        ax.set_ylabel(unit, fontsize=8, color=TEXT_2)
    ax.tick_params(colors=TEXT_2, labelsize=7, length=0)
    ax.grid(axis='y', color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for s in ax.spines.values():
        s.set_visible(False)


def _vals(rows, col):
    out = []
    for r in rows:
        v = r.get(col)
        if v not in (None, ''):
            try:
                out.append((r, float(v)))
            except (TypeError, ValueError):
                pass
    return out


def _box(ax, x, pos, color=BOX):
    ax.boxplot([x], positions=[pos], widths=0.5, showfliers=False, patch_artist=True,
               boxprops=dict(facecolor='none', edgecolor=color, linewidth=1),
               whiskerprops=dict(color=color, linewidth=1), capprops=dict(color=color, linewidth=1),
               medianprops=dict(color=TEXT, linewidth=1.5))


def distribution_figure(rows, columns, units, flags, out_path, id_key='case', ncols=4, title=None):
    flagged = {(f[id_key], f['measurement']) for f in flags}
    cols = [c for c in columns if _vals(rows, c)]
    if not cols:
        return None
    nrows = int(np.ceil(len(cols) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3.4, nrows * 3.0), squeeze=False)
    fig.patch.set_facecolor(SURFACE)
    rng = np.random.default_rng(0)
    for i, col in enumerate(cols):
        ax = axes[i // ncols, i % ncols]
        pairs = _vals(rows, col)
        x = np.array([v for _, v in pairs])
        _box(ax, x, 0)
        jitter = rng.uniform(-0.18, 0.18, len(pairs))
        is_flag = np.array([(r.get(id_key), col) in flagged for r, _ in pairs])
        ax.scatter(jitter[~is_flag], x[~is_flag], s=22, color=SERIES[0], edgecolor=SURFACE,
                   linewidth=1, zorder=3)
        if is_flag.any():
            ax.scatter(jitter[is_flag], x[is_flag], s=30, color=CRITICAL, marker='D',
                       edgecolor=SURFACE, linewidth=1, zorder=4)
            for (r, v), j, f in zip(pairs, jitter, is_flag):
                if f:
                    ax.annotate(str(r.get(id_key)), (j, v), xytext=(6, 0), textcoords='offset points',
                                fontsize=6.5, color=TEXT, va='center')
        ax.set_xlim(-0.6, 0.9)
        ax.set_xticks([])
        _style(ax, f'{col}  (n={len(x)})', units.get(col, ''))
    for i in range(len(cols), nrows * ncols):
        axes[i // ncols, i % ncols].axis('off')
    fig.suptitle(title or 'Measurement distributions  (red diamond = flagged outlier, |robust z| > 3.5)',
                 fontsize=11, color=TEXT, x=0.01, ha='left')
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out_path, dpi=130, facecolor=SURFACE)
    plt.close(fig)
    return out_path


def comparison_figure(groups, columns, units, out_path, ncols=4):
    names = list(groups)
    cols = [c for c in columns if any(_vals(groups[g], c) for g in names)]
    if not cols:
        return None
    nrows = int(np.ceil(len(cols) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3.6, nrows * 3.1), squeeze=False)
    fig.patch.set_facecolor(SURFACE)
    rng = np.random.default_rng(0)
    for i, col in enumerate(cols):
        ax = axes[i // ncols, i % ncols]
        for k, g in enumerate(names):
            x = np.array([v for _, v in _vals(groups[g], col)])
            if len(x) == 0:
                continue
            color = SERIES[k % len(SERIES)]
            _box(ax, x, k)
            ax.scatter(k + rng.uniform(-0.17, 0.17, len(x)), x, s=16, color=color,
                       edgecolor=SURFACE, linewidth=0.8, zorder=3)
        ax.set_xticks(range(len(names)))
        ax.set_xticklabels(names, rotation=20, ha='right')
        _style(ax, col, units.get(col, ''))
    for i in range(len(cols), nrows * ncols):
        axes[i // ncols, i % ncols].axis('off')
    handles = [plt.Line2D([], [], marker='o', linestyle='', color=SERIES[k % len(SERIES)], label=g)
               for k, g in enumerate(names)]
    fig.legend(handles=handles, loc='upper right', ncol=len(names), frameon=False, fontsize=9,
               labelcolor=TEXT)
    fig.suptitle('Measurement comparison by cohort', fontsize=11, color=TEXT, x=0.01, ha='left')
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_path, dpi=130, facecolor=SURFACE)
    plt.close(fig)
    return out_path


def agreement_figure(pairs, out_path, ncols=4):
    """pairs: {label: [(auto, manual, case), ...]} -> for each pair a
    scatter (auto vs manual, identity line) above a Bland-Altman plot."""
    labels = [k for k, v in pairs.items() if len(v) >= 3]
    if not labels:
        return None
    ncols = min(ncols, len(labels))
    nblocks = int(np.ceil(len(labels) / ncols))
    fig, axes = plt.subplots(nblocks * 2, ncols, figsize=(ncols * 3.6, nblocks * 6.4), squeeze=False)
    fig.patch.set_facecolor(SURFACE)
    for i, lab in enumerate(labels):
        b, c = divmod(i, ncols)
        ax, bx = axes[2 * b, c], axes[2 * b + 1, c]
        x = np.array([p[0] for p in pairs[lab]])
        y = np.array([p[1] for p in pairs[lab]])
        lo, hi = float(min(x.min(), y.min())), float(max(x.max(), y.max()))
        pad = 0.05 * (hi - lo or 1)
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color=BOX, linewidth=1, linestyle='--')
        ax.scatter(y, x, s=18, color=SERIES[0], edgecolor=SURFACE, linewidth=0.8, zorder=3)
        ax.set_xlabel('manual', fontsize=7, color=TEXT_2)
        _style(ax, lab, 'automated')
        mean, diff = (x + y) / 2, x - y
        bias, sd = diff.mean(), diff.std(ddof=1)
        for v, ls in ((bias, '-'), (bias - 1.96 * sd, ':'), (bias + 1.96 * sd, ':')):
            bx.axhline(v, color=TEXT_2, linewidth=1, linestyle=ls)
        bx.axhline(0, color=GRID, linewidth=1)
        bx.scatter(mean, diff, s=18, color=SERIES[0], edgecolor=SURFACE, linewidth=0.8, zorder=3)
        bx.set_xlabel('mean of auto and manual', fontsize=7, color=TEXT_2)
        _style(bx, f'Bland-Altman: bias {bias:+.2f}, LoA {bias - 1.96 * sd:+.2f} / {bias + 1.96 * sd:+.2f}',
               'auto - manual')
    for i in range(len(labels), nblocks * ncols):
        b, c = divmod(i, ncols)
        axes[2 * b, c].axis('off')
        axes[2 * b + 1, c].axis('off')
    fig.suptitle('Automated vs manual measurements', fontsize=11, color=TEXT, x=0.01, ha='left')
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out_path, dpi=130, facecolor=SURFACE)
    plt.close(fig)
    return out_path
