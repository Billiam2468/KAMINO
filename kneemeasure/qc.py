"""
QC image assembly from the per-case payloads captured by pipeline.measure_case:
one paginated grid per measurement (plus a segmentation overview grid),
and optional one-page-per-case composites.
"""
import hashlib
import os
import pickle
import textwrap

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from . import validation
from .pipeline import MEASURES

# key -> (panel width in, panel height in, annotation fontsize)
_PANELS = {
    'overview': (3.2, 3.6, 7),
    'kra': (4.0, 4.2, 8),
    'tt_tg': (6.0, 6.3, 6),
    'ptltr': (6.0, 6.3, 9),
    'ttl': (4.0, 4.2, 9),
    'trochlear': (4.0, 4.2, 8),
    'bump_height': (4.0, 4.2, 8),
    'cdi': (4.0, 4.5, 7),
}
GRID_TITLES = {'overview': 'Segmentation overview (lateral projection)'}
GRID_TITLES.update({m.key: m.title for m in MEASURES})


def _renderer(key):
    if key == 'overview':
        return validation.render_overview, lambda p: ''
    mod = next(m.module for m in MEASURES if m.key == key)
    return mod.render_panel, mod.panel_title


def _case_label(result):
    side = result.get('knee_side') or '?'
    return f"{result.get('name') or os.path.basename(result['path'])} ({side})"


def _failure_text(result, key):
    if result.get('error'):
        return result['error']
    if key == 'overview':
        return 'no overview'
    note = result['notes'].get(key)
    return note or 'no QC image'


def _wrap(text, width=38):
    return '\n'.join(textwrap.wrap(str(text), width)) or str(text)


def grid_inputs(result, key):
    """The part of one case's result that its panel in the `key` grid is
    drawn from (all a worker process needs, and what a page's fingerprint
    covers)."""
    return {k: result.get(k) for k in ('name', 'path', 'knee_side', 'warnings', 'cautions', 'error')} | {
        'notes': {key: (result.get('notes') or {}).get(key)},
        'payloads': {key: (result.get('payloads') or {}).get(key)}}


def build_grid(results, key, out_dir, ncols=5, per_page=30, dpi=120):
    """Writes <out_dir>/<key>_grid[_pN].png and returns the paths. A page
    whose inputs (and this code) are unchanged since it was last written is
    left as it is: its fingerprint is kept next to it as .<name>.key."""
    from .cache import code_fingerprint
    render, title_fn = _renderer(key)
    w, h, fs = _PANELS[key]
    paths = []
    pages = [results[i:i + per_page] for i in range(0, len(results), per_page)] or [[]]
    for page_no, page in enumerate(pages, start=1):
        name = f'{key}_grid.png' if len(pages) == 1 else f'{key}_grid_p{page_no}.png'
        path, key_path = os.path.join(out_dir, name), os.path.join(out_dir, f'.{name}.key')
        fingerprint = hashlib.sha1(pickle.dumps(
            (code_fingerprint(), ncols, per_page, dpi, len(pages), [grid_inputs(r, key) for r in page]),
            protocol=4)).hexdigest()
        paths.append(path)
        if os.path.exists(path) and os.path.exists(key_path) and open(key_path).read() == fingerprint:
            continue
        n = max(len(page), 1)
        cols = min(ncols, n)
        rows = int(np.ceil(n / cols))
        fig, axes = plt.subplots(rows, cols, figsize=(cols * w, rows * h + 0.6), squeeze=False)
        for i, r in enumerate(page):
            ax = axes[i // cols, i % cols]
            payload = r.get('payloads', {}).get(key)
            label = _case_label(r)
            if payload is None:
                ax.text(0.5, 0.5, f"{label}\n\nFAILED\n{_wrap(_failure_text(r, key))}",
                        ha='center', va='center', fontsize=8, color='red', transform=ax.transAxes)
                ax.axis('off')
                continue
            try:
                render(ax, payload, fontsize=fs)
                sub = title_fn(payload)
                warn = ' [!]' if (r.get('warnings') or r.get('cautions')) else ''
                ax.set_title(f"{label}{warn}" + (f"\n{sub}" if sub else ''), fontsize=8)
            except Exception as e:
                ax.clear()
                ax.text(0.5, 0.5, f"{label}\nrender error: {e}", ha='center', va='center',
                        fontsize=8, color='red', transform=ax.transAxes)
                ax.axis('off')
        for i in range(len(page), rows * cols):
            axes[i // cols, i % cols].axis('off')
        suffix = f' (page {page_no}/{len(pages)})' if len(pages) > 1 else ''
        fig.suptitle(GRID_TITLES[key] + suffix, fontsize=12)
        fig.tight_layout(rect=(0, 0, 1, 1 - 0.5 / (rows * h + 0.6)))
        fig.savefig(path, dpi=dpi)
        plt.close(fig)
        with open(key_path, 'w') as f:
            f.write(fingerprint)
    return paths


def build_case_sheet(result, out_path, keys, dpi=110):
    """One page per case: every measurement's QC panel side by side, plus
    the case's notes and segmentation warnings underneath."""
    cols = 4
    rows = int(np.ceil(len(keys) / cols))
    lines = [f"{_case_label(result)}   {result['path']}"]
    lines += [f"{k}: {v}" for k, v in result['notes'].items() if v]
    lines += [f"WARNING: {w}" for w in result.get('warnings', [])]
    lines += [f"CAUTION: {w}" for w in result.get('cautions', [])]
    text = '\n'.join(_wrap(l, 170) for l in lines)
    text_h = 0.25 + 0.17 * (text.count('\n') + 1)
    fig = plt.figure(figsize=(cols * 4.5, rows * 4.5 + text_h))
    grid = fig.add_gridspec(rows + 1, cols, height_ratios=[4.5] * rows + [text_h])
    for i, key in enumerate(keys):
        ax = fig.add_subplot(grid[i // cols, i % cols])
        payload = result.get('payloads', {}).get(key)
        if payload is None:
            ax.text(0.5, 0.5, f"{GRID_TITLES[key]}\nFAILED\n{_wrap(_failure_text(result, key))}",
                    ha='center', va='center', fontsize=8, color='red', transform=ax.transAxes)
            ax.axis('off')
            continue
        render, title_fn = _renderer(key)
        render(ax, payload, fontsize=_PANELS[key][2])
        ax.set_title(f"{GRID_TITLES[key]}\n{title_fn(payload)}", fontsize=8)
    text_ax = fig.add_subplot(grid[rows, :])
    text_ax.axis('off')
    text_ax.text(0, 1, text, va='top', fontsize=9, family='monospace', transform=text_ax.transAxes)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)
