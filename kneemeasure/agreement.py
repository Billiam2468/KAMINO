"""
Automated-vs-manual agreement for a measurement run, against a manual
measurement sheet (one row per knee).

For each (automated column, manual column) pair:
  n, means, bias (auto - manual), SD of differences, 95% limits of
  agreement, mean absolute error, Pearson r and Spearman rho, ICC(2,1)
  (two-way random, absolute agreement, single measure) and ICC(3,1)
  (two-way mixed, consistency) -- the Shrout & Fleiss forms -- and a
  sign check: for a signed measurement, whether auto correlates better
  with -manual than with +manual (a flipped sign convention shows up as a
  strongly negative r).

Matching: a case is matched to a sheet row by ID (exact, or with a trailing
L/R knee suffix removed) and, when the sheet has several rows per ID
(bilateral patients), by knee side.
"""
import csv
import math
import re

import numpy as np
from scipy import stats as sps

from .volume import parse_side


def _num(v):
    if v is None:
        return None
    try:
        f = float(str(v).strip())
    except ValueError:
        return None
    return f if math.isfinite(f) else None


def icc(x, y):
    """(ICC(2,1), ICC(3,1)) for two raters x, y over n subjects."""
    data = np.column_stack([x, y]).astype(float)
    n, k = data.shape
    if n < 3:
        return None, None
    grand = data.mean()
    ss_rows = k * ((data.mean(axis=1) - grand) ** 2).sum()
    ss_cols = n * ((data.mean(axis=0) - grand) ** 2).sum()
    ss_total = ((data - grand) ** 2).sum()
    ss_err = ss_total - ss_rows - ss_cols
    msr = ss_rows / (n - 1)
    msc = ss_cols / (k - 1)
    mse = ss_err / ((n - 1) * (k - 1))
    denom2 = msr + (k - 1) * mse + k * (msc - mse) / n
    denom3 = msr + (k - 1) * mse
    icc2 = (msr - mse) / denom2 if denom2 else None
    icc3 = (msr - mse) / denom3 if denom3 else None
    return (None if icc2 is None else float(icc2)), (None if icc3 is None else float(icc3))


def load_manual(path, id_column, side_column=None):
    with open(path, newline='', encoding='utf-8-sig', errors='replace') as f:
        rows = list(csv.DictReader(f))
    if rows and id_column not in rows[0]:
        raise ValueError(f'{path}: no column {id_column!r}')
    if side_column and rows and side_column not in rows[0]:
        raise ValueError(f'{path}: no column {side_column!r}')
    index = {}
    for r in rows:
        key = (r.get(id_column) or '').strip().lower()
        if key:
            index.setdefault(key, []).append(r)
    return index


def match(index, case, side, side_column):
    """(row or None, how)."""
    candidates = [case]
    m = re.fullmatch(r'(.+?)[_-]?([LR])', case, flags=re.I)
    if m:
        candidates.append(m.group(1))
    for key in candidates:
        rows = index.get(key.strip().lower())
        if not rows:
            continue
        if len(rows) == 1:
            r = rows[0]
            if side_column and side and parse_side(r.get(side_column, '')) not in (None, side):
                return None, f'sheet row is for the {parse_side(r.get(side_column))} knee'
            return r, 'id'
        if not side_column:
            return None, f'{len(rows)} sheet rows for this ID and no side column to choose'
        same = [r for r in rows if parse_side((r.get(side_column) or '').strip(' ,')) == side]
        if len(same) == 1:
            return same[0], 'id + side'
        return None, f'{len(rows)} sheet rows for this ID, {len(same)} for side {side}'
    return None, 'not in sheet'


def agreement(pairs_xy):
    x = np.array([p[0] for p in pairs_xy], dtype=float)
    y = np.array([p[1] for p in pairs_xy], dtype=float)
    n = len(x)
    out = {'n': n}
    if n < 3:
        return out
    d = x - y
    sd = float(d.std(ddof=1))
    r = sps.pearsonr(x, y)
    rho = sps.spearmanr(x, y)
    icc2, icc3 = icc(x, y)
    out.update({
        'auto_mean': float(x.mean()), 'manual_mean': float(y.mean()),
        'bias': float(d.mean()), 'sd_diff': sd,
        'loa_low': float(d.mean() - 1.96 * sd), 'loa_high': float(d.mean() + 1.96 * sd),
        'mae': float(np.abs(d).mean()),
        'pearson_r': float(r.statistic), 'pearson_p': float(r.pvalue),
        'spearman_rho': float(rho.statistic),
        'icc_2_1': icc2, 'icc_3_1': icc3,
    })
    neg_r = float(sps.pearsonr(x, -y).statistic)
    flipped_icc, _ = icc(x, -y)
    out['icc_2_1_if_sign_flipped'] = flipped_icc
    out['sign_check'] = ('auto agrees better with -manual: sign convention may be flipped'
                         if (flipped_icc or -1) > (icc2 or -1) + 0.1 and neg_r > 0.3 else '')
    return out
