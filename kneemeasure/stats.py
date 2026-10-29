"""
Statistics over a cohort's measurement table (no pandas dependency):

- describe: n / mean / SD / median / IQR / range per measurement, for the
  whole cohort and by knee side.
- flag_outliers: per-measurement robust z-scores (median / MAD); a value
  with |robust z| > threshold is flagged for QC review -- the cases that
  turned out to be bad segmentations in earlier cohorts all showed up this
  way (e.g. KRA -40.8 deg, TT-TG 0.05 mm).
- compare_groups: between-cohort comparison per measurement -- Welch t,
  Mann-Whitney U, Cohen's d / Hedges' g for two groups; one-way ANOVA and
  Kruskal-Wallis (+ pairwise Mann-Whitney) for more. p-values are also
  Holm-adjusted across the measurements in the table.
"""
import math

import numpy as np
from scipy import stats as sps

MAD_TO_SD = 1.4826


def _vals(rows, col):
    out = []
    for r in rows:
        v = r.get(col)
        if v is None or v == '':
            continue
        try:
            f = float(v)
        except (TypeError, ValueError):
            continue
        if math.isfinite(f):
            out.append(f)
    return np.array(out, dtype=float)


def _summary(x):
    n = len(x)
    if n == 0:
        return {'n': 0, 'mean': None, 'sd': None, 'median': None, 'q1': None, 'q3': None,
                'iqr': None, 'min': None, 'max': None}
    q1, med, q3 = np.percentile(x, [25, 50, 75])
    return {'n': n, 'mean': float(x.mean()), 'sd': float(x.std(ddof=1)) if n > 1 else None,
            'median': float(med), 'q1': float(q1), 'q3': float(q3), 'iqr': float(q3 - q1),
            'min': float(x.min()), 'max': float(x.max())}


def describe(rows, columns, units=None, split_by='knee_side'):
    """One row per (measurement, subgroup): subgroup 'all' plus each value
    of split_by (e.g. L / R)."""
    units = units or {}
    groups = [('all', rows)]
    if split_by:
        for g in sorted({r.get(split_by) for r in rows if r.get(split_by)}):
            groups.append((f'{split_by}={g}', [r for r in rows if r.get(split_by) == g]))
    out = []
    for col in columns:
        for gname, grows in groups:
            x = _vals(grows, col)
            out.append({'measurement': col, 'unit': units.get(col, ''), 'subgroup': gname,
                        'n_cases': len(grows), 'n_missing': len(grows) - len(x), **_summary(x)})
    return out


def flag_outliers(rows, columns, threshold=3.5, min_n=8, id_key='case'):
    """[{case, measurement, value, median, robust_z}] for values far from
    the cohort median. Skipped for measurements with < min_n values."""
    flags = []
    for col in columns:
        x = _vals(rows, col)
        if len(x) < min_n:
            continue
        med = float(np.median(x))
        mad = float(np.median(np.abs(x - med))) * MAD_TO_SD
        if mad == 0:
            q1, q3 = np.percentile(x, [25, 75])
            mad = float(q3 - q1) / 1.349
        if mad == 0:
            continue
        for r in rows:
            v = r.get(col)
            if v is None or v == '':
                continue
            z = (float(v) - med) / mad
            if abs(z) > threshold:
                flags.append({id_key: r.get(id_key), 'measurement': col, 'value': float(v),
                              'cohort_median': med, 'robust_z': z})
    return flags


def holm(pvals):
    """Holm-Bonferroni adjusted p-values (None entries passed through)."""
    idx = [i for i, p in enumerate(pvals) if p is not None and not math.isnan(p)]
    order = sorted(idx, key=lambda i: pvals[i])
    m = len(order)
    adj = [None] * len(pvals)
    running = 0.0
    for rank, i in enumerate(order):
        running = max(running, min(1.0, (m - rank) * pvals[i]))
        adj[i] = running
    return adj


def _cohens_d(a, b):
    na, nb = len(a), len(b)
    if na < 2 or nb < 2:
        return None, None
    sp = math.sqrt(((na - 1) * a.var(ddof=1) + (nb - 1) * b.var(ddof=1)) / (na + nb - 2))
    if sp == 0:
        return None, None
    d = (a.mean() - b.mean()) / sp
    g = d * (1 - 3 / (4 * (na + nb) - 9))
    return float(d), float(g)


def _safe(fn, *args, **kw):
    try:
        return fn(*args, **kw)
    except Exception:
        return None


def compare_groups(groups, columns, units=None, min_n=3):
    """groups: {group_name: rows}. Returns (per_measurement_rows,
    pairwise_rows). With exactly two groups the per-measurement row holds
    the two-sample tests; with more, the omnibus tests, and pairwise_rows
    holds every pair."""
    units = units or {}
    names = list(groups)
    table, pairwise = [], []
    for col in columns:
        data = {g: _vals(groups[g], col) for g in names}
        row = {'measurement': col, 'unit': units.get(col, '')}
        for g in names:
            s = _summary(data[g])
            row[f'{g}_n'] = s['n']
            row[f'{g}_mean'] = s['mean']
            row[f'{g}_sd'] = s['sd']
            row[f'{g}_median'] = s['median']
            row[f'{g}_iqr'] = s['iqr']
        usable = [g for g in names if len(data[g]) >= min_n]
        if len(names) == 2 and len(usable) == 2:
            a, b = data[names[0]], data[names[1]]
            t = _safe(sps.ttest_ind, a, b, equal_var=False)
            u = _safe(sps.mannwhitneyu, a, b, alternative='two-sided')
            d, g = _cohens_d(a, b)
            row.update({'mean_difference': float(a.mean() - b.mean()),
                        'welch_t_p': None if t is None else float(t.pvalue),
                        'mann_whitney_p': None if u is None else float(u.pvalue),
                        'cohens_d': d, 'hedges_g': g})
        elif len(usable) >= 3:
            f = _safe(sps.f_oneway, *[data[g] for g in usable])
            k = _safe(sps.kruskal, *[data[g] for g in usable])
            row.update({'anova_p': None if f is None else float(f.pvalue),
                        'kruskal_p': None if k is None else float(k.pvalue)})
        table.append(row)

        if len(names) > 2:
            for i in range(len(names)):
                for j in range(i + 1, len(names)):
                    a, b = data[names[i]], data[names[j]]
                    if len(a) < min_n or len(b) < min_n:
                        continue
                    u = _safe(sps.mannwhitneyu, a, b, alternative='two-sided')
                    t = _safe(sps.ttest_ind, a, b, equal_var=False)
                    d, g = _cohens_d(a, b)
                    pairwise.append({'measurement': col, 'group_a': names[i], 'group_b': names[j],
                                     'mean_a': float(a.mean()), 'mean_b': float(b.mean()),
                                     'welch_t_p': None if t is None else float(t.pvalue),
                                     'mann_whitney_p': None if u is None else float(u.pvalue),
                                     'cohens_d': d, 'hedges_g': g})

    for key in ('welch_t_p', 'mann_whitney_p', 'anova_p', 'kruskal_p'):
        if any(key in r for r in table):
            adj = holm([r.get(key) for r in table])
            for r, a in zip(table, adj):
                if key in r:
                    r[f'{key}_holm'] = a
    if pairwise:
        for key in ('welch_t_p', 'mann_whitney_p'):
            for col in columns:
                rows_c = [r for r in pairwise if r['measurement'] == col]
                for r, a in zip(rows_c, holm([r[key] for r in rows_c])):
                    r[f'{key}_holm_pairs'] = a
    return table, pairwise
