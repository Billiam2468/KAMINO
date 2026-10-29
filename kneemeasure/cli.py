"""
Command line:

    python -m kneemeasure pull "Scan Pulls/Controls"          # one pull, its own layout
    python -m kneemeasure pulls "Scan Pulls"                  # every pull + status + comparison
    python -m kneemeasure app                                 # the app: pick folders, measure, review
    python -m kneemeasure review "Scan Pulls/Controls"        # open a pull straight into review
    python -m kneemeasure run SEGS... --laterality MAP.csv --out OUT [options]
    python -m kneemeasure check SEGS... --out OUT
    python -m kneemeasure compare NAME=RUN_DIR NAME=RUN_DIR ... --out OUT
    python -m kneemeasure validate RUN_DIR --manual SHEET.csv --id-column MRN --pair AUTO=MANUAL ... --out OUT
    python -m kneemeasure convert SEGS... --out OUT
    python -m kneemeasure stamp-laterality SEGS... --laterality MAP.csv (--out OUT | --in-place --backup DIR)

Run `python -m kneemeasure <command> -h` for each command's options.
"""
import argparse
import os
import shutil
import sys

import nrrd

from . import __version__
from .laterality import LateralityTable, resolve
from .pipeline import MEASURES
from .pulls import SOURCES, inspect_pull, run_pull, run_all
from .runner import run_cohort, compare_runs, validate_run
from .volume import find_segmentations, case_name_from_path, normalize_slicer_layers, \
    header_laterality, write_with_laterality


EXTRAS_HELP = ('also write statistics, QC flags, segmentation checks, landmarks, Slicer markups, a summary '
               'and change history (default: just measurements.csv and the QC grids)')


def _add_laterality_args(p):
    p.add_argument('--laterality', metavar='CSV',
                   help='CSV/TSV mapping case IDs to knee side (e.g. mrn_code_mapping.csv). '
                        "Optional if the NRRD headers carry a 'laterality' field.")
    p.add_argument('--id-column', action='append', dest='id_columns', metavar='COL',
                   help='table column(s) holding the case ID (default: code / folder_name / id / name ...)')
    p.add_argument('--side-column', metavar='COL', help='table column holding L/R (default: laterality / side)')
    p.add_argument('--laterality-priority', choices=['header', 'table'],
                   help='which source wins if header and table disagree (default: skip the case)')


def _exclusions(values):
    """--exclude takes case names and/or text files listing one per line."""
    out = []
    for v in values or []:
        if os.path.isfile(v):
            with open(v) as f:
                out += [l.strip() for l in f if l.strip() and not l.startswith('#')]
        else:
            out.append(v)
    return out


def cmd_run(a):
    keys = [k.strip() for k in a.measures.split(',')] if a.measures else None
    if keys:
        unknown = set(keys) - {m.key for m in MEASURES}
        if unknown:
            sys.exit(f"unknown measurement(s) {sorted(unknown)}; choose from {[m.key for m in MEASURES]}")
    run_cohort(a.segs, a.out, laterality=a.laterality, id_columns=a.id_columns, side_column=a.side_column,
               priority=a.laterality_priority, carry=a.carry or [], measure_keys=keys, jobs=a.jobs,
               qc_grids=not a.no_qc, case_sheets=a.case_sheets, slicer_export=not a.no_slicer,
               outlier_z=a.outlier_z, cohort_name=a.name or os.path.basename(os.path.abspath(a.out)),
               quiet=a.quiet, exclude=_exclusions(a.exclude), rerun=a.rerun, extras=a.extras)


def cmd_pull(a):
    spec = inspect_pull(a.pull, a.source)
    if not spec['eligible']:
        sys.exit(f"{spec['name']}: can't be measured -- {spec['reason']}")
    print(f"{spec['name']}: {len(spec['files'])} segmentation(s) "
          f"({', '.join(f'{v} {k}' for k, v in spec['counts'].items())}), {len(spec['exclude'])} excluded "
          f"-> {a.out or spec['out']}")
    if a.dry_run:
        return
    run_pull(a.pull, out=a.out, source=a.source, jobs=a.jobs, rerun=a.rerun, quiet=a.quiet,
             case_sheets=a.case_sheets, extras=a.extras)


def cmd_pulls(a):
    specs = run_all(a.root, report=a.report, jobs=a.jobs, rerun=a.rerun, only=a.only, quiet=a.quiet,
                    extras=a.extras)
    print(f"{len(specs)} pull(s) measured" + (f"; status and comparison in {a.report}" if a.report else ''))


def cmd_review(a):
    from .project import Project
    from .webapp import serve
    spec = inspect_pull(a.pull)
    if not spec['eligible']:
        sys.exit(f"{spec['name']}: can't be reviewed -- {spec['reason']}")
    files = {os.path.basename(f)[:-5]: f for f in spec['files']}
    for code in spec['exclude']:
        files.pop(code, None)
    project = Project(files=files, root=spec['path'], name=spec['name'], scans_dir=os.path.join(spec['path'], 'sorted'),
                      sheet=spec['laterality'])
    serve(project=project, port=a.port, reviewer=a.reviewer or '', open_browser=not a.no_browser,
          case_info={c: v for c, v in spec['case_info'].items() if c in files})


def cmd_app(a):
    from .webapp import serve
    serve(segs_dir=a.segmentations, scans_dir=a.scans, port=a.port, reviewer=a.reviewer or '',
          open_browser=not a.no_browser)


def cmd_check(a):
    run_cohort(a.segs, a.out, laterality=a.laterality, id_columns=a.id_columns, side_column=a.side_column,
               priority=a.laterality_priority, measure_keys=[], jobs=a.jobs, case_sheets='none',
               slicer_export=False, cohort_name=a.name or 'segmentation check', quiet=a.quiet)


def cmd_compare(a):
    groups = {}
    for g in a.groups:
        if '=' not in g:
            sys.exit(f'expected NAME=RUN_DIR, got {g!r}')
        name, path = g.split('=', 1)
        groups[name] = path
    table, _ = compare_runs(groups, a.out)
    print(f"wrote {os.path.join(a.out, 'comparison.csv')} ({len(table)} measurements, {len(groups)} groups)")


def cmd_validate(a):
    pairs = []
    for p in a.pair:
        if '=' not in p:
            sys.exit(f'expected AUTO_COLUMN=MANUAL_COLUMN, got {p!r}')
        pairs.append(tuple(s.strip() for s in p.split('=', 1)))
    results, per_case = validate_run(a.run, a.manual, a.id_column, pairs, a.out, side_column=a.side_column,
                                     exclude=_exclusions(a.exclude))
    n_matched = sum(1 for c in per_case if c['matched'] in ('id', 'id + side'))
    print(f"{n_matched}/{len(per_case)} cases matched to the manual sheet")
    for r in results:
        if r.get('icc_2_1') is not None:
            print(f"{r['automated_column']:28s} {r['subset']:13s} n={r['n']:3d}  bias={r['bias']:+7.2f}  "
                  f"ICC(2,1)={r['icc_2_1']:.2f}  r={r['pearson_r']:.2f}  {r['sign_check']}")
    print(f"wrote {os.path.join(a.out, 'agreement.csv')}")


def _out_path(src, out_dir):
    return os.path.join(out_dir, os.path.basename(src))


def cmd_convert(a):
    os.makedirs(a.out, exist_ok=True)
    for path in find_segmentations(a.segs):
        data, header = nrrd.read(path)
        new, new_header, notes = normalize_slicer_layers(data, header)
        dst = _out_path(path, a.out)
        if notes:
            nrrd.write(dst, new, new_header)
            print(f"{case_name_from_path(path)}: {'; '.join(notes)} -> {dst}")
        elif a.copy_unchanged:
            shutil.copy2(path, dst)
            print(f"{case_name_from_path(path)}: already standard, copied")
        else:
            print(f"{case_name_from_path(path)}: already standard, skipped")


def cmd_stamp(a):
    if not a.laterality:
        sys.exit('--laterality is required')
    if bool(a.out) == bool(a.in_place):
        sys.exit('give exactly one of --out DIR or --in-place')
    if a.in_place and not a.backup:
        sys.exit('--in-place requires --backup DIR (originals are copied there first)')
    table = LateralityTable(a.laterality, id_columns=a.id_columns, side_column=a.side_column)
    for d in (a.out, a.backup):
        if d:
            os.makedirs(d, exist_ok=True)
    for path in find_segmentations(a.segs):
        name = case_name_from_path(path)
        header_side = header_laterality(nrrd.read_header(path))
        side, source, note = resolve(name, None, table)
        if side is None:
            print(f"{name}: SKIPPED -- {note}")
            continue
        if header_side and header_side != side and not a.overwrite:
            print(f"{name}: SKIPPED -- header already says {header_side}, table says {side} (use --overwrite)")
            continue
        if a.in_place:
            shutil.copy2(path, _out_path(path, a.backup))
            dst = path
        else:
            dst = _out_path(path, a.out)
        write_with_laterality(path, dst, side)
        print(f"{name}: laterality:={side} -> {dst}")


def main(argv=None):
    ap = argparse.ArgumentParser(prog='kneemeasure', description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--version', action='version', version=f'kneemeasure {__version__}')
    sub = ap.add_subparsers(dest='command', required=True)

    p = sub.add_parser('pull', help='measure one scan pull laid out by Knee-Segmentation-Pipeline')
    p.add_argument('pull', help='the pull folder (contains sorted/mrn_code_mapping.csv)')
    p.add_argument('--source', choices=SOURCES,
                   help='segmentations to use (default from kneemeasure.toml, else best: cleaned, then Dixon)')
    p.add_argument('--out', help='results folder (default: <pull>/kneemeasure)')
    p.add_argument('--jobs', type=int, default=4)
    p.add_argument('--rerun', action='store_true', help='re-measure everything, ignoring cached results')
    p.add_argument('--case-sheets', choices=['flagged', 'all', 'none'], default='flagged',
                   help='per-case QC pages (with --extras)')
    p.add_argument('--extras', action='store_true', help=EXTRAS_HELP)
    p.add_argument('--dry-run', action='store_true', help='show what would be measured, then stop')
    p.add_argument('--quiet', action='store_true')
    p.set_defaults(func=cmd_pull)

    p = sub.add_parser('pulls', help='measure every pull under a folder and print a status table')
    p.add_argument('root', help='folder holding the pulls (e.g. "Scan Pulls")')
    p.add_argument('--report', metavar='DIR', help='also write status.csv/.md and a cross-pull comparison here')
    p.add_argument('--extras', action='store_true', help=EXTRAS_HELP)
    p.add_argument('--only', action='append', metavar='PULL', help='only this pull folder (repeatable)')
    p.add_argument('--jobs', type=int, default=4)
    p.add_argument('--rerun', action='store_true')
    p.add_argument('--quiet', action='store_true')
    p.set_defaults(func=cmd_pulls)

    p = sub.add_parser('app', help='the browser app: pick folders, measure, review the tuberosity')
    p.add_argument('segmentations', nargs='?', help='optional: open this segmentations folder right away')
    p.add_argument('--scans', help='optional scans folder (e.g. sorted/)')
    p.add_argument('--reviewer', help='name recorded with review decisions (default: login name)')
    p.add_argument('--port', type=int, default=8765, help='local port (default 8765; 0 = any free port)')
    p.add_argument('--no-browser', action='store_true', help="don't open the browser automatically")
    p.set_defaults(func=cmd_app)

    p = sub.add_parser('review', help='open a scan pull straight into the tuberosity review')
    p.add_argument('pull', help='the pull folder')
    p.add_argument('--reviewer', help='default reviewer initials')
    p.add_argument('--port', type=int, default=8765, help='local port (default 8765; 0 = any free port)')
    p.add_argument('--no-browser', action='store_true', help="don't open the browser automatically")
    p.set_defaults(func=cmd_review)

    p = sub.add_parser('run', help='measure a cohort: CSVs, statistics, QC grids, Slicer markups')
    p.add_argument('segs', nargs='+', help='segmentation .nrrd files and/or folders')
    p.add_argument('--out', required=True, help='output folder')
    _add_laterality_args(p)
    p.add_argument('--carry', action='append', metavar='COL',
                   help='laterality-table column to copy into measurements.csv (e.g. --carry mrn)')
    p.add_argument('--measures', help=f"comma-separated subset of: {','.join(m.key for m in MEASURES)}")
    p.add_argument('--jobs', type=int, default=1, help='parallel worker processes (default 1)')
    p.add_argument('--name', help='cohort name for titles and summary.md')
    p.add_argument('--no-qc', action='store_true', help='skip QC images')
    p.add_argument('--case-sheets', choices=['flagged', 'all', 'none'], default='flagged',
                   help='one-page QC sheets per case (default: only flagged cases)')
    p.add_argument('--no-slicer', action='store_true', help='skip the 3D Slicer markups export')
    p.add_argument('--outlier-z', type=float, default=3.5, help='robust z threshold for outlier flags')
    p.add_argument('--rerun', action='store_true', help='re-measure everything, ignoring cached results')
    p.add_argument('--extras', action='store_true', help=EXTRAS_HELP)
    p.add_argument('--exclude', action='append', metavar='CASE_OR_FILE',
                   help='case name to skip, or a text file of names (repeatable)')
    p.add_argument('--quiet', action='store_true')
    p.set_defaults(func=cmd_run)

    p = sub.add_parser('check', help='segmentation sanity checks + overview grid only (no measurements)')
    p.add_argument('segs', nargs='+')
    p.add_argument('--out', required=True)
    _add_laterality_args(p)
    p.add_argument('--jobs', type=int, default=1)
    p.add_argument('--name')
    p.add_argument('--quiet', action='store_true')
    p.set_defaults(func=cmd_check)

    p = sub.add_parser('compare', help='compare measurement runs between cohorts')
    p.add_argument('groups', nargs='+', metavar='NAME=RUN_DIR',
                   help='a run output folder (or its measurements.csv) per group')
    p.add_argument('--out', required=True)
    p.set_defaults(func=cmd_compare)

    p = sub.add_parser('validate', help='agreement of a run with a manual measurement sheet')
    p.add_argument('run', help='run output folder (or its measurements.csv)')
    p.add_argument('--manual', required=True, help='manual measurement sheet (CSV, one row per knee)')
    p.add_argument('--id-column', required=True, help='sheet column matching the case name (e.g. MRN)')
    p.add_argument('--side-column', help='sheet column with the knee side (needed for bilateral patients)')
    p.add_argument('--pair', action='append', required=True, metavar='AUTO=MANUAL',
                   help='automated column = manual column, e.g. --pair "KRA_deg=Knee Rotation Angle (MRI)"')
    p.add_argument('--exclude', action='append', metavar='CASE', help='case to leave out (repeatable)')
    p.add_argument('--out', required=True)
    p.set_defaults(func=cmd_validate)

    p = sub.add_parser('convert', help='write 3D standard-label copies of 3D Slicer layered (4D) exports')
    p.add_argument('segs', nargs='+')
    p.add_argument('--out', required=True)
    p.add_argument('--copy-unchanged', action='store_true', help='also copy files that need no conversion')
    p.set_defaults(func=cmd_convert)

    p = sub.add_parser('stamp-laterality', help="write a 'laterality' field into NRRD headers from a table")
    p.add_argument('segs', nargs='+')
    _add_laterality_args(p)
    p.add_argument('--out', help='write stamped copies here')
    p.add_argument('--in-place', action='store_true', help='modify the files themselves (needs --backup)')
    p.add_argument('--backup', help='folder to copy originals into before an in-place stamp')
    p.add_argument('--overwrite', action='store_true', help='replace an existing, different header laterality')
    p.set_defaults(func=cmd_stamp)

    a = ap.parse_args(argv)
    a.func(a)


if __name__ == '__main__':
    main()
