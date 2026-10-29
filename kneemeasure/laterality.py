"""
Knee laterality (L/R) for each case. Most measurements need it: medial vs
lateral is decided from the knee side, and a wrong side silently flips the
sign of KRA and swaps medial/lateral landmarks.

Sources, per case:
  1. a 'laterality' field in the segmentation's NRRD header
     (laterality:=L), written by `kneemeasure stamp-laterality`;
  2. a mapping table (CSV/TSV) with a case-ID column and a laterality
     column -- e.g. the mrn_code_mapping.csv sheets in each cohort folder.

If both exist and disagree the case is NOT measured (unless a priority is
given), since the earlier cohorts showed mapping sheets do contain wrong
sides.
"""
import csv

from .volume import parse_side

ID_COLUMN_CANDIDATES = ('code', 'case', 'case_id', 'id', 'name', 'filename', 'file', 'folder_name')
SIDE_COLUMN_CANDIDATES = ('laterality', 'side', 'knee_side', 'knee', 'leg')


class LateralityTable:
    """A case-ID -> side mapping loaded from a CSV/TSV file."""

    def __init__(self, path, id_columns=None, side_column=None):
        self.path = path
        with open(path, newline='', encoding='utf-8-sig') as f:
            sample = f.read(4096)
            f.seek(0)
            try:
                dialect = csv.Sniffer().sniff(sample, delimiters=',\t;')
            except csv.Error:
                dialect = csv.excel
            reader = csv.DictReader(f, dialect=dialect)
            self.fieldnames = reader.fieldnames or []
            self.rows = [r for r in reader if any((v or '').strip() for v in r.values())]

        lower = {c.lower().strip(): c for c in self.fieldnames}
        if side_column:
            if side_column not in self.fieldnames:
                raise ValueError(f'{path}: no column {side_column!r} (columns: {self.fieldnames})')
            self.side_column = side_column
        else:
            self.side_column = next((lower[c] for c in SIDE_COLUMN_CANDIDATES if c in lower), None)
            if self.side_column is None:
                raise ValueError(f'{path}: no laterality column found (looked for {SIDE_COLUMN_CANDIDATES})')
        if id_columns:
            missing = [c for c in id_columns if c not in self.fieldnames]
            if missing:
                raise ValueError(f'{path}: no column(s) {missing} (columns: {self.fieldnames})')
            self.id_columns = list(id_columns)
        else:
            self.id_columns = [lower[c] for c in ID_COLUMN_CANDIDATES if c in lower]
            if not self.id_columns:
                self.id_columns = [self.fieldnames[0]]

        self._index = {}
        for row in self.rows:
            for col in self.id_columns:
                key = (row.get(col) or '').strip()
                if key:
                    self._index.setdefault(key.lower(), []).append(row)

    def lookup(self, case_name):
        """(side or None, row or None, how). Tries an exact (case-
        insensitive) ID match, then the name with a trailing _suffix
        (e.g. _20240523 scan dates) removed."""
        for key, how in ((case_name, 'id'), (case_name.rsplit('_', 1)[0], 'id prefix')):
            rows = self._index.get(key.strip().lower())
            if not rows:
                continue
            unique = {id(r): r for r in rows}.values()
            sides = {parse_side(r.get(self.side_column)) for r in unique}
            if len(unique) > 1 and len(sides) > 1:
                return None, None, f'ambiguous: {len(unique)} table rows with different sides'
            row = next(iter(unique))
            side = parse_side(row.get(self.side_column))
            if side is None:
                return None, row, f"unreadable laterality value {row.get(self.side_column)!r}"
            return side, row, how
        return None, None, 'not in table'


def resolve(case_name, header_side, table, priority=None):
    """(side, source, note). priority: None (refuse on conflict),
    'header' or 'table'."""
    table_side, row, how = (None, None, 'no table')
    if table is not None:
        table_side, row, how = table.lookup(case_name)

    if header_side and table_side and header_side != table_side:
        if priority == 'header':
            return header_side, 'header', f'conflict: table says {table_side}; header used'
        if priority == 'table':
            return table_side, 'table', f'conflict: header says {header_side}; table used'
        return None, 'conflict', f'laterality conflict: header {header_side} vs table {table_side}'
    if header_side:
        return header_side, 'header', None
    if table_side:
        return table_side, 'table' if how == 'id' else f'table ({how})', None
    return None, 'missing', f'no laterality ({how})'


def table_row(table, case_name):
    if table is None:
        return None
    return table.lookup(case_name)[1]
