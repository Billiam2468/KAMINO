# Command line

For measuring whole scan pulls or other folders of segmentations without the app, and for
cohort statistics. The Python package and command are called `kneemeasure` (the project's
original name), as are the files it writes.

## Install

Python 3.11+:

```bash
pip install .            # or: pip install -r requirements.txt, then run from this folder
kneemeasure app          # the same app as the download
kneemeasure --help       # every command
```

## Scan pulls

The normal workflow is per scan pull. It expects the folder layout that
`Knee-Segmentation-Pipeline/run_pipeline.py` writes:
`sorted/mrn_code_mapping.csv`, `cleaned_segs/`, `dixon_segmentations/segmentations/`.

```bash
python -m kneemeasure pulls "/path/to/Scan Pulls" --jobs 8
python -m kneemeasure pull "/path/to/Scan Pulls/Controls" --jobs 8
```

The first command measures every pull. The second measures one; add `--dry-run` to see
what it would use first.

- **Segmentation source:** the cleaned segmentation if there is one, otherwise the Dixon
  model output. It's recorded per scan in the `segmentation_source` column; change it with
  `--source cleaned|dixon|mapss`.
- **Laterality:** `sorted/mrn_code_mapping.csv`.
- **Results:** `<pull>/kneemeasure/`.
- **Pulls that can't be measured** (no laterality sheet, no segmentations) are skipped,
  and `pulls` prints why. It prints one status line per pull; `--report DIR` also writes a
  status table and a cross-pull comparison to DIR.
- **Reruns are incremental.** A scan is re-measured only if its segmentation file's
  contents, knee side, review, options or the kneemeasure code changed; everything else comes
  from the hidden cache in `kneemeasure/.kneemeasure/`, so renaming files costs nothing. Use
  `--rerun` to force a full recompute.
- **Change tracking** (with `--extras`): `changes.csv` lists added and removed scans,
  changed segmentation files and values that moved by more than 0.01, and the previous tables
  are kept in `history/`.

Per-pull settings go in an optional `<pull>/kneemeasure.toml`:

```toml
name = "Controls"            # display name (default: the folder name)
source = "best"              # best | cleaned | dixon | mapss
[exclude]                    # scans to leave out, with the reason (logged in the run)
6YGMS = "tibia mostly outside the scan"
```

Five pulls (151 scans) take about 4 minutes the first time with 8 workers. A rerun with
nothing changed only redraws the QC images.

## Opening the app on a folder

- `python -m kneemeasure app <segmentations folder> --scans <scans folder>` opens a folder
  directly.
- `python -m kneemeasure review <pull>` opens a scan pull straight into review, using the
  pull's cleaned/model selection and its `sorted/` scans.

## Any folder of segmentations

```bash
python -m kneemeasure run "/path/to/segmentations" \
    --laterality "/path/to/mrn_code_mapping.csv" --carry mrn \
    --name "Controls" --out "/path/to/results/Controls" --jobs 8
```

A 30-case cohort takes about 40 seconds with 8 workers.

Compare cohorts, and check against manual measurements:

```bash
python -m kneemeasure compare Controls=results/Controls "Cartilage damage=results/CartDmg" --out results/comparison

python -m kneemeasure validate results/ORS --manual "MPFL+ Master.csv" --id-column MRN \
    --side-column "Laterality - use this unless says bilateral, then use laterality column 'R'" \
    --pair "KRA_deg=Knee Rotation Angle (MRI) (+ve = tibial ER relative to femur" \
    --pair "TT_TG_mm=TT-TG (MRI)" --out results/ORS_validation
```

## Commands

| command | what it does |
|---|---|
| `pull` | one scan pull, using its own layout and settings; results in `<pull>/kneemeasure/` |
| `pulls` | every pull under a folder, with a printed status line per pull (`--report DIR` for files) |
| `app` | the browser app: pick folders, set knee sides, measure, prescreen and review the tuberosity |
| `review` | open a scan pull straight into the app's tuberosity review |
| `run` | the full suite on a cohort, with all outputs below |
| `check` | segmentation checks and the overview grid only, no measurements (fast triage of a new batch) |
| `compare` | between-cohort statistics and figure from two or more `run` folders |
| `validate` | agreement of a `run` with a manual measurement sheet |
| `convert` | 3D standard-label copies of Slicer layered exports |
| `stamp-laterality` | write `laterality:=L/R` into NRRD headers from a sheet |

`run` options worth knowing:

- `--jobs N` runs N workers in parallel.
- `--measures kra,tt_tg` runs only some measurements.
- `--exclude CASE_OR_FILE` skips known-bad cases; give a name or a text file of names.
- `--case-sheets all|flagged|none` controls the per-case QC pages.
- `--carry COL` copies sheet columns (e.g. `mrn`) into the output.
- `--no-slicer` and `--no-qc` skip the Slicer markups and the QC images.
- `--rerun` ignores the cache.
- `--extras` writes the extra files listed below.

## Inputs

**Segmentations.** NRRD label maps (`.nrrd` / `.seg.nrrd`) using the standard labels:

| label | structure |
|---|---|
| 1 | femur bone |
| 2 | femur cartilage |
| 3 | patella bone |
| 4 | patella cartilage |
| 5 | tibia bone |

3D Slicer "layered" exports (4D files, produced when segments overlap) are converted
automatically on load. `Segment_k` maps to label k, or segments are mapped by anatomical
name if they have one. The conversion is noted in the case's warnings. `kneemeasure
convert` writes converted 3D copies if you want to fix the files themselves.

**Laterality.** Most measurements need the knee side: it decides medial vs lateral, and a
wrong side flips the sign of KRA. It comes from either:

- **A mapping sheet** (`--laterality`), a CSV/TSV with a case-ID column (`code`,
  `folder_name`, `id`, `name`, ... are detected automatically, or use `--id-column`) and a
  `laterality`/`side` column (L/R/Left/Right). A filename like `11U4C_20190219` also
  matches ID `11U4C`.
- **The NRRD header**, a `laterality:=L` field. `kneemeasure stamp-laterality` writes it
  from a sheet, either to copies (`--out`) or in place (`--in-place --backup DIR`). Once
  stamped, no sheet is needed.

If both sources exist and disagree, the case is not measured and is flagged, unless
`--laterality-priority header|table` is given. Cases with no side still get the
side-independent measurements (bump height, CDI) and the segmentation checks.

## Outputs

By default a run writes only:

| file | contents |
|---|---|
| `measurements.csv` | one row per scan: every value, knee side and where it came from, segmentation source, `TT_review` and the out-of-window flag, tibia length, and `notes` (why a measurement failed, segmentation warnings, cautions) |
| `qc/<measurement>_grid.png` | QC overlays for every scan (paginated at 30), plus a segmentation overview grid |
| `.kneemeasure/` | hidden: the rerun cache, `run_info.json` (inputs, parameters, versions) and `log.txt` |

With `--extras` it also writes the following. A later run without `--extras` removes them
again.

| file | contents |
|---|---|
| `summary.md` | success counts, statistics table, every QC flag |
| `statistics.csv` | n / mean / SD / median / IQR / range for all cases and by knee side; TT-dependent measurements repeated with likely out-of-window cases excluded |
| `qc_flags.csv` | failures, segmentation warnings, laterality problems, cautions, robust-z outliers |
| `segmentation_check.csv` | per-case label volumes, cartilage-on-bone fraction, tibia length, patella position |
| `landmarks.csv` | every landmark in the reoriented voxel frame and in LPS mm |
| `measurement_definitions.csv` | column, unit, definition |
| `plots/distributions.png` | every measurement's distribution, with outliers labeled |
| `qc/<measurement>_grid.png` | QC overlays for every case (paginated at 30), plus a segmentation overview grid |
| `qc/cases/<case>.png` | one page per flagged case with all of its panels and notes |
| `slicer/<case>.mrk.json` | landmarks as a 3D Slicer point list; drag it into Slicer with the segmentation |
| `changes.csv`, `history/` | what changed since the previous run |

Cases marked `[!]` in the QC grids have a segmentation warning or caution.
