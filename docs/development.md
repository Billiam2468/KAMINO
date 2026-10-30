# Working on the code

## Public repository

This repository is **public**: never commit scans, segmentations, results, mapping
sheets, review files or MRNs. After cloning, turn on the check that stops such commits:

```bash
git config core.hooksPath .githooks
```

`.githooks/pre-commit` refuses a commit with data-type files (`.nrrd`, `.nii.gz`, `.csv`,
...) or numbers shaped like MRNs; `git commit --no-verify` overrides a false alarm.

## Code map

```
kneemeasure/
  volume.py              I/O, LPS reorientation, true-physical transforms, Slicer layer conversion
  labels.py              label constants
  trochlea.py            cartilage facet/sulcus detection, reference slice, deepest groove
  pbcl.py, tibial_pcl.py posterior condylar lines (femur, tibia)
  tibial_tuberosity.py   tuberosity point
  landmarks.py           per-case landmark cache shared by all measurements
  kra.py, tt_tg.py, ptltr.py, ttl.py, trochlear_morphology.py, bump_height.py, cdi.py
                         one module per measurement: compute / qc_payload / render_panel / landmark_points
  pipeline.py            measure_case(): one file in, plain-dict result out (no volumes kept)
  validation.py          segmentation checks + overview panel
  runner.py              cohort run, compare, validate; all file outputs
  pulls.py               scan-pull layout: source selection, kneemeasure.toml, run-all + status
  cache.py               per-case result cache keyed on file content + code fingerprint
  review.py              landmark review decisions (kneemeasure_review.json)
  project.py             app projects: folders, matching scans, knee-side rules, folder browsing
  scans.py               NIfTI reader + display alignment of scans to segmentations
  webapp.py, webapp/     the local web app (standard-library server + one HTML page)
  desktop.py             what the downloadable app runs: opens the browser, quits when unused
  stats.py, agreement.py descriptive stats, outlier flags, group tests, ICC / Bland-Altman
  qc.py, qc_draw.py, plots.py   QC grids, case sheets, statistics figures
  slicer.py              3D Slicer markups export
  laterality.py          mapping sheets, header field, conflict rules
  cli.py                 command line
packaging/               PyInstaller recipe for the downloadable app
.github/workflows/       release.yml: builds the Windows and Mac apps for each version tag
```

**Adding a measurement:**

1. Write a module with `COLUMNS`, `compute(case, lm)`, `qc_payload`, `render_panel`,
   `panel_title` and `landmark_points`, following `kra.py`.
2. Add it to `MEASURES` in `pipeline.py`.
3. Add a panel size in `qc.py`.

The CSVs, statistics, grids and Slicer export pick it up automatically.

**Toward Slicer / a GUI.** Every measurement is a pure function of the loaded volume.
Landmarks are exported in LPS mm as Slicer point lists. A Slicer module could therefore
call `measure_case`, show the markups, and let a user drag a point. Re-measuring from
edited landmarks would need each `compute` to accept landmark overrides, which isn't built
yet.

## Releasing a new version

1. Set the new version in `kneemeasure/__init__.py` (`__version__ = '0.2.0'`) and commit.
2. Tag it and push the tag:

   ```bash
   git tag v0.2.0
   git push origin main v0.2.0
   ```
3. GitHub Actions (`.github/workflows/release.yml`) runs the tests, builds the Windows and
   Mac apps, and attaches `KAMINO-Windows.zip` and `KAMINO-Mac.zip` to a release
   for that tag (about 10 minutes; see the Actions tab). Add release notes on the release
   page. Users' apps show the update link from then on.

To build without releasing, use **Run workflow** on the Actions tab and download the zips
from the run, or build locally on the platform you're on:

```bash
pip install . pyinstaller
pyinstaller --noconfirm packaging/KAMINO.spec    # -> dist/
```
