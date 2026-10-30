# Using the app: full guide

Every step, button and shortcut of the KAMINO app. The [README](../README.md) has the short
version.

Start it by double-clicking the downloaded app, or from the command line with
`kneemeasure app` (see [command-line.md](command-line.md)). It starts a small web app on
this computer and opens it in the browser. No data leaves the machine.

The app has five steps, shown in order along the top: **Setup**, then the three landmarks a
person should check (the tibial tuberosity, the tibial posterior condylar line (TPCL), and
the trochlear sulcus), then **Results**.

- **Each landmark step** starts with a **prescreen** and goes on to a **review** of the scans
  picked there. Its link in the step bar shows how many scans still have no decision, or a
  tick once they all do.
- **Moving on:** the prescreen's button opens the review, or goes straight to the next step
  when no scan is picked. In the review, **Go on to …** (under the decision buttons) leads to
  the next step, and turns blue once every scan in the list has a decision.

## 1. Setup & measure

- **Pick folders** with the folder browser (or type a path and press Enter). Picking the
  segmentations folder opens it; there's no separate Open step.
  - **Segmentations** (`.nrrd`): required.
  - **Scans:** optional. Each scan is shown under its segmentation; scans are matched by
    name, `<code>/water.nii.gz` (the pipeline's `sorted/` layout) or `<code>.nii.gz`.
  - **Suggestion:** inside a scan pull, picking the segmentations folder fills in the
    pull's `sorted/` folder automatically.
- **Knee side** for each scan comes from, first match wins:
  1. what was entered in the app;
  2. the segmentation's header (`laterality:=L`);
  3. the scan's NIfTI description (`laterality=L`);
  4. a `mrn_code_mapping.csv` next to the folders, if there is one.

  Scans with none of these get L/R buttons; each answer is saved, so it's asked once.
- **One button** shows the next step:
  1. **Set the knee side for N scans first** (disabled) while any side is missing; only
     those scans are listed.
  2. **Measure N scans** runs everything, with a progress bar, then goes on to the
     prescreen. Later, when scans are added or a segmentation changes, it measures just
     those.
  3. **Start review** once the results are up to date.

  The app writes nothing next to the scans except the decisions file (below). Its results
  are kept in your user cache folder (`~/.cache/kneemeasure/` on Linux,
  `~/Library/Caches/kneemeasure/` on a Mac, `%LOCALAPPDATA%\kneemeasure\` on Windows), so
  reopening a folder is instant. Download what you need from the Results page.
- **Start over** (on the folder's card, when there are decisions) removes every
  hand-checked landmark for the folder, so each scan goes back to its automatic landmarks.
  Knee sides are kept, and the decisions file is first copied to
  `kneemeasure_review.backup-<date>.json` next to it.

## 2. Prescreen tuberosities

A grid with one tile per scan.

- **Each tile:** the axial and sagittal slices through the automatic tuberosity point, over
  the MRI, with the scan's TT-TG and tibia length. A red dashed line marks the bottom of
  the scan when the tibia reaches it.
- **Selecting:** nothing starts out selected: click the tiles that need a closer look.
  "Select flagged" picks the scans flagged as maybe out of view or without an automatic
  point.
- **"Accept N · review M":** accepts the automatic point for every unselected scan
  (recorded as "accepted at prescreen") and opens the review with only the selected scans.
  "show all" in the review list brings the rest back. With no scan selected, the button
  accepts the rest and goes on to the next step. If an unselected scan has no automatic point, the page
  warns you first, since that scan would be left without a decision.
- **After a review:** a tile shows the decided point (green) instead of the automatic one,
  on the axial scan if it was marked there.
- **Loading:** tiles appear as the scans load in the background, a few seconds for a
  typical pull.

## 3. Tuberosity review

The tuberosity is often below the scan's field of view, and then the automatic point isn't
the tuberosity. TT-TG, sTT-TG and TTL depend on it, so every scan should be reviewed. For
each scan:

- **Views:** the axial and sagittal views show the MRI (when a scan was given) with the
  segmentation over it.
  - **Segmentation opacity:** 20% to start; set it with the slider, or press `S` to toggle it.
    Cartilage is drawn in darker colours at 1.75 times that opacity, since it's the
    brightest tissue on the Dixon.
  - **Zoom:** Ctrl+wheel or `+` / `−`; drag to pan. Each scan opens showing the whole knee.
  - **Side view:** the sagittal view shows the whole height of the scan, so you can see
    whether the tibia runs off the bottom.
- **Scrolling:** the mouse wheel over a view scrolls that view's slices, or use the arrow
  keys: up/down for axial, left/right for sagittal. White lines show where each view cuts
  the other.
- **Placing a point:** click in either view.
- **The axial scan** (`V`, or the Dixon / Axial buttons): when the scans folder is a
  pipeline `sorted/` folder, each scan also has `axial.nii.gz` (usually the PD axial),
  which reaches 20–30 mm further down the tibia than the Dixon.
  - **What you see:** the axial scan positioned on the segmentation, with slices continuing
    below the Dixon's bottom (negative z; the sagittal view marks where the Dixon stops). It
    is shown at its own resolution (0.44–0.66 mm pixels, finer than the Dixon's 0.9 mm),
    and a point clicked on it is kept to a fraction of a voxel. Its slices are 4 mm thick and
    about 5 mm apart, so the sagittal view of it is coarse from top to bottom, and each
    scroll step moves one axial slice. Scans flagged
    as out of view open on it, and scrolling down past the Dixon's bottom slice switches to
    it.
  - **Check the alignment first:** where both scans show the knee, the segmentation should
    sit on the bones of the axial scan. Both series come from the same session and share
    scanner coordinates, so it does unless the patient moved in between.
  - **When it lines up:** click the tuberosity and use your point. TT-TG, sTT-TG and TTL
    combine it with the groove and the posterior lines from the segmentation; they only
    need its position across the slice, so a point below the segmentation is fine.
- **Decide:**
  - **Accept** (`A`).
  - **Use my clicked point** (`Enter`).
  - **Tuberosity not in view** (`N`): TT-TG, sTT-TG and TTL are left blank instead of
    reported wrong.
  - **Axial doesn't line up** (`M`, on the axial scan): the tuberosity is visible there,
    but the segmentation doesn't match the axial scan, so TT-TG, sTT-TG and TTL are left
    blank, to be measured manually. This is rare: across the 151 scans in the pulls, the
    Dixon and axial lined up within 1 mm for 132 and within 2 mm for 150 (in-plane, 2026-10-29).
- **Order:** scans with less than 37 mm of tibia in the scan come first.
- **After each decision:** the scan's values are recalculated on the spot, and the app
  moves to the next unreviewed scan. Click a scan in the list, or press `J` / `K`, to move
  through the list yourself.

Everything a person enters (decisions and knee sides) is saved immediately to
`<root>/kneemeasure_review.json`, with the login name and time. "Undo decision" removes a
scan's decision.

- **Where it lives:** outside the results folder, so rebuilding never loses it.
- **How it's applied:** every later run (`kneemeasure pull`, or **Update results** on the
  app's Results page) uses it, and only re-measures the scans that changed.
- **In the output:** `measurements.csv` has a `TT_review` column (auto / accepted /
  replaced / replaced on axial / not in view / misaligned).

## 4. Tibial posterior line: prescreen and review

The tibial posterior condylar line (TPCL) sets the rotation reference for **KRA** and
**TTL**.

- **Prescreen:** a tile per scan shows the tibial slice the automatic line was taken from,
  with the line, KRA and TTL. "Select flagged" picks the scans with no automatic line, or
  a KRA beyond ±15°. Unselected scans get the automatic line accepted.
- **Review:** for each selected scan:
  - **Accept the automatic line** (`A`).
  - **Use the line on this slice** (`Enter`): scroll to a slice where both posterior
    condyles are clear. The bitangent across the notch on that slice is previewed in cyan.
  - **Use my two points** (`P`): click the medial and lateral condyle points yourself
    (`Esc` clears them), for when the slice method can't find both condyles.

  KRA and TTL are recalculated on each decision.
- **Storage:** decisions go in the same `kneemeasure_review.json`, under `tpcl`.
  `measurements.csv` has a `TPCL_review` column (auto / accepted / replaced).

## 5. Cartilage sulcus: prescreen and review

The cartilage sulcus and the two cartilage facet peaks on the trochlear reference slice
(the cartilage surface, not the bone under it). This is the first slice, from the top,
whose cartilage groove is at least 2 mm deep. They set the **sulcus
angle**, **LTI**, the **facet depths and lengths** and **trochlear depth**, and the sulcus
is the trochlear groove point for **TT-TG**.

- **Prescreen:** a tile per scan shows the front of the femur on the automatic slice, with
  the facet peaks and sulcus, the sulcus angle and TT-TG. "Select flagged" picks the scans
  with no automatic sulcus, or a sulcus angle outside 140–170°. Unselected scans get the
  automatic sulcus accepted.
- **Review:** for each selected scan:
  - **Accept the automatic sulcus** (`A`).
  - **Use the sulcus on this slice** (`Enter`): scroll to the slice you want. The sulcus
    and facet peaks found on it (same method, without the 2 mm rule) are previewed in
    cyan.
  - **Use my three points** (`P`): click the two facet peaks and the deepest point of the
    groove, in any order (`Esc` clears them). The middle point from side to side is the
    sulcus; which outer point is medial or lateral follows the knee side.

  The trochlear values and TT-TG are recalculated on each decision.
- **Storage:** decisions go in `kneemeasure_review.json` under `sulcus`.
  `measurements.csv` has a `Sulcus_review` column (auto / accepted / replaced).
- **Not covered:** bump height uses its own slice (the deepest groove), which this step
  doesn't change.

## 6. Results

- **Decisions:** how many scans have a decision for each landmark, and of which kind.
  Scans without one use the automatic landmark.
- **Update results:** re-runs the measurements with every decision (only the scans whose
  segmentation or decisions changed), then redoes the QC grids. The page says
  when decisions have changed since the results were made.
- **Download measurements.csv:** one row per scan with every measurement, the knee side and
  each landmark's decision (columns in [measurements.md](measurements.md)). Your browser saves it wherever you
  choose.
- **QC grids:** each measurement's grid, with the measurement drawn on every scan. Click a
  page to open it full size. **Download all grids (.zip)** saves every page as PNG.

## Closing, updates and the log

- **Closing it:** just close the browser tab. KAMINO stops by itself about 10
  minutes later (never while it's measuring); opening it again starts it again.
- **Updates:** when a newer version is out, a "Version ... is available" link appears at
  the top of the page. Download it the same way; your review decisions are kept in your
  data folders, and the app re-measures with the new version when you update the results.
- **Log:** messages and errors go to `.config/kneemeasure/app.log` in your home folder.
