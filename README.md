<h1 align="center">KAMINO</h1>

<p align="center"><b>K</b>nee <b>A</b>natomy <b>M</b>easurements from <b>I</b>maging and <b>N</b>etwork <b>O</b>utput</p>

<div align="center">
<table align="center">
  <tr align="center"><th>K</th><th>A</th><th>M</th><th>I</th><th>N</th><th>O</th></tr>
  <tr align="center"><td>Knee</td><td>Anatomy</td><td>Measurements</td><td>Imaging</td><td>Network</td><td>Output</td></tr>
</table>
</div>

Automated patellofemoral and tibiofemoral measurements from knee MRI segmentations, with
QC images of every measurement and an app for checking the landmarks a person should see.
Everything runs on your own computer; no scan or result leaves it.

![Every KAMINO measurement drawn on one knee's segmentation](docs/images/measurements.png)

*One knee's segmentation with each measurement drawn on it, as in KAMINO's QC grids
(green: femur, blue: femoral cartilage, purple: patella, red: patellar cartilage, orange: tibia).*

## What it measures

Each measurement is a column in `measurements.csv`, one row per scan.

| Measurement | Column | What it is |
|---|---|---|
| **KRA**, knee rotation angle | `KRA_deg` | Tibial vs femoral posterior condylar line; + = tibial external rotation |
| **TT-TG** | `TT_TG_mm` | Tibial tuberosity to trochlear groove, along the femoral posterior condylar line |
| **sTT-TG**, sagittal TT-TG | `sTT_TG_mm` | The same, front to back; − = tuberosity behind the groove |
| **PTLTR** | `PTLTR_mm` | Patella's lateral edge to the lateral trochlear ridge |
| **TTL**, tibial tubercle lateralization | `TTL` | Tuberosity's position across the tibia's width, from the medial edge (a fraction) |
| **LTI**, lateral trochlear inclination | `LTI_deg` | Lateral facet vs the femoral posterior condylar line |
| **Sulcus angle** | `sulcus_angle_deg` | Angle between the trochlear cartilage facets |
| **Trochlear depth** | `trochlear_depth_mm` | Groove depth below the facets; plus each facet's depth |
| **Facet lengths** | `medial_/lateral_facet_length_mm` | Groove to each facet peak; plus `facet_asymmetry_ratio` |
| **Bump height** | `bone_bump_mm`, `cartilage_bump_mm` | Trochlear bump above the anterior femoral cortex |
| **CDI** | `CDI` | Caton-Deschamps index |

The file also has, per scan:

- **Knee side.**
- **Review columns:** whether each landmark was checked in the app (`TT_review`,
  `TPCL_review`, `Sulcus_review`).
- **`notes`:** why a measurement failed, and segmentation warnings.

Exact definitions, validation against manual measurements and limitations:
[docs/measurements.md](docs/measurements.md) and [VALIDATION.md](VALIDATION.md).

## Download and run

1. **Download** the zip for your computer from the
   [latest release](https://github.com/Billiam2468/KAMINO/releases/latest):
   `KAMINO-Windows.zip`, or `KAMINO-Mac.zip` (Apple Silicon, M1 or later).
2. **Unzip it.** On Windows, right-click it and choose **Extract All**; KAMINO doesn't run
   from inside the zip.
3. **Double-click KAMINO** (`KAMINO.exe` on Windows). Your browser opens on the setup page.

The app isn't signed with a paid certificate, so the first time:

- **Windows:** at "Windows protected your PC", click **More info**, then **Run anyway**.
- **Mac:** if it can't be opened, go to **System Settings > Privacy & Security** and click
  **Open Anyway**. If it says the app is "damaged", run
  `xattr -dr com.apple.quarantine /Applications/KAMINO.app` in Terminal.

To quit, close the browser tab; KAMINO stops by itself about 10 minutes later. When a new
version is out, a link appears at the top of the page.

### What you need

- **Segmentations** (required): one `.nrrd` label map per knee, labelled 1 femur,
  2 femoral cartilage, 3 patella, 4 patellar cartilage, 5 tibia.
- **Scans** (optional, recommended): one subfolder per segmentation, with the same name,
  holding `water.nii.gz` (the Dixon) and `axial.nii.gz` (the PD axial). This is the
  `sorted/` folder the Knee-Segmentation-Pipeline writes.

## Using the app

The steps run along the top of the page. Each landmark step starts with a **prescreen**,
a grid with one tile per scan. Click the tiles that look wrong; the rest get the
automatic landmark. Then **review** the ones you picked.

### 1. Setup

- Pick your segmentations folder, and the scans folder if you have one. Inside a scan
  pull, the scans folder is filled in for you.
- The knee side is read from the files when the pipeline wrote it. Otherwise, click L or R
  once per scan.
- Click **Measure**. Later, only new or changed scans are measured again.

### 2. Tuberosity

TT-TG, sTT-TG and TTL depend on this point. The tuberosity is often below the Dixon's
field of view, so check every scan.

- Scroll with the mouse wheel or the arrow keys. Click either view to place a point.
- **Axial scan** (`V`): the PD axial reaches 20–30 mm further down the tibia. First check
  the segmentation sits on its bones, then mark the tuberosity there.
- Decide:
  - **Accept** (`A`) or **Use my clicked point** (`Enter`).
  - **Tuberosity not in view** (`N`) or **Axial doesn't line up** (`M`): TT-TG, sTT-TG
    and TTL are left blank, to be measured manually.

### 3. Tibial line

The tibial posterior condylar line, used for KRA and TTL. Choose one:

- **Accept** (`A`).
- **Use the line on this slice** (`Enter`): scroll to a slice where both condyles are clear.
- **Click the two condyle points** (`P`).

### 4. Cartilage sulcus

The groove and the two facet peaks on the trochlear reference slice, used for the sulcus
angle, LTI, trochlear depth, facet lengths and TT-TG. Choose one:

- **Accept** (`A`).
- **Use the sulcus on this slice** (`Enter`).
- **Click the two peaks and the groove** (`P`).

### 5. Results

- **Update results** after changing decisions.
- Download `measurements.csv` and the **QC grids**: each measurement drawn on every scan.
  Look through them before trusting a new cohort's numbers.

Your decisions and knee sides are saved as you go, in `kneemeasure_review.json` in the
scan pull folder (or the folder above your segmentations folder). **Start over** clears the
decisions, after backing the file up.

Every button and shortcut: [docs/app.md](docs/app.md).

## More

- [docs/app.md](docs/app.md): the full app guide.
- [docs/measurements.md](docs/measurements.md): definitions, landmark rules, validation
  and limitations.
- [docs/command-line.md](docs/command-line.md): measuring scan pulls or any folder from the
  command line, and every output file.
- [docs/development.md](docs/development.md): the code, the public-repo data check, and
  releasing a version.
