# Measurements: methods, validation and limitations

How each value in `measurements.csv` is calculated is in `Measurement methods.xlsx`; this
page has the column definitions, the landmark rules and how far the values can be trusted.

All are computed in 3D from the segmentation. Distances that end at a structure's edge run
edge to edge, as measured by hand: to the outer edge of its outermost voxel, not the voxel's
centre (TTL's T and t, PTLTR's patella end, CDI's cartilage length and numerator, both
patellar widths).

- **Axial measurements** use the reoriented voxel frame scaled to mm. Relative
  angles and distances are unaffected by any in-plane rotation in the scan header.
- **The tibial tuberosity** search runs in true physical coordinates, because "most
  anterior" is an absolute direction.


| column | definition |
|---|---|
| `KRA_deg` | Knee rotation angle: tibial posterior condylar line (TPCL) vs femoral posterior bicondylar line (PBCL). **+ = tibial external rotation** relative to the femur. |
| `TT_TG_mm` | TT-TG along the PBCL; the groove is the sulcus of the Dejour cartilage reference slice |
| `sTT_TG_mm` | sagittal TT-TG, perpendicular to the PBCL; **negative = TT posterior to TG** |
| `PTLTR_mm` | patella's lateral border (the lateral end of the PW_max line, on its widest slice) to the lateral trochlear bony ridge on the tendon cut (the first slice below the patella, or the nearest slice with a usable ridge), along the posterior condylar line on that slice (unsigned; Mistovich 2018, with the patella's edge for the tendon's) |
| `TTL`, `TTL_t_mm`, `TTL_T_mm` | tubercle lateralization t/T along the TPCL (t = medial tibial edge to TT, T = tibial width) |
| `LTI_deg` | lateral trochlear inclination at the Dejour reference slice (+ = lateral facet anterior to PBCL) |
| `sulcus_angle_deg` | cartilage sulcus angle at the reference slice |
| `medial_/lateral_facet_depth_mm` | facet height above the sulcus, perpendicular to the PBCL |
| `trochlear_depth_mm` | sulcus depth below the facet chord, perpendicular to the PBCL |
| `medial_/lateral_facet_length_mm`, `facet_asymmetry_ratio` | sulcus to facet-peak distances; asymmetry = medial / lateral |
| `bone_bump_mm`, `cartilage_bump_mm` | trochlear bump height above the anterior femoral cortex line, on the sagittal slice through the deepest groove |
| `CDI`, `CDI_numerator_mm`, `patella_cartilage_length_mm` | Caton-Deschamps index on the sagittal slice with the longest patellar cartilage |
| `PW_effective_mm` | effective patellar width (Hevesi 2023): widest medial-lateral extent of the patellar cartilage, parallel to the PBCL ([below](#patellar-width-and-tilt)) |
| `PW_max_mm` | maximum patellar width (Chen 2024 / Staubli): bony patella, tip to tip along its own axis |
| `patellar_tilt_deg` | lateral patellar tilt (Chen 2024): PW_max line vs PBCL; **+ = lateral border posterior** (the lines converge laterally) |
| `PW_effective_slice`, `PW_max_slice` | the axial slices (z, as the app numbers them) the widths were measured on; information only, not in the statistics |
| `PW_effective_offset_mm`, `PW_max_offset_mm` | those slices' distance from the PBCL slice, + = proximal; information only |

## Landmark rules

- **Dejour reference slice:** the first slice, cranial to caudal, with a cartilage sulcus
  of at least 2 mm.
- **PBCL:** the bitangent of the posterior femoral condyles, on the slice where it lies
  most posterior.
- **TPCL:** one slice, the first below the tibial plateau where the posterior notch is clearly
  present between two convex condylar lobes; the line is the bitangent of the posterior
  contour across the notch.
- **Tuberosity:** a fixed shaft plane through the tibia, then the most anterior and
  lateral point.

Each module's docstring gives the full method and the reasons behind each rule.

## Patellar width and tilt

Two published widths, plus tilt, so a small patella can be told apart from a tilted one.
Like the literature and the manual measurements, all three are 2D, in-plane on axial slices:
the segmentation's own slices, which are parallel to the PD axial the raters read (median
0.0 deg apart; 4 of 215 scans over 5 deg).

- **`PW_effective_mm`**, Hevesi et al. (Am J Sports Med 2023), our manual "Effective
  Patellar Width": the widest distance between the medial and lateral borders of the
  patellar articular cartilage (label 4), measured parallel to the PBCL (moved, not
  rotated, to each slice). Every slice with cartilage is tried; the widest wins. A tilted
  patella reads narrower, since this is a projection onto the femur's axis.
- **`PW_max_mm`**, Chen et al. (Orthop J Sports Med 2024), after Staubli: the maximum width
  of the bony patella (label 3) along the patella's own axis, the line through its medial
  and lateral tips. Every slice is tried; the widest wins.
  - **Axis:** it starts as the slice's longest chord. Near the long axis a chord's length
    hardly changes with its angle, so on a voxel grid the longest chord can sit a voxel off
    the axis and tip it by several degrees; the axis is redrawn through the centres of the
    two tips (the voxels within 2 mm of each end), roughly where a reader clicks the borders.
    On a synthetic patella this recovers tilt to within 1 deg.
- **`patellar_tilt_deg`**, Chen 2024: on the PW_max slice, the angle between the PW_max line
  and the PBCL. **Positive = lateral tilt:** the patella's lateral border is more posterior
  than the PBCL would put it, so the lines converge laterally (Chen reports positive lateral
  tilt in both controls and patients). It is side-aware, like KRA, so left and right knees
  read alike.
- **Slices:** the two widths are taken where each is largest, so they can come from
  different slices. `PW_effective_slice`, `PW_max_slice` and their offsets from the PBCL
  slice say which; the run log and `summary.md` say how many knees differ and by how much.
- **Edge to edge:** each width runs from the outer corner of its outermost voxels, as for
  TTL's tibial width. Along a direction at angle a to the grid that adds |cos a| + |sin a|
  voxels: 1 along a grid axis, up to 1.41 at 45 deg. On 85 knees it added a median 1.19
  voxels to PW_effective (90% of knees under 1.33) and 1.23 to PW_max (under 1.39), whose
  axis can sit at a steeper angle to the grid.
- **Structures:** each is its largest 3D piece, so stray voxels can't widen it.
- **Flags in `notes`:**
  - label 3 or 4 missing;
  - a second bony piece over 5% of the largest (a bipartite fragment, which isn't measured);
  - either label in more than one piece on its chosen slice;
  - the patella reaching the edge of the scan;
  - PW_effective larger than PW_max. It happens on about 1 knee in 10, by up to about 2 mm,
    where the segmented cartilage wraps a little past the end of the bone;
  - a tilt over 45 deg: worth a look in the QC. On the J-Sign knees these were real, severely
    tilted patellae.
- **QC grid** ("Patellar width and tilt"): the PW_max slice with the bone (purple) and
  cartilage (red), the PW_effective slice's cartilage see-through (orange) when it's another
  slice, PW_max (blue; diamond = lateral end), PW_effective parallel to the PBCL with
  perpendicular reference lines at its ends (gold), and the PBCL's direction (red).

## Tuberosity out of the field of view

Many knee MRIs cut the tibia off above the
tuberosity. The detector still returns the most anterior point, but it isn't the
tuberosity, so TT-TG / sTT-TG / TTL come out wrong, often implausibly small. On the ORS
cohort, checked against the "out of window" column in the tracking table:

- every in-window case had at least 38 mm of tibia segmented;
- a threshold of under 37 mm caught 12 of 15 out-of-window cases with 0 of 30 false flags.

These cases are measured but flagged (a caution in `notes`). `statistics.csv` (with `--extras`), `compare` and `validate` also
report TT-dependent results with them excluded. In the app, the tuberosity of a flagged scan
can be marked on the axial scan instead, which usually reaches it (see the tuberosity review
in [app.md](app.md#3-tuberosity-review)); a reviewed scan is no longer flagged.

## Verification

- **Exact reproduction.** Re-measuring the Controls (30 cases) and Isolated Cartilage
  Damage (16 cases) cohorts reproduced every value from the original scripts (510/510 and
  272/272 values), which were run on the same files on 2026-10-14 to 2026-10-20.
- **Tests.** `python -m unittest discover -s tests -v` runs the unit tests, plus a
  regression test on 5 real cases (normal, converted 4D, short tibia, 0.6 mm and 0.4 mm
  scans) when the T7 drive is mounted.
- **Agreement with manual measurements** (ORS validation cohort, 45 knees, 2 known
  segmentation failures excluded): see [`VALIDATION.md`](../VALIDATION.md).

## Limitations

- The automated values have only been compared with manual measurements on the ORS cohort
  ([Verification](#verification)). Check the QC grids before trusting a new cohort's numbers.
- Some measurements use a different definition from the manual sheet, so agreement is
  expected to be limited:
  - LTI uses the cartilage reference slice.
  - TTL is a ratio from the medial tibial edge.
  - Sulcus angle is cartilage-based.
- Robust-z outlier flags are relative to the cohort itself. A uniformly bad cohort won't
  trigger them, but the segmentation checks will.
- Segmentations must be NRRD. NIfTI is read only for display scans, so the ACL Revisions
  `.nii.gz` segmentations would need converting first, and that cohort has no laterality
  sheet.
