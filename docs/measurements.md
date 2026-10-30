# Measurements: methods, validation and limitations

How each value in `measurements.csv` is calculated is in `Measurement methods.xlsx`; this
page has the column definitions, the landmark rules and how far the values can be trusted.

All are computed in 3D from the segmentation.

- **Axial measurements** use the reoriented voxel frame scaled to mm. Relative
  angles and distances are unaffected by any in-plane rotation in the scan header.
- **The tibial tuberosity** search runs in true physical coordinates, because "most
  anterior" is an absolute direction.


| column | definition |
|---|---|
| `KRA_deg` | Knee rotation angle: tibial posterior condylar line (TPCL) vs femoral posterior bicondylar line (PBCL). **+ = tibial external rotation** relative to the femur. |
| `TT_TG_mm` | TT-TG along the PBCL; the groove is the sulcus of the Dejour cartilage reference slice |
| `sTT_TG_mm` | sagittal TT-TG, perpendicular to the PBCL; **negative = TT posterior to TG** |
| `PTLTR_mm` | patella lateral margin to the lateral trochlear bony ridge on the tendon cut (the first slice below the patella, or the nearest slice with a usable ridge), along the posterior condylar line on that slice (unsigned; Mistovich 2018, with the patella's edge for the tendon's) |
| `TTL`, `TTL_t_mm`, `TTL_T_mm` | tubercle lateralization t/T along the TPCL (t = medial tibial edge to TT, T = tibial width) |
| `LTI_deg` | lateral trochlear inclination at the Dejour reference slice (+ = lateral facet anterior to PBCL) |
| `sulcus_angle_deg` | cartilage sulcus angle at the reference slice |
| `medial_/lateral_facet_depth_mm` | facet height above the sulcus, perpendicular to the PBCL |
| `trochlear_depth_mm` | sulcus depth below the facet chord, perpendicular to the PBCL |
| `medial_/lateral_facet_length_mm`, `facet_asymmetry_ratio` | sulcus to facet-peak distances; asymmetry = medial / lateral |
| `bone_bump_mm`, `cartilage_bump_mm` | trochlear bump height above the anterior femoral cortex line, on the sagittal slice through the deepest groove |
| `CDI`, `CDI_numerator_mm`, `patella_cartilage_length_mm` | Caton-Deschamps index on the sagittal slice with the longest patellar cartilage |

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
