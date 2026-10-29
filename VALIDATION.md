# Agreement with manual measurements: ORS validation cohort

**Data**

- Cohort: 45 knees from the ORS validation cohort (model segmentations, knee side from the
  cohort's own list). Two knees with known segmentation failures are excluded.
- Manual reference: `MPFL+ Master(ISO COMBINED NEW MASTER SHEET).csv`, matched by MRN and,
  for bilateral patients, by knee side. All 45 knees matched.
- Run date: 2026-10-29 (re-run for the single-slice TPCL and the bitangent PBCL; first run
  2026-10-21).

**Reproduce**

```bash
python -m kneemeasure validate "<run>/PFI ORS Validation" --manual "<sheet>" --id-column MRN \
    --side-column "Laterality - use this unless says bilateral, then use laterality column 'R'" \
    --pair "KRA_deg=Knee Rotation Angle (MRI) (+ve = tibial ER relative to femur" ...
```

**Statistics**

- Bias = automated - manual. LoA = 95% limits of agreement.
- ICC(2,1) = two-way random, absolute agreement.
- "TT in window" excludes cases with under 37 mm of tibia in the scan (see the README).

| measurement | manual column | knees | bias | LoA | MAE | r | ICC(2,1) |
|---|---|---|---|---|---|---|---|
| KRA (deg) | Knee Rotation Angle (MRI) | 45 | -1.3 | -9.9 / +7.4 | 3.4 | 0.83 | **0.82** |
| sTT-TG (mm), TT in window | sTT-TG (MRI) | 32 | +0.5 | -4.4 / +5.5 | 2.0 | 0.90 | **0.90** |
| sTT-TG (mm), all | | 43 | -0.3 | -6.1 / +5.6 | 2.3 | 0.85 | 0.86 |
| PTLTR (mm) | PTLTR | 43 | -2.8 | -13.0 / +7.5 | 4.7 | 0.80 | **0.73** |
| Bone bump (mm) | Bump Height_Bone (MRI) | 45 | -0.7 | -3.9 / +2.4 | 1.4 | 0.65 | 0.60 |
| CDI | CDI (MRI) | 25 | -0.01 | -0.34 / +0.33 | 0.14 | 0.58 | 0.58 |
| Cartilage bump (mm) | Bump Height_Cartilage (MRI) | 45 | -0.7 | -5.2 / +3.7 | 1.8 | 0.46 | 0.44 |
| TT-TG (mm), TT in window | TT-TG (MRI) | 32 | -0.3 | -11.9 / +11.3 | 4.0 | 0.41 | 0.41 |
| TT-TG (mm), all | | 43 | -3.2 | -20.8 / +14.3 | 6.2 | -0.05 | -0.04 |
| TTL t (mm), TT in window | Tubercle Lateralization (MRI) | 33 | +1.4 | -9.7 / +12.4 | 4.0 | 0.40 | 0.39 |
| TTL ratio, TT in window | Tubercle Lateralization % (MRI) | 33 | +0.04 | -0.11 / +0.18 | 0.06 | 0.06 | 0.05 |
| LTI (deg) | Lateral Trochlear Inclination (Cade) | 43 | +8.2 | -1.3 / +17.7 | 8.6 | 0.58 | 0.26 |
| Sulcus angle (deg) | Sulcus Angle (MRI) BONE | 43 | +5.2 | -20.9 / +31.3 | 11.7 | 0.39 | 0.17 |

Per-knee values and figures are in `<run>/PFI ORS Validation/manual_agreement/`:
`agreement.csv`, `agreement_cases.csv` and `agreement.png` (scatter plus Bland-Altman for
each measurement).

## What this says

- **The KRA sign convention is correct.** Automated KRA tracks the manual "+ve = tibial ER"
  column at r = 0.83, and the sign-flip check finds no reversal. There's a -1.3 deg bias.
- **sTT-TG, KRA and PTLTR agree well.**
- **TT-TG has a specific problem: the mediolateral position of a landmark.**
  - The anteroposterior component (sTT-TG) agrees at ICC 0.90, while TT-TG reaches only
    0.41. The two manual raters agree with each other at ICC 0.89 (from the IRR sheet), so
    this isn't noise in the reference.
  - Varying the tuberosity tie-break (scratch analysis, no code changed), with ICC(2,1)
    against manual TT-TG on in-window knees:

    | rule | ICC |
    |---|---|
    | current, 1 mm tie tolerance with the most lateral point chosen | 0.41 |
    | 2 mm tolerance | 0.51 |
    | 0 mm tolerance (pure most-anterior) | 0.12 |
    | the deepest groove instead of the reference-slice groove | 0.28 |

  - So the tuberosity's mediolateral pick is the sensitive part. 32 knees is too few to
    retune on safely, so the default is unchanged.
- **Out-of-window tuberosities ruin TT-TG/TTL.** Excluding them takes TT-TG agreement from
  -0.05 to 0.41.
- **LTI, sulcus angle and TTL differ in definition from the manual sheet**, so low
  agreement is expected rather than a malfunction:
  - LTI: the automated value uses the cartilage reference slice and reads about 8 deg
    higher, with r = 0.57.
  - Sulcus angle: the manual value is bony, the automated one cartilage-based.
  - TTL: the automated value measures t from the medial tibial edge along the TPCL.

  The two manual LTI raters (Cade, Divesh) give essentially the same agreement as each
  other.
- **The tibial posterior condylar line (TPCL) is now taken from one slice** (2026-10-29):
  the first slice below the tibial plateau with a clear posterior notch between two
  condylar lobes, with the bitangent of the posterior contour as the line. The plateau is
  where the tibia's cross-section first reaches 90% of its largest.

  | TPCL method | KRA ICC(2,1) |
  |---|---|
  | previous: the 5 deepest-notch slices, condyle peaks averaged | 0.79 |
  | the deepest-notch slice, bitangent, no averaging | 0.78 |
  | first notch counting from the top of the segmented tibia | 0.56 |
  | **first notch below the plateau (area >= 90% of max), bitangent** | **0.82** |
  | same, plateau at 85% / 95% of max area | 0.81 / 0.84 |

  - **Why the plateau matters:** the top ~10 mm of the segmented tibia cuts through the
    tibial spines and articular surfaces, which show false notches.
  - **Why it's more robust:** the deepest-notch rule sometimes drifted far down the shaft.
    For example, Y3082 took its line from about 13 mm further down than the new slice,
    which is 11 mm below the top of the tibia.
  - **Effect on TTL:** it also uses this line, and its agreement barely changed (t: 0.41 to
    0.39; the ratio, already poor, 0.17 to 0.05).
- **The femoral PBCL is now the true bitangent** (2026-10-29): the line tangent to the backs
  of both condyles, on the slice where it lies most posterior. It previously joined the two
  condyles' local maxima (peak tops).
  - **Line change:** median 0.1 deg, at most 2.2 deg.
  - **Agreement:** KRA, TT-TG, sTT-TG and LTI are unchanged.
  - **PTLTR fix:** PTLTR first dropped from 0.76 to 0.67. Its search for the lateral
    trochlear ridge divided lateral from medial at the midpoint of the PBCL's two points,
    and a shift of a few voxels let it reach onto the medial facet (two knees went from
    about 20 to 46 mm). It now divides at each slice's own sulcus, which gives 0.75. With the
    old PBCL that split changes nothing (0.76 either way).
  - **Bitangent on the old method's slice:** worse on every measurement (KRA 0.45),
    because on that slice the deepest concavity isn't always the intercondylar notch.
- **PTLTR now follows Mistovich et al. 2018 for the ridge and the line** (2026-10-29): the
  lateral trochlear ridge is taken on the tendon cut (the first slice below the patella),
  or the nearest slice with a usable ridge, and PTLTR is measured along the posterior
  condylar line on that same slice, as in the paper's Figure 2B. The knee's usual PBCL is
  unchanged and still used by every other measurement; PTLTR falls back to it when the
  ridge's slice has no notch (4 of 43 knees). The patella's lateral edge still stands in
  for the tendon's, which isn't segmented.

  | PTLTR ridge / line | ICC(2,1) |
  |---|---|
  | previous: most anterior lateral point over the trochlea / usual PBCL | 0.75 |
  | tendon cut / usual PBCL | 0.73 |
  | most anterior / line on the ridge's slice | 0.74 |
  | **tendon cut / line on the ridge's slice (now)** | **0.73** |

  - **Why it barely moves:** in 16 of the 43 knees the tendon cut is above the trochlea
    (the patella rides high with the knee straight), so the ridge comes from the top of the
    trochlea, where the old method usually found it too; the ridge stayed on the same slice
    in 19 knees and moved a median 1.8 mm. On the ridge's slice the line is within a median
    0.6 deg (at most 3.7 deg) of the usual PBCL.
  - **Kept for fidelity, not agreement:** the change follows the published method; agreement
    with the manual sheet is within noise of before. The manual values (mean 16.8 mm) are
    larger than tendon-based values would be, which suggests the raters measured to the
    patella's edge; the remaining bias (-2.8 mm) may come from how they measured.
