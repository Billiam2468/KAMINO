"""
KAMINO (Knee Anatomy Measurements from Imaging and Network Output); Python
package kneemeasure -- standalone automated knee MRI measurements from labeled
segmentations (1 femur bone, 2 femur cartilage, 3 patella bone,
4 patella cartilage, 5 tibia bone).

Python use:

    from kneemeasure.pipeline import measure_case
    r = measure_case('case.nrrd', 'L')
    r['values']          # {'KRA_deg': ..., 'TT_TG_mm': ..., ...}

    from kneemeasure.runner import run_cohort
    run_cohort(['segs/'], 'results/', laterality='mapping.csv')
"""
__version__ = '0.1.0'
GITHUB_REPO = 'Billiam2468/KAMINO'     # where releases (the downloadable app) are published
