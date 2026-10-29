"""
Standard segmentation label convention used throughout this package:

    1 = femur bone        2 = femur cartilage
    3 = patella bone      4 = patella cartilage
    5 = tibia bone

There is no separate tibia cartilage label.
"""

FEM_B = 1
FEM_C = 2
PAT_B = 3
PAT_C = 4
TIB_B = 5

STANDARD_LABELS = (FEM_B, FEM_C, PAT_B, PAT_C, TIB_B)

LABEL_NAMES = {
    FEM_B: 'femur bone',
    FEM_C: 'femur cartilage',
    PAT_B: 'patella bone',
    PAT_C: 'patella cartilage',
    TIB_B: 'tibia bone',
}

CARTILAGE_LABELS = (FEM_C, PAT_C)
