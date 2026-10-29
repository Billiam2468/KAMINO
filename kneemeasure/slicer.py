"""
3D Slicer export of the automated landmarks, so they can be reviewed (and,
later, corrected) on top of the original scan + segmentation.

Each case gets one markups point list, <case>.mrk.json, in true physical
LPS mm (Slicer's native file coordinate system) -- drag it into Slicer
alongside the segmentation. Point labels say which measurement and
landmark each point belongs to.
"""
import json

SCHEMA = ('https://raw.githubusercontent.com/slicer/slicer/main/Modules/Loadable/Markups/'
          'Resources/Schema/markups-schema-v1.0.3.json#')


def write_markups(path, case_name, landmarks):
    """landmarks: [(measure, label, (x, y, z) index, (L, P, S) mm)]."""
    points = []
    for i, (measure, label, _idx, lps) in enumerate(landmarks, start=1):
        points.append({
            'id': str(i),
            'label': f'{label}',
            'description': f'{measure}',
            'position': [round(c, 3) for c in lps],
            'orientation': [-1.0, -0.0, -0.0, -0.0, -1.0, -0.0, 0.0, 0.0, 1.0],
            'selected': True,
            'locked': True,
            'visibility': True,
            'positionStatus': 'defined',
        })
    doc = {
        '@schema': SCHEMA,
        'markups': [{
            'type': 'Fiducial',
            'coordinateSystem': 'LPS',
            'coordinateUnits': 'mm',
            'locked': True,
            'labelFormat': '%N-%d',
            'controlPoints': points,
            'display': {'visibility': True, 'textScale': 2.0, 'glyphScale': 2.0,
                        'glyphType': 'Sphere3D', 'pointLabelsVisibility': True},
        }],
    }
    with open(path, 'w') as f:
        json.dump(doc, f, indent=2)
