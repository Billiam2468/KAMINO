"""
Per-case result cache, so re-running a cohort only re-measures what changed.

A cached result is reused only if ALL of these match what produced it:
the segmentation file's content (SHA-1, not its name or date), the knee
side, the measurement set and QC options, and the kneemeasure code itself
(a fingerprint of every .py file in the package -- any code change
invalidates the whole cache, so a cached number can never come from older
code). Renamed files are fine; re-segmented or edited files are re-measured.
"""
import glob
import hashlib
import json
import os
import pickle

_CODE_FINGERPRINT = None


def file_sha1(path, chunk=1 << 20):
    h = hashlib.sha1()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(chunk), b''):
            h.update(block)
    return h.hexdigest()


def code_fingerprint():
    global _CODE_FINGERPRINT
    if _CODE_FINGERPRINT is None:
        h = hashlib.sha1()
        here = os.path.dirname(os.path.abspath(__file__))
        for p in sorted(glob.glob(os.path.join(here, '*.py'))):
            h.update(os.path.basename(p).encode())
            with open(p, 'rb') as f:
                h.update(f.read())
        _CODE_FINGERPRINT = h.hexdigest()[:12]
    return _CODE_FINGERPRINT


def cache_key(sha1, knee_side, measure_keys, qc, landmarks, overrides=None):
    review = hashlib.sha1(json.dumps(overrides or {}, sort_keys=True).encode()).hexdigest()[:12]
    return '|'.join([code_fingerprint(), sha1, knee_side or '-', review,
                     ','.join(measure_keys) if measure_keys is not None else 'all',
                     f'qc={int(bool(qc))}', f'lm={int(bool(landmarks))}'])


class ResultCache:
    def __init__(self, directory):
        self.dir = directory
        os.makedirs(directory, exist_ok=True)

    def _path(self, name):
        return os.path.join(self.dir, f'{name}.pkl')

    def get(self, name, key):
        p = self._path(name)
        if not os.path.exists(p):
            return None
        try:
            with open(p, 'rb') as f:
                stored_key, result = pickle.load(f)
        except Exception:
            return None
        return result if stored_key == key else None

    def put(self, name, key, result):
        tmp = self._path(name) + '.tmp'
        with open(tmp, 'wb') as f:
            pickle.dump((key, result), f, protocol=pickle.HIGHEST_PROTOCOL)
        os.replace(tmp, self._path(name))

    def prune(self, keep_names):
        """Drop cache entries for cases no longer in the cohort."""
        keep = set(keep_names)
        for p in glob.glob(os.path.join(self.dir, '*.pkl')):
            if os.path.basename(p)[:-4] not in keep:
                os.remove(p)
