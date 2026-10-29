# PyInstaller recipe for the downloadable app, run from the repository folder:
#
#     pyinstaller packaging/KAMINO.spec
#
# Gives dist/KAMINO/ (Windows: KAMINO.exe inside) or, on a Mac,
# dist/KAMINO.app. GitHub Actions builds both for each release
# (.github/workflows/release.yml).
import os
import re
import sys

ROOT = os.path.abspath(os.path.join(SPECPATH, '..'))
PKG = os.path.join(ROOT, 'kneemeasure')
VERSION = re.search(r"__version__ = '([^']+)'", open(os.path.join(PKG, '__init__.py')).read()).group(1)

datas = [(os.path.join(PKG, 'webapp', 'index.html'), 'kneemeasure/webapp')]
# the package's .py files as well: the result cache is keyed on their contents (cache.code_fingerprint),
# so a new version never reuses results measured by an older one
datas += [(os.path.join(PKG, f), 'kneemeasure') for f in sorted(os.listdir(PKG)) if f.endswith('.py')]

a = Analysis([os.path.join(SPECPATH, 'KAMINO.py')], pathex=[ROOT], datas=datas,
             hiddenimports=['kneemeasure.webapp', 'kneemeasure.cli'],
             excludes=['tkinter', 'PyQt5', 'PyQt6', 'PySide2', 'PySide6', 'IPython', 'pytest'])
pyz = PYZ(a.pure)
exe = EXE(pyz, a.scripts, [], exclude_binaries=True, name='KAMINO', console=False, upx=False)
coll = COLLECT(exe, a.binaries, a.datas, name='KAMINO', upx=False)
if sys.platform == 'darwin':
    app = BUNDLE(coll, name='KAMINO.app', bundle_identifier='org.kamino.app',
                 info_plist={'CFBundleShortVersionString': VERSION, 'CFBundleVersion': VERSION,
                             'NSHighResolutionCapable': True})
