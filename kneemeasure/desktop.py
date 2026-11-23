"""
What the downloadable app (KAMINO.exe / KAMINO.app) runs on a
double-click: the same app as `python -m kneemeasure app`, opened in the
browser.

- The same version already running? Just opens it again in the browser.
- An older (or newer) version running? Asks it to stop, then starts this
  one in its place -- otherwise a new download would keep showing the old
  version. If it won't stop (it's measuring, or it's too old to know how),
  this one starts on a free port instead.
- Port 8765 taken by something else? Uses a free one.
- There's no window to close, so it stops by itself once no KAMINO page
  has been open for IDLE_MINUTES (never while measuring). Opening it again
  starts it again.
- Its messages go to ~/.config/kneemeasure/app.log.
"""
import json
import os
import sys
import time
import urllib.request
import webbrowser

from . import __version__

PORT = 8765
IDLE_MINUTES = 10
LOG = os.path.join(os.path.expanduser('~'), '.config', 'kneemeasure', 'app.log')


def _running_version(port=None):
    """The version of the KAMINO running on the port, or None if none is."""
    try:
        with urllib.request.urlopen(f'http://127.0.0.1:{port or PORT}/api/version', timeout=1) as r:
            j = json.load(r)
        return j.get('version') if j.get('app') == 'kneemeasure' else None
    except Exception:
        return None


def _stop_running(port=None, wait_s=5.0):
    """Ask the KAMINO on the port to stop; True once it has."""
    port = port or PORT
    req = urllib.request.Request(f'http://127.0.0.1:{port}/api/quit', data=b'{}', method='POST',
                                 headers={'Content-Type': 'application/json'})
    try:
        with urllib.request.urlopen(req, timeout=2):
            pass
    except Exception:
        return False                     # refused (measuring) or too old to have /api/quit
    end = time.time() + wait_s
    while time.time() < end:
        if _running_version(port) is None:
            return True
        time.sleep(0.2)
    return False


def _log_to_file():
    """No terminal: keep what would have been printed (errors included)."""
    os.makedirs(os.path.dirname(LOG), exist_ok=True)
    if os.path.exists(LOG) and os.path.getsize(LOG) > 1_000_000:
        os.replace(LOG, LOG + '.old')
    f = open(LOG, 'a', buffering=1, encoding='utf-8')
    sys.stdout = sys.stderr = f


def main():
    if sys.stdout is None or getattr(sys, 'frozen', False):
        _log_to_file()
    running = _running_version()
    if running == __version__:
        webbrowser.open(f'http://127.0.0.1:{PORT}/#setup')
        return
    if running is not None:
        stopped = _stop_running()
        print(f'KAMINO {running} was running: ' + ('stopped it, starting ' if stopped else
              "it didn't stop, so starting on another port: ") + __version__)
    from .webapp import serve
    serve(port=PORT, any_port_if_busy=True, idle_exit_minutes=IDLE_MINUTES)
