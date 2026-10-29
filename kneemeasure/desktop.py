"""
What the downloadable app (KAMINO.exe / KAMINO.app) runs on a
double-click: the same app as `python -m kneemeasure app`, opened in the
browser.

- Already running? Just opens it again in the browser.
- Port 8765 taken by something else? Uses a free one.
- There's no window to close, so it stops by itself once no KAMINO page
  has been open for IDLE_MINUTES (never while measuring). Opening it again
  starts it again.
- Its messages go to ~/.config/kneemeasure/app.log.
"""
import json
import os
import sys
import urllib.request
import webbrowser

PORT = 8765
IDLE_MINUTES = 10
LOG = os.path.join(os.path.expanduser('~'), '.config', 'kneemeasure', 'app.log')


def _already_running():
    try:
        with urllib.request.urlopen(f'http://127.0.0.1:{PORT}/api/version', timeout=1) as r:
            return json.load(r).get('app') == 'kneemeasure'
    except Exception:
        return False


def _log_to_file():
    """No terminal: keep what would have been printed (errors included)."""
    os.makedirs(os.path.dirname(LOG), exist_ok=True)
    if os.path.exists(LOG) and os.path.getsize(LOG) > 1_000_000:
        os.replace(LOG, LOG + '.old')
    f = open(LOG, 'a', buffering=1, encoding='utf-8')
    sys.stdout = sys.stderr = f


def main():
    if _already_running():
        webbrowser.open(f'http://127.0.0.1:{PORT}/#setup')
        return
    if sys.stdout is None or getattr(sys, 'frozen', False):
        _log_to_file()
    from .webapp import serve
    serve(port=PORT, any_port_if_busy=True, idle_exit_minutes=IDLE_MINUTES)
