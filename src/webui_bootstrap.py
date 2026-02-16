from __future__ import annotations

import os
import threading
import time
import webbrowser
from pathlib import Path
import sys

from streamlit.web import cli as stcli


def _resolve_webui_path() -> Path:
    if getattr(sys, "frozen", False):
        base = Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parent))
    else:
        base = Path(__file__).resolve().parent

    candidate = base / "src" / "webui.py"
    if candidate.exists():
        return candidate

    # development fallback
    project_root = base.parent
    fallback = project_root / "src" / "webui.py"
    if fallback.exists():
        return fallback

    raise FileNotFoundError(f"Could not locate webui.py. searched: {candidate} and {fallback}")


def _open_browser_later(url: str, delay: int = 2) -> None:
    def _launcher() -> None:
        time.sleep(delay)
        try:
            webbrowser.open_new_tab(url)
        except Exception:
            pass

    threading.Thread(target=_launcher, daemon=True).start()


def main() -> None:
    webui_script = _resolve_webui_path()
    host = os.environ.get("DH5A_WEBUI_HOST", "127.0.0.1")
    port = os.environ.get("DH5A_WEBUI_PORT", "8501")
    url = f"http://{host}:{port}"

    os.environ.setdefault("STREAMLIT_BROWSER_GATHER_USAGE_STATS", "false")
    os.environ.setdefault("PYTHONUTF8", "1")

    _open_browser_later(url, delay=2)
    sys.argv = [
        "streamlit",
        "run",
        str(webui_script),
        "--server.address",
        host,
        "--server.port",
        port,
        "--server.headless",
        "true",
    ]
    stcli.main()


if __name__ == "__main__":
    main()

