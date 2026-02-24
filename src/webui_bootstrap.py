from __future__ import annotations

import os
import traceback
import threading
import time
import webbrowser
from pathlib import Path
import sys
from typing import Optional

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


def _write_error_log(exc: Exception, webui_script: Optional[Path] = None) -> Path:
    log_dir = Path(os.environ.get("APPDATA", Path.home() / "AppData" / "Local")).expanduser()
    log_dir = log_dir / "DH5a-UTG"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / "webui_bootstrap_error.log"
    with log_path.open("a", encoding="utf-8") as f:
        f.write("=" * 80 + "\n")
        f.write(f"Error: {exc}\n")
        f.write(f"Script: {webui_script}\n")
        f.write(f"Executable: {Path(sys.executable)}\n")
        f.write(f"WorkingDir: {Path.cwd()}\n")
        traceback.print_exc(file=f)
    return log_path


def _hold_and_exit(message: str, detail: Optional[Path] = None) -> None:
    print(message)
    if detail:
        print(f"Log: {detail}")
    print("Window will stay for 15 seconds for inspection.")
    try:
        time.sleep(15)
    except Exception:
        pass


def main() -> None:
    try:
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
    except Exception as exc:
        log_path = _write_error_log(exc, webui_script if "webui_script" in locals() else None)
        _hold_and_exit(
            "DH5a-UTG WebUI startup failed. See message above for details.",
            detail=log_path,
        )
        raise


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"Unhandled exception: {exc}")
        # Keep a little time to make short-lived console close events visible.
        time.sleep(15)
        raise
