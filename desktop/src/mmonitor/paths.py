"""
Path resolution for MMonitor.

Works in three modes:
1. pip-installed or editable install  – paths relative to __file__
2. PyInstaller frozen bundle          – paths via sys.executable / _internal
3. Legacy dev checkout fallback       – walk up to find project root
"""
import os
import sys
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

# ── Package directory (where this file lives) ──────────────────────
_PACKAGE_DIR = Path(__file__).resolve().parent  # .../mmonitor/

# ── Resources & workflows inside the package ───────────────────────
RESOURCES_DIR = str(_PACKAGE_DIR / 'resources')
IMAGES_DIR = str(_PACKAGE_DIR / 'resources' / 'images')
WORKFLOWS_DIR = str(_PACKAGE_DIR / 'workflows')


def _find_project_root() -> str:
    """Find the MMonitor project root (contains server/, desktop/, etc.).

    Used for locating the Django server directory and legacy resource paths.
    Returns the package dir itself as fallback if no project root is found
    (e.g. when pip-installed outside a checkout).
    """
    if getattr(sys, 'frozen', False):
        return _find_frozen_root()

    # Walk up from package dir looking for the project root
    current = _PACKAGE_DIR
    for parent in [current] + list(current.parents):
        if (parent / 'server' / 'manage.py').exists():
            return str(parent)
        if (parent / 'desktop' / 'setup.py').exists():
            return str(parent)
    # Fallback: package dir itself
    return str(_PACKAGE_DIR)


def _find_frozen_root() -> str:
    """Resolve root when running as a PyInstaller bundle."""
    app_path = os.path.dirname(sys.executable)

    if sys.platform == 'darwin' and 'Contents/MacOS' in app_path:
        app_path = os.path.abspath(os.path.join(app_path, '..', '..', '..'))
        internal_dir = os.path.join(os.path.dirname(app_path), "MMonitor", "_internal")
        if os.path.exists(internal_dir):
            return os.path.dirname(internal_dir)
        return app_path

    if '_internal' in app_path:
        root_dir = app_path
        while os.path.basename(root_dir) != '_internal' and os.path.dirname(root_dir) != root_dir:
            root_dir = os.path.dirname(root_dir)
        if os.path.basename(root_dir) == '_internal':
            return os.path.dirname(root_dir)

    return app_path


# ── Legacy compatibility names ─────────────────────────────────────
ROOT = _find_project_root()
PROJECT_ROOT = ROOT
SERVER_DIR = os.path.join(ROOT, "server")
SRC_DIR = os.path.join(ROOT, "desktop", "src")
LIB_DIR = str(_PACKAGE_DIR / 'lib')

# Frozen-app overrides for RESOURCES/IMAGES when running under PyInstaller
if getattr(sys, 'frozen', False):
    _internal = os.path.join(PROJECT_ROOT, "_internal")
    if os.path.exists(_internal):
        RESOURCES_DIR = os.path.join(_internal, "mmonitor", "src", "resources")
        IMAGES_DIR = os.path.join(RESOURCES_DIR, "images")
        LIB_DIR = os.path.join(_internal, "lib")
    CENTRIFUGER_DIR = os.path.join(LIB_DIR, "centrifuger-realtime", "centrifuger")
else:
    CENTRIFUGER_DIR = os.path.join(ROOT, "lib", "centrifuger-realtime", "centrifuger")

logger.debug("RESOURCES_DIR: %s", RESOURCES_DIR)
logger.debug("IMAGES_DIR: %s", IMAGES_DIR)
logger.debug("ROOT: %s", ROOT)
