"""Shared helpers for the Telos NGS RNA-seq workflow."""
from __future__ import annotations

import logging
from pathlib import Path
import sys
from datetime import datetime, timezone


def repo_root() -> Path:
    """Return Telos_NGS_Project root for scripts stored in 2_scripts/2_RNA."""
    return Path(__file__).resolve().parents[2]


def resolve_from_repo(path_text: str, root: Path | None = None) -> Path:
    """Resolve absolute paths unchanged and relative paths from the repo root."""
    path = Path(path_text).expanduser()
    if path.is_absolute():
        return path
    return (root or repo_root()) / path


def setup_logger(name: str, log_dir: Path) -> tuple[logging.Logger, Path]:
    log_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path = log_dir / f"{name}_{timestamp}.log"

    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    logger.propagate = False

    formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
    file_handler = logging.FileHandler(log_path, encoding="utf-8")
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    return logger, log_path


def begin_run(done_path: Path) -> None:
    """Remove a previous success marker so it cannot be mistaken for this run."""
    done_path.parent.mkdir(parents=True, exist_ok=True)
    done_path.unlink(missing_ok=True)


def write_done_marker(done_path: Path, text: str) -> None:
    done_path.parent.mkdir(parents=True, exist_ok=True)
    done_path.write_text(text.rstrip() + "\n", encoding="utf-8")
