"""Shared helpers for the cell-type annotation pipeline (paths, logging, LLM env)."""

from __future__ import annotations

import json
import logging
import os
import re
import sys
import warnings
from pathlib import Path
from typing import Any

warnings.filterwarnings("ignore")
os.environ.setdefault("PYTHONWARNINGS", "ignore")

# This file lives in <repo>/pipeline/ — project root is one level up.
_PIPELINE_DIR = Path(__file__).resolve().parent
REPO_ROOT = _PIPELINE_DIR.parent

# Default: <repo>/outputs/Exchange/
DEFAULT_OUT_SUB = "outputs"


def setup_logging(step: str) -> logging.Logger:
    logging.basicConfig(
        level=logging.INFO,
        format=f"%(asctime)s [%(levelname)s] [{step}] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stdout,
        force=True,
    )
    for name in (
        "scanpy",
        "anndata",
        "h5py",
        "matplotlib",
        "numba",
        "joblib",
        "omicverse",
        "openai",
        "httpx",
        "httpcore",
    ):
        logging.getLogger(name).setLevel(logging.ERROR)
    log = logging.getLogger(step)

    def _filter(record: logging.LogRecord) -> bool:
        return record.levelno != logging.WARNING

    for handler in logging.getLogger().handlers:
        handler.addFilter(_filter)
    for handler in log.handlers:
        handler.addFilter(_filter)
    return log


def get_out_dir() -> Path:
    root = os.environ.get("CELLTYPE_PIPELINE_OUT", "").strip()
    if root:
        return Path(root).resolve()
    return (REPO_ROOT / DEFAULT_OUT_SUB).resolve()


def exchange_dir() -> Path:
    d = get_out_dir() / "Exchange"
    d.mkdir(parents=True, exist_ok=True)
    return d


def figures_dir() -> Path:
    d = get_out_dir() / "figures"
    d.mkdir(parents=True, exist_ok=True)
    return d


def working_adata_path() -> Path:
    """Intermediate AnnData (large); removed by run_workflow.sh after success."""
    return get_out_dir() / ".workflow_adata.h5ad"


def final_adata_path() -> Path:
    """Single user-facing output h5ad (annotations + AUCell)."""
    name = os.environ.get("OUTPUT_H5AD_NAME", "SC_annotated.h5ad").strip() or "SC_annotated.h5ad"
    return get_out_dir() / name


def get_species() -> str:
    s = os.environ.get("SPECIES", "human").strip().lower()
    if s in ("mouse", "mm", "mus"):
        return "mouse"
    return "human"


def get_llm_config() -> dict[str, str]:
    base = os.environ.get("LLM_BASE_URL", "").strip()
    model = os.environ.get("LLM_MODEL", "").strip()
    key = os.environ.get("LLM_API_KEY", "").strip()
    if not base or not model:
        raise RuntimeError(
            "LLM_BASE_URL and LLM_MODEL must be set (use .env and/or run_workflow.sh prompts)."
        )
    return {"base_url": base, "api_key": key, "model": model}


def slug_geneset_name(name: str, used: set[str] | None = None) -> str:
    """Safe obs column prefix for Omicverse geneset_aucell (alphanumeric + underscore)."""
    s = re.sub(r"[^0-9a-zA-Z]+", "_", name.strip())
    s = re.sub(r"_+", "_", s).strip("_") or "celltype"
    s = s[:80]
    if used is not None:
        base = s
        i = 0
        while s in used:
            i += 1
            suffix = f"_{i}"
            s = (base[: max(0, 80 - len(suffix))] + suffix)[:80]
        used.add(s)
    return s


def save_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False, indent=2)


def merge_progress(updates: dict[str, Any]) -> None:
    p = exchange_dir() / "progress.json"
    prev: dict[str, Any] = {}
    if p.exists():
        with open(p, encoding="utf-8") as f:
            prev = json.load(f)
    prev.update(updates)
    save_json(p, prev)
