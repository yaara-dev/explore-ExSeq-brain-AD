"""General utilities used across the pipeline."""

from __future__ import annotations

import json
import os
import random
import re
from pathlib import Path
from typing import Any

import numpy as np
import torch


def limit_cpu_threads(n_threads: int = 1) -> None:
    """Set common numerical-library thread limits before heavy computation."""
    value = str(n_threads)
    for variable in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ.setdefault(variable, value)


def set_seed(seed: int) -> None:
    """Set deterministic seeds for Python, NumPy, and PyTorch."""
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def normalize_cell_id(value: Any) -> str:
    """Normalize identifiers such as C4_123 and C04-123 to C04_123."""
    text = str(value).strip()
    matches = re.findall(r"([A-Za-z])0?(\d{1,2})[_-](\d+)", text)
    if not matches:
        return text
    letter, fov_number, cell_number = matches[-1]
    return f"{letter.upper()}{int(fov_number):02d}_{int(cell_number)}"


def is_background_cell(cell_id: Any) -> bool:
    """Return True for segmentation-background IDs ending in '_0'."""
    return str(cell_id).strip().endswith("_0")


def normalize_log1p(matrix: np.ndarray, target_sum: float) -> np.ndarray:
    """Library-size normalize each cell and apply log(1+x)."""
    x = np.asarray(matrix, dtype=np.float32)
    row_sums = x.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    x = x / row_sums * float(target_sum)
    return np.log1p(x).astype(np.float32)


def dump_json(path: Path, payload: Any) -> None:
    """Write JSON with stable formatting and NumPy-safe fallbacks."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True, default=str)
