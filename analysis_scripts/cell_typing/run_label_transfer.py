#!/usr/bin/env python3
"""Command-line entry point for supervised STELLAR encoder label transfer."""

from __future__ import annotations

import argparse

from stellar_transfer.config import load_config
from stellar_transfer.pipeline import run_pipeline
from stellar_transfer.utils import limit_cpu_threads


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Train the published STELLAR graph encoder with supervised cross-entropy "
            "on an annotated spatial reference and transfer labels to a target tissue."
        )
    )
    parser.add_argument("--config", required=True, help="Path to a YAML configuration file.")
    parser.add_argument(
        "--threads", type=int, default=1,
        help="Thread limit for BLAS/OpenMP libraries (default: 1).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    limit_cpu_threads(args.threads)
    config = load_config(args.config)
    run_pipeline(config)


if __name__ == "__main__":
    main()
