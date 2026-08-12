#!/usr/bin/env python3
"""
Join region_genes + cell_types CSVs and write the 5-file normalized format
used by the ExSeq brain AD data explorer.

Sample names and input file paths come from data/sample_mapping.json.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_MAPPING = REPO_ROOT / "data" / "sample_mapping.json"
DEFAULT_OUTPUT = REPO_ROOT / "data" / "viewer_normalized"

REQUIRED_REGION_COLUMNS = (
    "region_name",
    "gene",
    "global_x",
    "global_y",
    "Z",
    "cell",
    "fov",
)


def load_mapping(path: Path) -> dict:
    with open(path) as f:
        return json.load(f)


def normalize_cell_id(cell_id: str) -> str:
    """
    Align region_genes cell IDs with cell_types_stellar cell_id values.

    Some samples encode FOV row in the prefix (e.g. 1H3_0 from FOV_31H3) while
    cell typing uses the short form (H03_0). Steps:
      1. Pad single-digit column index: 1H3_0 -> 1H03_0
      2. Drop leading row digit when present: 1H03_0 -> H03_0
    IDs already in short form (e.g. B04_0) pass through unchanged.
    """
    s = str(cell_id).strip()

    m = re.match(r"^(\d+)([A-Za-z])(\d)_(\d+)$", s)
    if m:
        s = f"{m.group(1)}{m.group(2)}0{m.group(3)}_{m.group(4)}"

    m2 = re.match(r"^\d+([A-Za-z]+\d+)_(\d+)$", s)
    if m2:
        return f"{m2.group(1)}_{m2.group(2)}"

    return s


def join_cell_types(regions_df: pd.DataFrame, cell_types_path: Path) -> pd.DataFrame:
    cell_types = pd.read_csv(cell_types_path)
    cell_types.columns = cell_types.columns.str.strip()

    if "cell_id" not in cell_types.columns or "predicted_cell_type" not in cell_types.columns:
        raise ValueError(
            f"{cell_types_path} must contain cell_id and predicted_cell_type columns"
        )

    cell_type_map = dict(
        zip(
            cell_types["cell_id"].astype(str).str.strip(),
            cell_types["predicted_cell_type"].astype(str).str.strip(),
        )
    )

    out = regions_df.copy()
    out["cell"] = out["cell"].astype(str).str.strip()
    out["cell_lookup"] = out["cell"].map(normalize_cell_id)
    out["cell_type"] = out["cell_lookup"].map(cell_type_map).fillna("unassigned")
    out = out.drop(columns=["cell_lookup"])
    return out


def apply_coordinate_scale(df: pd.DataFrame, scale: float) -> pd.DataFrame:
    out = df.copy()
    out["X"] = out["global_x"].astype(float) / scale
    out["Y"] = out["global_y"].astype(float) / scale
    out["Z"] = out["Z"].astype(float) / scale
    return out


def build_normalized_tables(df: pd.DataFrame, sample_name: str) -> dict[str, pd.DataFrame]:
    df = df.copy()
    df["region_name"] = df["region_name"].fillna("unassigned").astype(str)
    df["cell_type"] = df["cell_type"].fillna("unassigned").astype(str)
    df["gene"] = df["gene"].astype(str)
    df["fov"] = df["fov"].astype(str)
    df["cell"] = df["cell"].astype(str)

    points = pd.DataFrame(
        {
            "point_id": range(1, len(df) + 1),
            "gene": df["gene"].values,
            "Z": df["Z"].values,
            "X": df["X"].values,
            "Y": df["Y"].values,
            "fov": df["fov"].values,
        }
    )

    region_cols = ["region_name"]
    if "region_area" in df.columns:
        region_cols.append("region_area")
    if "region_proportion" in df.columns:
        region_cols.append("region_proportion")

    regions = (
        df[region_cols]
        .drop_duplicates(subset=["region_name"])
        .reset_index(drop=True)
    )
    regions.insert(0, "region_id", range(1, len(regions) + 1))
    region_id_map = dict(zip(regions["region_name"], regions["region_id"]))

    cells = df[["cell"]].drop_duplicates(subset=["cell"]).reset_index(drop=True)
    cells.insert(0, "cell_id", range(1, len(cells) + 1))
    cell_id_map = dict(zip(cells["cell"], cells["cell_id"]))

    points_regions = pd.DataFrame(
        {
            "point_id": points["point_id"],
            "region_id": df["region_name"].map(region_id_map),
        }
    )

    points_cells = pd.DataFrame(
        {
            "point_id": points["point_id"],
            "cell_id": df["cell"].map(cell_id_map),
            "cell_type": df["cell_type"].values,
        }
    )

    return {
        f"{sample_name}_points.csv": points,
        f"{sample_name}_regions.csv": regions,
        f"{sample_name}_points_regions.csv": points_regions,
        f"{sample_name}_cells.csv": cells,
        f"{sample_name}_points_cells.csv": points_cells,
    }


def convert_sample(
    sample_key: str,
    sample_cfg: dict,
    regions_genes_dir: Path,
    cell_types_dir: Path,
    output_dir: Path,
    scale: float,
) -> None:
    published_name = sample_cfg["published_name"]
    regions_path = regions_genes_dir / sample_cfg["regions_genes"]
    cell_types_path = cell_types_dir / sample_cfg["cell_types"]

    if not regions_path.exists():
        raise FileNotFoundError(f"Missing regions_genes file: {regions_path}")
    if not cell_types_path.exists():
        raise FileNotFoundError(f"Missing cell_types file: {cell_types_path}")

    regions_df = pd.read_csv(regions_path)
    regions_df.columns = regions_df.columns.str.strip()

    missing = [c for c in REQUIRED_REGION_COLUMNS if c not in regions_df.columns]
    if missing:
        raise ValueError(f"{regions_path} missing required columns: {missing}")

    merged = join_cell_types(regions_df, cell_types_path)
    merged = apply_coordinate_scale(merged, scale)
    tables = build_normalized_tables(merged, published_name)

    output_dir.mkdir(parents=True, exist_ok=True)
    for filename, table in tables.items():
        table.to_csv(output_dir / filename, index=False)

    assigned = (merged["cell_type"] != "unassigned").sum()
    print(
        f"{sample_key} -> {published_name}: "
        f"{len(merged):,} points, {len(tables[f'{published_name}_cells.csv']):,} cells, "
        f"{assigned:,} assigned cell types"
    )


def select_samples(mapping: dict, new_only: bool, sample_keys: list[str] | None) -> dict:
    samples = mapping["samples"]

    if sample_keys:
        unknown = [k for k in sample_keys if k not in samples]
        if unknown:
            raise ValueError(f"Unknown sample keys: {unknown}")
        return {k: samples[k] for k in sample_keys}

    if new_only:
        return {k: v for k, v in samples.items() if not v.get("existing_in_repo", False)}

    return samples


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert region_genes + cell_types CSVs into 5 normalized files per sample."
    )
    parser.add_argument(
        "--mapping",
        type=Path,
        default=DEFAULT_MAPPING,
        help=f"Path to sample_mapping.json (default: {DEFAULT_MAPPING})",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Output directory (default: {DEFAULT_OUTPUT})",
    )
    parser.add_argument(
        "--scale",
        type=float,
        default=3.5,
        help="Divide global_x/global_y/Z by this factor (default: 3.5)",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Process all samples in the mapping, including existing ones (default: new only)",
    )
    parser.add_argument(
        "--sample",
        action="append",
        dest="samples",
        metavar="KEY",
        help="Process one sample key from the mapping (repeatable)",
    )
    args = parser.parse_args()

    new_only = not args.all
    if args.samples:
        new_only = False

    mapping = load_mapping(args.mapping)
    regions_genes_dir = Path(mapping["regions_genes_dir"])
    cell_types_dir = Path(mapping["cell_types_dir"])
    selected = select_samples(mapping, new_only=new_only, sample_keys=args.samples)

    if not selected:
        print("No samples selected.")
        return

    print(f"Output directory: {args.output_dir}")
    print(f"Processing {len(selected)} sample(s)...")

    for sample_key, sample_cfg in selected.items():
        convert_sample(
            sample_key=sample_key,
            sample_cfg=sample_cfg,
            regions_genes_dir=regions_genes_dir,
            cell_types_dir=cell_types_dir,
            output_dir=args.output_dir,
            scale=args.scale,
        )

    print("Done.")


if __name__ == "__main__":
    main()
