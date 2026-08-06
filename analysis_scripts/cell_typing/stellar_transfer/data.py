"""Input parsing, filtering, label harmonization, and gene alignment."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .config import GraphConfig, ReferenceConfig, TargetConfig
from .graph import build_graph
from .utils import is_background_cell, normalize_cell_id, normalize_log1p


DEFAULT_NON_GENE_COLUMNS = {
    "cell", "cell_id", "cell_norm", "region_name", "region_area", "region_proportion",
    "global_x", "global_y", "Z", "fov", "total_counts", "total_counts_raw",
    "n_genes_detected", "assigned_label", "confidence", "supporting_genes",
    "candidate_labels", "auto_label", "auto_confidence", "auto_supporting_genes",
    "top_genes", "score_details", "candidate_labels_v2", "auto_label_v2",
    "auto_confidence_v2", "auto_supporting_genes_v2", "top_genes_v2",
    "score_details_v2", "candidate_labels_v3", "auto_label_v3",
    "auto_confidence_v3", "auto_supporting_genes_v3", "top_genes_v3",
    "score_details_v3", "n_astro_core_hits", "n_immune_core_hits",
    "n_oligo_core_hits", "auto_label_v5", "auto_confidence_v5",
    "original_label_before_v5", "reference_label_source_v5", "astro_rescued_v5",
    "astro_probability_clean_classifier", "cell_type_label",
}


@dataclass
class LoadedDataset:
    x: np.ndarray
    genes: np.ndarray
    cell_ids: np.ndarray
    positions: np.ndarray
    edges: np.ndarray
    table: pd.DataFrame
    labels: np.ndarray | None = None
    label_to_id: dict[str, int] | None = None
    id_to_label: dict[int, str] | None = None
    points_for_polygons: pd.DataFrame | None = None


def _prepare_reference_index(frame: pd.DataFrame) -> pd.DataFrame:
    if "cell_id" in frame.columns:
        result = frame.set_index("cell_id")
    elif frame.columns[0].startswith("Unnamed"):
        result = frame.rename(columns={frame.columns[0]: "cell_id"}).set_index("cell_id")
    elif not pd.api.types.is_numeric_dtype(frame[frame.columns[0]]):
        result = frame.rename(columns={frame.columns[0]: "cell_id"}).set_index("cell_id")
    else:
        raise ValueError("Reference CSV must contain a cell_id column or a non-numeric first ID column.")
    result.index = result.index.astype(str).str.strip()
    if result.index.duplicated().any():
        examples = result.index[result.index.duplicated()].unique()[:5].tolist()
        raise ValueError(f"Reference cell IDs are not unique. Examples: {examples}")
    return result


def _reference_gene_columns(frame: pd.DataFrame, config: ReferenceConfig) -> list[str]:
    if config.gene_columns is not None:
        missing = [column for column in config.gene_columns if column not in frame.columns]
        if missing:
            raise ValueError(f"Configured reference gene columns are missing: {missing[:10]}")
        return list(config.gene_columns)

    excluded = DEFAULT_NON_GENE_COLUMNS | set(config.non_gene_columns)
    columns = [
        column for column in frame.columns
        if column not in excluded and pd.api.types.is_numeric_dtype(frame[column])
    ]
    if not columns:
        raise ValueError(
            "No numeric reference gene columns were found. Supply reference.gene_columns "
            "or extend reference.non_gene_columns in the YAML configuration."
        )
    return columns


def _canonicalize_label(raw_label: object, config: ReferenceConfig) -> str | None:
    if pd.isna(raw_label):
        return None
    label = str(raw_label).strip()
    if label in set(config.exclude_labels):
        return None
    if label in config.label_map:
        return config.label_map[label]
    if label in set(config.label_map.values()):
        return label
    return None


def _collapse_gene_columns(expression: pd.DataFrame) -> pd.DataFrame:
    renamed = expression.copy()
    renamed.columns = [str(column).strip().upper() for column in renamed.columns]
    return renamed.T.groupby(level=0, sort=False).sum().T


def _encode_labels(labels: pd.Series) -> tuple[np.ndarray, dict[str, int], dict[int, str]]:
    unique_labels = sorted(labels.astype(str).unique())
    label_to_id = {label: index for index, label in enumerate(unique_labels)}
    id_to_label = {index: label for label, index in label_to_id.items()}
    encoded = labels.map(label_to_id).to_numpy(dtype=np.int64)
    return encoded, label_to_id, id_to_label


def load_reference(
    path: Path,
    config: ReferenceConfig,
    graph_config: GraphConfig,
    target_sum: float,
    seed: int,
    remove_background_cell_zero: bool,
) -> LoadedDataset:
    """Load and preprocess the annotated reference cell-by-gene CSV."""
    frame = _prepare_reference_index(pd.read_csv(path))

    if remove_background_cell_zero:
        frame = frame.loc[~frame.index.map(is_background_cell)].copy()
    if config.label_column not in frame.columns:
        raise ValueError(f"Reference label column '{config.label_column}' was not found.")

    if config.keep_confidence is not None:
        if config.confidence_column is None or config.confidence_column not in frame.columns:
            raise ValueError(
                "reference.keep_confidence was provided, but the configured confidence column is missing."
            )
        confidence = frame[config.confidence_column].astype(str).str.strip().str.lower()
        frame = frame.loc[confidence.isin(config.keep_confidence)].copy()

    frame["cell_type_label"] = frame[config.label_column].map(lambda value: _canonicalize_label(value, config))
    frame = frame.loc[frame["cell_type_label"].notna()].copy()

    gene_columns = _reference_gene_columns(frame, config)
    expression = _collapse_gene_columns(frame[gene_columns].apply(pd.to_numeric, errors="raise"))
    total_counts = (
        pd.to_numeric(frame["total_counts"], errors="coerce")
        if "total_counts" in frame.columns
        else expression.sum(axis=1)
    )
    keep = total_counts.fillna(0) >= config.min_total_counts
    frame = frame.loc[keep].copy()
    expression = expression.loc[keep].copy()

    label_counts = frame["cell_type_label"].value_counts()
    retained_labels = label_counts[label_counts >= config.min_cells_per_label].index
    frame = frame.loc[frame["cell_type_label"].isin(retained_labels)].copy()
    expression = expression.loc[frame.index].copy()
    if frame.empty:
        raise ValueError("No reference cells remain after label harmonization and filtering.")

    if config.balance:
        sampled_indices: list[str] = []
        for _, subset in frame.groupby("cell_type_label", sort=True):
            if len(subset) > config.max_cells_per_label:
                subset = subset.sample(config.max_cells_per_label, random_state=seed)
            sampled_indices.extend(subset.index.tolist())
        rng = np.random.default_rng(seed)
        sampled_indices = list(np.asarray(sampled_indices)[rng.permutation(len(sampled_indices))])
        frame = frame.loc[sampled_indices].copy()
        expression = expression.loc[sampled_indices].copy()

    required_coordinates = {"global_x", "global_y"}
    if not required_coordinates.issubset(frame.columns):
        raise ValueError("Reference CSV must contain global_x and global_y columns.")
    positions = frame[["global_x", "global_y"]].apply(pd.to_numeric, errors="raise").to_numpy(np.float32)
    edges = build_graph(
        positions,
        mode=graph_config.mode,
        k=graph_config.k,
        radius=graph_config.radius,
        make_undirected=graph_config.make_undirected,
    )
    labels, label_to_id, id_to_label = _encode_labels(frame["cell_type_label"])

    output_table = frame.copy()
    output_table.insert(0, "cell_id", output_table.index.astype(str))
    return LoadedDataset(
        x=normalize_log1p(expression.to_numpy(np.float32), target_sum),
        genes=expression.columns.to_numpy(dtype=object),
        cell_ids=frame.index.to_numpy(dtype=object),
        positions=positions,
        edges=edges,
        table=output_table.reset_index(drop=True),
        labels=labels,
        label_to_id=label_to_id,
        id_to_label=id_to_label,
    )


def load_target(
    count_path: Path,
    points_path: Path,
    config: TargetConfig,
    graph_config: GraphConfig,
    target_sum: float,
) -> LoadedDataset:
    """Load target genes-by-cells counts and transcript/point metadata."""
    counts = pd.read_csv(count_path)
    points = pd.read_csv(points_path)
    if counts.empty or counts.columns[0].strip().lower() != "gene":
        raise ValueError("Target count table must have a first column named 'gene'.")
    if "cell" not in points.columns:
        raise ValueError("Target points CSV must contain a 'cell' column.")
    for coordinate in ("global_x", "global_y"):
        if coordinate not in points.columns:
            raise ValueError(f"Target points CSV must contain '{coordinate}'.")

    genes = counts.iloc[:, 0].astype(str).str.strip().str.upper()
    normalized_ids = [normalize_cell_id(column) for column in counts.columns[1:]]
    expression = pd.DataFrame(
        counts.iloc[:, 1:].T.to_numpy(dtype=np.float32),
        index=normalized_ids,
        columns=genes,
    )
    expression = expression.groupby(level=0, sort=False).sum()
    expression = expression.T.groupby(level=0, sort=False).sum().T

    points = points.copy()
    points["cell_norm"] = points["cell"].map(normalize_cell_id)
    if config.remove_background_cell_zero:
        expression = expression.loc[~expression.index.map(is_background_cell)].copy()
        points = points.loc[~points["cell_norm"].map(is_background_cell)].copy()

    aggregations: dict[str, tuple[str, str]] = {
        "global_x": ("global_x", "mean"),
        "global_y": ("global_y", "mean"),
    }
    for optional_column in ("region_name", "Z", "fov"):
        if optional_column in points.columns:
            aggregations[optional_column] = (
                optional_column,
                "mean" if optional_column == "Z" else "first",
            )
    metadata = points.groupby("cell_norm", sort=False).agg(**aggregations)
    metadata.index = metadata.index.astype(str)

    total_counts = expression.sum(axis=1)
    expression = expression.loc[total_counts >= config.min_total_counts].copy()
    total_counts = total_counts.loc[expression.index]
    common_cells = sorted(set(expression.index) & set(metadata.index))
    if not common_cells:
        raise ValueError("No overlapping cell IDs were found between target counts and points metadata.")

    expression = expression.loc[common_cells].copy()
    metadata = metadata.loc[common_cells].copy()
    total_counts = total_counts.loc[common_cells]
    positions = metadata[["global_x", "global_y"]].apply(pd.to_numeric, errors="raise").to_numpy(np.float32)
    edges = build_graph(
        positions,
        mode=graph_config.mode,
        k=graph_config.k,
        radius=graph_config.radius,
        make_undirected=graph_config.make_undirected,
    )

    table = metadata.copy()
    table.insert(0, "cell_id", expression.index.astype(str))
    table["total_counts_raw"] = total_counts.to_numpy()
    for column, default in (("region_name", "NA"), ("Z", np.nan), ("fov", "NA")):
        if column not in table.columns:
            table[column] = default

    polygon_columns = ["cell_norm", "global_x", "global_y"]
    if "Z" in points.columns:
        polygon_columns.append("Z")
    points_for_polygons = points.loc[
        points["cell_norm"].isin(common_cells), polygon_columns
    ].copy()

    return LoadedDataset(
        x=normalize_log1p(expression.to_numpy(np.float32), target_sum),
        genes=expression.columns.to_numpy(dtype=object),
        cell_ids=expression.index.to_numpy(dtype=object),
        positions=positions,
        edges=edges,
        table=table.reset_index(drop=True),
        points_for_polygons=points_for_polygons,
    )


def align_genes(reference: LoadedDataset, target: LoadedDataset) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Return reference and target matrices restricted to shared genes in stable order."""
    reference_index = {str(gene).strip().upper(): i for i, gene in enumerate(reference.genes)}
    target_index = {str(gene).strip().upper(): i for i, gene in enumerate(target.genes)}
    shared = sorted(set(reference_index) & set(target_index))
    if not shared:
        raise ValueError("Reference and target have no shared genes after case normalization.")
    ref_columns = np.asarray([reference_index[gene] for gene in shared], dtype=np.int64)
    target_columns = np.asarray([target_index[gene] for gene in shared], dtype=np.int64)
    return reference.x[:, ref_columns], target.x[:, target_columns], shared
