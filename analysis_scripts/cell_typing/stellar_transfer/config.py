"""Configuration loading and validation."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


class ConfigError(ValueError):
    """Raised when a configuration file is incomplete or invalid."""


def _require(mapping: dict[str, Any], key: str, section: str) -> Any:
    if key not in mapping:
        raise ConfigError(f"Missing required key '{section}.{key}'.")
    return mapping[key]


def _resolve_path(value: str | Path, base_dir: Path) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = (base_dir / path).resolve()
    return path


@dataclass(frozen=True)
class PathsConfig:
    stellar_repo: Path
    reference_csv: Path
    target_count_csv: Path
    target_points_csv: Path
    output_dir: Path


@dataclass(frozen=True)
class ReferenceConfig:
    label_column: str
    confidence_column: str | None
    keep_confidence: tuple[str, ...] | None
    min_total_counts: float
    min_cells_per_label: int
    balance: bool
    max_cells_per_label: int
    exclude_labels: tuple[str, ...]
    label_map: dict[str, str]
    gene_columns: tuple[str, ...] | None
    non_gene_columns: tuple[str, ...]


@dataclass(frozen=True)
class TargetConfig:
    min_total_counts: float
    remove_background_cell_zero: bool


@dataclass(frozen=True)
class GraphConfig:
    mode: str
    k: int
    radius: float
    make_undirected: bool
    partition_parts: int


@dataclass(frozen=True)
class TrainingConfig:
    seed: int
    target_sum: float
    epochs: int
    learning_rate: float
    weight_decay: float


@dataclass(frozen=True)
class PlotConfig:
    enabled: bool
    invert_y: bool
    point_size: float
    dpi: int
    polygon_z_half_width: float
    split_excitatory_by_region: bool
    excitatory_label: str
    ca_plot_label: str
    dg_plot_label: str
    ca_like_regions: tuple[str, ...]
    colors: dict[str, str]


@dataclass(frozen=True)
class RunConfig:
    paths: PathsConfig
    reference: ReferenceConfig
    target: TargetConfig
    graph: GraphConfig
    training: TrainingConfig
    plots: PlotConfig
    config_path: Path


def load_config(path: str | Path) -> RunConfig:
    """Load a YAML configuration file into validated dataclasses."""
    config_path = Path(path).expanduser().resolve()
    with config_path.open("r", encoding="utf-8") as handle:
        raw = yaml.safe_load(handle) or {}

    base_dir = config_path.parent
    paths_raw = _require(raw, "paths", "root")
    ref_raw = _require(raw, "reference", "root")
    target_raw = raw.get("target", {})
    graph_raw = _require(raw, "graph", "root")
    train_raw = _require(raw, "training", "root")
    plot_raw = raw.get("plots", {})

    paths = PathsConfig(
        stellar_repo=_resolve_path(_require(paths_raw, "stellar_repo", "paths"), base_dir),
        reference_csv=_resolve_path(_require(paths_raw, "reference_csv", "paths"), base_dir),
        target_count_csv=_resolve_path(_require(paths_raw, "target_count_csv", "paths"), base_dir),
        target_points_csv=_resolve_path(_require(paths_raw, "target_points_csv", "paths"), base_dir),
        output_dir=_resolve_path(_require(paths_raw, "output_dir", "paths"), base_dir),
    )

    keep_conf = ref_raw.get("keep_confidence", ["high", "medium"])
    gene_columns = ref_raw.get("gene_columns")
    reference = ReferenceConfig(
        label_column=str(_require(ref_raw, "label_column", "reference")),
        confidence_column=(
            None if ref_raw.get("confidence_column") in (None, "")
            else str(ref_raw["confidence_column"])
        ),
        keep_confidence=(
            None if keep_conf is None
            else tuple(str(x).strip().lower() for x in keep_conf)
        ),
        min_total_counts=float(ref_raw.get("min_total_counts", 40)),
        min_cells_per_label=int(ref_raw.get("min_cells_per_label", 10)),
        balance=bool(ref_raw.get("balance", True)),
        max_cells_per_label=int(ref_raw.get("max_cells_per_label", 250)),
        exclude_labels=tuple(str(x).strip() for x in ref_raw.get("exclude_labels", [])),
        label_map={str(k).strip(): str(v).strip() for k, v in ref_raw.get("label_map", {}).items()},
        gene_columns=(None if gene_columns is None else tuple(str(x) for x in gene_columns)),
        non_gene_columns=tuple(str(x) for x in ref_raw.get("non_gene_columns", [])),
    )

    target = TargetConfig(
        min_total_counts=float(target_raw.get("min_total_counts", 40)),
        remove_background_cell_zero=bool(target_raw.get("remove_background_cell_zero", True)),
    )

    mode = str(graph_raw.get("mode", "knn")).strip().lower()
    if mode not in {"knn", "radius"}:
        raise ConfigError("graph.mode must be 'knn' or 'radius'.")
    graph = GraphConfig(
        mode=mode,
        k=int(graph_raw.get("k", 3)),
        radius=float(graph_raw.get("radius", 80.0)),
        make_undirected=bool(graph_raw.get("make_undirected", True)),
        partition_parts=int(graph_raw.get("partition_parts", 1)),
    )
    if graph.k < 1:
        raise ConfigError("graph.k must be at least 1.")
    if graph.partition_parts < 1:
        raise ConfigError("graph.partition_parts must be at least 1.")

    training = TrainingConfig(
        seed=int(train_raw.get("seed", 1)),
        target_sum=float(train_raw.get("target_sum", 1e4)),
        epochs=int(train_raw.get("epochs", 60)),
        learning_rate=float(train_raw.get("learning_rate", 1e-3)),
        weight_decay=float(train_raw.get("weight_decay", 0.0)),
    )
    if training.epochs < 1:
        raise ConfigError("training.epochs must be at least 1.")

    plots = PlotConfig(
        enabled=bool(plot_raw.get("enabled", True)),
        invert_y=bool(plot_raw.get("invert_y", False)),
        point_size=float(plot_raw.get("point_size", 7)),
        dpi=int(plot_raw.get("dpi", 300)),
        polygon_z_half_width=float(plot_raw.get("polygon_z_half_width", 2)),
        split_excitatory_by_region=bool(plot_raw.get("split_excitatory_by_region", True)),
        excitatory_label=str(plot_raw.get("excitatory_label", "Excitatory_Neuron")),
        ca_plot_label=str(plot_raw.get("ca_plot_label", "CA_GLUT")),
        dg_plot_label=str(plot_raw.get("dg_plot_label", "DG_GLUT")),
        ca_like_regions=tuple(str(x) for x in plot_raw.get("ca_like_regions", ["CA1", "CA3", "SO"])),
        colors={str(k): str(v) for k, v in plot_raw.get("colors", {}).items()},
    )

    return RunConfig(
        paths=paths,
        reference=reference,
        target=target,
        graph=graph,
        training=training,
        plots=plots,
        config_path=config_path,
    )
