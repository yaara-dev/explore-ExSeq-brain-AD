"""Publication-friendly diagnostic plots."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Polygon
from scipy.spatial import ConvexHull

from .config import PlotConfig


def add_plot_labels(predictions: pd.DataFrame, config: PlotConfig) -> pd.DataFrame:
    """Optionally split excitatory predictions by anatomical region for plotting only."""
    result = predictions.copy()
    result["plot_label"] = result["pred_label"].astype(str)
    if not config.split_excitatory_by_region or "region_name" not in result.columns:
        return result
    mask = result["pred_label"].astype(str).eq(config.excitatory_label)
    result.loc[mask, "plot_label"] = np.where(
        result.loc[mask, "region_name"].isin(config.ca_like_regions),
        config.ca_plot_label,
        config.dg_plot_label,
    )
    return result


def plot_prediction_counts(predictions: pd.DataFrame, output_dir: Path, config: PlotConfig) -> None:
    counts = predictions["pred_label"].value_counts()
    plt.figure(figsize=(7, 4.5))
    counts.plot(kind="bar")
    plt.ylabel("Number of cells")
    plt.xlabel("Predicted cell type")
    plt.title("Supervised STELLAR encoder predictions")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(output_dir / "prediction_counts_barplot.png", dpi=config.dpi)
    plt.close()


def plot_spatial_predictions(predictions: pd.DataFrame, output_dir: Path, config: PlotConfig) -> None:
    label_column = "plot_label" if "plot_label" in predictions.columns else "pred_label"
    plt.figure(figsize=(10, 9))
    for label in sorted(predictions[label_column].dropna().unique()):
        subset = predictions.loc[predictions[label_column] == label]
        y = -subset["global_y"] if config.invert_y else subset["global_y"]
        plt.scatter(
            subset["global_x"], y, s=config.point_size, alpha=0.85, linewidths=0,
            label=f"{label} (n={len(subset)})", color=config.colors.get(label, "gray"),
        )
    plt.xlabel("global_x")
    plt.ylabel("-global_y" if config.invert_y else "global_y")
    plt.title("Target tissue predictions")
    plt.gca().set_aspect("equal", adjustable="box")
    plt.legend(markerscale=2, bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    plt.tight_layout()
    plt.savefig(output_dir / "spatial_predictions.png", dpi=config.dpi)
    plt.savefig(output_dir / "spatial_predictions.pdf")
    plt.close()


def plot_spatial_confidence(predictions: pd.DataFrame, output_dir: Path, config: PlotConfig) -> None:
    y = -predictions["global_y"] if config.invert_y else predictions["global_y"]
    plt.figure(figsize=(9, 8))
    scatter = plt.scatter(
        predictions["global_x"], y, c=predictions["pred_confidence"],
        s=config.point_size, alpha=0.9, linewidths=0,
    )
    plt.xlabel("global_x")
    plt.ylabel("-global_y" if config.invert_y else "global_y")
    plt.title("Prediction confidence")
    plt.gca().set_aspect("equal", adjustable="box")
    plt.colorbar(scatter, label="maximum class probability")
    plt.tight_layout()
    plt.savefig(output_dir / "spatial_prediction_confidence.png", dpi=config.dpi)
    plt.close()


def plot_spatial_polygons(
    predictions: pd.DataFrame,
    points: pd.DataFrame,
    output_dir: Path,
    config: PlotConfig,
) -> None:
    """Plot RNA-position convex hulls; these are not segmentation-derived cell boundaries."""
    label_column = "plot_label" if "plot_label" in predictions.columns else "pred_label"
    cell_to_label = predictions.set_index("cell_id")[label_column].astype(str).to_dict()
    plt.figure(figsize=(10, 9))
    axis = plt.gca()
    plotted: set[str] = set()

    for cell_id, subset in points.groupby("cell_norm", sort=False):
        if cell_id not in cell_to_label:
            continue
        if "Z" in subset.columns and subset["Z"].notna().any():
            middle_z = subset["Z"].median()
            central = subset.loc[
                subset["Z"].between(
                    middle_z - config.polygon_z_half_width,
                    middle_z + config.polygon_z_half_width,
                )
            ]
            if len(central) >= 3:
                subset = central

        coordinates = subset[["global_x", "global_y"]].dropna().to_numpy(float)
        if coordinates.shape[0] < 3:
            continue
        try:
            hull = ConvexHull(coordinates)
        except Exception:
            continue
        hull_coordinates = coordinates[hull.vertices].copy()
        if config.invert_y:
            hull_coordinates[:, 1] *= -1

        label = cell_to_label[cell_id]
        color = config.colors.get(label, "gray")
        axis.add_patch(
            Polygon(
                hull_coordinates, closed=True, facecolor=color, edgecolor=color,
                alpha=0.55, linewidth=0.25,
                label=label if label not in plotted else None,
            )
        )
        plotted.add(label)

    plt.xlabel("global_x")
    plt.ylabel("-global_y" if config.invert_y else "global_y")
    plt.title("RNA-position convex hulls by predicted cell type")
    axis.set_aspect("equal", adjustable="box")
    axis.autoscale_view()
    if plotted:
        plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    plt.tight_layout()
    plt.savefig(output_dir / "spatial_predictions_rna_hulls.png", dpi=config.dpi)
    plt.savefig(output_dir / "spatial_predictions_rna_hulls.pdf")
    plt.close()


def make_all_plots(
    predictions: pd.DataFrame,
    points: pd.DataFrame | None,
    output_dir: Path,
    config: PlotConfig,
) -> None:
    if not config.enabled:
        return
    plot_prediction_counts(predictions, output_dir, config)
    plot_spatial_predictions(predictions, output_dir, config)
    plot_spatial_confidence(predictions, output_dir, config)
    if points is not None and not points.empty:
        plot_spatial_polygons(predictions, points, output_dir, config)
