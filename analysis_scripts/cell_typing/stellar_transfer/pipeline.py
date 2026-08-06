"""End-to-end supervised STELLAR encoder label-transfer pipeline."""

from __future__ import annotations

import platform
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import yaml

from .config import RunConfig
from .data import align_genes, load_reference, load_target
from .graph import graph_stats
from .model import (
    build_graph_dataset,
    import_upstream_components,
    predict,
    train_supervised_encoder,
)
from .plotting import add_plot_labels, make_all_plots
from .utils import dump_json, set_seed


def _package_versions() -> dict[str, str]:
    versions = {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "torch": torch.__version__,
    }
    try:
        import sklearn
        versions["scikit-learn"] = sklearn.__version__
    except Exception:
        pass
    try:
        import scipy
        versions["scipy"] = scipy.__version__
    except Exception:
        pass
    try:
        import torch_geometric
        versions["torch-geometric"] = torch_geometric.__version__
    except Exception:
        versions["torch-geometric"] = "not reported"
    return versions


def run_pipeline(config: RunConfig) -> Path:
    """Run label transfer and return the prediction-table path."""
    output_dir = config.paths.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    set_seed(config.training.seed)

    print("Loading reference...")
    reference = load_reference(
        config.paths.reference_csv,
        config.reference,
        config.graph,
        config.training.target_sum,
        config.training.seed,
        config.target.remove_background_cell_zero,
    )
    print("Loading target...")
    target = load_target(
        config.paths.target_count_csv,
        config.paths.target_points_csv,
        config.target,
        config.graph,
        config.training.target_sum,
    )
    reference_x, target_x, shared_genes = align_genes(reference, target)

    print(f"Reference cells: {reference_x.shape[0]}")
    print(f"Target cells:    {target_x.shape[0]}")
    print(f"Shared genes:    {len(shared_genes)}")
    reference_counts = pd.Series(reference.labels).map(reference.id_to_label).value_counts()
    print("Reference label counts:")
    print(reference_counts.to_string())

    GraphDataset, Encoder = import_upstream_components(config.paths.stellar_repo)
    dataset = build_graph_dataset(
        GraphDataset,
        reference_x,
        reference.labels,
        target_x,
        reference.edges,
        target.edges,
    )
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")
    print("Training supervised STELLAR encoder...")
    model, history = train_supervised_encoder(
        Encoder,
        dataset,
        n_classes=len(reference.id_to_label),
        epochs=config.training.epochs,
        learning_rate=config.training.learning_rate,
        weight_decay=config.training.weight_decay,
        partition_parts=config.graph.partition_parts,
        device=device,
    )

    ref_ids, ref_labels, ref_confidence, _ = predict(
        model, dataset.labeled_data, reference.id_to_label, device
    )
    target_ids, target_labels, target_confidence, probabilities = predict(
        model, dataset.unlabeled_data, reference.id_to_label, device
    )
    true_reference_labels = pd.Series(reference.labels).map(reference.id_to_label).astype(str).to_numpy()
    training_accuracy = float((ref_labels == true_reference_labels).mean())

    predictions = target.table.copy().reset_index(drop=True)
    predictions["pred_label"] = target_labels
    predictions["pred_id"] = target_ids
    predictions["pred_confidence"] = target_confidence
    for class_id, class_name in reference.id_to_label.items():
        predictions[f"prob_{class_name}"] = probabilities[:, class_id]
    predictions = add_plot_labels(predictions, config.plots)

    prediction_path = output_dir / "target_predictions_supervised_stellar_encoder.csv"
    predictions.to_csv(prediction_path, index=False)
    predictions["pred_label"].value_counts().rename_axis("pred_label").reset_index(
        name="n_cells"
    ).to_csv(output_dir / "target_prediction_counts.csv", index=False)
    if "plot_label" in predictions.columns:
        predictions["plot_label"].value_counts().rename_axis("plot_label").reset_index(
            name="n_cells"
        ).to_csv(output_dir / "target_prediction_counts_plot_labels.csv", index=False)
    if "region_name" in predictions.columns:
        pd.crosstab(predictions["region_name"], predictions["pred_label"]).to_csv(
            output_dir / "prediction_counts_by_region_diagnostic_only.csv"
        )

    pd.DataFrame(
        {
            "cell_id": reference.cell_ids,
            "true_label": true_reference_labels,
            "pred_label": ref_labels,
            "pred_id": ref_ids,
            "pred_confidence": ref_confidence,
        }
    ).to_csv(output_dir / "reference_training_check.csv", index=False)
    pd.DataFrame(history).to_csv(output_dir / "training_history.csv", index=False)
    (output_dir / "shared_genes.txt").write_text("\n".join(shared_genes) + "\n", encoding="utf-8")
    dump_json(output_dir / "label_mapping.json", reference.label_to_id)
    torch.save(model.state_dict(), output_dir / "model_state_dict.pt")
    shutil.copy2(config.config_path, output_dir / "config_used.yaml")

    summary = {
        "method": "supervised_STELLAR_encoder",
        "method_note": (
            "The published STELLAR GraphSAGE encoder was trained using cross-entropy "
            "on reference labels only. The original open-world/semi-supervised STELLAR "
            "objective and novel-class discovery procedure were not used."
        ),
        "reference_cells": int(reference_x.shape[0]),
        "target_cells": int(target_x.shape[0]),
        "shared_genes": int(len(shared_genes)),
        "reference_training_accuracy_sanity_check": training_accuracy,
        "reference_label_counts": reference_counts.to_dict(),
        "target_prediction_counts": pd.Series(target_labels).value_counts().to_dict(),
        "reference_graph": graph_stats(reference_x.shape[0], reference.edges),
        "target_graph": graph_stats(target_x.shape[0], target.edges),
        "graph_mode": config.graph.mode,
        "knn_k": config.graph.k if config.graph.mode == "knn" else None,
        "radius": config.graph.radius if config.graph.mode == "radius" else None,
        "make_undirected": config.graph.make_undirected,
        "partition_parts": config.graph.partition_parts,
        "epochs": config.training.epochs,
        "learning_rate": config.training.learning_rate,
        "weight_decay": config.training.weight_decay,
        "seed": config.training.seed,
        "target_sum": config.training.target_sum,
        "device": str(device),
        "versions": _package_versions(),
    }
    dump_json(output_dir / "run_summary.json", summary)
    make_all_plots(predictions, target.points_for_polygons, output_dir, config.plots)

    print("Done.")
    print(f"Predictions: {prediction_path}")
    return prediction_path
