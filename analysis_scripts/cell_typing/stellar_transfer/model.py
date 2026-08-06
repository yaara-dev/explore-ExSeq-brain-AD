"""Compatibility layer and supervised training for the published STELLAR encoder."""

from __future__ import annotations

import copy
import importlib
import sys
from pathlib import Path
from typing import Any

import numpy as np
import torch


class StellarImportError(ImportError):
    """Raised when the upstream STELLAR checkout cannot be imported safely."""


def import_upstream_components(stellar_repo: Path) -> tuple[type, type]:
    """Import GraphDataset and Encoder from a local checkout of snap-stanford/stellar."""
    repo = stellar_repo.resolve()
    required = [repo / "datasets.py", repo / "models" / "__init__.py"]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise StellarImportError(
            "The configured STELLAR repository is incomplete. Missing: " + ", ".join(missing)
        )

    sys.path.insert(0, str(repo))
    try:
        datasets_module = importlib.import_module("datasets")
        models_module = importlib.import_module("models")
    except Exception as exc:
        raise StellarImportError(
            f"Could not import upstream STELLAR from {repo}: {exc}"
        ) from exc

    datasets_file = Path(getattr(datasets_module, "__file__", "")).resolve()
    models_file = Path(getattr(models_module, "__file__", "")).resolve()
    if repo not in datasets_file.parents or repo not in models_file.parents:
        raise StellarImportError(
            "Python imported a different package named 'datasets' or 'models'. "
            "Run in a clean environment and ensure paths.stellar_repo points to the "
            "snap-stanford/stellar checkout."
        )
    try:
        return datasets_module.GraphDataset, models_module.Encoder
    except AttributeError as exc:
        raise StellarImportError(
            "The upstream checkout does not expose GraphDataset and Encoder as expected."
        ) from exc


def build_graph_dataset(
    graph_dataset_class: type,
    reference_x: np.ndarray,
    reference_y: np.ndarray,
    target_x: np.ndarray,
    reference_edges: np.ndarray,
    target_edges: np.ndarray,
) -> Any:
    """Instantiate the upstream GraphDataset using edge lists shaped [n_edges, 2]."""
    return graph_dataset_class(
        reference_x.astype(np.float32),
        reference_y.astype(np.int64),
        target_x.astype(np.float32),
        reference_edges.astype(np.int64),
        target_edges.astype(np.int64),
    )


def _partitioned_loader(data: Any, partition_parts: int) -> Any:
    try:
        from torch_geometric.loader import ClusterData, ClusterLoader
    except ImportError:
        from torch_geometric.data import ClusterData, ClusterLoader

    try:
        clustered = ClusterData(data, num_parts=partition_parts, recursive=False)
        return ClusterLoader(clustered, batch_size=1, shuffle=True, num_workers=0)
    except Exception as exc:
        raise RuntimeError(
            "Partitioned graph training could not be initialized. This usually means "
            "that the installed PyTorch Geometric build lacks its METIS/partitioning "
            "backend. Set graph.partition_parts: 1 for full-graph training, or install "
            "the matching PyG optional dependencies."
        ) from exc


def train_supervised_encoder(
    encoder_class: type,
    dataset: Any,
    n_classes: int,
    epochs: int,
    learning_rate: float,
    weight_decay: float,
    partition_parts: int,
    device: torch.device,
) -> tuple[torch.nn.Module, list[dict[str, float | int]]]:
    """Train the STELLAR GraphSAGE encoder with cross-entropy on reference labels only."""
    input_dim = int(dataset.labeled_data.x.shape[1])
    model = encoder_class(input_dim, int(n_classes)).to(device)
    optimizer = torch.optim.Adam(
        model.parameters(), lr=float(learning_rate), weight_decay=float(weight_decay)
    )
    criterion = torch.nn.CrossEntropyLoss()

    if partition_parts == 1:
        full_reference_graph = copy.deepcopy(dataset.labeled_data).to(device)
        partitioned_batches = None
    else:
        full_reference_graph = None
        partitioned_batches = _partitioned_loader(dataset.labeled_data, partition_parts)

    history: list[dict[str, float | int]] = []
    for epoch in range(1, int(epochs) + 1):
        batches = [full_reference_graph] if full_reference_graph is not None else partitioned_batches
        model.train()
        total_loss = 0.0
        n_batches = 0
        for batch in batches:
            if full_reference_graph is None:
                batch = batch.to(device)
            optimizer.zero_grad(set_to_none=True)
            logits, _, _ = model(batch)
            loss = criterion(logits, batch.y)
            loss.backward()
            optimizer.step()
            total_loss += float(loss.item())
            n_batches += 1

        mean_loss = total_loss / max(n_batches, 1)
        reference_accuracy = evaluate_accuracy(model, dataset.labeled_data, device)
        history.append(
            {
                "epoch": epoch,
                "loss": mean_loss,
                "reference_training_accuracy": reference_accuracy,
            }
        )
        if epoch == 1 or epoch % 20 == 0 or epoch == epochs:
            print(
                f"epoch {epoch:03d} | loss={mean_loss:.4f} | "
                f"reference_training_accuracy={reference_accuracy:.3f}"
            )
    return model, history


def predict(
    model: torch.nn.Module,
    data: Any,
    id_to_label: dict[int, str],
    device: torch.device,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Predict class IDs, labels, maximum probability, and full probabilities."""
    model.eval()
    with torch.no_grad():
        graph = copy.deepcopy(data).to(device)
        logits, _, _ = model(graph)
        probabilities = torch.softmax(logits, dim=1).cpu().numpy()
    predicted_ids = probabilities.argmax(axis=1).astype(np.int64)
    confidence = probabilities.max(axis=1).astype(float)
    predicted_labels = np.asarray(
        [id_to_label.get(int(index), f"UNKNOWN_{index}") for index in predicted_ids],
        dtype=object,
    )
    return predicted_ids, predicted_labels, confidence, probabilities


def evaluate_accuracy(model: torch.nn.Module, data: Any, device: torch.device) -> float:
    """Compute training-graph accuracy as a sanity check, not held-out validation."""
    model.eval()
    with torch.no_grad():
        graph = copy.deepcopy(data).to(device)
        logits, _, _ = model(graph)
        predictions = logits.argmax(dim=1)
        return float((predictions == graph.y).float().mean().item())
