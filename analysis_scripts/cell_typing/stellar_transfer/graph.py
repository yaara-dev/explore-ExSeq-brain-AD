"""Spatial graph construction."""

from __future__ import annotations

import numpy as np
from sklearn.neighbors import NearestNeighbors, radius_neighbors_graph


def build_graph(
    positions: np.ndarray,
    mode: str = "knn",
    k: int = 3,
    radius: float = 80.0,
    make_undirected: bool = True,
) -> np.ndarray:
    """Build an edge list with shape [n_edges, 2] from XY coordinates."""
    xy = np.asarray(positions, dtype=np.float32)
    if xy.ndim != 2 or xy.shape[1] != 2:
        raise ValueError(f"positions must have shape [n_cells, 2]; received {xy.shape}.")
    if xy.shape[0] < 2:
        raise ValueError("At least two cells are required to build a graph.")
    if not np.isfinite(xy).all():
        raise ValueError("Spatial coordinates contain NaN or infinite values.")

    edge_set: set[tuple[int, int]] = set()
    mode = mode.strip().lower()

    if mode == "knn":
        effective_k = min(int(k) + 1, xy.shape[0])
        neighbors = NearestNeighbors(n_neighbors=effective_k, algorithm="auto").fit(xy)
        indices = neighbors.kneighbors(xy, return_distance=False)
        for source, row in enumerate(indices):
            for destination in row[1:]:
                edge_set.add((int(source), int(destination)))
                if make_undirected:
                    edge_set.add((int(destination), int(source)))
    elif mode == "radius":
        matrix = radius_neighbors_graph(
            xy,
            radius=float(radius),
            mode="connectivity",
            include_self=False,
        ).tocoo()
        for source, destination in zip(matrix.row, matrix.col):
            edge_set.add((int(source), int(destination)))
            if make_undirected:
                edge_set.add((int(destination), int(source)))
    else:
        raise ValueError("Graph mode must be either 'knn' or 'radius'.")

    if not edge_set:
        raise ValueError("No graph edges were created. Increase radius or use KNN mode.")
    return np.asarray(sorted(edge_set), dtype=np.int64)


def graph_stats(n_nodes: int, edges: np.ndarray) -> dict[str, float | int]:
    """Summarize node degree and graph size."""
    edge_array = np.asarray(edges, dtype=np.int64)
    degree = np.zeros(int(n_nodes), dtype=np.int64)
    np.add.at(degree, edge_array[:, 0], 1)
    np.add.at(degree, edge_array[:, 1], 1)
    isolates = int((degree == 0).sum())
    return {
        "num_nodes": int(n_nodes),
        "num_directed_edges": int(edge_array.shape[0]),
        "mean_incident_degree": float(degree.mean()),
        "median_incident_degree": float(np.median(degree)),
        "max_incident_degree": int(degree.max()),
        "isolates": isolates,
        "isolate_fraction": float(isolates / max(int(n_nodes), 1)),
    }
