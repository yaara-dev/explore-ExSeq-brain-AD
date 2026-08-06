from pathlib import Path

import numpy as np

from stellar_transfer.config import load_config
from stellar_transfer.data import align_genes, load_reference, load_target
from stellar_transfer.graph import build_graph, graph_stats
from stellar_transfer.utils import normalize_cell_id, normalize_log1p


ROOT = Path(__file__).resolve().parents[1]


def test_cell_id_normalization():
    assert normalize_cell_id("C4_123") == "C04_123"
    assert normalize_cell_id("prefix_X49G09_123") == "G09_123"
    assert normalize_cell_id("unmatched") == "unmatched"


def test_normalization_is_finite():
    matrix = np.array([[0, 0], [1, 3]], dtype=float)
    result = normalize_log1p(matrix, 10000)
    assert result.shape == matrix.shape
    assert np.isfinite(result).all()


def test_knn_graph_has_edges():
    positions = np.array([[0, 0], [1, 0], [4, 0]], dtype=float)
    edges = build_graph(positions, mode="knn", k=1, make_undirected=True)
    stats = graph_stats(3, edges)
    assert edges.shape[1] == 2
    assert stats["num_directed_edges"] >= 4


def test_toy_inputs_load_and_align():
    config = load_config(ROOT / "examples" / "toy_config.yaml")
    reference = load_reference(
        config.paths.reference_csv,
        config.reference,
        config.graph,
        config.training.target_sum,
        config.training.seed,
        config.target.remove_background_cell_zero,
    )
    target = load_target(
        config.paths.target_count_csv,
        config.paths.target_points_csv,
        config.target,
        config.graph,
        config.training.target_sum,
    )
    ref_x, target_x, genes = align_genes(reference, target)
    assert ref_x.shape[1] == target_x.shape[1] == len(genes) == 12
    assert len(reference.id_to_label) == 6
    assert target.points_for_polygons is not None
