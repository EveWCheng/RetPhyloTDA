import networkx as nx
import numpy as np
import pytest
from itertools import combinations

from find_cycles import CycleFinder


def _small_reticulate_graph():
    """A tiny sim_bdh-shaped network: undirected, int node ids (deliberately NOT
    0..n-1 or in insertion order, like real post-filter_nodes ids), is_leaf flags,
    and both 'length'/'time_length' edge weights that disagree on which path is
    shortest -- so a test using them can tell which attribute dressing actually used.

        root(100) --- h(250) --- leafA(3)
          |             |  \
          |             |   leafC(77)
          +--- leafB(410)+

    'length' shortest paths:   leafA-leafB=2 (via h), leafA-leafC=3, leafB-leafC=3 (via h)
    'time_length' shortest:    leafA-leafB=3 (via root!), leafA-leafC=51, leafB-leafC=52 (via root!)
    """
    G = nx.Graph()
    G.add_node(100, is_leaf=False)   # root
    G.add_node(250, is_leaf=False)   # h (internal / hybrid junction)
    G.add_node(3, is_leaf=True)      # leafA
    G.add_node(410, is_leaf=True)    # leafB
    G.add_node(77, is_leaf=True)     # leafC

    G.add_edge(100, 250, length=1.0, time_length=1.0)
    G.add_edge(100, 410, length=10.0, time_length=1.0)
    G.add_edge(250, 3, length=1.0, time_length=1.0)
    G.add_edge(250, 410, length=1.0, time_length=20.0)
    G.add_edge(250, 77, length=2.0, time_length=50.0)
    return G


class _RecordingCycleFinder(CycleFinder):
    """Captures the dist_matrix immediately before/after dressing, and the
    original_G snapshot used, without changing any pipeline behavior."""

    def dress_distance_matrix_with_choice(self, dist_matrix, original_G):
        self.pre_dress_matrix = dist_matrix.copy()
        self.dress_original_G = original_G.copy()
        result = super().dress_distance_matrix_with_choice(dist_matrix, original_G)
        self.post_dress_matrix = result.copy()
        return result


def _run(tmp_path, dress_distance_matrix):
    G = _small_reticulate_graph()
    cf = _RecordingCycleFinder(
        G,
        threshold_mode=["cyclelength"],
        cycle_qualify_mode=[],
        output_dir=str(tmp_path),
        weight_attr="time_length",
        vis=False,
        use_data_prep=False,
        dress_distance_matrix=dress_distance_matrix,
    )
    cf.find_cycles()
    return cf


def test_dress_distance_matrix_writes_correct_cells(tmp_path):
    cf = _run(tmp_path, "true_distance_between_tips")

    assert cf.G.number_of_nodes() == 5
    assert cf.pre_dress_matrix.shape == (5, 5)

    node_to_index = {node: idx for idx, node in cf.index_to_node.items()}
    # sanity: node ids keep their identity but do NOT land at matrix positions
    # equal to their own id (100 != 0, etc) -- the mismatch a raw
    # dist_matrix[node, node] bug would hit.
    assert node_to_index[100] == 0   # root, first node inserted
    assert node_to_index[250] == 1   # h
    assert node_to_index[3] == 2     # leafA
    assert node_to_index[410] == 3   # leafB
    assert node_to_index[77] == 4    # leafC
    assert set(cf.index_to_node.values()) == {100, 250, 3, 410, 77}

    # independently-computed true ('length'-weighted) tip distances on a
    # freshly built copy of the pristine graph -- not derived from the code
    # under test, so this isn't circular.
    expected_G = _small_reticulate_graph()
    expected = {}
    for u, v in combinations([3, 410, 77], 2):
        expected[frozenset((u, v))] = nx.shortest_path_length(expected_G, u, v, weight="length")
    assert expected[frozenset((3, 410))] == 2
    assert expected[frozenset((3, 77))] == 3
    assert expected[frozenset((410, 77))] == 3

    for u, v in combinations([3, 410, 77], 2):
        iu, iv = node_to_index[u], node_to_index[v]
        want = expected[frozenset((u, v))]
        assert cf.post_dress_matrix[iu, iv] == pytest.approx(want)
        assert cf.post_dress_matrix[iv, iu] == pytest.approx(want)  # symmetry preserved

    # and dressing actually changed something -- these pairs' true 'length'
    # distance differs from the 'time_length' distance that filled the matrix
    # first, so a wrong-attribute or wrong-cell bug would leave the pre-dress value.
    assert cf.pre_dress_matrix[node_to_index[3], node_to_index[410]] == pytest.approx(3.0)
    assert cf.pre_dress_matrix[node_to_index[410], node_to_index[77]] == pytest.approx(52.0)


def test_dress_distance_matrix_leaves_other_cells_untouched(tmp_path):
    cf = _run(tmp_path, "true_distance_between_tips")
    node_to_index = {node: idx for idx, node in cf.index_to_node.items()}
    leaf_indices = {node_to_index[n] for n in (3, 410, 77)}

    n = cf.pre_dress_matrix.shape[0]
    for i in range(n):
        for j in range(n):
            if i in leaf_indices and j in leaf_indices and i != j:
                continue  # these are the cells dressing is allowed to change
            assert cf.pre_dress_matrix[i, j] == pytest.approx(cf.post_dress_matrix[i, j]), (
                f"cell ({i}, {j}) changed but isn't a leaf<->leaf pair"
            )


def test_dress_distance_matrix_off_by_default_is_a_no_op(tmp_path):
    cf = _run(tmp_path, None)
    assert cf.pre_dress_matrix == pytest.approx(cf.post_dress_matrix)


def test_dress_distance_matrix_single_leaf_pair_is_a_no_op(tmp_path):
    # only one leaf -> combinations(leaf_nodes, 2) is empty -> nothing to write,
    # must not raise (e.g. from an empty node_to_index lookup or similar)
    G = nx.Graph()
    G.add_node("root", is_leaf=False)
    G.add_node("leaf1", is_leaf=True)
    G.add_edge("root", "leaf1", length=1.0, time_length=1.0)

    cf = _RecordingCycleFinder(
        G,
        threshold_mode=["cyclelength"],
        cycle_qualify_mode=[],
        output_dir=str(tmp_path),
        weight_attr="time_length",
        vis=False,
        use_data_prep=True,
        dress_distance_matrix="true_distance_between_tips",
    )
    cf.find_cycles()
    assert cf.pre_dress_matrix == pytest.approx(cf.post_dress_matrix)


# ── dress log: ratio + reticulation_edges_on_path ────────────────────────────
# call dress_distance_matrix_with_choice() directly with a hand-built dist_matrix
# and original_G, isolating the logging logic from Populate_Edge / harmonics.

def _dress_cf(tmp_path, dress="true_distance_between_tips"):
    return CycleFinder(
        nx.Graph(), threshold_mode=[], cycle_qualify_mode=[],
        output_dir=str(tmp_path), vis=False, dress_distance_matrix=dress,
    )


def _reticulate_original_G():
    """m--A is a reticulation shortcut; the tree-only A<->B route goes A--r--m--B.

    length shortest paths: A-B = 2 (A-m-B, via reticulation), B-C = 3 (B-m-r-C, tree).
    """
    G = nx.Graph()
    G.add_node("r", is_leaf=False, label="r")
    G.add_node("m", is_leaf=False, label="m")
    G.add_node("A", is_leaf=True, label="A")
    G.add_node("B", is_leaf=True, label="B")
    G.add_node("C", is_leaf=True, label="C")
    G.add_edge("r", "A", length=10.0, edge_type="tree")
    G.add_edge("r", "m", length=1.0, edge_type="tree")
    G.add_edge("m", "B", length=1.0, edge_type="tree")
    G.add_edge("m", "A", length=1.0, edge_type="reticulation", inher_weight=0.3)
    G.add_edge("r", "C", length=1.0, edge_type="tree")
    return G


def test_dress_log_records_ratio_and_reticulation_edges_on_path(tmp_path):
    og = _reticulate_original_G()
    cf = _dress_cf(tmp_path)
    cf.index_to_node = {0: "r", 1: "m", 2: "A", 3: "B", 4: "C"}
    dist = np.zeros((5, 5))
    for (i, j), d in {(2, 3): 12.0, (2, 4): 11.0, (3, 4): 3.0}.items():
        dist[i, j] = dist[j, i] = d

    out = cf.dress_distance_matrix_with_choice(dist, og)
    log = {(e["i"], e["j"]): e for e in cf.log["dress_distance_matrix"]}

    ab = log[("A", "B")]
    assert ab["old_dist"] == 12.0
    assert ab["new_dist"] == 2.0
    assert ab["changed"] is True
    assert ab["ratio"] == pytest.approx((2.0 - 12.0) / 12.0)
    assert ab["reticulation_edges_on_path"] == [{"edge": ["A", "m"], "inher_weight": 0.3}]
    assert out[2, 3] == pytest.approx(2.0)
    assert out[3, 2] == pytest.approx(2.0)

    bc = log[("B", "C")]
    assert bc["reticulation_edges_on_path"] == []
    assert bc["ratio"] == pytest.approx(0.0)
    assert bc["changed"] is False


def test_dress_log_ratio_is_none_when_old_dist_zero(tmp_path):
    og = nx.Graph()
    og.add_node("r", is_leaf=False, label="r")
    og.add_node("A", is_leaf=True, label="A")
    og.add_node("B", is_leaf=True, label="B")
    og.add_edge("r", "A", length=0.0, edge_type="tree")
    og.add_edge("r", "B", length=0.0, edge_type="tree")

    cf = _dress_cf(tmp_path)
    cf.index_to_node = {0: "r", 1: "A", 2: "B"}
    cf.dress_distance_matrix_with_choice(np.zeros((3, 3)), og)

    entry = cf.log["dress_distance_matrix"][0]
    assert (entry["i"], entry["j"]) == ("A", "B")
    assert entry["old_dist"] == 0.0
    assert entry["ratio"] is None


def test_dress_log_reticulation_edge_uses_str_node_when_label_missing(tmp_path):
    og = nx.Graph()
    og.add_node(1, is_leaf=False)  # no label attr
    og.add_node(2, is_leaf=True, label="A")
    og.add_node(3, is_leaf=True, label="B")
    og.add_edge(1, 2, length=1.0, edge_type="reticulation", inher_weight=0.25)
    og.add_edge(1, 3, length=1.0, edge_type="tree")

    cf = _dress_cf(tmp_path)
    cf.index_to_node = {0: 1, 1: 2, 2: 3}
    dist = np.zeros((3, 3))
    dist[1, 2] = dist[2, 1] = 5.0
    cf.dress_distance_matrix_with_choice(dist, og)

    entry = cf.log["dress_distance_matrix"][0]
    assert entry["reticulation_edges_on_path"] == [{"edge": ["A", "1"], "inher_weight": 0.25}]


def test_dress_log_not_created_when_dressing_off(tmp_path):
    cf = _dress_cf(tmp_path, dress=None)
    cf.index_to_node = {0: "r", 1: "m", 2: "A", 3: "B", 4: "C"}
    cf.dress_distance_matrix_with_choice(np.zeros((5, 5)), _reticulate_original_G())
    assert "dress_distance_matrix" not in cf.log
