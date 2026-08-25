import networkx as nx
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

    h-leafB and h-leafC carry large time_length values so Populate_Edge (run with
    weight_attr="time_length") inserts a phantom node into each of them, shifting
    the distance-matrix row/col order away from the original node labels.
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
        use_data_prep=True,
        dress_distance_matrix=dress_distance_matrix,
    )
    cf.find_cycles()
    return cf


def test_dress_distance_matrix_writes_correct_cells(tmp_path):
    cf = _run(tmp_path, "true_distance_between_tips")

    # Populate_Edge (weight_attr="time_length") should have split the two
    # expensive h-leafB / h-leafC edges, appending exactly 2 phantom nodes.
    assert cf.G.number_of_nodes() == 5 + 2
    assert cf.pre_dress_matrix.shape == (7, 7)

    node_to_index = {node: idx for idx, node in cf.index_to_node.items()}
    # sanity: original node labels keep their identity across the mutation,
    # but do NOT land at matrix positions equal to their own label (100 != 0, etc)
    # -- this is exactly the mismatch a raw dist_matrix[node, node] bug would hit.
    assert node_to_index[100] == 0   # root, first node inserted
    assert node_to_index[250] == 1   # h
    assert node_to_index[3] == 2     # leafA
    assert node_to_index[410] == 3   # leafB
    assert node_to_index[77] == 4    # leafC
    assert set(cf.index_to_node.values()) - {100, 250, 3, 410, 77} == {
        "added_node_5", "added_node_6",
    }

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
    # distance differs from the 'time_length' distance Populate_Edge originally
    # filled in, so if dressing used the wrong attribute (or hit the wrong
    # cells) this would still equal the pre-dress value.
    assert cf.pre_dress_matrix[node_to_index[3], node_to_index[410]] == pytest.approx(3.0)
    assert cf.pre_dress_matrix[node_to_index[410], node_to_index[77]] == pytest.approx(52.0)


def test_dress_distance_matrix_never_selects_phantom_nodes_as_leaves(tmp_path):
    cf = _run(tmp_path, "true_distance_between_tips")
    node_to_index = {node: idx for idx, node in cf.index_to_node.items()}
    phantom_nodes = {"added_node_5", "added_node_6"}

    # cf.dress_original_G is the pre-Populate_Edge snapshot the method actually
    # scanned for leaves -- confirm the phantom nodes never even existed in it,
    # so they structurally cannot have been picked up as tips.
    assert phantom_nodes.isdisjoint(cf.dress_original_G.nodes)
    assert phantom_nodes <= set(cf.G.nodes)  # they DO exist in the post-populate graph

    # and confirm those two phantom nodes' rows/cols in the distance matrix,
    # which is indexed over the post-populate graph, were left untouched
    phantom_indices = {node_to_index[n] for n in phantom_nodes}
    for idx in phantom_indices:
        assert cf.pre_dress_matrix[idx, :] == pytest.approx(cf.post_dress_matrix[idx, :])
        assert cf.pre_dress_matrix[:, idx] == pytest.approx(cf.post_dress_matrix[:, idx])


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
