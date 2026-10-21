import json

import pytest
import networkx as nx
import numpy as np

from find_cycles import CycleFinder
from main_result_analysis import MainResultAnalysis


NODES_CSV = """id,label,type,is_leaf
-1,-1,internal,False
-2,-2,internal,False
-3,-3,internal,False
sp4,sp4,leaf,True
hyb5,hyb5,hyb_leaf,True
"""

EDGES_CSV = """from,to,edge_type,length,time_length,inher_weight
-1,-2,tree,1.0,1.0,
-1,-3,tree,2.0,2.0,0.7
-2,-3,reticulation,0.5,0.0,0.3
-2,sp4,tree,3.0,3.0,
-3,hyb5,tree,1.5,1.5,
"""


def _circle_dist_matrix(n_points):
    """Return the Euclidean distance matrix of n points evenly spaced on a unit circle."""
    angles = np.linspace(0, 2 * np.pi, n_points, endpoint=False)
    points = np.column_stack([np.cos(angles), np.sin(angles)])
    return np.linalg.norm(points[:, None, :] - points[None, :, :], axis=-1)


# every edge in a cycle's weight history should be sampled at the same filtration times
def test_weight_history_times_line_up_across_edges_in_each_cycle(tmp_path):
    n_points = 8
    G = nx.cycle_graph(n_points)
    index_to_name = {i: f"sp{i}" for i in range(n_points)}
    cf = CycleFinder(G, threshold_mode=[], output_dir=str(tmp_path), dist_matrix=_circle_dist_matrix(n_points), index_to_name=index_to_name, sim_label="sim0", harmonic_histories=True)
    cf.find_cycles()

    histories = MainResultAnalysis(str(tmp_path)).read_weight_histories(0)
    assert len(histories) > 0

    for cycle_index, edges in histories.items():
        expected_times = None
        for simplex, points in edges.items():
            times = [t for t, _ in points]
            if expected_times is None:
                expected_times = times
            assert times == expected_times, f"cycle {cycle_index} edge {simplex} has different history times"


# the CSVs should come back as an undirected label-keyed graph with typed node and edge attributes
def test_read_filtered_graph_builds_undirected_typed_graph(tmp_path):
    sim_csv_dir = tmp_path / "phylo_csv" / "sim0"
    sim_csv_dir.mkdir(parents=True)
    (sim_csv_dir / "sim0_filtered_nodes.csv").write_text(NODES_CSV)
    (sim_csv_dir / "sim0_filtered_edges.csv").write_text(EDGES_CSV)

    G, D = MainResultAnalysis(str(tmp_path)).read_filtered_graph(0)

    assert isinstance(G, nx.Graph) and not G.is_directed()
    assert set(G.nodes()) == {"-1", "-2", "-3", "sp4", "hyb5"}
    assert G.nodes["sp4"] == {"label": "sp4", "type": "leaf", "is_leaf": True}
    assert G.nodes["hyb5"] == {"label": "hyb5", "type": "hyb_leaf", "is_leaf": True}
    assert G.nodes["-1"] == {"label": "-1", "type": "internal", "is_leaf": False}

    assert G.number_of_edges() == 5
    assert G.has_edge("-3", "-1")
    assert G.edges["-2", "-3"] == {"edge_type": "reticulation", "length": 0.5, "time_length": 0.0, "inher_weight": 0.3}
    assert G.edges["-1", "-3"] == {"edge_type": "tree", "length": 2.0, "time_length": 2.0, "inher_weight": 0.7}
    assert G.edges["-2", "sp4"] == {"edge_type": "tree", "length": 3.0, "time_length": 3.0}

    assert D.is_directed()
    assert D.has_edge("-1", "-3") and not D.has_edge("-3", "-1")
    assert D.edges["-2", "-3"] == G.edges["-2", "-3"]


# midpoints split into their parents; the leader is kept only when exactly two distinct nodes remain
@pytest.mark.parametrize("u, v, expected", [
    ("sp12", "hyb18", ("hyb18", "sp12")),
    ("hyb18", "sp12", ("hyb18", "sp12")),
    ("sp12_hyb18", "sp12", ("hyb18", "sp12")),
    ("sp12", "sp12_hyb18", ("hyb18", "sp12")),
    ("sp12_hyb18", "sp4", None),
    ("sp12_hyb18", "sp4_sp5", None),
])
def test_traced_pair(u, v, expected):
    assert MainResultAnalysis().traced_pair(u, v) == expected


# history leaders use |weight|, keep ties, and skip time points where every weight is zero
def test_history_leaders_uses_abs_weight_keeps_ties_and_skips_zero_points(tmp_path):
    rip = {"harmonic_cycles": [{"cycle_index": 0, "weight_history": [
        {"simplex": [0, 1], "history": [[1.0, 0.0], [2.0, -0.8], [3.0, 0.5]]},
        {"simplex": [1, 2], "history": [[1.0, 0.0], [2.0, 0.8], [3.0, 0.2]]},
        {"simplex": [0, 2], "history": [[1.0, 0.0], [2.0, 0.3], [3.0, -0.6]]},
    ]}]}
    sim_dir = tmp_path / "cycle_outputs" / "sim0"
    sim_dir.mkdir(parents=True)
    (sim_dir / "rip.json").write_text(json.dumps(rip))

    leaders = MainResultAnalysis(str(tmp_path)).history_leaders(0)

    assert leaders == {0: [
        (2.0, 0.8, [(0, 1), (1, 2)]),
        (3.0, 0.6, [(0, 2)]),
    ]}


# weight leaders use |weight_before_projection|, keep ties, and skip cycles where every weight is zero
def test_weight_leaders_uses_abs_weight_keeps_ties_and_skips_zero_cycles(tmp_path):
    rip = {"harmonic_cycles": [
        {"cycle_index": 0, "edges": [
            {"simplex": [0, 1], "weight_before_projection": -0.8},
            {"simplex": [1, 2], "weight_before_projection": 0.8},
            {"simplex": [0, 2], "weight_before_projection": 0.3},
        ]},
        {"cycle_index": 1, "edges": [
            {"simplex": [0, 1], "weight_before_projection": 0.0},
            {"simplex": [1, 2], "weight_before_projection": 0.0},
            {"simplex": [0, 2], "weight_before_projection": 0.0},
        ]},
    ]}
    sim_dir = tmp_path / "cycle_outputs" / "sim0"
    sim_dir.mkdir(parents=True)
    (sim_dir / "rip.json").write_text(json.dumps(rip))

    leaders = MainResultAnalysis(str(tmp_path)).weight_leaders(0)

    assert leaders == {0: (0.8, [(0, 1), (1, 2)])}
