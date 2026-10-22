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


# the CSVs should come back as a directed label-keyed graph with typed node and edge attributes
def test_read_filtered_graph_builds_directed_typed_graph(tmp_path):
    sim_csv_dir = tmp_path / "phylo_csv" / "sim0"
    sim_csv_dir.mkdir(parents=True)
    (sim_csv_dir / "sim0_filtered_nodes.csv").write_text(NODES_CSV)
    (sim_csv_dir / "sim0_filtered_edges.csv").write_text(EDGES_CSV)

    D = MainResultAnalysis(str(tmp_path)).read_filtered_graph(0)

    assert D.is_directed()
    assert set(D.nodes()) == {"-1", "-2", "-3", "sp4", "hyb5"}
    assert D.nodes["sp4"] == {"label": "sp4", "type": "leaf", "is_leaf": True}
    assert D.nodes["hyb5"] == {"label": "hyb5", "type": "hyb_leaf", "is_leaf": True}
    assert D.nodes["-1"] == {"label": "-1", "type": "internal", "is_leaf": False}

    assert D.number_of_edges() == 5
    assert D.has_edge("-1", "-3") and not D.has_edge("-3", "-1")
    assert D.edges["-2", "-3"] == {"edge_type": "reticulation", "length": 0.5, "time_length": 0.0, "inher_weight": 0.3}
    assert D.edges["-1", "-3"] == {"edge_type": "tree", "length": 2.0, "time_length": 2.0, "inher_weight": 0.7}
    assert D.edges["-2", "sp4"] == {"edge_type": "tree", "length": 3.0, "time_length": 3.0}


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

    leaders = MainResultAnalysis(str(tmp_path))._history_leaders_(0)

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

    assert leaders == {0: [(0.8, [(0, 1), (1, 2)])]}


def _ancestry_network():
    """Directed test network: sp1 and sp2 meet at -2, with three reticulations joining their lineages and three that should not count."""
    D = nx.DiGraph()
    tree_edges = [
        ("-1", "-2"), ("-1", "-3"),
        ("-2", "-4"), ("-2", "-5"),
        ("-4", "-6"), ("-4", "sp7"),
        ("-6", "sp1"),
        ("-5", "-8"), ("-5", "sp3"),
        ("-8", "sp2"),
        ("-3", "sp9"), ("-3", "sp10"),
    ]
    for a, b in tree_edges:
        D.add_edge(a, b, edge_type="tree")
    reticulations = [
        ("-6", "-8"),
        ("-4", "sp2"),
        ("-5", "sp1"),
        ("-3", "-6"),
        ("-2", "-8"),
        ("sp7", "-6"),
    ]
    for a, b in reticulations:
        D.add_edge(a, b, edge_type="reticulation")
    return D


# only reticulations joining sp1's lineage to sp2's lineage count, however far up each end is
def test_ancestor_reticulations_finds_every_reticulation_joining_the_two_lineages():
    D = _ancestry_network()
    analysis = MainResultAnalysis()

    expected = {("-6", "-8"), ("-4", "sp2"), ("-5", "sp1")}
    assert set(analysis.ancestor_reticulations(D, "sp1", "sp2")) == expected
    assert set(analysis.ancestor_reticulations(D, "sp2", "sp1")) == expected


NETWORK_NODES_CSV = """id,label,type,is_leaf
-1,-1,internal,False
-2,-2,internal,False
-3,-3,internal,False
sp1,sp1,leaf,True
sp2,sp2,leaf,True
sp4,sp4,leaf,True
sp5,sp5,leaf,True
"""

NETWORK_EDGES_CSV = """from,to,edge_type,length,time_length,inher_weight
-1,-2,tree,1.0,1.0,
-1,-3,tree,1.0,1.0,
-2,sp1,tree,1.0,1.0,
-2,sp4,tree,1.0,1.0,
-3,sp2,tree,1.0,1.0,
-3,sp5,tree,1.0,1.0,
-2,-3,reticulation,0.5,0.0,0.3
sp4,sp5,reticulation,0.5,0.0,0.3
"""


# leaders trace to node pairs, midpoint leaders that give three nodes are skipped, and pairs are grouped by reticulation
@pytest.mark.parametrize("sel_op, expected, expected_missed", [
    ("only_at_birth", {("-2", "-3"): [("sp1", "sp2")]}, [("sp4", "sp5")]),
    ("all", {("-2", "-3"): [("sp1", "sp2"), ("sp4", "sp5")], ("sp4", "sp5"): [("sp4", "sp5")]}, []),
])
def test_from_pairs_to_ret_each_sim_collects_reticulations_from_leader_pairs(tmp_path, sel_op, expected, expected_missed):
    sim_csv_dir = tmp_path / "phylo_csv" / "sim0"
    sim_csv_dir.mkdir(parents=True)
    (sim_csv_dir / "sim0_filtered_nodes.csv").write_text(NETWORK_NODES_CSV)
    (sim_csv_dir / "sim0_filtered_edges.csv").write_text(NETWORK_EDGES_CSV)

    proc_dir = tmp_path / "proc_phylo_outputs" / "sim0"
    proc_dir.mkdir(parents=True)
    (proc_dir / "populated_headers.txt").write_text("sp1\nsp2\nsp4\nsp1_-2\nsp5\n")

    rip = {"harmonic_cycles": [{"cycle_index": 0, "weight_history": [
        {"simplex": [0, 1], "history": [[1.0, 1.0], [2.0, 0.2]]},
        {"simplex": [2, 4], "history": [[1.0, 0.5], [2.0, 0.7]]},
        {"simplex": [1, 3], "history": [[1.0, 0.0], [2.0, -0.7]]},
    ]}]}
    sim_dir = tmp_path / "cycle_outputs" / "sim0"
    sim_dir.mkdir(parents=True)
    (sim_dir / "rip.json").write_text(json.dumps(rip))

    reticulations, missed = MainResultAnalysis(str(tmp_path), sel_op=sel_op).from_pairs_to_ret_each_sim(0)

    assert reticulations == expected
    assert missed == expected_missed
