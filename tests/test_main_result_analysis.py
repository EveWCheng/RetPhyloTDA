import json

import pytest

import main_result_analysis as mra

FILTERED_EDGES_HEADER = "from,to,edge_type,length,time_length,inher_weight\n"
NODES_HEADER = "id,label,type,is_leaf\n"


@pytest.fixture
def sim_dirs(tmp_path, monkeypatch):
    phylo = tmp_path / "phylo_csv"
    proc = tmp_path / "proc"
    phylo.mkdir()
    proc.mkdir()
    monkeypatch.setattr(mra, "PHYLO_CSV_DIR", str(phylo))
    monkeypatch.setattr(mra, "PROC_OUTPUTS_DIR", str(proc))
    return phylo, proc


def _write_filtered_edges(phylo, sim_id, rows):
    """rows: (from, to, edge_type, length, inher_weight-or-'')"""
    body = "".join(f"{u},{v},{et},{length},{length},{inher}\n" for u, v, et, length, inher in rows)
    (phylo / f"sim{sim_id}_filtered_edges.csv").write_text(FILTERED_EDGES_HEADER + body)


def _write_nodes(phylo, sim_id, rows):
    """rows: (id, label, type, is_leaf-bool)"""
    body = "".join(f"{i},{label},{t},{is_leaf}\n" for i, label, t, is_leaf in rows)
    (phylo / f"sim{sim_id}_nodes.csv").write_text(NODES_HEADER + body)


def _write_log(proc, sim_id, dress_entries):
    d = proc / f"sim{sim_id}"
    d.mkdir()
    (d / "find_cycles_log.json").write_text(json.dumps({"dress_distance_matrix": dress_entries}))


# ── read_leaf_labels ────────────────────────────────────────────────────────

def test_read_leaf_labels_uses_is_leaf_not_type(sim_dirs):
    phylo, _ = sim_dirs
    _write_nodes(phylo, 0, [
        (1, "sp1", "leaf", True),
        (2, "hyb2", "hyb_leaf", True),
        (3, "hyb3", "hyb_leaf", False),   # former hybrid tip that later speciated
        (4, "-4", "internal", False),
    ])
    assert mra.read_leaf_labels(0) == {"sp1", "hyb2"}


# ── read_network_graph ──────────────────────────────────────────────────────

def test_read_network_graph_builds_weighted_typed_leaf_tagged_graph(sim_dirs):
    phylo, _ = sim_dirs
    _write_filtered_edges(phylo, 0, [
        ("r", "A", "tree", 1.0, ""),
        ("m", "A", "reticulation", 2.0, "0.3"),
    ])
    _write_nodes(phylo, 0, [
        (1, "r", "internal", False),
        (2, "A", "leaf", True),
        (3, "m", "internal", False),
    ])
    G = mra.read_network_graph(0)
    assert not G.is_directed()
    assert G.edges["r", "A"] == {"length": 1.0, "edge_type": "tree"}
    assert G.edges["m", "A"]["edge_type"] == "reticulation"
    assert G.nodes["A"]["is_leaf"] is True
    assert G.nodes["r"]["is_leaf"] is False
    assert G.nodes["m"]["is_leaf"] is False


def test_read_network_graph_keeps_shorter_of_parallel_edges(sim_dirs):
    phylo, _ = sim_dirs
    _write_filtered_edges(phylo, 0, [
        ("x", "y", "tree", 5.0, ""),
        ("x", "y", "reticulation", 2.0, "0.4"),  # shorter -> wins
        ("x", "y", "tree", 9.0, ""),             # longer -> ignored
    ])
    _write_nodes(phylo, 0, [(1, "x", "leaf", True), (2, "y", "leaf", True)])
    G = mra.read_network_graph(0)
    assert G.edges["x", "y"]["length"] == 2.0
    assert G.edges["x", "y"]["edge_type"] == "reticulation"


def test_read_network_graph_missing_file_returns_none(sim_dirs):
    assert mra.read_network_graph(999) is None


# ── reticulate_edge_shortening ──────────────────────────────────────────────

def test_reticulate_edge_shortening_prepopulates_every_edge(sim_dirs):
    phylo, proc = sim_dirs
    _write_filtered_edges(phylo, 0, [
        ("a", "H1", "reticulation", 1.0, "0.3"),
        ("b", "H2", "reticulation", 1.0, "0.1"),
    ])
    _write_log(proc, 0, [])
    assert mra.reticulate_edge_shortening(0) == {
        frozenset({"a", "H1"}): (0, 0.0, None),
        frozenset({"b", "H2"}): (0, 0.0, None),
    }


def test_reticulate_edge_shortening_counts_pairs_and_takes_min_ratio(sim_dirs):
    phylo, proc = sim_dirs
    _write_filtered_edges(phylo, 0, [("a", "H1", "reticulation", 1.0, "0.3")])
    _write_log(proc, 0, [
        {"i": "L1", "j": "L2", "ratio": -0.2, "reticulation_edges_on_path": [{"edge": ["a", "H1"], "inher_weight": 0.3}]},
        {"i": "L3", "j": "L4", "ratio": -0.5, "reticulation_edges_on_path": [{"edge": ["H1", "a"], "inher_weight": 0.3}]},  # reversed, same edge
        {"i": "L5", "j": "L6", "ratio": -0.1, "reticulation_edges_on_path": []},
        {"i": "L7", "j": "L8", "ratio": None, "reticulation_edges_on_path": [{"edge": ["a", "H1"], "inher_weight": 0.3}]},   # None skipped
    ])
    # max_shortening -0.5 came from the (L3, L4) pair
    assert mra.reticulate_edge_shortening(0) == {frozenset({"a", "H1"}): (2, -0.5, ("L3", "L4"))}


def test_reticulate_edge_shortening_missing_log_returns_none(sim_dirs):
    phylo, _ = sim_dirs
    _write_filtered_edges(phylo, 0, [("a", "H1", "reticulation", 1.0, "0.3")])
    assert mra.reticulate_edge_shortening(0) is None


# ── reticulate_edge_table_lines ─────────────────────────────────────────────

def test_reticulate_edge_table_lines_no_log():
    assert mra.reticulate_edge_table_lines([], None, set()) == ["reticulate-edge table: no find_cycles log"]


def test_reticulate_edge_table_lines_sorted_by_shortening_and_flagged():
    retic_edges = [("a", "H1", 0.3), ("b", "H2", None)]
    shortening = {
        frozenset({"a", "H1"}): (4, -0.8, ("L3", "L4")),
        frozenset({"b", "H2"}): (0, 0.0, None),
    }
    on_path = {frozenset({"a", "H1"})}

    body = mra.reticulate_edge_table_lines(retic_edges, shortening, on_path)[2:]

    assert body[0].split() == ["a", "--", "H1", "4", "-0.800000", "L3", "--", "L4", "0.3000", "True"]
    assert body[1].split() == ["b", "--", "H2", "0", "0.000000", "-", "None", "False"]


# ── top_cycle_edge_paths ────────────────────────────────────────────────────

def _cycle(*edges):
    return {"edges": [{"nodes": [u, v], "weight": w} for u, v, w in edges]}


def test_top_cycle_edge_paths_ranks_by_abs_weight_and_traces_path(sim_dirs):
    phylo, _ = sim_dirs
    _write_filtered_edges(phylo, 0, [
        ("A", "m", "tree", 1.0, ""),
        ("m", "B", "reticulation", 1.0, "0.3"),
        ("A", "C", "tree", 4.0, ""),
    ])
    _write_nodes(phylo, 0, [
        (1, "A", "leaf", True), (2, "B", "leaf", True),
        (3, "C", "leaf", True), (4, "m", "internal", False),
    ])
    cycles = [_cycle(("A", "B", 1.0), ("A", "C", -0.5), ("A", "m", 0.1))]

    out = mra.top_cycle_edge_paths(0, cycles, top_k=2, tip_to_tip_only=False)

    assert [e["edge"] for e in out[0]] == [("A", "B"), ("A", "C")]  # by |weight|
    assert out[0][0]["path"] == [
        ("A", "m", "tree", 1.0),
        ("m", "B", "reticulation", 1.0),
    ]


def test_top_cycle_edge_paths_includes_edges_tied_with_the_kth(sim_dirs):
    phylo, _ = sim_dirs
    _write_filtered_edges(phylo, 0, [
        ("A", "B", "tree", 1.0, ""),
        ("C", "D", "tree", 1.0, ""),
        ("A", "C", "tree", 1.0, ""),
    ])
    _write_nodes(phylo, 0, [
        (1, "A", "leaf", True), (2, "B", "leaf", True),
        (3, "C", "leaf", True), (4, "D", "leaf", True),
    ])
    # two edges tied at |weight| 1.0, one at 0.3; top_k=1 -> both 1.0 edges
    cycles = [_cycle(("A", "B", 1.0), ("C", "D", -1.0), ("A", "C", 0.3))]

    out = mra.top_cycle_edge_paths(0, cycles, top_k=1, tip_to_tip_only=False)
    assert {e["edge"] for e in out[0]} == {("A", "B"), ("C", "D")}


def test_top_cycle_edge_paths_none_path_when_endpoint_off_graph(sim_dirs):
    phylo, _ = sim_dirs
    _write_filtered_edges(phylo, 0, [("A", "B", "tree", 1.0, "")])
    _write_nodes(phylo, 0, [(1, "A", "leaf", True), (2, "B", "leaf", True)])
    cycles = [_cycle(("A", "added_node_9", 1.0))]

    out = mra.top_cycle_edge_paths(0, cycles, top_k=1, tip_to_tip_only=False)
    assert out[0][0]["path"] is None


def test_top_cycle_edge_paths_tip_to_tip_only_skips_internal_endpoint_edges(sim_dirs):
    phylo, _ = sim_dirs
    _write_filtered_edges(phylo, 0, [
        ("A", "m", "tree", 1.0, ""),
        ("m", "B", "tree", 1.0, ""),
    ])
    _write_nodes(phylo, 0, [
        (1, "A", "leaf", True), (2, "B", "leaf", True), (3, "m", "internal", False),
    ])
    # highest |weight| edge touches internal node m; only A--B is leaf-to-leaf
    cycles = [_cycle(("A", "m", 1.0), ("A", "B", 0.2))]

    kept = mra.top_cycle_edge_paths(0, cycles, top_k=1, tip_to_tip_only=True)
    assert kept[0][0]["edge"] == ("A", "B")

    unfiltered = mra.top_cycle_edge_paths(0, cycles, top_k=1, tip_to_tip_only=False)
    assert unfiltered[0][0]["edge"] == ("A", "m")


# ── read_dress_change_ratios ────────────────────────────────────────────────

def test_read_dress_change_ratios_reads_ratio_field_and_skips_none(sim_dirs):
    _, proc = sim_dirs
    _write_log(proc, 0, [{"ratio": -0.5}, {"ratio": None}, {"ratio": 0.0}])
    assert mra.read_dress_change_ratios(0) == [-0.5, 0.0]


def test_read_dress_change_ratios_missing_file_returns_none(sim_dirs):
    assert mra.read_dress_change_ratios(999) is None
