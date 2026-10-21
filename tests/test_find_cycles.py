import networkx as nx
from find_cycles import CycleFinder
from filter_cycle import FilterCycle


def _fc(min_cycle_length=0):
    cf = CycleFinder(G=nx.Graph(), threshold_mode=["cyclelength", "marker"], output_dir="/tmp/find_cycles_test", dist_matrix=None, index_to_name=None)
    return FilterCycle(cf, cycle_qualify_mode=[], min_cycle_length=min_cycle_length)


def _cycle(birth, edges):
    return {"birth": birth, "edges": [{"simplex": s, "weight": w} for s, w in edges]}


# only edges with |weight| above WEIGHT_ZERO_TOL count towards the cycle length
def test_select_threshold_cyclelength_counts_edges_above_weight_tol():
    cycle = _cycle(1.0, [
        ([0, 1], 0.5),
        ([1, 2], -0.3),
        ([2, 3], 0.01),
        ([3, 4], 0.005),
        ([4, 0], 0.0),
    ])
    expected_edge_count = 2

    fc_at_count = _fc(min_cycle_length=expected_edge_count)
    fc_at_count.WEIGHT_ZERO_TOL = 0.01
    fc_below_count = _fc(min_cycle_length=expected_edge_count - 1)
    fc_below_count.WEIGHT_ZERO_TOL = 0.01

    assert fc_below_count._select_threshold_cyclelength(cycle) is True
    assert fc_at_count._select_threshold_cyclelength(cycle) is False
