import networkx as nx
from find_cycles import CycleFinder


def _cf(min_cycle_length=0):
    return CycleFinder(G=None, threshold_mode=["cyclelength", "marker"], cycle_qualify_mode=["marker"], output_dir="/tmp/find_cycles_test", min_cycle_length=min_cycle_length)


def _cycle(birth, edges):
    return {"birth": birth, "edges": [{"simplex": s, "weight": w} for s, w in edges]}


def _cf_with_graph(G, index_to_node):
    cf = CycleFinder(G=G, threshold_mode=[], cycle_qualify_mode=["crossover"], output_dir="/tmp/find_cycles_test")
    cf.index_to_node = index_to_node
    return cf


# only edges above the weight tolerance should contribute, and only names
# containing "hyb" or "sp" (skipping non-string labels) count as markers
def test_marker_nodes_filters_by_weight_and_name_pattern():
    cf = _cf()
    cf.index_to_name = {0: "sp0", 1: "internal", 2: "hyb1", 3: 42}
    cycle = {"edges": [
        {"simplex": [0, 1], "weight": 0.5},
        {"simplex": [1, 2], "weight": 0.005},
        {"simplex": [2, 3], "weight": -0.3},
    ]}

    assert cf._marker_nodes(cycle) == {"sp0", "hyb1"}


# a marker can only trigger qualification for the first cycle it appears in;
# a cycle skipped for being too small (<= min_cycle_length edges) must NOT
# mark its markers as seen, so a later cycle can still qualify with them
def test_generate_threshold_cycle_keys_marker_and_min_length_rules():
    cf = _cf(min_cycle_length=1)
    cf.index_to_name = {0: "root", 1: "sp1", 2: "sp2", 3: "hyb1", 4: "other"}
    cf.cycle_log = {"harmonic_cycles": [
        _cycle(1.0, [([0, 4], 0.5), ([4, 0], 0.5)]),  # no markers -> doesn't qualify
        _cycle(2.0, [([1, 0], 0.5), ([0, 4], 0.5)]),  # new marker sp1 -> qualifies
        _cycle(3.0, [([1, 4], 0.5), ([4, 0], 0.5)]),  # sp1 already seen -> doesn't qualify
        _cycle(4.0, [([3, 0], 0.5)]),                  # new marker hyb1, but only 1 edge
                                                         # (<= min_cycle_length=1) -> skipped,
                                                         # hyb1 NOT marked seen
        _cycle(5.0, [([3, 4], 0.5), ([4, 0], 0.5)]),   # hyb1 "new" again -> qualifies
    ]}

    cf.generate_threshold_cycle_keys()

    assert cf.thresholds == [2.0, 5.0]
    assert cf.qualifying_cycle_keys == [
        cf.cycle_log["harmonic_cycles"][1],
        cf.cycle_log["harmonic_cycles"][4],
    ]


# a cycle disqualifies as soon as any two of its nodes have disjoint sources
def test_qualifying_cycle_crossover_disjoint_sources_disqualifies():
    G = nx.Graph()
    G.add_node("n0", sources={1, 2})
    G.add_node("n1", sources={2, 3})
    G.add_node("n2", sources={3, 4})  # disjoint from n0's {1, 2}
    cf = _cf_with_graph(G, {0: "n0", 1: "n1", 2: "n2"})
    cycle = _cycle(1.0, [([0, 1], 0.5), ([1, 2], 0.5), ([2, 0], 0.5)])

    assert cf.qualifying_cycle_crossover(cycle) is False


# if every pair of nodes in the cycle shares at least one source, it qualifies
def test_qualifying_cycle_crossover_shared_source_qualifies():
    G = nx.Graph()
    G.add_node("n0", sources={1, 2})
    G.add_node("n1", sources={2, 3})
    G.add_node("n2", sources={2, 4})  # all three share source 2
    cf = _cf_with_graph(G, {0: "n0", 1: "n1", 2: "n2"})
    cycle = _cycle(1.0, [([0, 1], 0.5), ([1, 2], 0.5), ([2, 0], 0.5)])

    assert cf.qualifying_cycle_crossover(cycle) is True


# "fixed" mode should gate cycle_qualify by a pre-set thresholds list rather
# than deriving thresholds from cycle data, and must not grow that list
def test_generate_threshold_cycle_keys_fixed_mode_does_not_grow_thresholds():
    cf = CycleFinder(G=None, threshold_mode=["fixed"], cycle_qualify_mode=[], output_dir="/tmp/find_cycles_test", thresholds=[1])
    cf.cycle_log = {"harmonic_cycles": [
        _cycle(0.5, [([0, 1], 0.5)]),  # born before threshold 1 -> qualifies
        _cycle(2.0, [([0, 1], 0.5)]),  # born after threshold 1 -> doesn't qualify
    ]}

    cf.generate_threshold_cycle_keys()

    assert cf.thresholds == [1]
    assert cf.qualifying_cycle_keys == [cf.cycle_log["harmonic_cycles"][0]]
