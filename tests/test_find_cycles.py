import networkx as nx
from find_cycles import CycleFinder
from filter_cycle import FilterCycle


def _fc(min_cycle_length=0):
    cf = CycleFinder(G=nx.Graph(), threshold_mode=["cyclelength", "marker"], output_dir="/tmp/find_cycles_test", dist_matrix=None, index_to_name=None)
    return FilterCycle(cf, cycle_qualify_mode=[], min_cycle_length=min_cycle_length)


def _cycle(birth, edges):
    return {"birth": birth, "edges": [{"simplex": s, "weight": w} for s, w in edges]}


# only edges above the weight tolerance should contribute, and only names
# containing "hyb" or "sp" (skipping non-string labels) count as markers
def test_marker_nodes_filters_by_weight_and_name_pattern():
    fc = _fc()
    fc.index_to_name = {0: "sp0", 1: "internal", 2: "hyb1", 3: 42}
    cycle = {"edges": [
        {"simplex": [0, 1], "weight": 0.5},
        {"simplex": [1, 2], "weight": 0.005},
        {"simplex": [2, 3], "weight": -0.3},
    ]}

    assert fc._marker_nodes(cycle) == {"sp0", "hyb1"}


# a marker can only trigger qualification for the first cycle it appears in;
# a cycle skipped for being too small (<= min_cycle_length edges) must NOT
# mark its markers as seen, so a later cycle can still qualify with them
def test_generate_threshold_cycle_keys_marker_and_min_length_rules():
    fc = _fc(min_cycle_length=1)
    fc.index_to_name = {0: "root", 1: "sp1", 2: "sp2", 3: "hyb1", 4: "other"}
    fc.cycle_log = {"harmonic_cycles": [
        _cycle(1.0, [([0, 4], 0.5), ([4, 0], 0.5)]),  # no markers -> doesn't qualify
        _cycle(2.0, [([1, 0], 0.5), ([0, 4], 0.5)]),  # new marker sp1 -> qualifies
        _cycle(3.0, [([1, 4], 0.5), ([4, 0], 0.5)]),  # sp1 already seen -> doesn't qualify
        _cycle(4.0, [([3, 0], 0.5)]),                  # new marker hyb1, but only 1 edge
                                                         # (<= min_cycle_length=1) -> skipped,
                                                         # hyb1 NOT marked seen
        _cycle(5.0, [([3, 4], 0.5), ([4, 0], 0.5)]),   # hyb1 "new" again -> qualifies
    ]}

    fc.generate_threshold_cycle_keys()

    assert fc.thresholds == [2.0, 5.0]
    assert fc.qualifying_cycle_keys == [
        fc.cycle_log["harmonic_cycles"][1],
        fc.cycle_log["harmonic_cycles"][4],
    ]


# "fixed" mode should gate cycle_qualify by a pre-set thresholds list rather
# than deriving thresholds from cycle data, and must not grow that list
def test_generate_threshold_cycle_keys_fixed_mode_does_not_grow_thresholds():
    cf = CycleFinder(G=nx.Graph(), threshold_mode=["fixed"], output_dir="/tmp/find_cycles_test", dist_matrix=None, index_to_name=None, thresholds=[1])
    fc = FilterCycle(cf, cycle_qualify_mode=[])
    fc.cycle_log = {"harmonic_cycles": [
        _cycle(0.5, [([0, 1], 0.5)]),  # born before threshold 1 -> qualifies
        _cycle(2.0, [([0, 1], 0.5)]),  # born after threshold 1 -> doesn't qualify
    ]}

    fc.generate_threshold_cycle_keys()

    assert fc.thresholds == [1]
    assert fc.qualifying_cycle_keys == [fc.cycle_log["harmonic_cycles"][0]]
