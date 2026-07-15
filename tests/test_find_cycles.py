from find_cycles import _marker_nodes, _cycle_key, _select_thresholds


# only edges above the weight tolerance should contribute, and only names
# containing "hyb" or "sp" (skipping non-string labels) count as markers
def test_marker_nodes_filters_by_weight_and_name_pattern():
    index_to_name = {0: "sp0", 1: "internal", 2: "hyb1", 3: 42}
    cycle = {"edges": [
        {"simplex": [0, 1], "weight": 0.5},
        {"simplex": [1, 2], "weight": 0.005},  # below WEIGHT_ZERO_TOL, ignored
        {"simplex": [2, 3], "weight": -0.3},
    ]}

    assert _marker_nodes(cycle, index_to_name) == {"sp0", "hyb1"}


# the key should be a canonical, order-independent representation of a
# cycle's edges, so two equivalent edge lists hash/compare equal
def test_cycle_key_is_order_independent():
    drawn_edges = [([0, 1], 0.5), ([2, 3], -0.3)]

    key = _cycle_key(drawn_edges)
    reversed_key = _cycle_key(list(reversed(drawn_edges)))

    assert key == (((0, 1), 0.5), ((2, 3), -0.3))
    assert key == reversed_key


def _cycle(birth, edges):
    return {"birth": birth, "edges": [{"simplex": s, "weight": w} for s, w in edges]}


# a marker can only trigger qualification for the first cycle it appears in;
# a cycle skipped for being too small (<= min_cycle_length edges) must NOT
# mark its markers as seen, so a later cycle can still qualify with them
def test_select_thresholds_marker_and_min_length_rules():
    index_to_name = {0: "root", 1: "sp1", 2: "sp2", 3: "hyb1", 4: "other"}
    cycle_log = {"harmonic_cycles": [
        _cycle(1.0, [([0, 4], 0.5), ([4, 0], 0.5)]),  # no markers -> doesn't qualify
        _cycle(2.0, [([1, 0], 0.5), ([0, 4], 0.5)]),  # new marker sp1 -> qualifies
        _cycle(3.0, [([1, 4], 0.5), ([4, 0], 0.5)]),  # sp1 already seen -> doesn't qualify
        _cycle(4.0, [([3, 0], 0.5)]),                  # new marker hyb1, but only 1 edge
                                                         # (<= min_cycle_length=1) -> skipped,
                                                         # hyb1 NOT marked seen
        _cycle(5.0, [([3, 4], 0.5), ([4, 0], 0.5)]),   # hyb1 "new" again -> qualifies
    ]}

    thresholds, keys = _select_thresholds(cycle_log, index_to_name, min_cycle_length=1)

    assert thresholds == [2.0, 5.0]
    assert keys == {
        (((0, 4), 0.5), ((1, 0), 0.5)),
        (((3, 4), 0.5), ((4, 0), 0.5)),
    }
