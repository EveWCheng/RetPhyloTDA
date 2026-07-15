import pytest
import networkx as nx

from sim_bdh import (
    PhyloNetwork,
    SimState,
    _nchoose2,
    _collapse_hyb_nodes,
    _suppress_unary_nodes,
    _assign_labels,
    enumerate_gene_trees,
)


def test_nchoose2_counts_unordered_pairs():
    assert _nchoose2(0) == 0
    assert _nchoose2(1) == 0
    assert _nchoose2(2) == 1
    assert _nchoose2(4) == 6


def test_phylo_network_properties():
    G = nx.DiGraph()
    G.add_node(1, is_leaf=False)
    G.add_node(2, is_leaf=True)
    G.add_node(3, is_leaf=True, extinct=True)
    G.add_node(4, is_leaf=True, is_hyb_leaf=True)
    G.add_edges_from([(1, 2), (1, 3), (1, 4)])
    phy = PhyloNetwork(G=G, nleaves=2)

    assert sorted(phy.leaves) == [2, 3, 4]
    assert phy.extinct == [3]
    assert phy.hyb_tips == [4]
    assert phy.Nnode == 1


def test_gene_tree_filters_edges_by_gene_membership():
    G = nx.DiGraph()
    G.add_node(1, is_leaf=False)
    G.add_node(2, is_leaf=True)
    G.add_node(3, is_leaf=True)
    G.add_edge(1, 2, genes={0, 1})
    G.add_edge(1, 3, genes={1})
    phy = PhyloNetwork(G=G, nleaves=2)

    assert set(phy.gene_tree(0).edges()) == {(1, 2)}
    assert set(phy.gene_tree(1).edges()) == {(1, 2), (1, 3)}


def _hyb_test_graph():
    # root(1) -> sp2, ext3(extinct), hyb_node(4); hyb_node(4) also has a
    # second parent sp6, and a successor sp5 -- lets filter_nodes exercise
    # both extinct-removal and hyb-node collapsing in the same graph.
    G = nx.DiGraph()
    G.add_node(1, label="root", is_leaf=False)
    G.add_node(2, label="sp2", is_leaf=True)
    G.add_node(3, label="ext3", is_leaf=True, extinct=True)
    G.add_node(4, label="-4", is_leaf=False, is_hyb_node=True)
    G.add_node(5, label="sp5", is_leaf=True)
    G.add_node(6, label="sp6", is_leaf=True)
    G.add_edge(1, 2)
    G.add_edge(1, 3)
    G.add_edge(1, 4, length=1.0, time_length=1.0)
    G.add_edge(6, 4, length=2.0, time_length=2.0)
    G.add_edge(4, 5, length=3.0, time_length=3.0)
    return PhyloNetwork(G=G, nleaves=3)


# all_nodes mode should drop extinct nodes but leave hyb_node junctions intact
def test_filter_nodes_all_nodes_drops_extinct_but_keeps_hyb_node():
    phy = _hyb_test_graph()

    result = phy.filter_nodes(which_nodes="all_nodes")

    assert set(result.nodes()) == {1, 2, 4, 5, 6}
    assert set(result.edges()) == {(1, 2), (1, 4), (6, 4), (4, 5)}


# no_hyb_nodes mode should drop extinct nodes AND collapse hyb_node junctions,
# summing branch lengths across the splice like _collapse_hyb_nodes does
def test_filter_nodes_no_hyb_nodes_drops_extinct_and_collapses_hyb_node():
    phy = _hyb_test_graph()

    result = phy.filter_nodes(which_nodes="no_hyb_nodes")

    assert set(result.nodes()) == {1, 2, 5, 6}
    assert set(result.edges()) == {(1, 2), (1, 5), (6, 5)}
    assert result[1][5]["length"] == 4.0
    assert result[6][5]["length"] == 5.0


def test_collapse_hyb_nodes_splices_predecessors_to_successors():
    G = nx.DiGraph()
    G.add_node("p1")
    G.add_node("p2")
    G.add_node("h", is_hyb_node=True)
    G.add_node("c")
    G.add_edge("p1", "h", length=1.0, time_length=1.0)
    G.add_edge("p2", "h", length=2.0, time_length=2.0)
    G.add_edge("h", "c", length=3.0, time_length=3.0)

    result = _collapse_hyb_nodes(G)

    assert "h" not in result.nodes()
    assert result["p1"]["c"]["length"] == 4.0
    assert result["p2"]["c"]["length"] == 5.0


def test_suppress_unary_nodes_splices_out_middle_node():
    tree = nx.DiGraph()
    tree.add_node("root", is_leaf=False)
    tree.add_node("A", is_leaf=True)
    tree.add_node("mid", is_leaf=False)
    tree.add_node("leafB", is_leaf=True)
    tree.add_edge("root", "A")
    tree.add_edge("root", "mid")
    tree.add_edge("mid", "leafB")

    result = _suppress_unary_nodes(tree)

    assert "mid" not in result.nodes()
    assert set(result.edges()) == {("root", "A"), ("root", "leafB")}


def test_suppress_unary_nodes_removes_non_leaf_dead_ends():
    tree = nx.DiGraph()
    tree.add_node("root", is_leaf=False)
    tree.add_node("leafA", is_leaf=True)
    tree.add_node("leafB", is_leaf=True)
    tree.add_node("dangling", is_leaf=False)
    tree.add_edge("root", "leafA")
    tree.add_edge("root", "leafB")
    tree.add_edge("root", "dangling")

    result = _suppress_unary_nodes(tree)

    assert "dangling" not in result.nodes()
    assert set(result.edges()) == {("root", "leafA"), ("root", "leafB")}


def test_suppress_unary_nodes_root_chain_collapses_with_no_reattachment():
    # sharp edge: when the root itself ends up unary (in-degree 0,
    # out-degree 1), it is deleted with no reattachment step -- unlike an
    # internal unary node, which gets spliced out via add_edge(pred, succ).
    # A root -> mid -> leaf chain therefore collapses all the way down to a
    # single isolated leaf node with zero edges, rather than leaving 'leaf'
    # standing alone as the new root. This only shows up when the root has
    # no surviving sibling branch to keep its out-degree >= 2 (see the
    # splice-out-middle-node test above, where 'A' keeps root branching).
    tree = nx.DiGraph()
    tree.add_node("root", is_leaf=False)
    tree.add_node("mid", is_leaf=False)
    tree.add_node("leaf", is_leaf=True)
    tree.add_edge("root", "mid")
    tree.add_edge("mid", "leaf")

    result = _suppress_unary_nodes(tree)

    assert set(result.nodes()) == {"leaf"}
    assert set(result.edges()) == set()


def test_assign_labels_uses_correct_prefix_per_node_kind():
    G = nx.DiGraph()
    G.add_node(1, is_leaf=False)
    G.add_node(2, is_leaf=True)
    G.add_node(3, is_leaf=True, extinct=True)
    G.add_node(4, is_leaf=True, is_hyb_leaf=True)

    _assign_labels(G)

    assert G.nodes[1]["label"] == "-1"
    assert G.nodes[2]["label"] == "sp2"
    assert G.nodes[3]["label"] == "ext3"
    assert G.nodes[4]["label"] == "hyb4"



def test_enumerate_gene_trees_degenerative_absorbs_secondary_with_no_trace():
    G = nx.DiGraph()
    G.add_node(1, is_leaf=False)
    G.add_node(2, is_leaf=False)
    G.add_node(3, is_leaf=True)
    G.add_node(4, is_leaf=False, is_hyb_node=True)
    G.add_node(5, is_leaf=False)
    G.add_node(6, is_leaf=True, is_hyb_leaf=True)
    G.add_edge(1, 2)
    G.add_edge(1, 3)
    G.add_edge(2, 4, inher_weight=0.7)
    G.add_edge(2, 5)
    G.add_edge(4, 6)
    G.add_edge(5, 4, inher_weight=0.3)

    results = enumerate_gene_trees(G, n_samples="all")

    assert len(results) == 2
    by_weight = {round(w, 5): tree for tree, w in results}
    assert set(by_weight.keys()) == {0.3, 0.7}

    # both resolutions collapse to the identical tree -- node 5 leaves no
    # trace either way, since it never had its own continuation leaf
    for tree in by_weight.values():
        assert set(tree.nodes()) == {1, 3, 6}
        assert set(tree.edges()) == {(1, 3), (1, 6)}


def test_enumerate_gene_trees_neutral_preserves_both_parents_as_distinct_leaves():
    G = nx.DiGraph()
    G.add_node(1, is_leaf=False)
    G.add_node(2, is_leaf=False)
    G.add_node(3, is_leaf=True)
    G.add_node(4, is_leaf=False, is_hyb_node=True)
    G.add_node(5, is_leaf=False)
    G.add_node(6, is_leaf=True, is_hyb_leaf=True)
    G.add_node(7, is_leaf=True)
    G.add_edge(1, 2)
    G.add_edge(1, 3)
    G.add_edge(2, 4, inher_weight=0.7)
    G.add_edge(2, 5)
    G.add_edge(4, 6)
    G.add_edge(5, 4, inher_weight=0.3)
    G.add_edge(5, 7)

    results = enumerate_gene_trees(G, n_samples="all")

    assert len(results) == 2
    by_weight = {round(w, 5): tree for tree, w in results}
    assert set(by_weight.keys()) == {0.3, 0.7}

    tree_a = by_weight[0.7]  # edge (2,4) kept: node 2 survives as hub, node 5 pruned
    assert set(tree_a.nodes()) == {1, 2, 3, 6, 7}
    assert set(tree_a.edges()) == {(1, 2), (1, 3), (2, 6), (2, 7)}

    tree_b = by_weight[0.3]  # edge (5,4) kept: node 5 survives as hub, node 2 pruned
    assert set(tree_b.nodes()) == {1, 3, 5, 6, 7}
    assert set(tree_b.edges()) == {(1, 3), (1, 5), (5, 6), (5, 7)}


def test_enumerate_gene_trees_generating_keeps_both_parents_and_a_new_hyb_node():
    G = nx.DiGraph()
    G.add_node(1, is_leaf=False)
    G.add_node(2, is_leaf=False)
    G.add_node(3, is_leaf=True)
    G.add_node(4, is_leaf=False)
    G.add_node(5, is_leaf=False)
    G.add_node(6, is_leaf=True)
    G.add_node(7, is_leaf=True)
    G.add_node(8, is_leaf=False, is_hyb_node=True)
    G.add_node(9, is_leaf=True, is_hyb_leaf=True)
    G.add_edge(1, 2)
    G.add_edge(1, 3)
    G.add_edge(2, 4)
    G.add_edge(2, 5)
    G.add_edge(4, 6)
    G.add_edge(4, 8, inher_weight=0.7)
    G.add_edge(5, 7)
    G.add_edge(5, 8, inher_weight=0.3)
    G.add_edge(8, 9)

    results = enumerate_gene_trees(G, n_samples="all")

    assert len(results) == 2
    by_weight = {round(w, 5): tree for tree, w in results}
    assert set(by_weight.keys()) == {0.3, 0.7}

    tree_a = by_weight[0.7]  # edge (4,8) kept: node 4 survives as hub, node 5 pruned
    assert set(tree_a.nodes()) == {1, 2, 3, 4, 6, 7, 9}
    assert set(tree_a.edges()) == {(1, 2), (1, 3), (2, 4), (2, 7), (4, 6), (4, 9)}

    tree_b = by_weight[0.3]  # edge (5,8) kept: node 5 survives as hub, node 4 pruned
    assert set(tree_b.nodes()) == {1, 2, 3, 5, 6, 7, 9}
    assert set(tree_b.edges()) == {(1, 2), (1, 3), (2, 5), (2, 6), (5, 7), (5, 9)}


# mrca=True should start the sim with a root and two active leaf lineages
def test_init_mrca_true_creates_root_with_two_leaves():
    state = SimState(mrca=True, Ngene=2)

    assert state.leaves == {2, 3}
    assert state.G.nodes[1]["is_leaf"] is False
    assert state.G.nodes[2]["is_leaf"] is True
    assert state.G.nodes[3]["is_leaf"] is True
    assert set(state.G.edges()) == {(1, 2), (1, 3)}
    assert state.G[1][2]["genes"] == {0, 1}
    assert state.G[1][3]["genes"] == {0, 1}


# mrca=False should start the sim with just a root and a single leaf lineage
def test_init_mrca_false_creates_root_with_one_leaf():
    state = SimState(mrca=False, Ngene=0)

    assert state.leaves == {2}
    assert state.G.number_of_nodes() == 2
    assert set(state.G.edges()) == {(1, 2)}


# a trait_model missing any of the required callback/initial keys should be rejected up front
def test_init_trait_model_requires_all_keys():
    with pytest.raises(ValueError):
        SimState(mrca=True, trait_model={"initial": [0, 0]})


# trait_model['initial'] must supply one value per starting lineage (2 for mrca=True)
def test_init_trait_model_initial_length_must_match_mrca():
    trait_model = {
        "initial": [0],  # mrca=True needs 2 starting values, not 1
        "time_fxn": lambda x, t: x,
        "spec_fxn": lambda x: (x, x),
        "hyb_event_fxn": lambda *a: 0,
        "hyb_compatibility_fxn": lambda *a: True,
    }
    with pytest.raises(ValueError):
        SimState(mrca=True, trait_model=trait_model)


# speciation should seal the parent's incoming edge and replace it with two active leaf children
def test_speciation_seals_parent_and_adds_two_leaf_children():
    state = SimState(mrca=False, Ngene=2)

    state.speciation(2)

    assert state.leaves == {3, 4}
    assert state.G.nodes[2]["is_leaf"] is False
    assert set(state.G.successors(2)) == {3, 4}
    assert state.G[2][3]["genes"] == {0, 1}
    assert state.G[2][4]["genes"] == {0, 1}


# extinction should seal the parent edge, flag the node extinct, and drop it from active leaves
def test_extinction_seals_and_marks_extinct():
    state = SimState(mrca=False, Ngene=0)
    state.time = 3.0

    state.extinction(2)

    assert state.leaves == set()
    assert state.G.nodes[2]["is_leaf"] is False
    assert state.G.nodes[2]["extinct"] is True
    assert state.G[1][2]["length"] == 3.0


# for two still-active (unsealed) leaves, tip_distance should sum the elapsed time on each branch
def test_tip_distance_sums_elapsed_time_for_active_leaves():
    state = SimState(mrca=True, Ngene=0)
    state.time = 5.0  # neither leaf's incoming edge has been sealed yet

    assert state.tip_distance(2, 3) == 10.0  # 5.0 elapsed on each branch


# once an edge is sealed, tip_distance must use its frozen length, not further elapsed time
def test_tip_distance_uses_sealed_length_not_elapsed_time():
    state = SimState(mrca=True, Ngene=0)
    state.time = 3.0
    state.speciation(2)  # seals edge (1,2) at length 3.0; new leaves 4,5 replace it
    state.time = 10.0    # time keeps moving, but (1,2)'s length must stay frozen at 3.0

    # path 4 -> 2 -> 1 -> 3: (2,4) elapsed since creation (10-3=7), (1,2) sealed
    # at 3.0 regardless of the extra elapsed time, (1,3) elapsed since root (10-0=10)
    assert state.tip_distance(4, 3) == 7.0 + 3.0 + 10.0


def _speciated_state():
    # root(1) -> leaf1(2), leaf2(3); leaf1 speciates into 4,5.
    # leaf2(3) is left untouched so root keeps two children throughout.
    state = SimState(mrca=True, Ngene=0)
    state.speciation(2)
    return state


# hyb_generating should create a brand-new hyb_node fed by both parents, each parent
# also keeping its own separate continuation leaf -- the shape assumed by the
# "generating" enumerate_gene_trees test above
def test_hyb_generating_matches_shape_used_in_enumerate_gene_trees_tests():
    state = _speciated_state()

    state.hyb_generating(4, 5, inher=0.3)

    assert state.leaves == {3, 6, 7, 9}
    assert set(state.G.successors(4)) == {6, 8}  # own leaf (6) + new hyb_node (8)
    assert set(state.G.successors(5)) == {7, 8}  # own leaf (7) + new hyb_node (8)
    assert set(state.G.successors(8)) == {9}
    assert state.G[4][8]["inher_weight"] == pytest.approx(0.7)
    assert state.G[5][8]["inher_weight"] == pytest.approx(0.3)
    assert state.G[5][8]["edge_type"] == "reticulation"
    assert state.G.nodes[8]["is_hyb_node"] is True
    assert state.G.nodes[9]["is_hyb_leaf"] is True


# hyb_degenerative should make the primary parent itself the hyb_node, with the
# secondary parent fully absorbed (no continuation leaf) -- the shape assumed by
# the "degenerative" enumerate_gene_trees test above
def test_hyb_degenerative_matches_shape_used_in_enumerate_gene_trees_tests():
    state = _speciated_state()

    state.hyb_degenerative(4, 5, inher=0.3)

    assert state.leaves == {3, 6}  # secondary (5) gets no continuation leaf
    assert set(state.G.successors(4)) == {6}
    assert set(state.G.successors(5)) == {4}  # secondary feeds straight into primary
    assert state.G[2][4]["inher_weight"] == pytest.approx(0.7)
    assert state.G[5][4]["inher_weight"] == pytest.approx(0.3)
    assert state.G[5][4]["edge_type"] == "reticulation"
    assert state.G.nodes[4]["is_hyb_node"] is True  # primary itself becomes the hyb_node
    assert state.G.nodes[6]["is_hyb_leaf"] is True


# hyb_neutral should make the primary parent the hyb_node like hyb_degenerative,
# but give the secondary parent its own donor leaf -- the shape assumed by the
# "neutral" enumerate_gene_trees test above
def test_hyb_neutral_matches_shape_used_in_enumerate_gene_trees_tests():
    state = _speciated_state()

    state.hyb_neutral(4, 5, inher=0.3)

    assert state.leaves == {3, 6, 7}  # secondary (5) keeps its own donor leaf (7)
    assert set(state.G.successors(4)) == {6}
    assert set(state.G.successors(5)) == {4, 7}  # donor leaf + reticulation into primary
    assert state.G[2][4]["inher_weight"] == pytest.approx(0.7)
    assert state.G[5][4]["inher_weight"] == pytest.approx(0.3)
    assert state.G[5][4]["edge_type"] == "reticulation"
    assert state.G.nodes[4]["is_hyb_node"] is True
    assert state.G.nodes[6]["is_hyb_leaf"] is True
