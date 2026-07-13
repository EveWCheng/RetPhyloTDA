Gene tree resolution / distance matrix fixes (2026-07-13)
============================================================

Context
-------
tree_main.py's process_gene_trees() pipeline is:

    phy (network, sim_bdh.py)
      -> filter_nodes()            drop extinct tips, collapse is_hyb_node junctions
      -> enumerate_gene_trees()    resolve reticulation into single-topology gene trees
      -> networkx_to_tree_json()   compute each node's descendant-leaf-set "label"
      -> merge_trees()             build one grouping graph (merged_G) across all gene trees
      -> floyd_warshall_numpy()    distance matrix over merged_G, fed into harmonic_cycle/TDA

This session found and fixed several bugs along that chain, all stemming from
enumerate_gene_trees() not producing a clean, fully-resolved single tree.


1. index_to_name lost track of node identity (tree_main.py)
-------------------------------------------------------------
merged_G's nodes were bare tuples with no "label" attribute, so
merged_G.nodes(data="label") returned None for everything.

Fix: merge_trees() (tree_addition.py) now returns
nx.convert_node_labels_to_integers(builder.G, label_attribute="label")
so merged_G's nodes are integers 0..n-1 (matching floyd_warshall_numpy's
matrix indices directly), each carrying its original tuple identity under
"label". tree_main.py builds index_to_name = dict(merged_G.nodes(data="label")).


2. visualize() ignored the "label" attribute (tree_addition.py)
-------------------------------------------------------------
visualize() always displayed a node's raw internal id instead of its
computed "label" (e.g. the leaf-set list networkx_to_tree_groups just
assigned). Fixed to prefer attrs["label"] when present, falling back to
the node id otherwise (needed for the tuple-keyed builder.G call site,
which never carries a "label" attribute).


3. Unary nodes after reticulation resolution caused duplicate labels
-------------------------------------------------------------
(sim_bdh.py: _suppress_unary_nodes, called from enumerate_gene_trees)

enumerate_gene_trees() resolves each reticulation node (in_degree > 1) by
keeping exactly one incoming edge and dropping the rest. This can strand an
ancestor with only one surviving child (in_degree<=1, out_degree==1) -- not
a genuine branch point, just bookkeeping debris from the edge-drop. Since
networkx_to_tree_groups() labels a node by the union of its children's
labels, such a node ends up with the exact same label as its lone child
(e.g. two distinct nodes both labeled [1]), which later collided when
TreeBuilder used that label as node identity.

Fix: _suppress_unary_nodes() splices these out (reconnects predecessor
directly to successor) after each resolved tree is built, including the
unary-root case (in_degree==0, out_degree==1, e.g. a root whose only
child already covers the full leaf set).

Splice edges get weight=1 (the hop count is intentionally *not* preserved
across the collapsed edge -- an explicit, deliberate choice, not an
oversight).


4. dedupe compared the wrong representation (sim_bdh.py)
-------------------------------------------------------------
enumerate_gene_trees(..., dedupe=True) used to fingerprint each resolved
combo by its raw pre-splice edge set. Two different reticulation combos
can differ only in edges that get spliced away by _suppress_unary_nodes,
collapsing into the *same* simplified tree post-splice while still being
treated as distinct pre-splice -- so true duplicates went undetected.
Fixed by computing the dedupe key from the tree's edges *after*
_suppress_unary_nodes() runs.

Note: dedupe is still off by default, and even with dedupe=True the
function does not resample to replace a discarded duplicate, so
len(results) can come out smaller than n_samples.


5. enumerate_gene_trees() discarded all node attributes
-------------------------------------------------------------
Each resolved tree used to be built as:
    tree = nx.DiGraph(); tree.add_edges_from(kept)
add_edges_from() with bare (u, v) tuples creates brand new nodes with
empty attribute dicts -- is_leaf, timecreation, label, etc. from the
simulation were silently lost for every node in every resolved tree.

Fix: tree = G.edge_subgraph(kept).copy(), which keeps each node's and
edge's original attribute dict from the source graph.


6. Fake leaves from fully-dropped subtrees (sim_bdh.py)
-------------------------------------------------------------
This was only detectable once fix #5 preserved is_leaf (see below).

_suppress_unary_nodes() only removed nodes left with out_degree==1. A
node can also lose *all* of its children in one combo -- e.g. a node
whose two children are themselves both reticulation nodes, and both of
this node's edges into them happen to lose their respective reticulation
draws. That node ends up with out_degree==0 despite is_leaf=False, and
was being silently treated as a genuine tip by every downstream
leaf-counting step (out_degree==0 was used as the leaf test), inflating
the leaf count and fabricating a tip that was never actually sampled.

Confirmed empirically: from the same filtered_G (8 real leaves), some
resolved combos produced 9 "leaves" -- an internal node (is_leaf=False)
counted as a tip purely because it lost every child in that combo.

Fix: _suppress_unary_nodes() now also prunes any node with
out_degree==0 that is not flagged is_leaf, before re-scanning for unary
nodes each pass -- removing one dead end can itself turn its parent into
a new unary node or a new dead end, which the loop's next pass catches.


7. SimParams gained a second stopping criterion
-------------------------------------------------------------
SimParams.stopping_num_leaves (default None) stops _sim_one()'s
Gillespie loop as soon as len(state.leaves) reaches that count, in
addition to the existing age-based stop. Checked at the top of the loop,
before rates/timestep are drawn for that iteration -- mirrors how the
age check stops before firing the next event, just without capping
state.time to anything (that's specific to the age criterion).
tree_main.py sets STOPPING_NUM_LEAVES = 8 and passes it through.


Still open
-------------------------------------------------------------
- networkx_to_tree_groups() (tree_addition.py) writes its own
  leaf-descendant-set label into the same "label" attribute the
  simulation already uses for a descriptive string (e.g. 'sp7', 'hyb8').
  Its "if 'label' already present, skip" guard assumes label can only
  have been set by its own earlier recursion -- true before fix #5, no
  longer true now that nodes arrive with a pre-existing simulation
  label. This causes internal nodes to keep their stale string label
  instead of the intended leaf-index list, which crashes TreeBuilder
  downstream ('str' object has no attribute 'remove'). A fix was drafted
  (clear any inherited "label" at the top of networkx_to_tree_groups
  before computing its own) but was reverted pending a decision on
  whether to do that, or give the two schemes separate attribute names.
