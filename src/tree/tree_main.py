import os
import re
import shutil
import sys
import numpy as np
import networkx as nx
from tqdm.std import TqdmDefaultWriteLock

import cProfile, pstats
from pstats import SortKey

HERE = os.path.dirname(os.path.abspath(__file__))
SHARED_DIR = os.path.join(HERE, os.pardir, "shared")
if SHARED_DIR not in sys.path:
    sys.path.insert(0, SHARED_DIR)

TqdmDefaultWriteLock.mp_lock = None

from sim_bdh import SimState, SimParams, _sim_one, enumerate_gene_trees
from export import export_csv, export_filtered
from network_lab_tda.tree_edit.tree_addition import networkx_to_tree_json, merge_trees, visualize
from find_cycles import CycleFinder
from filter_cycle import FilterCycle
from tree_filter_cycle import TreeFilterCycle
from shared_node_utils import filter_shared_nodes_by_spread
from polymorphic_edges import write_polymorphic_edges, write_reticulate_edges

TREE_GROUP_OUTPUTS_DIR = os.path.join(HERE, os.pardir, "outputs", "tree_group_outputs")
TREE_GROUPS_DIR = os.path.join(TREE_GROUP_OUTPUTS_DIR, "tree_groups")
MERGED_TREE_DIR = os.path.join(TREE_GROUP_OUTPUTS_DIR, "merged_tree")
PHYLO_CSV_DIR = os.path.join(TREE_GROUP_OUTPUTS_DIR, "phylo_csv")

# ── Parameters ────────────────────────────────────────────────────────────────

# simulation time horizon (tree age)
AGE      = 4
# whether to condition on/track most recent common ancestor
MRCA     = True
# speciation rate
LAMBDA   = 0.5
# extinction rate
MU       = 0.08
# hybridization rate
NU       = 0.01
HYBPROPS = [1, 0,0]   # [lineage generating, degenerative, neutral]
# stop simulation once tree reaches this many leaves
STOPPING_NUM_LEAVES = 40
# minimum cycle length filter for TDA cycle detection
MIN_CYCLE_LENGTH = 4
# cap on how many qualifying cycles get rendered as HTML; None means render all of them
MAX_PLOT_CYCLES = 30
# number of gene trees to simulate
Ngene = 200
# trait evolution model (none used)
TRAIT_MODEL = None
# number of gene-tree samples to enumerate, "all" or an integer -- "all" enumerates
# the full cartesian product of reticulate-node choices (2^(number of reticulation
# events)), which grows exponentially with tree size and can blow up long before
# STOPPING_NUM_LEAVES gets large
N_SAMPLES = Ngene
# fixed thresholds to plot cycles at; set to None to derive thresholds dynamically via THRESHOLD_MODE instead
THRESHOLDS = [1]
# threshold-selection strategies CycleFinder runs per cycle
# available options: "cyclelength", "marker", "fixed"; [] means every cycle passes (no filtering)
THRESHOLD_MODE = ["fixed"]
# cycle-qualification strategies CycleFinder runs per cycle
# available options: "marker", "crossover"; [] means every selected cycle qualifies (no filtering)
CYCLE_QUALIFY_MODE = []
# list of units print_most_shared_units reports on, one output file per entry
# available options: "edge", "node"
SHARING_UNIT = ["edge", "node"]
# if True, sharing_edge_frequency/sharing_nodes_frequency skip edges G marks as "true_edge"
# (edges present in every merged input tree, per tree_addition.add_G_edge)
DELETE_TRUE_EDGES = False
# shared_nodes_all.txt lines whose max_tip_spread (hop count on the reticulation-free
# tree backbone) is smaller than this are copied into shared_nodes_all_filtered.txt
MAX_SHARED_NODE_SPREAD = 3


# draws hybrid inheritance probability
hyb_inher_fxn = lambda: np.random.uniform(0, 1)
# function for hybridization rate (unset, uses default NU)
hyb_rate_fxn  = None
filter_option = ["tree_by_tree_deletion",0]


# ── Run ───────────────────────────────────────────────────────────────────────

def process_gene_trees(phy, which_nodes: str = "no_hyb_nodes"):
    filtered_G = phy.filter_nodes(which_nodes=which_nodes)
#    print(filtered_G.nodes(data="is_hyb_node"))
    os.makedirs(TREE_GROUP_OUTPUTS_DIR, exist_ok=True)
    visualize(filtered_G, output=os.path.join(TREE_GROUP_OUTPUTS_DIR, "filtered_G.html"))
    enumerated_trees = enumerate_gene_trees(filtered_G,n_samples=N_SAMPLES)
    os.makedirs(TREE_GROUPS_DIR, exist_ok=True)
    for i, (tree, count) in enumerate(enumerated_trees.items()):
        networkx_to_tree_json(tree, os.path.join(TREE_GROUPS_DIR, f"gene_tree_{i}.json"), numeric_labels=True)
        visualize(tree, output=os.path.join(TREE_GROUPS_DIR, f"gene_tree_{i}.html"))

    os.makedirs(MERGED_TREE_DIR, exist_ok=True)
    merged_G = merge_trees(input_dir=TREE_GROUPS_DIR, output_dir=MERGED_TREE_DIR)
    return merged_G,enumerated_trees 


def leaf_labels_by_number(filtered_G):
    """{leaf number: original label (e.g. 'hyb53', 'sp51')} -- merging strips
    labels down to bare numbers, so this is captured before that happens."""
    labels = {}
    for n, attrs in filtered_G.nodes(data=True):
        if not attrs.get("is_leaf"):
            continue
        match = re.search(r"\d+", str(attrs.get("label", "")))
        if match:
            labels[int(match.group())] = attrs["label"]
    return labels


def build_dist_matrix(G, weight_attr, output_path, populated_header_fn="populated_headers.txt"):
    # graphs whose edges lack the weight_attr (e.g. merged_G) fall back to weight 1
    G_undirected = G.to_undirected() if G.is_directed() else G
    dist_matrix = nx.floyd_warshall_numpy(G_undirected, weight=weight_attr)
    index_to_name = {
        i: attrs.get("label", node)
        for i, (node, attrs) in enumerate(G_undirected.nodes(data=True))
    }
    # main_result_analysis.read_index_to_name() reads this back to name simplices
    with open(os.path.join(output_path, populated_header_fn), "w") as f:
        f.write("\n".join(str(v) for v in index_to_name.values()))
    np.savetxt(os.path.join(output_path, "populated_distance_matrix.txt"), dist_matrix)
    return G_undirected, dist_matrix, index_to_name


def find_cycles_in_merged_tree(merged_G,enumerated_trees, filter_option, leaf_labels):
    output_path = os.path.join(TREE_GROUP_OUTPUTS_DIR, "proc_phylo_outputs")
    os.makedirs(output_path, exist_ok=True)
    G_undirected, dist_matrix, index_to_name = build_dist_matrix(merged_G, "length", output_path)
    cf = CycleFinder(G_undirected, threshold_mode=THRESHOLD_MODE, output_dir=TREE_GROUP_OUTPUTS_DIR, dist_matrix=dist_matrix, index_to_name=index_to_name, thresholds=THRESHOLDS)
    cf.find_cycles()
    fc = FilterCycle(cf, cycle_qualify_mode=CYCLE_QUALIFY_MODE, min_cycle_length=MIN_CYCLE_LENGTH, vis=True, max_plot_cycles=MAX_PLOT_CYCLES)
    fc.visualize()
    tfc = TreeFilterCycle(fc, sharing_unit=SHARING_UNIT, delete_true_edges=DELETE_TRUE_EDGES)
    if filter_option[0] == "nodes_only":
        tfc.print_most_shared_units()
    elif filter_option[0]== "tree_by_tree_deletion":
        tfc.tree_by_tree_delete(enumerated_trees,filter_option[1], leaf_labels)
    return cf.cycle_output_path


def main(seed=43, which_nodes: str = "no_hyb_nodes"):
    if seed is not None:
        np.random.seed(seed)

    params = SimParams(age=AGE, lambda_=LAMBDA, mu=MU, nu=NU, hybprops=HYBPROPS,hyb_inher_fxn=hyb_inher_fxn,hyb_rate_fxn=hyb_rate_fxn,stopping_num_leaves=STOPPING_NUM_LEAVES)
    state = SimState(mrca=MRCA, Ngene=Ngene, trait_model=TRAIT_MODEL)
    phy = _sim_one(state,params)['phy']
    print("original tree is calcualted")
    if phy != 0:
        if os.path.exists(TREE_GROUP_OUTPUTS_DIR):
            shutil.rmtree(TREE_GROUP_OUTPUTS_DIR)
        os.makedirs(TREE_GROUP_OUTPUTS_DIR, exist_ok=True)
        export_csv(phy, PHYLO_CSV_DIR, prefix="sim0_")
        filtered_G = phy.filter_nodes(which_nodes=which_nodes)
        export_filtered(filtered_G, PHYLO_CSV_DIR, prefix="sim0_")
        leaf_labels = leaf_labels_by_number(filtered_G)
        merged_G,enumerated_trees  = process_gene_trees(phy, which_nodes=which_nodes)
        print("the gene trees have been merged")
        cycle_output_path = find_cycles_in_merged_tree(merged_G,enumerated_trees, filter_option, leaf_labels)
        write_polymorphic_edges(enumerated_trees, cycle_output_path, leaf_labels, prefix="sim0_")
        write_reticulate_edges(enumerated_trees, cycle_output_path, leaf_labels, prefix="sim0_")
        if filter_option[0] == "nodes_only":
            filter_shared_nodes_by_spread(
                phy.G,
                os.path.join(cycle_output_path, "shared_nodes_all.txt"),
                os.path.join(cycle_output_path, "shared_nodes_all_filtered.txt"),
                MAX_SHARED_NODE_SPREAD,
            )
    else:
        print("Tree died")

if __name__ == "__main__":
    pr = cProfile.Profile()
    pr.enable()
    main()
    pr.disable()
    sortby = SortKey.CUMULATIVE
    with open("profile", "w") as f:
        ps = pstats.Stats(pr, stream=f).sort_stats(sortby)
        ps.print_stats()

