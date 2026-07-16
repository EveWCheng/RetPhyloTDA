import os
import shutil
import numpy as np
from tqdm.std import TqdmDefaultWriteLock

import cProfile, pstats
from pstats import SortKey

TqdmDefaultWriteLock.mp_lock = None

from sim_bdh import SimState, SimParams, _sim_one, enumerate_gene_trees
from export import export_csv
from network_lab_tda.tree_edit.tree_addition import networkx_to_tree_json, merge_trees, visualize
from find_cycles import CycleFinder

HERE = os.path.dirname(os.path.abspath(__file__))
TREE_GROUP_OUTPUTS_DIR = os.path.join(HERE, os.pardir, "Outputs", "tree_group_outputs")
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
MU       = 0.1
# hybridization rate
NU       = 0.01
HYBPROPS = [1, 0,0]   # [lineage generating, degenerative, neutral]
# stop simulation once tree reaches this many leaves
STOPPING_NUM_LEAVES = 30
# minimum cycle length filter for TDA cycle detection
MIN_CYCLE_LENGTH = 4
# number of gene trees to simulate
Ngene = 2
# trait evolution model (none used)
TRAIT_MODEL = None
# number of gene-tree samples to enumerate
N_SAMPLES = 2
# fixed thresholds to plot cycles at; set to None to derive thresholds dynamically via THRESHOLD_MODE instead
THRESHOLDS = [1]
# threshold-selection strategies CycleFinder runs per cycle
THRESHOLD_MODE = ["fixed"]
# cycle-qualification strategies CycleFinder runs per cycle
CYCLE_QUALIFY_MODE = ["crossover"]

# draws hybrid inheritance probability
hyb_inher_fxn = lambda: np.random.uniform(0, 1)
# function for hybridization rate (unset, uses default NU)
hyb_rate_fxn  = None


# ── Run ───────────────────────────────────────────────────────────────────────

def process_gene_trees(phy, which_nodes: str = "no_hyb_nodes"):
    filtered_G = phy.filter_nodes(which_nodes=which_nodes)
#    print(filtered_G.nodes(data="is_hyb_node"))
    os.makedirs(TREE_GROUP_OUTPUTS_DIR, exist_ok=True)
    visualize(filtered_G, output=os.path.join(TREE_GROUP_OUTPUTS_DIR, "filtered_G.html"))
    enumerated_trees = enumerate_gene_trees(filtered_G,n_samples=N_SAMPLES)
    os.makedirs(TREE_GROUPS_DIR, exist_ok=True)
    for i, (tree, weight) in enumerate(enumerated_trees):
        networkx_to_tree_json(tree, os.path.join(TREE_GROUPS_DIR, f"gene_tree_{i}.json"), numeric_labels=True)
        visualize(tree, output=os.path.join(TREE_GROUPS_DIR, f"gene_tree_{i}.html"))

    os.makedirs(MERGED_TREE_DIR, exist_ok=True)
    merged_G = merge_trees(input_dir=TREE_GROUPS_DIR, output_dir=MERGED_TREE_DIR)

    cf = CycleFinder(merged_G, threshold_mode=THRESHOLD_MODE, cycle_qualify_mode=CYCLE_QUALIFY_MODE, output_dir=TREE_GROUP_OUTPUTS_DIR, thresholds=THRESHOLDS, min_cycle_length=MIN_CYCLE_LENGTH, use_data_prep=False, vis=True)
    return cf.find_cycles()


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
        process_gene_trees(phy, which_nodes=which_nodes)
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

