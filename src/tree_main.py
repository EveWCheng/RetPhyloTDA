import os
import shutil
from typing import Optional, Callable
import numpy as np
import networkx as nx
from tqdm.std import TqdmDefaultWriteLock

import cProfile, pstats
from pstats import SortKey

TqdmDefaultWriteLock.mp_lock = None

from sim_bdh import SimState, SimParams, _sim_one, enumerate_gene_trees
from export import export_csv
from find_cycles import find_cycles
from network_lab_tda.tree_edit.tree_addition import networkx_to_tree_json, merge_trees, visualize
from network_lab_tda.tda_analysis import harmonic_cycle
from network_lab_tda.tda_visualisation.tda_visual import tda_visual_from_jason

HERE = os.path.dirname(os.path.abspath(__file__))
TREE_GROUP_OUTPUTS_DIR = os.path.join(HERE, os.pardir, "Outputs", "tree_group_outputs")
TREE_GROUPS_DIR = os.path.join(TREE_GROUP_OUTPUTS_DIR, "tree_groups")
MERGED_TREE_DIR = os.path.join(TREE_GROUP_OUTPUTS_DIR, "merged_tree")
PHYLO_CSV_DIR = os.path.join(TREE_GROUP_OUTPUTS_DIR, "phylo_csv")

# ── Parameters ────────────────────────────────────────────────────────────────

AGE      = 4
MRCA     = True
LAMBDA   = 0.5
MU       = 0.1
NU       = 0.2
HYBPROPS = [0, 0, 1]   # [lineage generating, degenerative, neutral]
STOPPING_NUM_LEAVES = 8
MIN_CYCLE_LENGTH = 4
Ngene = 2
TRAIT_MODEL = None

hyb_inher_fxn = lambda: np.random.uniform(0, 1)
hyb_rate_fxn  = None


# ── Run ───────────────────────────────────────────────────────────────────────

def process_gene_trees(phy, which_nodes: str = "no_hyb_nodes"):
    filtered_G = phy.filter_nodes(which_nodes=which_nodes)
#    print(filtered_G.nodes(data="is_hyb_node"))
    os.makedirs(TREE_GROUP_OUTPUTS_DIR, exist_ok=True)
    visualize(filtered_G, output=os.path.join(TREE_GROUP_OUTPUTS_DIR, "filtered_G.html"))
    enumerated_trees = enumerate_gene_trees(filtered_G,n_samples=4)
    os.makedirs(TREE_GROUPS_DIR, exist_ok=True)
    for entry in os.scandir(TREE_GROUPS_DIR):
        if entry.is_file():
            os.remove(entry.path)
        else:
            shutil.rmtree(entry.path)
    for i, (tree, weight) in enumerate(enumerated_trees):
        networkx_to_tree_json(tree, os.path.join(TREE_GROUPS_DIR, f"gene_tree_{i}.json"), numeric_labels=True)
        visualize(tree, output=os.path.join(TREE_GROUPS_DIR, f"gene_tree_{i}.html"))

    os.makedirs(MERGED_TREE_DIR, exist_ok=True)
    merged_G = merge_trees(input_dir=TREE_GROUPS_DIR, output_dir=MERGED_TREE_DIR)

    index_to_name = dict(merged_G.nodes(data="label"))

    dist_matrix = nx.floyd_warshall_numpy(merged_G)
    hc = harmonic_cycle(dist_matrix, cycle_dim=1, sim_log=True, log_path=os.path.join(TREE_GROUP_OUTPUTS_DIR, "rip.json"))
    hc.run_harmonics()

    vis = tda_visual_from_jason(data=hc.log, log_path=os.path.join(TREE_GROUP_OUTPUTS_DIR, "harmonic_cycle_plots"),index_to_name=index_to_name)
    vis.cycle_plot()

    return hc


def main(seed=43, which_nodes: str = "no_hyb_nodes"):
    if seed is not None:
        np.random.seed(seed)

    params = SimParams(age=AGE, lambda_=LAMBDA, mu=MU, nu=NU, hybprops=HYBPROPS,hyb_inher_fxn=hyb_inher_fxn,hyb_rate_fxn=hyb_rate_fxn,stopping_num_leaves=STOPPING_NUM_LEAVES)
    state = SimState(mrca=MRCA, Ngene=Ngene, trait_model=TRAIT_MODEL)
    phy = _sim_one(state,params)['phy']
    print("original tree is calcualted")
    if phy != 0:
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

