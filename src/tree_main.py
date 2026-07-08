import os
from typing import Optional, Callable
import numpy as np
import networkx as nx
from tqdm.std import TqdmDefaultWriteLock

TqdmDefaultWriteLock.mp_lock = None

from sim_bdh import SimState, SimParams, _sim_one, enumerate_gene_trees
from export import export_csv
from find_cycles import find_cycles
from network_lab_tda.tree_edit.tree_addition import networkx_to_tree_json, merge_trees
from network_lab_tda.tda_analysis import harmonic_cycle

HERE = os.path.dirname(os.path.abspath(__file__))
TREE_GROUP_OUTPUTS_DIR = os.path.join(HERE, os.pardir, "Outputs", "tree_group_outputs")
TREE_GROUPS_DIR = os.path.join(TREE_GROUP_OUTPUTS_DIR, "tree_groups")

# ── Parameters ────────────────────────────────────────────────────────────────

AGE      = 4
MRCA     = True
LAMBDA   = 0.5
MU       = 0.1
NU       = 0.5
HYBPROPS = [1, 1, 1]   # [lineage generating, degenerative, neutral]
MIN_CYCLE_LENGTH = 4
Ngene = 2
TRAIT_MODEL = None

hyb_inher_fxn = lambda: np.random.uniform(0, 1)
hyb_rate_fxn  = None


# ── Run ───────────────────────────────────────────────────────────────────────

def process_gene_trees(phy, which_nodes: str = "no_hyb_nodes"):
    filtered_G = phy.filter_nodes(which_nodes=which_nodes)
    enumerated_trees = enumerate_gene_trees(filtered_G)

    os.makedirs(TREE_GROUPS_DIR, exist_ok=True)
    for i, (tree, weight) in enumerate(enumerated_trees):
        networkx_to_tree_json(tree, os.path.join(TREE_GROUPS_DIR, f"gene_tree_{i}.json"))

    merged_G = merge_trees(input_dir=TREE_GROUPS_DIR, output_dir=TREE_GROUPS_DIR)

    print("Merged")

    dist_matrix = nx.floyd_warshall_numpy(merged_G)
    hc = harmonic_cycle(dist_matrix, cycle_dim=1, sim_log=True, log_path=os.path.join(TREE_GROUP_OUTPUTS_DIR, "rip.json"))
    hc.run_harmonics()
    return hc


def main(seed=None, which_nodes: str = "no_hyb_nodes"):
    if seed is not None:
        np.random.seed(seed)

    params = SimParams(age=AGE, lambda_=LAMBDA, mu=MU, nu=NU, hybprops=HYBPROPS,hyb_inher_fxn=hyb_inher_fxn,hyb_rate_fxn=hyb_rate_fxn)
    state = SimState(mrca=MRCA, Ngene=Ngene, trait_model=TRAIT_MODEL)
    phy = _sim_one(state,params)['phy']
    if phy != 0:
        process_gene_trees(phy, which_nodes=which_nodes)
    else:
        print("Tree died")



if __name__ == "__main__":
    main(seed=42)
