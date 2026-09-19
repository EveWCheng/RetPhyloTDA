import json
import os
import shutil
import sys
import tomllib
from typing import Optional, Callable

import numpy as np
from tqdm.std import TqdmDefaultWriteLock
from network_lab_tda.data_prep.Data_Prep import Data_Prep
from network_lab_tda.data_prep.Populate_Edge import Populate_Edge

HERE = os.path.dirname(os.path.abspath(__file__))
SHARED_DIR = os.path.join(HERE, os.pardir, "shared")
if SHARED_DIR not in sys.path:
    sys.path.insert(0, SHARED_DIR)

# leaving it unset causes a "leakedsemaphore" warning from multiprocessing.resource_tracker at interpreter shutdown.
TqdmDefaultWriteLock.mp_lock = None

from sim_bdh import SimState, SimParams, _sim_one
from export import export_filtered
from find_cycles import CycleFinder
from filter_cycle import FilterCycle

SIM_OUTPUTS_DIR = os.path.join(HERE, os.pardir, os.pardir, "outputs", "sim_phylo_outputs")
PHYLO_CSV_DIR = os.path.join(SIM_OUTPUTS_DIR, "phylo_csv")
CONFIG_PATH = os.path.join(HERE, "sim_config.toml")


def load_config(path=CONFIG_PATH):
    with open(path, "rb") as f:
        return tomllib.load(f)["sim"]


def sim_bdh_age(age: float, numbsim: int,
                 lambda_: float, mu: float, nu: float,
                 hybprops: list[float], hyb_inher_fxn: Callable,
                 mrca: bool = False,
                 hyb_rate_fxn: Optional[Callable] = None,
                 Ngene: int = 0,
                 trait_model: Optional[dict] = None,
                 stopping_num_leaves: Optional[int] = None,
                 which_nodes: str = "no_hyb_nodes") -> list[dict]:
    params = SimParams(age=age, lambda_=lambda_, mu=mu, nu=nu, hybprops=hybprops, hyb_inher_fxn=hyb_inher_fxn, hyb_rate_fxn=hyb_rate_fxn, stopping_num_leaves=stopping_num_leaves)
    results = []
    for i in range(numbsim):
        state = SimState(mrca=mrca, Ngene=Ngene, trait_model=trait_model)
        result = _sim_one(state, params, which_nodes=which_nodes)
        if result["phy"] != 0:
            size = result['phy'].G.number_of_nodes()
            print(f"size for {i}: {size}")
            results.append(result)
    return results


hyb_inher_fxn = lambda: np.random.uniform(0, 1)
hyb_rate_fxn  = None


def should_populate_fxn(u, v, attrs):
    """Don't add phantom nodes to reticulation edges."""
    return attrs.get("edge_type") != "reticulation"


def _dist_matrix_from_dict(tree_only_G_undirected, distance):
    node_to_index = {node: idx for idx, node in enumerate(tree_only_G_undirected.nodes())}
    n = len(node_to_index)
    dist_matrix = np.zeros((n, n))
    for u, row in distance.items():
        for v, d in row.items():
            dist_matrix[node_to_index[u], node_to_index[v]] = d
    return dist_matrix


def write_distance_include_tips_changes(filtered_G, distance_include_tips_changes, output_path):
    label_of = dict(filtered_G.nodes(data="label"))
    by_label = {
        f"{label_of[u]}--{label_of[v]}": change
        for (u, v), change in distance_include_tips_changes.items()
        if u in label_of and v in label_of
    }
    with open(os.path.join(output_path, "distance_include_tips_changes.json"), "w") as f:
        json.dump(by_label, f, indent=2)


def build_dist_matrix(G, distance_include_tips, distance_only_tips, use_data_prep, weight_attr, output_path, populated_header_fn="populated_headers.txt", should_populate_fxn=None):
    if use_data_prep:
        dp = Data_Prep(G=G, log_path=output_path, headers=False, weight_attr=weight_attr)
        pe = Populate_Edge(G=dp.G, log_path=output_path, headers=False, populated_header_fn=populated_header_fn, max_node_per_edge=0, weight_attr=weight_attr, should_populate_fxn=should_populate_fxn)
        dist_matrix = pe.populate_edges()
        index_to_name = pe.index_to_name
        node_to_index = {node: idx for idx, node in enumerate(pe.G.nodes())}
        for u, row in distance_only_tips.items():
            for v, d in row.items():
                dist_matrix[node_to_index[u], node_to_index[v]] = d
    else:
        dist_matrix = _dist_matrix_from_dict(G, distance_include_tips)
        index_to_name = {
            i: attrs.get("label", node)
            for i, (node, attrs) in enumerate(G.nodes(data=True))
        }
    with open(os.path.join(output_path, populated_header_fn), "w") as f:
        f.write("\n".join(str(v) for v in index_to_name.values()))
        np.savetxt(os.path.join(output_path, "populated_distance_matrix.txt"), dist_matrix)
    return dist_matrix, index_to_name


# ── Run ───────────────────────────────────────────────────────────────────────

def main(config_path=CONFIG_PATH, gene_index: Optional[int] = None):
    config = load_config(config_path)
    seed = config["seed"]
    which_nodes = config["which_nodes"]
    weight_attr = config["weight_attr"]
    min_cycle_length = config["min_cycle_length"]

    if seed is not None:
        np.random.seed(seed)

    if os.path.exists(SIM_OUTPUTS_DIR):
        shutil.rmtree(SIM_OUTPUTS_DIR)
    os.makedirs(PHYLO_CSV_DIR)

    results = sim_bdh_age(
        age=config["age"],
        numbsim=config["numbsim"],
        lambda_=config["lambda_"],
        mu=config["mu"],
        nu=config["nu"],
        hybprops=config["hybprops"],
        hyb_inher_fxn=hyb_inher_fxn,
        mrca=False,
        hyb_rate_fxn=hyb_rate_fxn,  # None
        Ngene=0,
        trait_model=None,
        stopping_num_leaves=config["stopping_num_leaves"],
        which_nodes=which_nodes,
    )
    print(f"length of results:{len(results)}")

    for i, r in enumerate(results):
        print(f"sim{i}: not extinct")

        filtered_G = r['filtered_G']

        export_filtered(filtered_G, PHYLO_CSV_DIR, prefix=f"sim{i}_")
        max_edge_length = max(d for _, _, d in filtered_G.edges(data=weight_attr))
        # snapshot the network at each reticulation edge's own length, in addition to cycle-birth
        # thresholds (CycleFinder appends those automatically since "fixed" isn't in threshold_mode)
        retic_edge_lengths = [
            attrs[weight_attr]+1e-2 for _, _, attrs in filtered_G.edges(data=True) if attrs.get("edge_type") == "reticulation"
        ]
        print(retic_edge_lengths)
        undirected_filtered_G_tree = r['tree_only_G'].to_undirected()

        output_path = os.path.join(SIM_OUTPUTS_DIR, "proc_phylo_outputs", f"sim{i}")
        os.makedirs(output_path, exist_ok=True)

        write_distance_include_tips_changes(filtered_G, r['distance_include_tips_changes'], output_path)

        dist_matrix, index_to_name = build_dist_matrix(
            undirected_filtered_G_tree,
            r['distance_include_tips'],
            r['distance_only_tips'],
            use_data_prep=False,
            weight_attr=weight_attr,
            output_path=output_path,
            should_populate_fxn=should_populate_fxn,
        )
        cf = CycleFinder(undirected_filtered_G_tree, threshold_mode=["cyclelength", "marker"], output_dir=SIM_OUTPUTS_DIR, dist_matrix=dist_matrix, index_to_name=index_to_name, sim_label=f"sim{i}", rips_threshold=float('inf'), thresholds=retic_edge_lengths)
        cf.find_cycles()
        FilterCycle(cf, cycle_qualify_mode=[], which_nodes=which_nodes, min_cycle_length=min_cycle_length).visualize()


if __name__ == "__main__":
    main()
