import os
import shutil
from typing import Optional, Callable

import numpy as np
from tqdm.std import TqdmDefaultWriteLock

# Prevent tqdm from creating a multiprocessing.RLock() for its default write
# lock. This program is single-process, so the lock is never needed, and
# leaving it unset causes a "leaked semaphore" warning from
# multiprocessing.resource_tracker at interpreter shutdown.
TqdmDefaultWriteLock.mp_lock = None

from sim_bdh import SimState, SimParams, _sim_one
from export import export_csv, export_filtered_edges_csv
from find_cycles import CycleFinder

HERE = os.path.dirname(os.path.abspath(__file__))
SIM_OUTPUTS_DIR = os.path.join(HERE, os.pardir, "outputs", "sim_phylo_outputs")
PHYLO_CSV_DIR = os.path.join(SIM_OUTPUTS_DIR, "phylo_csv")


def sim_bdh_age(age: float, numbsim: int,
                 lambda_: float, mu: float, nu: float,
                 hybprops: list[float], hyb_inher_fxn: Callable,
                 mrca: bool = False,
                 hyb_rate_fxn: Optional[Callable] = None,
                 Ngene: int = 0,
                 trait_model: Optional[dict] = None,
                 stopping_num_leaves: Optional[int] = None) -> list[dict]:
    """Run numbsim independent BDH simulations and return a list of result dicts.

    trait_model, if given, must supply callables under these keys:
        initial                  -- sequence of starting trait value(s), one per root leaf
        time_fxn(trait, dt)      -- evolve one lineage's trait over elapsed time dt
        spec_fxn(trait)          -- trait -> (child1_trait, child2_trait) at speciation
        hyb_event_fxn(t1, t2, inher)         -- parents' traits + inheritance -> hybrid trait
        hyb_compatibility_fxn(t1, t2, hyb_trait) -- bool: whether the hybridization can occur
    """
    params = SimParams(age=age, lambda_=lambda_, mu=mu, nu=nu, hybprops=hybprops, hyb_inher_fxn=hyb_inher_fxn, hyb_rate_fxn=hyb_rate_fxn, stopping_num_leaves=STOPPING_NUM_LEAVES)
    results = []
    for i in range(numbsim):
        state = SimState(mrca=mrca, Ngene=Ngene, trait_model=trait_model)
        result = _sim_one(state, params)
        if result["phy"] != 0:
            size = result['phy'].G.number_of_nodes()
            print(f"size for {i}: {size}")
            results.append(result)
    return results


# ── Parameters ────────────────────────────────────────────────────────────────

AGE      = 10
NUMBSIM  = 30
LAMBDA   = 0.5
MU       = 0.1
NU       = 0.5
HYBPROPS = [1, 0, 0]   # [lineage generating, degenerative, neutral]
STOPPING_NUM_LEAVES = 10  # each sim also stops once it reaches this many leaves
MIN_CYCLE_LENGTH = 0
WEIGHT_ATTR = "length"  # edge attribute CycleFinder measures distance with: "length" (genetic) or "time_length" (time). Should always be time_length
# options: "true_distance_between_tips" (overwrite leaf-pair distances with true shortest-path tip distances),
#          None (no dressing, use dist_matrix as computed)
DRESS_DISTANCE_MATRIX = "true_distance_between_tips"

hyb_inher_fxn = lambda: np.random.uniform(0, 1)
hyb_rate_fxn  = None


def should_populate_fxn(u, v, attrs):
    """Don't add phantom nodes to reticulation edges."""
    return attrs.get("edge_type") != "reticulation"


# ── Run ───────────────────────────────────────────────────────────────────────

def main(seed=42, gene_index: Optional[int] = None, which_nodes: str = "no_hyb_nodes"):
    if seed is not None:
        np.random.seed(seed)

    if os.path.exists(SIM_OUTPUTS_DIR):
        shutil.rmtree(SIM_OUTPUTS_DIR)
    os.makedirs(PHYLO_CSV_DIR)

    results = sim_bdh_age(
        age=AGE,
        numbsim=NUMBSIM,
        lambda_=LAMBDA,
        mu=MU,
        nu=NU,
        hybprops=HYBPROPS,
        hyb_inher_fxn=hyb_inher_fxn,
        mrca=False,
        hyb_rate_fxn=hyb_rate_fxn,  # None
        Ngene=0,
        trait_model=None,
        stopping_num_leaves=STOPPING_NUM_LEAVES,
    )
    print(f"length of results:{len(results)}")

    for i, r in enumerate(results):
        phy = r['phy']
        if phy == 0:
            print(f"sim{i}: extinct")
            continue
        else:
            print(f"sim{i}: not extinct")

        export_csv(phy, PHYLO_CSV_DIR, prefix=f"sim{i}_")
        filtered_G = phy.filter_nodes(which_nodes=which_nodes)

        export_filtered_edges_csv(filtered_G, PHYLO_CSV_DIR, prefix=f"sim{i}_")
        max_edge_length = max(d for _, _, d in filtered_G.edges(data=WEIGHT_ATTR))
        # snapshot the network at each reticulation edge's own length, in addition to cycle-birth
        # thresholds (CycleFinder appends those automatically since "fixed" isn't in threshold_mode)
        retic_edge_lengths = [
            attrs[WEIGHT_ATTR] for _, _, attrs in filtered_G.edges(data=True) if attrs.get("edge_type") == "reticulation"
        ]
        CycleFinder(filtered_G, threshold_mode=["cyclelength", "marker"], cycle_qualify_mode=["marker"], output_dir=SIM_OUTPUTS_DIR, which_nodes=which_nodes, sim_label=f"sim{i}", min_cycle_length=MIN_CYCLE_LENGTH, weight_attr=WEIGHT_ATTR, rips_threshold=max_edge_length, dress_distance_matrix=DRESS_DISTANCE_MATRIX, should_populate_fxn=should_populate_fxn, thresholds=retic_edge_lengths).find_cycles()


if __name__ == "__main__":
    main()
