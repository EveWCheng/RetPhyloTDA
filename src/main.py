import os
from typing import Optional, Callable

import numpy as np
from tqdm.std import TqdmDefaultWriteLock

# Prevent tqdm from creating a multiprocessing.RLock() for its default write
# lock. This program is single-process, so the lock is never needed, and
# leaving it unset causes a "leaked semaphore" warning from
# multiprocessing.resource_tracker at interpreter shutdown.
TqdmDefaultWriteLock.mp_lock = None

from sim_bdh import SimState, SimParams, _sim_one
from export import export_csv
from find_cycles import CycleFinder
import cProfile, pstats
from pstats import SortKey

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, os.pardir, "Outputs", "phylo_outputs")


def sim_bdh_age(age: float, numbsim: int,
                 lambda_: float, mu: float, nu: float,
                 hybprops: list[float], hyb_inher_fxn: Callable,
                 mrca: bool = False,
                 hyb_rate_fxn: Optional[Callable] = None,
                 Ngene: int = 0,
                 trait_model: Optional[dict] = None) -> list[dict]:
    """Run numbsim independent BDH simulations and return a list of result dicts.

    trait_model, if given, must supply callables under these keys:
        initial                  -- sequence of starting trait value(s), one per root leaf
        time_fxn(trait, dt)      -- evolve one lineage's trait over elapsed time dt
        spec_fxn(trait)          -- trait -> (child1_trait, child2_trait) at speciation
        hyb_event_fxn(t1, t2, inher)         -- parents' traits + inheritance -> hybrid trait
        hyb_compatibility_fxn(t1, t2, hyb_trait) -- bool: whether the hybridization can occur
    """
    params = SimParams(age=age, lambda_=lambda_, mu=mu, nu=nu, hybprops=hybprops, hyb_inher_fxn=hyb_inher_fxn, hyb_rate_fxn=hyb_rate_fxn)
    results = []
    for i in range(numbsim):
        state = SimState(mrca=mrca, Ngene=Ngene, trait_model=trait_model)
        result = _sim_one(state, params)
        if result["phy"] != 0:
            size = result['phy'].G.number_of_nodes()
            print(f"size for {i}: {size}")
            if size < 100:
                results.append(result)    
    return results


# ── Parameters ────────────────────────────────────────────────────────────────

AGE      = 4
NUMBSIM  = 30
LAMBDA   = 0.5
MU       = 0.1
NU       = 0.5
HYBPROPS = [1, 0, 1]   # [lineage generating, degenerative, neutral]
MIN_CYCLE_LENGTH = 4

hyb_inher_fxn = lambda: np.random.uniform(0, 1)
hyb_rate_fxn  = None


# ── Run ───────────────────────────────────────────────────────────────────────

def main(seed=42, gene_index: Optional[int] = None, which_nodes: str = "no_hyb_nodes"):
    if seed is not None:
        np.random.seed(seed)

    if not os.path.exists(OUT_DIR):
        os.makedirs(OUT_DIR)

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
    )
    print(f"length of results:{len(results)}")

    for i, r in enumerate(results):
        phy = r['phy']
        if phy == 0:
            print(f"sim{i}: extinct")
            continue
        else:
            print(f"sim{i}: not extinct")

        export_csv(phy, OUT_DIR, prefix=f"sim{i}_")
        filtered_G = phy.filter_nodes(which_nodes=which_nodes)
        CycleFinder(filtered_G, threshold_mode=["cyclelength", "marker"], cycle_qualify_mode=["marker"], output_dir=os.path.join(HERE, os.pardir, "Outputs"), which_nodes=which_nodes, sim_label=f"sim{i}", min_cycle_length=MIN_CYCLE_LENGTH).find_cycles()

        if gene_index is not None:
            gtree = phy.gene_tree(gene_index)
#            print(f"  gene {gene_index}: {gtree.number_of_nodes()} nodes, ",f"{gtree.number_of_edges()} edges")


if __name__ == "__main__":
    pr = cProfile.Profile()
    pr.enable()
    main()
    pr.disable()
    sortby = SortKey.CUMULATIVE
    with open("profile", "w") as f:
        ps = pstats.Stats(pr, stream=f).sort_stats(sortby)
        ps.print_stats()
