import os
from typing import Optional, Callable

import numpy as np
from tqdm.std import TqdmDefaultWriteLock

# Prevent tqdm from creating a multiprocessing.RLock() for its default write
# lock. This program is single-process, so the lock is never needed, and
# leaving it unset causes a "leaked semaphore" warning from
# multiprocessing.resource_tracker at interpreter shutdown.
TqdmDefaultWriteLock.mp_lock = None

from sim_bdh import SimState, _sim_one
from export import export_csv
from find_cycles import find_cycles

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
    results = []
    for i in range(numbsim):
        state = SimState(mrca=mrca, Ngene=Ngene, trait_model=trait_model)
        result = _sim_one(state, age, lambda_, mu, nu, hybprops, hyb_inher_fxn, hyb_rate_fxn)
        if result["phy"] is not 0:
            size = result['phy'].G.number_of_nodes()
            print(f"size for {i}: {size}")
            if size < 90:
                results.append(result)    
    return results


def _filter_nodes(G, which_nodes="all_nodes"):
    """Return the induced subgraph for the requested node subset.

    all_nodes     -- every node except extinct tips (label starts with 'ext')
    species_only  -- only surviving species tips (label contains 'sp')
    """
    if which_nodes == "species_only":
        keep = [n for n, label in G.nodes(data="label") if "sp" in label]
    else:
        keep = [n for n, label in G.nodes(data="label") if not label.startswith("ext")]
    return G.subgraph(keep).copy()


# ── Parameters ────────────────────────────────────────────────────────────────

AGE      = 4.0
NUMBSIM  = 35
LAMBDA   = 0.5
MU       = 0.1
NU       = 0.5
HYBPROPS = [1, 1, 1]   # [lineage generating, degenerative, neutral]
MIN_CYCLE_LENGTH = 4

hyb_inher_fxn = lambda: np.random.uniform(0, 1)
hyb_rate_fxn  = None


# ── Run ───────────────────────────────────────────────────────────────────────

def main(seed=None, gene_index: Optional[int] = None, which_nodes: str = "all_nodes"):
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
            print(f"sim {i}: extinct")
            continue
        else:
            print(f"sim{i}: not extinct")

        export_csv(phy, OUT_DIR, prefix=f"sim{i}_")
        filtered_G = _filter_nodes(phy.G, which_nodes=which_nodes)
        find_cycles(filtered_G, which_nodes=which_nodes, sim_label=f"sim{i}", min_cycle_length=MIN_CYCLE_LENGTH)

        if gene_index is not None:
            gtree = phy.gene_tree(gene_index)
            print(f"  gene {gene_index}: {gtree.number_of_nodes()} nodes, "
                  f"{gtree.number_of_edges()} edges")


if __name__ == "__main__":
    main(seed=42)
