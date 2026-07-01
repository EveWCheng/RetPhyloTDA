import os
import numpy as np
from sim_bdh import sim_bdh_age
from export import export_csv

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, os.pardir, "Outputs", "phylo_outputs")


# ── Parameters ────────────────────────────────────────────────────────────────

AGE      = 4.0
NUMBSIM  = 10
LAMBDA   = 0.5
MU       = 0.1
NU       = 0.5
HYBPROPS = [1, 1, 1]   # [lineage generating, degenerative, neutral]

hyb_inher_fxn = lambda: np.random.uniform(0, 1)
hyb_rate_fxn  = None


# ── Run ───────────────────────────────────────────────────────────────────────

def main(seed=None):
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
        hyb_rate_fxn=hyb_rate_fxn,
        mrca=False,
    )

    for i, r in enumerate(results):
        phy = r['phy']
        if phy == 0:
            print(f"sim {i}: extinct")
            continue
        print(f"sim {i}: {phy.nleaves} leaves, "
              f"{len(phy.extinct)} extinct, "
              f"{len(phy.hyb_tips)} hyb_tips")
        export_csv(phy, OUT_DIR, prefix=f"sim{i}_")


if __name__ == "__main__":
    main(seed=42)
