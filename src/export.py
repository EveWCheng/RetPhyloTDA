"""
Visualization utilities for PhyloNetwork objects.
"""
from __future__ import annotations
import os
import csv
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from sim_bdh import PhyloNetwork


# ── CSV export ─────────────────────────────────────────────────────────────────

def export_csv(phy: PhyloNetwork, out_dir: str, prefix: str = ""):
    """Write nodes.csv and edges.csv for plotting in R with igraph."""
    os.makedirs(out_dir, exist_ok=True)
    nodes_path = os.path.join(out_dir, f"{prefix}nodes.csv")
    edges_path = os.path.join(out_dir, f"{prefix}edges.csv")

    with open(nodes_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id", "label", "type"])
        for n, attrs in phy.G.nodes(data=True):
            if attrs.get("extinct"):
                ntype = "extinct"
            elif attrs.get("is_hyb_leaf"):
                ntype = "hyb_leaf"
            elif attrs["is_leaf"]:
                ntype = "leaf"
            else:
                ntype = "internal"
            w.writerow([n, attrs['label'], ntype])

    with open(edges_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["from", "to", "edge_type", "length", "time_length"])
        for u, v, attrs in phy.G.edges(data=True):
            w.writerow([u, v, attrs["edge_type"], attrs["length"], attrs["time_length"]])

