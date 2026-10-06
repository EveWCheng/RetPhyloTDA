"""
Visualization utilities for PhyloNetwork objects.
"""
from __future__ import annotations
import os
import csv

# ── CSV export ─────────────────────────────────────────────────────────────────

def export_csv(G, out_dir: str, prefix: str = "", by_label: bool = False):
    """Write {prefix}nodes.csv and {prefix}edges.csv for plotting in R with igraph.

    by_label: use node labels as ids (and in from/to) instead of the graph's
    node keys -- for filtered graphs, whose ids differ from phy.G's numeric ones.
    """
    os.makedirs(out_dir, exist_ok=True)
    nodes_path = os.path.join(out_dir, f"{prefix}nodes.csv")
    edges_path = os.path.join(out_dir, f"{prefix}edges.csv")

    def node_id(n):
        return G.nodes[n]["label"] if by_label else n

    with open(nodes_path, "w", newline="") as f:
        w = csv.writer(f)
        # 'type' is a display category (hyb_leaf sticks even after the node speciates);
        # 'is_leaf' is the live tip status -- use it, not type, to decide what is a tip.
        w.writerow(["id", "label", "type", "is_leaf"])
        for n, attrs in G.nodes(data=True):
            if attrs.get("extinct"):
                ntype = "extinct"
            elif attrs.get("is_hyb_leaf"):
                ntype = "hyb_leaf"
            elif attrs["is_leaf"]:
                ntype = "leaf"
            else:
                ntype = "internal"
            w.writerow([node_id(n), attrs['label'], ntype, bool(attrs["is_leaf"])])

    with open(edges_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["from", "to", "edge_type", "length", "time_length", "inher_weight"])
        for u, v, attrs in G.edges(data=True):
            w.writerow([
                node_id(u),
                node_id(v),
                attrs["edge_type"],
                attrs["length"],
                attrs["time_length"],
                attrs.get("inher_weight", ""),
            ])
