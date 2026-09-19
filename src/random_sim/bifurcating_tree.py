import os
import random

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from network_lab_tda.tda_analysis import harmonic_cycle


def build_bifurcating_tree(n_leaves: int, seed: int | None = None) -> nx.DiGraph:
    """Random bifurcating tree grown by repeatedly splitting a random leaf into two children."""
    rng = random.Random(seed)
    G = nx.DiGraph()
    G.add_node(0)
    leaves = [0]
    next_id = 1
    while len(leaves) < n_leaves:
        parent = rng.choice(leaves)
        leaves.remove(parent)
        child1, child2 = next_id, next_id + 1
        next_id += 2
        G.add_edge(parent, child1, length=1)
        G.add_edge(parent, child2, length=1)
        leaves.extend([child1, child2])
    return G


def _hierarchy_pos(G: nx.DiGraph, root):
    """x from in-order leaf position, y from depth below root."""
    pos = {}

    def leaf_count(n):
        children = list(G.successors(n))
        return 1 if not children else sum(leaf_count(c) for c in children)

    def assign(n, x_left, x_right, depth):
        pos[n] = ((x_left + x_right) / 2, -depth)
        children = list(G.successors(n))
        total = sum(leaf_count(c) for c in children)
        cursor = x_left
        for c in children:
            frac = leaf_count(c) / total
            assign(c, cursor, cursor + frac * (x_right - x_left), depth + 1)
            cursor += frac * (x_right - x_left)

    assign(root, 0, 1, 0)
    return pos


def distance_matrix(G: nx.DiGraph):
    """All-pairs distance, weighted by each edge's 'length' attribute."""
    undirected = G.to_undirected()
    nodes = sorted(undirected.nodes())
    index = {node: i for i, node in enumerate(nodes)}
    dist = np.zeros((len(nodes), len(nodes)))
    for u, lengths in nx.all_pairs_dijkstra_path_length(undirected, weight="length"):
        for v, d in lengths.items():
            dist[index[u], index[v]] = d
    return dist, nodes


def find_harmonic_cycles(dist_matrix, nodes, log_path="bifurcating_tree_cycles/rip.json"):
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    hc = harmonic_cycle(dist_matrix, cycle_dim=1, sim_log=True, log_path=log_path)
    hc.run_harmonics(threshold=float('inf'), save=True)
    cycles = hc.log["harmonic_cycles"]
    print(f"{len(cycles)} harmonic cycle(s) found")
    for i, cycle in enumerate(cycles):
        cycle_nodes = sorted({nodes[idx] for edge in cycle["edges"] for idx in edge["simplex"]})
        print(f"  cycle {i}: birth={cycle['birth']:.4f} death={cycle['death']} n_edges={len(cycle['edges'])} nodes={cycle_nodes}")
    return cycles


def plot_tree(G: nx.DiGraph, root=0, out_path: str = "bifurcating_tree.png"):
    pos = _hierarchy_pos(G, root)
    node_colors = ["tab:orange" if G.out_degree(n) == 0 else "tab:blue" for n in G.nodes]

    plt.figure(figsize=(10, 6))
    nx.draw(G, pos, with_labels=True, node_color=node_colors, node_size=500, font_size=8, arrows=False)
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"saved to {out_path}")


def sweep_pair_distance(dist_matrix, nodes, a, b, values):
    """Override dist[a,b]/dist[b,a] with each value in turn, run harmonic_cycle, report cycle count."""
    ia, ib = nodes.index(a), nodes.index(b)
    for v in values:
        m = dist_matrix.copy()
        m[ia, ib] = v
        m[ib, ia] = v
        hc = harmonic_cycle(m, cycle_dim=1, sim_log=False)
        hc.run_harmonics(threshold=float('inf'), save=False)
        cycles = hc.log["harmonic_cycles"]
        print(f"dist({a},{b})={v}: {len(cycles)} cycle(s)" + (f"  birth={cycles[0]['birth']:.4f}" if cycles else ""))


def sweep_two_pairs_by_amount(dist_matrix, nodes, pair1, pair2, amounts):
    """Subtract each amount from both pairs' baseline distances simultaneously, run harmonic_cycle, report."""
    (a, b), (c, d) = pair1, pair2
    ia, ib, ic, id_ = nodes.index(a), nodes.index(b), nodes.index(c), nodes.index(d)
    base1, base2 = dist_matrix[ia, ib], dist_matrix[ic, id_]
    for s in amounts:
        v1, v2 = base1 - s, base2 - s
        m = dist_matrix.copy()
        m[ia, ib] = m[ib, ia] = v1
        m[ic, id_] = m[id_, ic] = v2
        hc = harmonic_cycle(m, cycle_dim=1, sim_log=False)
        hc.run_harmonics(threshold=float('inf'), save=False)
        cycles = hc.log["harmonic_cycles"]
        details = "; ".join(f"birth={c['birth']:.4f} death={c['death']}" for c in cycles)
        print(f"shorten={s}  dist({a},{b})={v1}  dist({c},{d})={v2}: {len(cycles)} cycle(s)  {details}")


if __name__ == "__main__":
    G = build_bifurcating_tree(n_leaves=10, seed=42)
    plot_tree(G)
    dist, nodes = distance_matrix(G)
    print("nodes:", nodes)
    print(dist)
    find_harmonic_cycles(dist, nodes)

    print()
    print("sweeping dist(13, 18) down from 1.0:")
    sweep_pair_distance(dist, nodes, 13, 18, [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0])
