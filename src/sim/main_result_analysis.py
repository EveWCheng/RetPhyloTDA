import csv
import json
import math
import os
import re

import networkx as nx

from sim_main import PHYLO_CSV_DIR, SIM_OUTPUTS_DIR

PROC_OUTPUTS_DIR = os.path.join(SIM_OUTPUTS_DIR, "proc_phylo_outputs")
CYCLE_OUTPUTS_DIR = os.path.join(SIM_OUTPUTS_DIR, "cycle_outputs")
ANALYSIS_REPORTS_DIR = os.path.join(SIM_OUTPUTS_DIR, "analysis_reports")

SIM_DIR_RE = re.compile(r"^sim(\d+)$")


def discover_sim_ids():
    ids = []
    for name in os.listdir(CYCLE_OUTPUTS_DIR):
        match = SIM_DIR_RE.match(name)
        if match:
            ids.append(int(match.group(1)))
    return sorted(ids)


def read_index_to_name(sim_id):
    path = os.path.join(PROC_OUTPUTS_DIR, f"sim{sim_id}", "populated_headers.txt")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return dict(enumerate(line.rstrip("\n") for line in f))


def read_edge_appears_at(sim_id):
    path = os.path.join(CYCLE_OUTPUTS_DIR, f"sim{sim_id}", "rip.json")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        data = json.load(f)
    appears_at = {}
    for simplex, t in zip(data["simplicies"], data["appears_at"]):
        if len(simplex) == 2:
            appears_at[frozenset(simplex)] = t
    return appears_at


def read_cycles(sim_id):
    path = os.path.join(CYCLE_OUTPUTS_DIR, f"sim{sim_id}", "rip.json")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        data = json.load(f)
    appears_at = read_edge_appears_at(sim_id)
    cycles = []
    for cycle in data["harmonic_cycles"]:
        edges = []
        for edge in cycle["edges"]:
            if math.isclose(edge["weight"], 0.0, abs_tol=1e-6):
                continue
            edges.append({
                "simplex": edge["simplex"],
                "weight": edge["weight"],
                "appears_at": appears_at[frozenset(edge["simplex"])],
            })
        cycles.append({"birth": cycle["birth"], "death": cycle["death"], "edges": edges})
    return cycles


def read_cycles_by_name(sim_id):
    index_to_name = read_index_to_name(sim_id)
    cycles = read_cycles(sim_id)
    if index_to_name is None or cycles is None:
        return None
    result = []
    for cycle in cycles:
        named_edges = []
        for edge in cycle["edges"]:
            nodes = [index_to_name[idx] for idx in edge["simplex"]]
            named_edges.append({"nodes": nodes, "weight": edge["weight"], "appears_at": edge["appears_at"]})
        result.append({"birth": cycle["birth"], "death": cycle["death"], "edges": named_edges})
    return result


def minimal_cycles_by_name(sim_id):
    cycles = read_cycles_by_name(sim_id)
    if cycles is None:
        return None
    result = []
    for cycle in cycles:
        graph = nx.Graph()
        for edge in cycle["edges"]:
            u, v = edge["nodes"]
            graph.add_edge(u, v)
        basis = nx.minimum_cycle_basis(graph)
        minimal_cycle = min(basis, key=len) if basis else None
        result.append(minimal_cycle)
    return result


def read_reticulate_edges_by_name(sim_id):
    path = os.path.join(PHYLO_CSV_DIR, f"sim{sim_id}_filtered_edges.csv")
    if not os.path.exists(path):
        return None
    with open(path, newline="") as f:
        return [
            (row["from"], row["to"], float(row["inher_weight"]) if row["inher_weight"] else None)
            for row in csv.DictReader(f)
            if row["edge_type"] == "reticulation"
        ]

 
def read_network_graph(sim_id):
    """Undirected graph from simN_filtered_edges.csv. Edges carry 'length' and 'edge_type'. Nodes carry 'is_leaf'. Dressed tip-tip cells are shortest paths on it."""
    path = os.path.join(PHYLO_CSV_DIR, f"sim{sim_id}_filtered_edges.csv")
    if not os.path.exists(path):
        return None
    G = nx.Graph()
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            u, v, length = row["from"], row["to"], float(row["length"])
            if G.has_edge(u, v) and G.edges[u, v]["length"] <= length:
                continue
            G.add_edge(u, v, length=length, edge_type=row["edge_type"])
    leaf_labels = read_leaf_labels(sim_id) or set()
    for node in G.nodes:
        G.nodes[node]["is_leaf"] = node in leaf_labels
    return G


def shortest_path_edges(net, u, v):
    """Edges on any shortest u<->v path, as (a, b, edge_type, length). Undirected, weight 'length'. Returns None if unreachable. Shortest paths are often non-unique here. That is because net has length=0.0 edges. A single nx.shortest_path would miss reticulations on tied routes."""
    if not nx.has_path(net, u, v):
        return None
    seen = set()
    out = []
    for path in nx.all_shortest_paths(net, u, v, weight="length"):
        for a, b in nx.utils.pairwise(path):
            key = frozenset((a, b))
            if key not in seen:
                seen.add(key)
                out.append((a, b, net.edges[a, b]["edge_type"], net.edges[a, b]["length"]))
    return out


def top_cycle_edge_paths(sim_id, cycles, top_k, tip_to_tip_only):
    """Per cycle, dicts {edge, weight, path} for top-ranked cycle edges. Ranks by |weight|, keeps top_k plus ties. Only leaf-to-leaf edges when tip_to_tip_only. path lists edges on any shortest network route. path is None if an endpoint is off-graph."""
    net = read_network_graph(sim_id)
    result = []
    for cycle in cycles:
        edges = cycle["edges"]
        if tip_to_tip_only and net is not None:
            edges = [e for e in edges if all(n in net and net.nodes[n]["is_leaf"] for n in e["nodes"])]
        ranked = sorted(edges, key=lambda e: abs(e["weight"]), reverse=True)
        cutoff = abs(ranked[top_k - 1]["weight"]) if len(ranked) >= top_k else 0.0
        entries = []
        for edge in [e for e in ranked if abs(e["weight"]) >= cutoff]:
            u, v = edge["nodes"]
            if net is None or u not in net or v not in net:
                path = None
            else:
                path = shortest_path_edges(net, u, v)
            entries.append({"edge": (u, v), "weight": edge["weight"], "path": path})
        result.append(entries)
    return result


def _rank_cycle_edges(cycle):
    """Rank and closing-flag lookup keyed by appears_at. Shared by reticulate and tip edge reporting."""
    distinct = sorted({edge["appears_at"] for edge in cycle["edges"]})
    rank = {t: r for r, t in enumerate(distinct, start=1)}
    last = distinct[-1] if distinct else None
    return rank, len(distinct), last


def reticulate_edges_per_cycle(cycles, retic_edges):
    if cycles is None or retic_edges is None:
        return None
    inher_weight_by_set = {frozenset(e[:2]): e[2] for e in retic_edges}
    result = []
    for cycle in cycles:
        rank, n, last = _rank_cycle_edges(cycle)
        found = []
        for edge in cycle["edges"]:
            edge_set = frozenset(edge["nodes"])
            if edge_set in inher_weight_by_set:
                found.append({
                    "edge": tuple(edge["nodes"]),
                    "weight": edge["weight"],
                    "inher_weight": inher_weight_by_set[edge_set],
                    "appears_at": edge["appears_at"],
                    "rank": rank[edge["appears_at"]],
                    "num_distinct_appears_at": n,
                    "is_closing": edge["appears_at"] == last,
                })
        result.append(found)
    return result


def read_leaf_labels(sim_id):
    """Labels of current tips from simN_filtered_nodes.csv. Uses live is_leaf, not the sticky 'type' category."""
    path = os.path.join(PHYLO_CSV_DIR, f"sim{sim_id}_filtered_nodes.csv")
    if not os.path.exists(path):
        return None
    with open(path, newline="") as f:
        return {row["label"] for row in csv.DictReader(f) if row["is_leaf"] == "True"}


def tip_edges_per_cycle(cycles, leaf_labels):
    """Edges within each cycle whose both endpoints are leaves. These are the tip-to-tip edges."""
    if cycles is None or leaf_labels is None:
        return None
    result = []
    for cycle in cycles:
        rank, n, last = _rank_cycle_edges(cycle)
        tip_edges = []
        for edge in cycle["edges"]:
            u, v = edge["nodes"]
            if u not in leaf_labels or v not in leaf_labels:
                continue
            tip_edges.append({
                "edge": (u, v),
                "weight": edge["weight"],
                "appears_at": edge["appears_at"],
                "rank": rank[edge["appears_at"]],
                "num_distinct_appears_at": n,
                "is_closing": edge["appears_at"] == last,
            })
        result.append(tip_edges)
    return result


def write_sim_report(sim_id, top_k, tip_to_tip_only):
    retic_edges = read_reticulate_edges_by_name(sim_id)
    cycles = read_cycles_by_name(sim_id)
    found_per_cycle = reticulate_edges_per_cycle(cycles, retic_edges)
    leaf_labels = read_leaf_labels(sim_id)
    tip_edges_per_cycle_result = tip_edges_per_cycle(cycles, leaf_labels)
    if retic_edges is None or cycles is None or found_per_cycle is None or tip_edges_per_cycle_result is None:
        return None

    lines = []
    lines.append(f"sim{sim_id} analysis")
    lines.append("")
    lines.append(f"reticulate edges ({len(retic_edges)}):")
    for edge in retic_edges:
        lines.append(f"  {edge[0]} -- {edge[1]}  inher_weight={edge[2]}")
    lines.append("")
    lines.append(f"cycles ({len(cycles)}):")
    for i in range(len(cycles)):
        cycle = cycles[i]
        found = found_per_cycle[i]
        tip_edges = tip_edges_per_cycle_result[i]
        lines.append(f"  cycle {i}: {len(cycle['edges'])} edges, birth={cycle['birth']}, death={cycle['death']}, {len(found)} reticulate edges found")
        for edge in found:
            closing_flag = "  CLOSING EDGE" if edge["is_closing"] else ""
            lines.append(
                f"    {edge['edge']}  weight={edge['weight']}  inher_weight={edge['inher_weight']}  "
                f"rank={edge['rank']}/{edge['num_distinct_appears_at']}{closing_flag}"
            )
        lines.append(f"    tip-to-tip edges ({len(tip_edges)}):")
        for edge in tip_edges:
            closing_flag = "  CLOSING EDGE" if edge["is_closing"] else ""
            lines.append(
                f"      {edge['edge']}  weight={edge['weight']}  "
                f"rank={edge['rank']}/{edge['num_distinct_appears_at']}{closing_flag}"
            )

    top_paths = top_cycle_edge_paths(sim_id, cycles, top_k, tip_to_tip_only)
    lines.append("")
    lines.append(f"top-{top_k} cycle-edge shortest paths (all edges on any co-shortest network path, incl. reticulation edges, weight=length):")
    for i, entries in enumerate(top_paths):
        if not entries:
            lines.append(f"  cycle {i}: no nonzero edges")
            continue
        lines.append(f"  cycle {i}:")
        for entry in entries:
            u, v = entry["edge"]
            lines.append(f"    {u} -- {v}  weight={entry['weight']}")
            if entry["path"] is None:
                lines.append("      no path (endpoint missing from network graph)")
                continue
            for a, b, edge_type, length in entry["path"]:
                marker = "  <-- reticulation" if edge_type == "reticulation" else ""
                lines.append(f"      {a} -- {b}  ({edge_type}, length={length}){marker}")

    reticulation_on_paths = {
        frozenset((a, b))
        for entries in top_paths
        for entry in entries if entry["path"]
        for a, b, edge_type, _ in entry["path"]
        if edge_type == "reticulation"
    }
    involved = [edge for edge in retic_edges if frozenset(edge[:2]) in reticulation_on_paths]
    missing = [edge for edge in retic_edges if frozenset(edge[:2]) not in reticulation_on_paths]
    lines.append("")
    lines.append(f"reticulate edges on at least one top-{top_k} cycle-edge path ({len(involved)}):")
    for edge in involved:
        lines.append(f"  {edge[0]} -- {edge[1]}  inher_weight={edge[2]}")
    lines.append(f"reticulate edges on no top-{top_k} cycle-edge path ({len(missing)}):")
    for edge in missing:
        lines.append(f"  {edge[0]} -- {edge[1]}  inher_weight={edge[2]}")

    os.makedirs(ANALYSIS_REPORTS_DIR, exist_ok=True)
    path = os.path.join(ANALYSIS_REPORTS_DIR, f"sim{sim_id}_analysis.txt")
    with open(path, "w") as f:
        for line in lines:
            f.write(line + "\n")
    return path


def main():
    # largest-|weight| cycle edges to trace paths for, per cycle.
    top_k = 3
    # rank only leaf-to-leaf cycle edges. Internal or phantom endpoints are possible.
    tip_to_tip_only = True
    sim_ids = discover_sim_ids()
    for sim_id in sim_ids:
        path = write_sim_report(sim_id, top_k, tip_to_tip_only)
        print(f"sim{sim_id}: wrote {path}")


if __name__ == "__main__":
    main()


