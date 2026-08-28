import csv
import json
import math
import os
import re
import statistics

import networkx as nx

from main import PHYLO_CSV_DIR, SIM_OUTPUTS_DIR

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
    """Undirected nx.Graph from simN_filtered_edges.csv (tree + reticulation edges), weight 'length', edges carry 'edge_type' -- the graph the dressed tip-tip cells are shortest paths on."""
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
    return G


def reticulate_edge_shortening(sim_id):
    """{frozenset({label_u, label_v}): (n_pairs, max_shortening)} over the dressed tip pairs whose shortest path crosses each reticulate edge."""
    log_path = os.path.join(PROC_OUTPUTS_DIR, f"sim{sim_id}", "find_cycles_log.json")
    retic_edges = read_reticulate_edges_by_name(sim_id)
    if not os.path.exists(log_path) or retic_edges is None:
        return None
    with open(log_path) as f:
        entries = json.load(f).get("dress_distance_matrix", [])
    result = {frozenset(edge[:2]): (0, 0.0) for edge in retic_edges}
    for entry in entries:
        ratio = entry["ratio"]
        if ratio is None:
            continue
        for retic in entry.get("reticulation_edges_on_path", []):
            key = frozenset(retic["edge"])
            n_pairs, max_shortening = result[key]
            result[key] = (n_pairs + 1, min(max_shortening, ratio))
    return result


def reticulate_edge_table_lines(retic_edges, shortening, reticulation_on_paths):
    """Report table: every reticulate edge with n_pairs, max_shortening, inher_weight, and whether it lies on a cycle closing-edge path."""
    if shortening is None:
        return ["reticulate-edge table: no find_cycles log"]
    lines = [
        "reticulate-edge table:",
        f"  {'edge':<22}{'n_pairs':>9}{'max_shortening':>16}{'inher_weight':>14}{'on_closing_path':>17}",
    ]
    rows = []
    for u, v, inher in retic_edges:
        key = frozenset((u, v))
        n_pairs, max_shortening = shortening.get(key, (0, 0.0))
        rows.append((f"{u} -- {v}", n_pairs, max_shortening, inher, key in reticulation_on_paths))
    for name, n_pairs, max_shortening, inher, on_path in sorted(rows, key=lambda r: r[2]):
        inher_str = f"{inher:.4f}" if inher is not None else "None"
        lines.append(f"  {name:<22}{n_pairs:>9}{max_shortening:>16.6f}{inher_str:>14}{str(on_path):>17}")
    return lines


def top_cycle_edge_paths(sim_id, cycles, top_k):
    """Per cycle, up to top_k dicts {edge, weight, path} for the largest-|weight| edges; path is the shortest network path's (a, b, edge_type, length) edges, or None if an endpoint is off-graph."""
    net = read_network_graph(sim_id)
    result = []
    for cycle in cycles:
        ranked = sorted(cycle["edges"], key=lambda e: abs(e["weight"]), reverse=True)
        entries = []
        for edge in ranked[:top_k]:
            u, v = edge["nodes"]
            if net is None or u not in net or v not in net:
                path = None
            else:
                nodes = nx.shortest_path(net, u, v, weight="length")
                path = [(a, b, net.edges[a, b]["edge_type"], net.edges[a, b]["length"]) for a, b in nx.utils.pairwise(nodes)]
            entries.append({"edge": (u, v), "weight": edge["weight"], "path": path})
        result.append(entries)
    return result


def _rank_cycle_edges(cycle):
    """rank/closing-flag lookup by appears_at, shared by reticulate/tip edge reporting."""
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
    path = os.path.join(PHYLO_CSV_DIR, f"sim{sim_id}_nodes.csv")
    if not os.path.exists(path):
        return None
    with open(path, newline="") as f:
        return {row["label"] for row in csv.DictReader(f) if row["type"] in ("leaf", "hyb_leaf")}


def tip_edges_per_cycle(cycles, leaf_labels):
    """edges within each cycle whose both endpoints are leaves (is_leaf), i.e. tip-to-tip edges."""
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


def write_sim_report(sim_id, top_k):
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

    dress_ratios = read_dress_change_ratios(sim_id)
    lines.append("")
    if dress_ratios is None:
        lines.append("dress_distance_matrix change ratios (new_dist - old_dist) / old_dist: no log")
    else:
        median = statistics.median(dress_ratios) if dress_ratios else None
        dress_max = max(dress_ratios) if dress_ratios else None
        dress_min = min(dress_ratios) if dress_ratios else None
        lines.append(f"dress_distance_matrix change ratios (new_dist - old_dist) / old_dist ({len(dress_ratios)}):")
        lines.append(f"  median={median}  max={dress_max}  min={dress_min}")

    top_paths = top_cycle_edge_paths(sim_id, cycles, top_k)
    lines.append("")
    lines.append(f"top-{top_k} cycle-edge shortest paths (network incl. reticulation edges, weight=length):")
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

    lines.append("")
    lines.extend(reticulate_edge_table_lines(retic_edges, reticulate_edge_shortening(sim_id), reticulation_on_paths))

    os.makedirs(ANALYSIS_REPORTS_DIR, exist_ok=True)
    path = os.path.join(ANALYSIS_REPORTS_DIR, f"sim{sim_id}_analysis.txt")
    with open(path, "w") as f:
        for line in lines:
            f.write(line + "\n")
    return path


def read_dress_change_ratios(sim_id):
    """Dressed (new_dist - old_dist) / old_dist ratios from simN find_cycles_log.json (old_dist == 0 cells are stored as None and skipped)."""
    path = os.path.join(PROC_OUTPUTS_DIR, f"sim{sim_id}", "find_cycles_log.json")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        data = json.load(f)
    ratios = []
    for entry in data.get("dress_distance_matrix", []):
        if entry["ratio"] is not None:
            ratios.append(entry["ratio"])
    return ratios


def summarize_dress_change_ratios(sim_ids=None):
    """Aggregate read_dress_change_ratios across sims into {ratios, median, max, min}."""
    if sim_ids is None:
        sim_ids = discover_sim_ids()
    ratios = []
    for sim_id in sim_ids:
        sim_ratios = read_dress_change_ratios(sim_id)
        if sim_ratios:
            ratios.extend(sim_ratios)
    return {
        "ratios": ratios,
        "median": statistics.median(ratios) if ratios else None,
        "max": max(ratios) if ratios else None,
        "min": min(ratios) if ratios else None,
    }


def main():
    # per cycle, how many largest-|harmonic weight| edges to trace shortest paths for (1 = closing edge only; 5 added no extra reticulations)
    top_k = 1
    sim_ids = discover_sim_ids()
    for sim_id in sim_ids:
        path = write_sim_report(sim_id, top_k)
        print(f"sim{sim_id}: wrote {path}")


if __name__ == "__main__":
    main()


