import csv
import itertools
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


def read_cycles(sim_id, weight_field="weight"):
    """Cycles from rip.json, each edge's 'weight' taken from weight_field."""
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
            if math.isclose(edge[weight_field], 0.0, abs_tol=1e-6):
                continue
            edges.append({
                "simplex": edge["simplex"],
                "weight": edge[weight_field],
                "appears_at": appears_at[frozenset(edge["simplex"])],
            })
        history = {tuple(e["simplex"]): e["history"] for e in cycle.get("weight_history", [])}
        cycles.append({"birth": cycle["birth"], "death": cycle["death"], "edges": edges, "history": history})
    return cycles


def read_cycles_by_name(sim_id, weight_field="weight"):
    index_to_name = read_index_to_name(sim_id)
    cycles = read_cycles(sim_id, weight_field)
    if index_to_name is None or cycles is None:
        return None
    result = []
    for cycle in cycles:
        named_edges = []
        for edge in cycle["edges"]:
            nodes = [index_to_name[idx] for idx in edge["simplex"]]
            named_edges.append({"nodes": nodes, "weight": edge["weight"], "appears_at": edge["appears_at"]})
        history = {tuple(index_to_name[idx] for idx in simplex): points for simplex, points in cycle["history"].items()}
        result.append({"birth": cycle["birth"], "death": cycle["death"], "edges": named_edges, "history": history})
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
    """Undirected network from simN_filtered_edges.csv with length, edge_type and is_leaf."""
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
    """Edges on any tied shortest u-v path as (a, b, edge_type, length); None if unreachable."""
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


def tree_only_graph(net):
    """Network with its reticulation edges removed."""
    tree = net.copy()
    tree.remove_edges_from([(a, b) for a, b, edge_type in net.edges(data="edge_type") if edge_type == "reticulation"])
    return tree


def route_reticulations(net, tree, u, v, max_routes):
    """Up to max_routes u-v routes, each removing the previous routes' reticulations, kept while shorter than tree-only."""
    tree_dist = nx.shortest_path_length(tree, u, v, weight="length") if nx.has_path(tree, u, v) else math.inf
    G = net.copy()
    routes = []
    while len(routes) < max_routes and nx.has_path(G, u, v):
        length = nx.shortest_path_length(G, u, v, weight="length")
        path = shortest_path_edges(G, u, v)
        retics = {frozenset((a, b)) for a, b, edge_type, _ in path if edge_type == "reticulation"}
        if routes and (not retics or length >= tree_dist - 1e-9):
            break
        routes.append({"path": path, "length": length, "shrink": length / tree_dist - 1 if tree_dist > 0 else 0.0, "retics": retics})
        if not retics:
            break
        G.remove_edges_from([tuple(r) for r in retics])
    return routes


def _edge_on_shortest_path(dist, u, v, a, b, length, tol):
    """True if d(u, a) + length + d(b, v) equals d(u, v), in either direction."""
    if v not in dist[u]:
        return False
    d_uv = dist[u][v]
    via_a_then_b = dist[u].get(a, math.inf) + length + dist[v].get(b, math.inf)
    via_b_then_a = dist[u].get(b, math.inf) + length + dist[v].get(a, math.inf)
    return (math.isclose(via_a_then_b, d_uv, rel_tol=tol, abs_tol=tol)
            or math.isclose(via_b_then_a, d_uv, rel_tol=tol, abs_tol=tol))


def reticulate_edges_off_tip_paths(net, retic_edges, tol=1e-9):
    """Reticulate edges on no shortest path between any two tips."""
    if net is None or retic_edges is None:
        return None
    leaves = [n for n in net if net.nodes[n]["is_leaf"]]
    dist = {u: nx.single_source_dijkstra_path_length(net, u, weight="length") for u in leaves}

    off = []
    for edge in retic_edges:
        a, b = edge[:2]
        if not net.has_edge(a, b) or net.edges[a, b]["edge_type"] != "reticulation":
            off.append(edge)
            continue
        length = net.edges[a, b]["length"]
        on_any_tip_path = any(
            _edge_on_shortest_path(dist, u, v, a, b, length, tol)
            for u, v in itertools.combinations(leaves, 2)
        )
        if not on_any_tip_path:
            off.append(edge)
    return off


def history_leaders(cycle, is_tip, tol=1e-9):
    """Tip-tip edges with top |weight| at any snapshot, with first lead time, lead share and weight."""
    history = {edge: points for edge, points in cycle["history"].items() if all(is_tip(n) for n in edge)}
    times = sorted({t for points in history.values() for t, _ in points})
    if not times:
        return []
    end = cycle["death"] if cycle["death"] is not None else times[-1]
    lifetime = end - cycle["birth"]
    leaders = {}
    for i, t in enumerate(times):
        state = {edge: abs([w for s, w in points if s <= t][-1]) for edge, points in history.items() if points[0][0] <= t}
        best = max(state.values(), default=0.0)
        if best <= tol:
            continue
        span = (times[i + 1] if i + 1 < len(times) else end) - t
        for edge, w in state.items():
            if w >= best - tol:
                lead = leaders.setdefault(edge, {"edge": edge, "first": t, "lead": 0.0, "weight": w})
                lead["lead"] += span / lifetime if lifetime > 0 else 1.0
    return sorted(leaders.values(), key=lambda lead: lead["first"])


def top_cycle_edge_paths(sim_id, cycles, top_k, tip_to_tip_only, max_routes=1, leaders="birth"):
    """Per cycle, top-k cycle edges by |weight| with their traced shortest paths and routes."""
    net = read_network_graph(sim_id)
    tree = tree_only_graph(net) if net is not None else None
    mid_tips = read_midpoints(sim_id)
    result = []
    for cycle in cycles:
        edges = cycle["edges"]
        if tip_to_tip_only and net is not None:
            edges = [e for e in edges if all(n in mid_tips or (n in net and net.nodes[n]["is_leaf"]) for n in e["nodes"])]
        if leaders == "history" and cycle["history"] and net is not None:
            is_tip = lambda n: n in mid_tips or (n in net and net.nodes[n]["is_leaf"])
            chosen = [{"nodes": list(lead["edge"]), "weight": lead["weight"], "first": lead["first"], "lead": lead["lead"]}
                      for lead in history_leaders(cycle, is_tip)]
        else:
            ranked = sorted(edges, key=lambda e: abs(e["weight"]), reverse=True)
            cutoff = abs(ranked[top_k - 1]["weight"]) if len(ranked) >= top_k else 0.0
            chosen = [e for e in ranked if abs(e["weight"]) >= cutoff]
        entries = []
        for edge in chosen:
            u, v = edge["nodes"]
            traced = [mid_tips[n] for n in (u, v) if n in mid_tips] or [(u, v)]
            parts = None if net is None or any(n not in net for pair in traced for n in pair) else [shortest_path_edges(net, a, b) for a, b in traced]
            if parts is None or None in parts:
                path = None
                routes = {}
            else:
                path = []
                for part in parts:
                    path += [p for p in part if p not in path]
                routes = {(a, b): route_reticulations(net, tree, a, b, max_routes) for a, b in traced}
            entries.append({"edge": (u, v), "weight": edge["weight"], "path": path, "traced": traced, "routes": routes,
                            "first": edge.get("first"), "lead": edge.get("lead")})
        result.append(entries)
    return result


def _rank_cycle_edges(cycle):
    """Rank and closing-flag lookup keyed by appears_at."""
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


def read_midpoints(sim_id):
    """Mid-tip names from midpoints.json mapped to their two parent tip labels."""
    path = os.path.join(PROC_OUTPUTS_DIR, f"sim{sim_id}", "midpoints.json")
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        return {m["name"]: tuple(m["parents"]) for m in json.load(f) if m["type"] == "mid-tip"}


def read_leaf_labels(sim_id):
    """Live tip labels from simN_filtered_nodes.csv, plus mid-tips."""
    path = os.path.join(PHYLO_CSV_DIR, f"sim{sim_id}_filtered_nodes.csv")
    if not os.path.exists(path):
        return None
    with open(path, newline="") as f:
        return {row["label"] for row in csv.DictReader(f) if row["is_leaf"] == "True"} | set(read_midpoints(sim_id))


def tip_edges_per_cycle(cycles, leaf_labels):
    """Cycle edges whose endpoints are both tips."""
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


def write_sim_report(sim_id, top_k, tip_to_tip_only, weight_field="weight", max_routes=1, leaders="birth"):
    retic_edges = read_reticulate_edges_by_name(sim_id)
    cycles = read_cycles_by_name(sim_id, weight_field)
    found_per_cycle = reticulate_edges_per_cycle(cycles, retic_edges)
    leaf_labels = read_leaf_labels(sim_id)
    tip_edges_per_cycle_result = tip_edges_per_cycle(cycles, leaf_labels)
    if retic_edges is None or cycles is None or found_per_cycle is None or tip_edges_per_cycle_result is None:
        return None
    mid_tips = read_midpoints(sim_id)

    lines = []
    lines.append(f"sim{sim_id} analysis")
    lines.append(f"weight = {weight_field}")
    lines.append(f"leaders = {leaders}, max_routes = {max_routes}")
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
            mid_flag = "  (mid-tip)" if any(n in mid_tips for n in edge["edge"]) else ""
            lines.append(
                f"      {edge['edge']}  weight={edge['weight']}  "
                f"rank={edge['rank']}/{edge['num_distinct_appears_at']}{closing_flag}{mid_flag}"
            )

    top_paths = top_cycle_edge_paths(sim_id, cycles, top_k, tip_to_tip_only, max_routes, leaders)
    lines.append("")
    lines.append(f"top-{top_k} cycle-edge shortest paths (all edges on any co-shortest network path, incl. reticulation edges, weight=length):")
    for i, entries in enumerate(top_paths):
        if not entries:
            lines.append(f"  cycle {i}: no nonzero edges")
            continue
        lines.append(f"  cycle {i}:")
        for entry in entries:
            u, v = entry["edge"]
            lead = f"  leads from t={entry['first']:.4f} for {entry['lead']:.0%} of lifetime" if entry["first"] is not None else ""
            lines.append(f"    {u} -- {v}  weight={entry['weight']}{lead}")
            if entry["traced"] != [(u, v)]:
                lines.append("      (mid-tip: traced parent path " + ", ".join(f"{a} -- {b}" for a, b in entry["traced"]) + ")")
            if entry["path"] is None:
                lines.append("      no path (endpoint missing from network graph)")
                continue
            for a, b, edge_type, length in entry["path"]:
                marker = "  <-- reticulation" if edge_type == "reticulation" else ""
                lines.append(f"      {a} -- {b}  ({edge_type}, length={length}){marker}")
            for (a, b), routes in entry["routes"].items():
                for n, route in enumerate(routes[1:], start=2):
                    retics = ", ".join(" -- ".join(sorted(r)) for r in route["retics"])
                    lines.append(f"      route {n} {a} -- {b}: length={route['length']:.4f}  shrink={route['shrink']:+.1%}  via {retics}")

    reticulation_on_paths = {
        frozenset((a, b))
        for entries in top_paths
        for entry in entries if entry["path"]
        for a, b, edge_type, _ in entry["path"]
        if edge_type == "reticulation"
    } | {
        r
        for entries in top_paths
        for entry in entries
        for routes in entry["routes"].values()
        for route in routes
        for r in route["retics"]
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

    off_tip_paths = reticulate_edges_off_tip_paths(read_network_graph(sim_id), retic_edges)
    if off_tip_paths is not None:
        lines.append(f"reticulate edges on no tip-to-tip shortest path ({len(off_tip_paths)}):")
        for edge in off_tip_paths:
            lines.append(f"  {edge[0]} -- {edge[1]}  inher_weight={edge[2]}")

    os.makedirs(ANALYSIS_REPORTS_DIR, exist_ok=True)
    path = os.path.join(ANALYSIS_REPORTS_DIR, f"sim{sim_id}_analysis_{weight_field}.txt")
    with open(path, "w") as f:
        for line in lines:
            f.write(line + "\n")
    return path


def main():
    top_k = 1
    tip_to_tip_only = True
    weight_field = "weight_before_projection"
    max_routes = 3
    leaders = "history"
    sim_ids = discover_sim_ids()
    for sim_id in sim_ids:
        path = write_sim_report(sim_id, top_k, tip_to_tip_only, weight_field, max_routes, leaders)
        print(f"sim{sim_id}: wrote {path}")


if __name__ == "__main__":
    main()


