import csv
import json
import math
import os
import re

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


def read_hyb_leaf_labels(sim_id):
    path = os.path.join(PHYLO_CSV_DIR, f"sim{sim_id}_nodes.csv")
    if not os.path.exists(path):
        return None
    with open(path, newline="") as f:
        return {row["label"] for row in csv.DictReader(f) if row["type"] == "hyb_leaf"}


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


def write_sim_report(sim_id):
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

    found_anywhere = set()
    for found in found_per_cycle:
        for edge in found:
            found_anywhere.add(frozenset(edge["edge"]))

    never_found = []
    for edge in retic_edges:
        if frozenset(edge[:2]) not in found_anywhere:
            never_found.append(edge)

    lines.append("")
    lines.append(f"reticulate edges never found in any cycle ({len(never_found)}):")
    for edge in never_found:
        lines.append(f"  {edge[0]} -- {edge[1]}  inher_weight={edge[2]}")

    tip_nodes = set()
    for tip_edges in tip_edges_per_cycle_result:
        for edge in tip_edges:
            u, v = edge["edge"]
            tip_nodes.add(u)
            tip_nodes.add(v)

    lines.append("")
    lines.append(f"nodes appearing in a tip-to-tip edge ({len(tip_nodes)}): {sorted(tip_nodes)}")

    hyb_leaf_labels = read_hyb_leaf_labels(sim_id)
    missing_hyb = sorted(hyb_leaf_labels - tip_nodes) if hyb_leaf_labels is not None else []
    lines.append(f"all hybrid species appear in a tip-to-tip edge: {not missing_hyb}")
    if missing_hyb:
        lines.append(f"  missing hybrid species ({len(missing_hyb)}): {missing_hyb}")

    os.makedirs(ANALYSIS_REPORTS_DIR, exist_ok=True)
    path = os.path.join(ANALYSIS_REPORTS_DIR, f"sim{sim_id}_analysis.txt")
    with open(path, "w") as f:
        for line in lines:
            f.write(line + "\n")
    return path


def main():
    sim_ids = discover_sim_ids()
    for sim_id in sim_ids:
        path = write_sim_report(sim_id)
        print(f"sim{sim_id}: wrote {path}")


if __name__ == "__main__":
    main()


