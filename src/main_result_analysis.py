import csv
import json
import math
import os
import re

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


def read_reticulate_edges_by_name(sim_id):
    path = os.path.join(PHYLO_CSV_DIR, f"sim{sim_id}_filtered_edges.csv")
    if not os.path.exists(path):
        return None
    with open(path, newline="") as f:
        return [
            (row["from"], row["to"])
            for row in csv.DictReader(f)
            if row["edge_type"] == "reticulation"
        ]


def reticulate_edges_per_cycle(sim_id):
    cycles = read_cycles_by_name(sim_id)
    retic_edges = read_reticulate_edges_by_name(sim_id)
    if cycles is None or retic_edges is None:
        return None
    retic_sets = [frozenset(e) for e in retic_edges]
    result = []
    for cycle in cycles:
        info_by_edge_set = {}
        for edge in cycle["edges"]:
            info_by_edge_set[frozenset(edge["nodes"])] = edge

        distinct_appears_at = sorted(set(edge["appears_at"] for edge in cycle["edges"]))
        rank_by_appears_at = {}
        for rank, value in enumerate(distinct_appears_at, start=1):
            rank_by_appears_at[value] = rank
        last_appears_at = distinct_appears_at[-1] if distinct_appears_at else None

        found = []
        for retic_set in retic_sets:
            if retic_set in info_by_edge_set:
                edge = info_by_edge_set[retic_set]
                found.append({
                    "edge": tuple(retic_set),
                    "weight": edge["weight"],
                    "appears_at": edge["appears_at"],
                    "rank": rank_by_appears_at[edge["appears_at"]],
                    "num_distinct_appears_at": len(distinct_appears_at),
                    "is_closing": edge["appears_at"] == last_appears_at,
                })
        result.append(found)
    return result


def write_sim_report(sim_id):
    retic_edges = read_reticulate_edges_by_name(sim_id)
    cycles = read_cycles_by_name(sim_id)
    found_per_cycle = reticulate_edges_per_cycle(sim_id)
    if retic_edges is None or cycles is None or found_per_cycle is None:
        return None

    lines = []
    lines.append(f"sim{sim_id} analysis")
    lines.append("")
    lines.append(f"reticulate edges ({len(retic_edges)}):")
    for edge in retic_edges:
        lines.append(f"  {edge[0]} -- {edge[1]}")
    lines.append("")
    lines.append(f"cycles ({len(cycles)}):")
    for i in range(len(cycles)):
        cycle = cycles[i]
        found = found_per_cycle[i]
        lines.append(f"  cycle {i}: {len(cycle['edges'])} edges, birth={cycle['birth']}, death={cycle['death']}, {len(found)} reticulate edges found")
        for edge in found:
            closing_flag = "  CLOSING EDGE" if edge["is_closing"] else ""
            lines.append(
                f"    {edge['edge']}  weight={edge['weight']}  "
                f"rank={edge['rank']}/{edge['num_distinct_appears_at']}{closing_flag}"
            )

    found_anywhere = set()
    for found in found_per_cycle:
        for edge in found:
            found_anywhere.add(frozenset(edge["edge"]))

    never_found = []
    for edge in retic_edges:
        if frozenset(edge) not in found_anywhere:
            never_found.append(edge)

    lines.append("")
    lines.append(f"reticulate edges never found in any cycle ({len(never_found)}):")
    for edge in never_found:
        lines.append(f"  {edge[0]} -- {edge[1]}")

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


