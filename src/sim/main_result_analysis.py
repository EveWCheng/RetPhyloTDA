import csv
import json
import os
import re

import networkx as nx

from sim_main import CONFIG_PATH, SIM_OUTPUTS_DIR, load_config


class MainResultAnalysis:

    def __init__(self, sim_outputs_dir=SIM_OUTPUTS_DIR, config_path=CONFIG_PATH, sel_op="all"):
        self.sim_outputs_dir = sim_outputs_dir
        self.config_path = config_path
        self.sel_op = sel_op
        self.phylo_csv_dir = os.path.join(sim_outputs_dir, "phylo_csv")
        self.proc_outputs_dir = os.path.join(sim_outputs_dir, "proc_phylo_outputs")
        self.cycle_outputs_dir = os.path.join(sim_outputs_dir, "cycle_outputs")
        self.analysis_reports_dir = os.path.join(sim_outputs_dir, "analysis_reports")

    def read_config(self):
        return load_config(self.config_path)

    def read_index_to_name(self, sim_id):
        """Return {matrix index: name} from sim{sim_id}'s populated_headers.txt."""
        path = os.path.join(self.proc_outputs_dir, f"sim{sim_id}", "populated_headers.txt")
        index_to_name = {}
        with open(path) as f:
            for index, line in enumerate(f):
                index_to_name[index] = line.rstrip("\n")
        return index_to_name

    def read_weight_histories(self, sim_id):
        """Return {cycle_index: {(i, j): [(t, weight), ...]}} from sim{sim_id}'s rip.json."""
        path = os.path.join(self.cycle_outputs_dir, f"sim{sim_id}", "rip.json")
        with open(path) as f:
            data = json.load(f)
        histories = {}
        for cycle in data["harmonic_cycles"]:
            if "weight_history" not in cycle:
                raise KeyError(f"sim{sim_id} has cycles but no weight_history; something weird has happened")
            edges = {}
            for edge in cycle["weight_history"]:
                i, j = edge["simplex"]
                points = []
                for t, weight in edge["history"]:
                    points.append((t, weight))
                edges[(i, j)] = points
            histories[cycle["cycle_index"]] = edges
        return histories

    def read_weights_before_projection(self, sim_id):
        """Return {cycle_index: {(i, j): weight_before_projection}} from sim{sim_id}'s rip.json."""
        path = os.path.join(self.cycle_outputs_dir, f"sim{sim_id}", "rip.json")
        with open(path) as f:
            data = json.load(f)
        weights = {}
        for cycle in data["harmonic_cycles"]:
            edges = {}
            for edge in cycle["edges"]:
                i, j = edge["simplex"]
                edges[(i, j)] = edge["weight_before_projection"]
            weights[cycle["cycle_index"]] = edges
        return weights

    def sim_ids(self):
        ids = []
        for name in os.listdir(self.cycle_outputs_dir):
            match = re.fullmatch(r"sim(\d+)", name)
            if match:
                ids.append(int(match.group(1)))
        return sorted(ids)

    def read_all_weight_histories(self):
        all_histories = {}
        for sim_id in self.sim_ids():
            all_histories[sim_id] = self.read_weight_histories(sim_id)
        return all_histories

    def read_all_weights(self):
        all_weights = {}
        for sim_id in self.sim_ids():
            all_weights[sim_id] = self.read_weights_before_projection(sim_id)
        return all_weights

    def read_filtered_graph(self, sim_id):
        sim_csv_dir = os.path.join(self.phylo_csv_dir, f"sim{sim_id}")
        D = nx.DiGraph()
        with open(os.path.join(sim_csv_dir, f"sim{sim_id}_filtered_nodes.csv"), newline="") as f:
            for row in csv.DictReader(f):
                attrs = {"label": row["label"], "type": row["type"], "is_leaf": row["is_leaf"] == "True"}
                D.add_node(row["id"], **attrs)
        with open(os.path.join(sim_csv_dir, f"sim{sim_id}_filtered_edges.csv"), newline="") as f:
            for row in csv.DictReader(f):
                attrs = {
                    "edge_type": row["edge_type"],
                    "length": float(row["length"]),
                    "time_length": float(row["time_length"]),
                }
                if row["inher_weight"] != "":
                    attrs["inher_weight"] = float(row["inher_weight"])
                D.add_edge(row["from"], row["to"], **attrs)
        return D

    def read_all_filtered_graphs(self):
        graphs = {}
        for sim_id in self.sim_ids():
            graphs[sim_id] = self.read_filtered_graph(sim_id)
        return graphs

    def history_leaders(self, sim_id, tol=1e-9):
        history_leaders = self._history_leaders_(sim_id, tol)
        leaders = {}
        for cycle in history_leaders.keys():
            cycle_leaders = history_leaders[cycle]
            if len(cycle_leaders) == 0:
                continue
            if self.sel_op == "only_at_birth":
                t, highest, tied = cycle_leaders[0]
                leaders[cycle] = [(highest, tied)]
            elif self.sel_op == "all":
                without_t = []
                for t, highest, tied in cycle_leaders:
                    without_t.append((highest, tied))
                leaders[cycle] = without_t
        return leaders

    def _history_leaders_(self, sim_id, tol=1e-9):
        """Return {cycle_index: [(t, max |weight|, [edges tied for it]), ...]} over each cycle's history."""
        leaders = {}
        for cycle_index, edges in self.read_weight_histories(sim_id).items():
            times = [t for t, _ in next(iter(edges.values()))]
            for edge, points in edges.items():
                if len(points) != len(times):
                    raise ValueError(f"sim{sim_id} cycle {cycle_index} edge {edge} history is not aligned")

            cycle_leaders = []
            for k, t in enumerate(times):
                highest = 0.0
                for points in edges.values():
                    if abs(points[k][1]) > highest:
                        highest = abs(points[k][1])
                if highest <= tol:
                    continue
                tied = []
                for edge, points in edges.items():
                    if abs(points[k][1]) >= highest - tol:
                        tied.append(edge)
                cycle_leaders.append((t, highest, tied))
            leaders[cycle_index] = cycle_leaders
        return leaders

    def weight_leaders(self, sim_id, tol=1e-9):
        leaders = {}
        for cycle_index, edges in self.read_weights_before_projection(sim_id).items():
            highest = 0.0
            for weight in edges.values():
                if abs(weight) > highest:
                    highest = abs(weight)
            if highest <= tol:
                continue
            tied = []
            for edge, weight in edges.items():
                if abs(weight) >= highest - tol:
                    tied.append(edge)
            leaders[cycle_index] = [(highest, tied)]
        return leaders

    def traced_pair(self, u, v):
        nodes = set()
        for name in (u, v):
            if "_" in name:
                parent_a, parent_b = name.split("_")
                nodes.add(parent_a)
                nodes.add(parent_b)
            else:
                nodes.add(name)
        if len(nodes) != 2:
            return None
        return tuple(sorted(nodes))

    def ancestor_reticulations(self, D, u, v):
        tree = D.copy()
        for a, b, edge_type in D.edges(data="edge_type"):
            if edge_type == "reticulation":
                tree.remove_edge(a, b)
        u_ancestors = nx.ancestors(tree, u)
        v_ancestors = nx.ancestors(tree, v)

        u_lineage = {u}
        for node in u_ancestors:
            if node not in v_ancestors:
                u_lineage.add(node)

        v_lineage = {v}
        for node in v_ancestors:
            if node not in u_ancestors:
                v_lineage.add(node)

        reticulations = []
        for a, b, edge_type in D.edges(data="edge_type"):
            if edge_type != "reticulation":
                continue
            if (a in u_lineage and b in v_lineage) or (a in v_lineage and b in u_lineage):
                reticulations.append((a, b))
        return reticulations

    def node_kind(self, D, node):
        """Return "tip" if the node is a leaf in D, otherwise "internal"."""
        if D.nodes[node]["is_leaf"]:
            return "tip"
        return "internal"

    def from_pairs_to_ret_each_sim(self, sim_id):
        D = self.read_filtered_graph(sim_id)
        index_to_name = self.read_index_to_name(sim_id)
        leaders = self.history_leaders(sim_id)

        leader_edges = []
        for cycle_leaders in leaders.values():
            for highest, tied in cycle_leaders:
                leader_edges.extend(tied)

        detected = {}
        for i, j in leader_edges:
            pair = self.traced_pair(index_to_name[i], index_to_name[j])
            if pair is None:
                continue
            for edge in self.ancestor_reticulations(D, pair[0], pair[1]):
                detected.setdefault(edge, []).append(pair)

        missed = []
        for a, b, edge_type in D.edges(data="edge_type"):
            if edge_type == "reticulation" and (a, b) not in detected:
                missed.append((a, b))
        return detected, missed


def main():
    analysis = MainResultAnalysis()
    lines = [f"sel_op: {analysis.sel_op}"]
    total_detected = 0
    total_missed = 0
    for sim_id in analysis.sim_ids():
        detected, missed = analysis.from_pairs_to_ret_each_sim(sim_id)
        D = analysis.read_filtered_graph(sim_id)
        total_detected += len(detected)
        total_missed += len(missed)
        lines.append(f"sim{sim_id}: {len(detected)} detected, {len(missed)} missed")
        for (a, b), pairs in detected.items():
            pair_texts = []
            for u, v in pairs:
                text = f"{u}({analysis.node_kind(D, u)})--{v}({analysis.node_kind(D, v)})"
                if text not in pair_texts:
                    pair_texts.append(text)
            lines.append(f"  detected {a}({analysis.node_kind(D, a)}) -- {b}({analysis.node_kind(D, b)})  by {', '.join(pair_texts)}")
        for a, b in missed:
            lines.append(f"  missed   {a} -- {b}")

    total = total_detected + total_missed
    lines.append("")
    if total == 0:
        lines.append("no reticulations in any simulation")
    else:
        lines.append(f"overall: {total_detected} of {total} reticulations detected ({total_detected / total:.1%})")

    os.makedirs(analysis.analysis_reports_dir, exist_ok=True)
    path = os.path.join(analysis.analysis_reports_dir, "detection_report.txt")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
