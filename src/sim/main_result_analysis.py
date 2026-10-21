import csv
import json
import os
import re

import networkx as nx

from sim_main import CONFIG_PATH, SIM_OUTPUTS_DIR, load_config


class MainResultAnalysis:

    def __init__(self, sim_outputs_dir=SIM_OUTPUTS_DIR, config_path=CONFIG_PATH):
        self.sim_outputs_dir = sim_outputs_dir
        self.config_path = config_path
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
        G = nx.Graph()
        D = nx.DiGraph()
        with open(os.path.join(sim_csv_dir, f"sim{sim_id}_filtered_nodes.csv"), newline="") as f:
            for row in csv.DictReader(f):
                attrs = {"label": row["label"], "type": row["type"], "is_leaf": row["is_leaf"] == "True"}
                G.add_node(row["id"], **attrs)
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
                G.add_edge(row["from"], row["to"], **attrs)
                D.add_edge(row["from"], row["to"], **attrs)
        return G, D

    def read_all_filtered_graphs(self):
        graphs = {}
        for sim_id in self.sim_ids():
            graphs[sim_id] = self.read_filtered_graph(sim_id)
        return graphs

    def history_leaders(self, sim_id, tol=1e-9):
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
        """Return {cycle_index: (max |weight_before_projection|, [edges tied for it])} for each cycle."""
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
            leaders[cycle_index] = (highest, tied)
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

    def short_paths(self, G, u, v, tol=1e-9):
        tree = G.copy()
        for a, b, edge_type in G.edges(data="edge_type"):
            if edge_type == "reticulation":
                tree.remove_edge(a, b)
        tree_length = nx.shortest_path_length(tree, u, v, weight="length")

        paths = []
        for path in nx.shortest_simple_paths(G, u, v, weight="length"):
            length = nx.path_weight(G, path, weight="length")
            if length >= tree_length - tol:
                break
            reticulations = []
            for a, b in nx.utils.pairwise(path):
                if G.edges[a, b]["edge_type"] == "reticulation":
                    reticulations.append((a, b))
            paths.append({"path": path, "length": length, "tree_length": tree_length, "reticulations": reticulations})
        return paths

    #need check
    def ancestor_reticulations(self, D, u, v):
        """Return reticulation edges joining u's lineage to v's lineage below their most recent common ancestor."""
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






def main():
    """Analyse simulation outputs and write reports."""
    analysis = MainResultAnalysis()


if __name__ == "__main__":
    main()
