import os
import shutil
import warnings
import copy
from collections import Counter
import networkx as nx
from network_lab_tda.data_prep.Data_Prep import Data_Prep
from network_lab_tda.data_prep.Populate_Edge import Populate_Edge
from network_lab_tda.tda_analysis import harmonic_cycle
from network_lab_tda.tda_visualisation.tda_visual import tda_visual_from_jason


class CycleFinder:
    WEIGHT_ZERO_TOL = 0.0

    def __init__(self, G, threshold_mode, cycle_qualify_mode, output_dir, populated_header_fn="populated_headers.txt", which_nodes="all_nodes", sim_label="", min_cycle_length=0, weight_attr="length", vis=False, use_data_prep=True, thresholds=None, sharing_unit="edge", sharing_which_cycles="all"):
        self.G = G
        self.populated_header_fn = populated_header_fn
        # "all_nodes": keep every non-extinct node
        # "no_hyb_nodes": all_nodes, with internal hybrid-junction nodes collapsed
        self.which_nodes = which_nodes
        self.sim_label = sim_label
        self.min_cycle_length = min_cycle_length
        self.weight_attr = weight_attr
        self.vis = vis
        self.use_data_prep = use_data_prep
        # list of threshold-selection strategies, ANDed together per cycle (getattr'd as _select_threshold_<mode>)
        # options: "cyclelength", "marker", "fixed"
        self.threshold_mode = threshold_mode
        # list of cycle-qualification strategies, ANDed together per cycle (getattr'd as qualifying_cycle_<mode>)
        # options: "marker" (requires "marker" in threshold_mode), "crossover"
        self.cycle_qualify_mode = cycle_qualify_mode
        self.thresholds = thresholds if thresholds is not None else []
        self.qualifying_cycle_keys = []
        self.seen_markers = set()
        # list of units print_most_shared_units reports on, one output file per entry
        # "edge": count each edge's simplex as-is (existing behavior)
        # "node": flatten each edge's simplex to its individual node indices before counting,
        # so edges that share the same underlying nodes are counted as the same unit
        self.sharing_unit = sharing_unit
        # "all": count over every detected cycle; "qualifying": count only over qualifying_cycle_keys
        self.sharing_which_cycles = sharing_which_cycles

        self.output_path = os.path.join(output_dir, "proc_phylo_outputs", sim_label)
        self.cycle_output_path = os.path.join(output_dir, "cycle_outputs", sim_label)
        vis_suffix = "all_nodes" if which_nodes == "all_nodes" else "leaf_nodes"
        self.vis_output_path = os.path.join(self.cycle_output_path, vis_suffix)

    def _prepare_dirs(self):
        if not os.path.exists(self.output_path):
            os.makedirs(self.output_path)
        if os.path.exists(self.vis_output_path):
            shutil.rmtree(self.vis_output_path)
        if not os.path.exists(self.cycle_output_path):
            os.makedirs(self.cycle_output_path)

    def _marker_nodes(self, cycle):
        names = set()
        for edge in cycle["edges"]:
            if abs(edge["weight"]) <= self.WEIGHT_ZERO_TOL:
                continue
            for idx in edge["simplex"]:
                name = self.index_to_name.get(idx)
                if isinstance(name, str) and ("hyb" in name or "sp" in name):
                    names.add(name)
        return names

    def _select_threshold_cyclelength(self,cycle):
            edges = [(edge["simplex"], edge["weight"]) for edge in cycle["edges"] if abs(edge["weight"]) > self.WEIGHT_ZERO_TOL]
            return len(edges) > self.min_cycle_length

    def _select_threshold_marker(self,cycle):
        marker_nodes = self._marker_nodes(cycle)
        new_markers = marker_nodes - self.seen_markers
        self.seen_markers |= marker_nodes
        return len(new_markers) != 0

    def _select_threshold_fixed(self,cycle):
        return any(cycle["birth"] <= t + 1e-5 for t in self.thresholds)

    def generate_threshold_cycle_keys(self):
       for cycle in self.cycle_log["harmonic_cycles"]:
           if self._resolve_threshold_selection(cycle):
               #if fixed: do not want more thresholds appended
               if "fixed" not in self.threshold_mode:
                   self.thresholds.append(cycle["birth"])
               if self._resolve_cycle_qualify(cycle):
                   self.qualifying_cycle_keys.append(cycle)

    def _resolve_threshold_selection(self,cycle):
        for mode in self.threshold_mode:
            if not getattr(self, f"_select_threshold_{mode}")(cycle):
                return False
        return True

    def qualifying_cycle_marker(self,cycle):
        if "marker" not in self.threshold_mode:
            raise ValueError("marker cycle_qualify used but marker is not in threshold_mode, so it was not properly run")
        return True

    def qualifying_cycle_crossover(self,cycle):
        nodes = {idx for edge in cycle["edges"] for idx in edge["simplex"]}
        node_sources = [self.G.nodes[self.index_to_node[idx]].get("sources", set()) for idx in nodes]
        for i in range(len(node_sources)):
            for j in range(i + 1, len(node_sources)):
                if not (node_sources[i] & node_sources[j]):
                    return False
        return True

    def _resolve_cycle_qualify(self,cycle):
        for mode in self.cycle_qualify_mode:
            if not getattr(self, f"qualifying_cycle_{mode}")(cycle):
                return False
        return True

    def _visualize(self):
        self.generate_threshold_cycle_keys()
        os.makedirs(self.vis_output_path)
        plotter = tda_visual_from_jason(
            plt_sep = True,
            data=self.cycle_log,
            thresholds=self.thresholds,
            index_to_name=self.index_to_name,
            log_path=self.vis_output_path,
            cycle_qualify=lambda cycle: cycle in self.qualifying_cycle_keys,
        )
        plotter.cycle_plot()

    def find_cycles(self):
        self._prepare_dirs()

        if self.use_data_prep:
            dp = Data_Prep(G=self.G, log_path=self.output_path, headers=False, weight_attr=self.weight_attr)
            pe = Populate_Edge(G=dp.G, log_path=self.output_path, headers=False, populated_header_fn=self.populated_header_fn, max_node_per_edge=1, weight_attr=self.weight_attr)
            dist_matrix = pe.populate_edges()
            self.index_to_name = pe.index_to_name
        else:
            dist_matrix = nx.floyd_warshall_numpy(self.G)
            self.index_to_name = dict(self.G.nodes(data="label"))
        self.index_to_node = dict(enumerate(self.G.nodes()))

        hc = harmonic_cycle(dist_matrix, cycle_dim=1, sim_log=True, log_path=os.path.join(self.cycle_output_path, "rip.json"))
        hc.run_harmonics(save=False)
        self.cycle_log = hc.log

        if self.vis:
            self._visualize()

        return self.cycle_log

    def _cycles_for(self):
        if self.sharing_which_cycles == "all":
            return self.cycle_log["harmonic_cycles"]
        elif self.sharing_which_cycles == "qualifying":
            return self.qualifying_cycle_keys
        raise ValueError(f"Unknown sharing_which_cycles option: {self.sharing_which_cycles}")

    def _named_point(self, idx):
        name = self.index_to_name[idx]
        if isinstance(name, (list, tuple, set, frozenset)):
            return tuple(name)
        return (name,)

    def sharing_edge_frequency(self):
        counts = Counter()
        for cycle in self._cycles_for():
            edge_keys = set()
            for edge in cycle["edges"]:
                if abs(edge["weight"]) <= self.WEIGHT_ZERO_TOL:
                    continue
                named_points = [self._named_point(idx) for idx in edge["simplex"]]
                edge_keys.add(frozenset(named_points))
            counts.update(edge_keys)
        return counts

    def _flatten_simplex(self, simplex):
        node_names = set()
        for elem in simplex:
            name = self.index_to_name[elem]
            if isinstance(name, (list, tuple, set, frozenset)):
                node_names.update(name)
            else:
                node_names.add(name)
        return frozenset(node_names)

    def sharing_nodes_frequency(self):
        counts = Counter()
        for cycle in self._cycles_for():
            node_set_keys = set()
            for edge in cycle["edges"]:
                if abs(edge["weight"]) <= self.WEIGHT_ZERO_TOL:
                    continue
                node_set_keys.add(self._flatten_simplex(edge["simplex"]))
            counts.update(node_set_keys)

        cumulative_counts = copy.deepcopy(counts)
        for node_set in counts:
            for other_node_set in counts:
                if node_set < other_node_set:
                    cumulative_counts[node_set] += counts[other_node_set]
        return cumulative_counts

    def _frequency_for_unit(self, unit):
        if unit == "edge":
            return self.sharing_edge_frequency()
        elif unit == "node":
            return self.sharing_nodes_frequency()
        raise ValueError(f"Unknown sharing_unit option: {unit}")

    def print_most_shared_units(self, top_n=None):
        for unit in self.sharing_unit:
            counts = self._frequency_for_unit(unit)
            log_fn = f"shared_{unit}s_{self.sharing_which_cycles}.txt"
            log_path = os.path.join(self.cycle_output_path, log_fn)
            with open(log_path, "w") as f:
                for key, count in counts.most_common(top_n):
                    names = tuple(sorted(key))
                    f.write(f"{names}: {count}\n")
