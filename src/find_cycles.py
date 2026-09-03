import os
import json
import shutil
import warnings
import copy
from collections import Counter
from itertools import combinations
import numpy as np
import networkx as nx
from network_lab_tda.data_prep.Data_Prep import Data_Prep
from network_lab_tda.data_prep.Populate_Edge import Populate_Edge
from network_lab_tda.tda_analysis import harmonic_cycle, harmonic_cycle_snapshot
from network_lab_tda.tda_visualisation.tda_visual import tda_visual_from_jason


def filter_edges(G: nx.DiGraph) -> nx.DiGraph:
    """Remove all reticulation edges from G."""
    G = G.copy()
    for u, v, edge_type in list(G.edges(data="edge_type")):
        if edge_type == "reticulation":
            G.remove_edge(u, v)
    return G


class CycleFinder:
    WEIGHT_ZERO_TOL = 0.0

    def __init__(self, G, threshold_mode, cycle_qualify_mode, output_dir, populated_header_fn="populated_headers.txt", which_nodes="all_nodes", sim_label="", min_cycle_length=0, weight_attr="length", vis=True, use_data_prep=True, thresholds=None, sharing_unit="edge", sharing_which_cycles="all", delete_true_edges=False, max_plot_cycles=None, rips_threshold=float('inf'), dress_distance_matrix=None, should_populate_fxn=None):
        self.G = G.to_undirected() if G.is_directed() else G
        self.populated_header_fn = populated_header_fn
        # cap on the Rips filtration passed to harmonic_cycle.run_harmonics() in the
        # non-"fixed" threshold_mode branch; float('inf') (default) reproduces the old
        # uncapped behavior (builds the complete simplex on every point). Callers can
        # pass e.g. the longest single edge length in G to prune the complex.
        self.rips_threshold = rips_threshold
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
        # running log accumulated throughout find_cycles(); flushed as JSON to
        # self.output_path/find_cycles_log.json at the end of find_cycles()
        self.log = {}
        # list of units print_most_shared_units reports on, one output file per entry
        # "edge": count each edge's simplex as-is (existing behavior)
        # "node": flatten each edge's simplex to its individual node indices before counting,
        # so edges that share the same underlying nodes are counted as the same unit
        self.sharing_unit = sharing_unit
        # "all": count over every detected cycle; "qualifying": count only over qualifying_cycle_keys
        self.sharing_which_cycles = sharing_which_cycles
        # if True, sharing_edge_frequency/sharing_nodes_frequency skip edges that G marks as "true_edge"
        # (edges present in every merged input tree, per tree_addition.add_G_edge)
        self.delete_true_edges = delete_true_edges
        # cap on how many qualifying cycles get rendered in _visualize; None means no cap
        self.max_plot_cycles = max_plot_cycles

        self.output_path = os.path.join(output_dir, "proc_phylo_outputs", sim_label)
        self.cycle_output_path = os.path.join(output_dir, "cycle_outputs", sim_label)
        vis_suffix = "all_nodes" if which_nodes == "all_nodes" else "leaf_nodes"
        self.vis_output_path = os.path.join(self.cycle_output_path, vis_suffix)
        # how the distance matrix is "dressed" (post-processed) after populate_edge
        # options: None (no-op), "true_distance_between_tips" (overwrite tip<->tip
        # cells with the true "length"-weighted shortest path, regardless of weight_attr)
        self.dress_distance_matrix = dress_distance_matrix
        # optional fn(u, v, attrs) -> bool; if given, only edges it approves get
        # phantom nodes added in Populate_Edge, others are left untouched
        self.should_populate_fxn = should_populate_fxn

    def _log(self, key, value):
        """Record key -> value in self.log (JSON-serialised on write) and echo it."""
        print(f"{key}: {value}")
        self.log[key] = value

    def _write_log(self):
        log_file = os.path.join(self.output_path, "find_cycles_log.json")
        with open(log_file, "w") as f:
            json.dump(self.log, f, indent=2, default=str)

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
        plotted_keys = self.qualifying_cycle_keys[:self.max_plot_cycles]
        plotter = tda_visual_from_jason(
            plt_sep = True,
            plt_together = False,
            data=self.cycle_log,
            thresholds=self.thresholds,
            index_to_name=self.index_to_name,
            log_path=self.vis_output_path,
            cycle_qualify=lambda cycle: cycle in plotted_keys,
        )
        plotter.cycle_plot()

    def dress_distance_matrix_with_choice(self, dist_matrix, original_G):
        if self.dress_distance_matrix == "true_distance_between_tips":
            node_to_index = {node: idx for idx, node in self.index_to_node.items()}
            leaf_nodes = [n for n, attrs in original_G.nodes(data=True) if attrs.get('is_leaf')]
            dress_log = self.log.setdefault("dress_distance_matrix", [])
            for i, j in combinations(leaf_nodes, 2):
                tip_dist = nx.shortest_path_length(original_G, i, j, weight='length')
                # length=0.0 edges make shortest paths non-unique; log every
                # reticulation on ANY co-shortest path, once per tip pair
                reticulation_edges_on_path = []
                seen = set()
                for path in nx.all_shortest_paths(original_G, i, j, weight='length'):
                    for u, v in nx.utils.pairwise(path):
                        attrs = original_G.edges[u, v]
                        if attrs.get("edge_type") == "reticulation" and frozenset((u, v)) not in seen:
                            seen.add(frozenset((u, v)))
                            reticulation_edges_on_path.append({
                                "edge": [original_G.nodes[u].get("label", str(u)), original_G.nodes[v].get("label", str(v))],
                                "inher_weight": attrs.get("inher_weight"),
                            })
                idx_i, idx_j = node_to_index[i], node_to_index[j]
                old_dist = float(dist_matrix[idx_i, idx_j])
                new_dist = float(tip_dist)
                dress_log.append({
                    "i": original_G.nodes[i].get("label", str(i)),
                    "j": original_G.nodes[j].get("label", str(j)),
                    "idx_i": int(idx_i),
                    "idx_j": int(idx_j),
                    "old_dist": old_dist,
                    "new_dist": new_dist,
                    "changed": old_dist != new_dist,
                    "ratio": (new_dist - old_dist) / old_dist if old_dist else None,
                    "reticulation_edges_on_path": reticulation_edges_on_path,
                })
                dist_matrix[idx_i, idx_j] = tip_dist
                dist_matrix[idx_j, idx_i] = tip_dist
        return dist_matrix

    def find_cycles(self):
        self._prepare_dirs()
        # snapshot before use_data_prep's Populate_Edge mutates self.G in place
        original_G = self.G.copy()
        
        print("preparing prep..G is getting filtered: reticulation edges deleted")
        self.G = filter_edges(self.G)
 
        if self.use_data_prep:
            dp = Data_Prep(G=self.G, log_path=self.output_path, headers=False, weight_attr=self.weight_attr)
            pe = Populate_Edge(G=dp.G, log_path=self.output_path, headers=False, populated_header_fn=self.populated_header_fn, max_node_per_edge=0, weight_attr=self.weight_attr, should_populate_fxn=self.should_populate_fxn)
            dist_matrix = pe.populate_edges()
            self.index_to_name = pe.index_to_name
        else:
            # weighted all-pairs distance keyed positionally (0..n-1), matching the
            # order gudhi/matilda index the Rips simplices by. weight=self.weight_attr
            # so genetic ("length") distances are used, not hop counts; graphs whose
            # edges lack that attr (e.g. tree_main's merged_G) fall back to weight 1.
            dist_matrix = nx.floyd_warshall_numpy(self.G, weight=self.weight_attr)
            self.index_to_name = {
                i: attrs.get("label", node)
                for i, (node, attrs) in enumerate(self.G.nodes(data=True))
            }
            # main_result_analysis.read_index_to_name() reads this back to name simplices
            with open(os.path.join(self.output_path, self.populated_header_fn), "w") as f:
                f.write("\n".join(str(v) for v in self.index_to_name.values()))
            np.savetxt(os.path.join(self.output_path, "populated_distance_matrix.txt"), dist_matrix)
        self.index_to_node = dict(enumerate(self.G.nodes()))
        dist_matrix = self.dress_distance_matrix_with_choice(dist_matrix, original_G)

        log_path = os.path.join(self.cycle_output_path, "rip.json")

        if self.thresholds and self.threshold_mode == ["fixed"]:
            hc = harmonic_cycle_snapshot(dist_matrix, thresholds=self.thresholds, cycle_dim=1, sim_log=True, log_path=log_path)
            hc.run_snapshot(save=False)
        else:
            hc = harmonic_cycle(dist_matrix, cycle_dim=1, sim_log=True, log_path=log_path)
            hc.run_harmonics(threshold=self.rips_threshold, save=True)
        self.cycle_log = hc.log

        if self.vis:
            self._visualize()

        self._write_log()
        return self.cycle_log

    def _cycles_for(self):
        if self.sharing_which_cycles == "all":
            return self.cycle_log["harmonic_cycles"]
        elif self.sharing_which_cycles == "qualifying":
            return self.qualifying_cycle_keys
        raise ValueError(f"Unknown sharing_which_cycles option: {self.sharing_which_cycles}")

    def _is_true_G_edge(self, edge):
        simplex = edge["simplex"]
        if len(simplex) != 2:
            return False
        n1 = self.index_to_node[simplex[0]]
        n2 = self.index_to_node[simplex[1]]
        if self.G.has_edge(n1, n2):
            return self.G.edges[n1, n2].get("label") == "true_edge"
        if self.G.has_edge(n2, n1):
            return self.G.edges[n2, n1].get("label") == "true_edge"
        return False

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
                if self.delete_true_edges and self._is_true_G_edge(edge):
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
                if self.delete_true_edges and self._is_true_G_edge(edge):
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

    def _tree_edge_set(self, tree):
        # networkx_to_tree_json already rewrote every node's "label" to the
        # leaf-label tuple for its clade, matching merged_G's node labels
        return {
            frozenset({tuple(tree.nodes[u]["label"]), tuple(tree.nodes[v]["label"])})
            for u, v in tree.edges()
        }

    def tree_by_tree_delete(self, enumerate_trees, number, leaf_labels=None):
        top_trees = sorted(enumerate_trees.items(), key=lambda kv: kv[1], reverse=True)[:number]

        edge_sets = [self._tree_edge_set(tree) for tree, _count in top_trees]
        delete_edges = set.intersection(*edge_sets) if edge_sets else set()

        counts = Counter()
        for cycle in self._cycles_for():
            edge_keys = set()
            for edge in cycle["edges"]:
                if abs(edge["weight"]) <= self.WEIGHT_ZERO_TOL:
                    continue
                key = frozenset(self._named_point(idx) for idx in edge["simplex"])
                if key in delete_edges:
                    continue
                edge_keys.add(key)
            counts.update(edge_keys)

        # ranked by inclusion
        cumulative_counts = copy.deepcopy(counts)
        to_delete = set()
        for edge in counts:
            a1, a2 = tuple(edge)
            for other_edge in counts:
                if edge == other_edge:
                    continue
                b1, b2 = tuple(other_edge)
                included = (set(a1) <= set(b1) and set(a2) <= set(b2)) or (set(a1) <= set(b2) and set(a2) <= set(b1))
                if included:
                    cumulative_counts[edge] += counts[other_edge]
                    to_delete.add(other_edge)
        counts = {e: c for e, c in cumulative_counts.items() if e not in to_delete}

        # complementary split -- disabled for now
        # seen = set()
        # to_delete = set()
        # cumulative_counts = copy.deepcopy(counts)
        # for edge in counts:
        #     if edge in seen:
        #         continue
        #     p1, p2 = tuple(edge)
        #     if len(p1) == len(p2):
        #         continue
        #     short, long_ = (p1, p2) if len(p1) < len(p2) else (p2, p1)
        #     complement = tuple(sorted(set(long_) - set(short)))
        #     candidate = frozenset({complement, long_})
        #     if candidate in counts and candidate not in seen:
        #         cumulative_counts[edge] += counts[candidate]
        #         to_delete.add(candidate)
        #         seen.add(edge)
        #         seen.add(candidate)
        # counts = {e: c for e, c in cumulative_counts.items() if e not in to_delete}

        # keep only edges that share a point with at least one other edge --
        # a point with no competing resolution is not a reticulation signal

        # step 1: group all edges by point. usually that's just the shorter of
        # an edge's two points, but if both points are the same length (no
        # single "short" one), index the edge under both of them
        edges_by_point = {}
        for edge in counts:
            point_a, point_b = tuple(edge)
            if len(point_a) < len(point_b):
                grouping_points = [point_a]
            elif len(point_b) < len(point_a):
                grouping_points = [point_b]
            else:
                grouping_points = [point_a, point_b]
            for point in grouping_points:
                if point not in edges_by_point:
                    edges_by_point[point] = []
                edges_by_point[point].append(edge)

        # step 2: an edge survives only if at least one of the points it was
        # grouped under has more than one edge in it (i.e. some other edge
        # shares that point)
        edges_to_keep = set()
        for point, edges in edges_by_point.items():
            if len(edges) > 1:
                for edge in edges:
                    edges_to_keep.add(edge)

        # step 3: filter counts down to just the survivors
        counts = {edge: count for edge, count in counts.items() if edge in edges_to_keep}

        leaf_labels = leaf_labels or {}

        def labeled(point):
            return tuple(leaf_labels.get(n, n) for n in sorted(point))

        log_fn = f"shared_edges_{self.sharing_which_cycles}_top{number}_deleted.txt"
        log_path = os.path.join(self.cycle_output_path, log_fn)
        with open(log_path, "w") as f:
            for key, count in sorted(counts.items(), key=lambda kv: kv[1], reverse=True):
                names = tuple(labeled(point) for point in sorted(key))
                f.write(f"{names}: {count}\n")

        # reticulate edges: for every point shared by 2+ edges, the "other"
        # point of each of those edges is a candidate parent -- report every
        # pair of parents (there can be more than 2 if the point is shared by
        # more than 2 edges) as a detected reticulation, annotated with
        # whether the two parent-edges are mutually exclusive across the
        # distinct resolved trees. A genuine reticulation choice never has
        # both edges (or neither) in the same tree -- exactly one always
        # wins. An ordinary structural coincidence (e.g. an unrelated sibling
        # relationship) can have both at once.
        tree_edge_sets = [self._tree_edge_set(tree) for tree in enumerate_trees]

        reticulate_fn = f"reticulate_edges_{self.sharing_which_cycles}_top{number}.txt"
        reticulate_path = os.path.join(self.cycle_output_path, reticulate_fn)
        with open(reticulate_path, "w") as f:
            for shared_point, edges in edges_by_point.items():
                if len(edges) <= 1:
                    continue
                parents = []
                for edge in edges:
                    point_a, point_b = tuple(edge)
                    parent = point_b if point_a == shared_point else point_a
                    parents.append(parent)
                for parent_a, parent_b in combinations(parents, 2):
                    edge_a = frozenset({shared_point, parent_a})
                    edge_b = frozenset({shared_point, parent_b})

                    only_a = 0
                    only_b = 0
                    neither = 0
                    for tree_edges in tree_edge_sets:
                        has_a = edge_a in tree_edges
                        has_b = edge_b in tree_edges
                        if has_a:
                            only_a += 1
                        elif has_b:
                            only_b += 1
                        else:
                            neither += 1
                    # a candidate pair only survives if the two sides split
                    # the trees evenly (only_a == only_b)
                    if only_a != only_b:
                        continue

                    f.write(
                        f"{labeled(shared_point)}: {labeled(parent_a)} -- {labeled(parent_b)} "
                        f"| only_a={only_a} only_b={only_b} neither={neither}\n"
                    )

    def print_most_shared_units(self, top_n=None):
        for unit in self.sharing_unit:
            counts = self._frequency_for_unit(unit)
            true_edges_suffix = "_no_true_edges" if self.delete_true_edges else ""
            log_fn = f"shared_{unit}s_{self.sharing_which_cycles}{true_edges_suffix}.txt"
            log_path = os.path.join(self.cycle_output_path, log_fn)
            with open(log_path, "w") as f:
                for key, count in counts.most_common(top_n):
                    names = tuple(sorted(key))
                    f.write(f"{names}: {count}\n")
