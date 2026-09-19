import os
import copy
from collections import Counter
from itertools import combinations


class TreeFilterCycle:
    """Reports on how qualifying cycles' edges/nodes are shared across a merged
    tree's gene trees, and derives reticulation edges from that sharing. Wraps
    a FilterCycle for the qualified cycles it operates on."""

    WEIGHT_ZERO_TOL = 0.0

    def __init__(self, filter_cycle, sharing_unit="edge", sharing_which_cycles="all", delete_true_edges=False):
        self.filter_cycle = filter_cycle
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

    # -- data pulled from the FilterCycle wrapping these cycles --

    @property
    def G(self):
        return self.filter_cycle.G

    @property
    def index_to_name(self):
        return self.filter_cycle.index_to_name

    @property
    def index_to_node(self):
        return self.filter_cycle.index_to_node

    @property
    def cycle_log(self):
        return self.filter_cycle.cycle_log

    @property
    def cycle_output_path(self):
        return self.filter_cycle.cycle_output_path

    @property
    def qualifying_cycle_keys(self):
        return self.filter_cycle.qualifying_cycle_keys

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

       # step 1: group all edges by point. usually that's just the shorter of an edge's two points, but if both points are the same length (no single "short" one), index the edge under both of them
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

                    f.write(
                        f"{labeled(shared_point)}: {labeled(parent_a)} -- {labeled(parent_b)} "
                        f"| only_a={only_a} only_b={only_b} neither={neither}\n"
                    )

    def _tree_edge_set(self, tree):
        # networkx_to_tree_json already rewrote every node's "label" to the
        # leaf-label tuple for its clade, matching merged_G's node labels
        return {
            frozenset({tuple(tree.nodes[u]["label"]), tuple(tree.nodes[v]["label"])})
            for u, v in tree.edges()
        }

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
