"""
Recover reticulation-driven edges directly from the set of distinct resolved
gene trees (enumerate_gene_trees' output), without going through cycle
detection at all.
"""
import os
from collections import Counter
from itertools import combinations


def _tree_edge_set(tree):
    return {
        frozenset({tuple(tree.nodes[u]["label"]), tuple(tree.nodes[v]["label"])})
        for u, v in tree.edges()
    }


def _ranked_by_inclusion(counts):
    """An edge picks up every other edge's count if that other edge's two
    endpoints each contain (superset) this edge's matching endpoint -- then
    delete the absorbed (larger/less specific) edge. Returns the surviving,
    boosted counts dict."""
    cumulative_counts = dict(counts)
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
    return {e: c for e, c in cumulative_counts.items() if e not in to_delete}


def polymorphic_edges(enumerate_trees, leaf_labels=None):
    total = len(enumerate_trees)
    edge_tree_counts = Counter()
    for tree in enumerate_trees:
        edge_tree_counts.update(_tree_edge_set(tree))

    counts = {edge: count for edge, count in edge_tree_counts.items() if count != total}

    counts = _ranked_by_inclusion(counts)

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

    edges_to_keep = set()
    for point, edges in edges_by_point.items():
        if len(edges) > 1:
            for edge in edges:
                edges_to_keep.add(edge)
    counts = {edge: count for edge, count in counts.items() if edge in edges_to_keep}

    leaf_labels = leaf_labels or {}
    results = []
    for edge, count in counts.items():
        names = tuple(
            tuple(leaf_labels.get(n, n) for n in point)
            for point in sorted(edge)
        )
        results.append((names, count, total))
    results.sort(key=lambda r: r[1], reverse=True)
    return results


def write_polymorphic_edges(enumerate_trees, out_dir, leaf_labels=None, prefix=""):
    results = polymorphic_edges(enumerate_trees, leaf_labels)
    os.makedirs(out_dir, exist_ok=True)
    log_path = os.path.join(out_dir, f"{prefix}polymorphic_edges.txt")
    with open(log_path, "w") as f:
        for names, count, total in results:
            f.write(f"{names}: {count}\n")
    return log_path


def reticulate_edges(enumerate_trees, leaf_labels=None):
    """For every point shared by 2+ candidate edges, report every pair of
    the edges' "other" points (parents) as a candidate reticulation, keeping
    only pairs where the parents split the distinct trees evenly (only_a ==
    only_b) -- the signature of a genuine binary reticulation choice, as
    opposed to an edge whose presence is entangled with some other choice.
    Same technique as TreeFilterCycle.tree_by_tree_delete's reticulate-edges
    step, but working entirely off enumerate_trees -- no cycle detection."""
    total = len(enumerate_trees)
    tree_edge_sets = [_tree_edge_set(tree) for tree in enumerate_trees]

    edge_tree_counts = Counter()
    for tree_edges in tree_edge_sets:
        edge_tree_counts.update(tree_edges)
    counts = {edge: count for edge, count in edge_tree_counts.items() if count != total}
    counts = _ranked_by_inclusion(counts)
    candidate_edges = list(counts)

    edges_by_point = {}
    for edge in candidate_edges:
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

    leaf_labels = leaf_labels or {}

    def labeled(point):
        return tuple(leaf_labels.get(n, n) for n in sorted(point))

    results = []
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

            results.append((labeled(shared_point), labeled(parent_a), labeled(parent_b), only_a, only_b, neither))
    return results


def write_reticulate_edges(enumerate_trees, out_dir, leaf_labels=None, prefix=""):
    results = reticulate_edges(enumerate_trees, leaf_labels)
    os.makedirs(out_dir, exist_ok=True)
    log_path = os.path.join(out_dir, f"{prefix}reticulate_edges.txt")
    with open(log_path, "w") as f:
        for shared_point, parent_a, parent_b, only_a, only_b, neither in results:
            f.write(f"{shared_point}: {parent_a} -- {parent_b} | only_a={only_a} only_b={only_b} neither={neither}\n")
    return log_path
