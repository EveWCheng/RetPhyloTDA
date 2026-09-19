import ast
import re
import networkx as nx


def max_tip_spread(G, labels):
    """
    Remove reticulation edges from G, then return the largest unweighted
    shortest-path distance (hop count) between any pair of nodes whose
    label's numeric part is in `labels`.
    """
    tree_only = G.copy()
    tree_only.remove_edges_from(
        [(u, v) for u, v, et in tree_only.edges(data="edge_type") if et == "reticulation"]
    )
    undirected = tree_only.to_undirected()

    label_to_node = {}
    for n, attrs in G.nodes(data=True):
        label = attrs.get("label")
        if label is None:
            continue
        match = re.search(r"\d+", str(label))
        if match:
            label_to_node[int(match.group())] = n

    nodes = [label_to_node[l] for l in labels if l in label_to_node]

    max_dist = 0
    for i, a in enumerate(nodes):
        for b in nodes[i + 1:]:
            d = nx.shortest_path_length(undirected, a, b)
            max_dist = max(max_dist, d)
    return max_dist


def filter_shared_nodes_by_spread(G, input_path, output_path, max_distance):
    """
    Read a shared_nodes-style file (lines like "(31, 35): 26") and copy only
    the lines whose max_tip_spread over the tuple's labels is smaller than
    max_distance.
    """
    with open(input_path) as f, open(output_path, "w") as out:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            labels_str, _, _count_part = line.rpartition(":")
            labels = ast.literal_eval(labels_str.strip())
            if max_tip_spread(G, labels) > max_distance:
                out.write(line + "\n")
