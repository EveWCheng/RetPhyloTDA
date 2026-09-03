import csv

import networkx as nx

from sim_bdh import PhyloNetwork
from export import export_csv


def _phy():
    """root -> sp2 (tip), root -> h (former hybrid tip that has since speciated)."""
    G = nx.DiGraph()
    G.add_node(1, label="root", is_leaf=False)
    G.add_node(2, label="sp2", is_leaf=True)
    G.add_node(3, label="hyb3", is_leaf=False, is_hyb_leaf=True)  # sticky flag, no longer a tip
    G.add_node(4, label="sp4", is_leaf=True)
    G.add_edge(1, 2, edge_type="tree", length=1.0, time_length=1.0)
    G.add_edge(1, 3, edge_type="tree", length=1.0, time_length=1.0)
    G.add_edge(3, 4, edge_type="tree", length=1.0, time_length=1.0)
    return PhyloNetwork(G=G, nleaves=2)


def test_export_csv_is_leaf_column_reflects_live_flag_not_type(tmp_path):
    export_csv(_phy(), str(tmp_path), prefix="sim0_")
    rows = {r["label"]: r for r in csv.DictReader(open(tmp_path / "sim0_nodes.csv"))}

    # the speciated former hybrid: display 'type' still hyb_leaf, but is_leaf is False
    assert rows["hyb3"]["type"] == "hyb_leaf"
    assert rows["hyb3"]["is_leaf"] == "False"

    assert rows["sp2"]["type"] == "leaf" and rows["sp2"]["is_leaf"] == "True"
    assert rows["root"]["type"] == "internal" and rows["root"]["is_leaf"] == "False"
