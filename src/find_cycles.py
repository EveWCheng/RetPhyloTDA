import os
import json
import shutil
import numpy as np
from network_lab_tda.data_prep.Data_Prep import Data_Prep
from network_lab_tda.data_prep.Populate_Edge import Populate_Edge
from network_lab_tda.tda_analysis import harmonic_cycle
from network_lab_tda.tda_visualisation.tda_visual import tda_visual_from_jason


HERE = os.path.dirname(os.path.abspath(__file__))


def _hyb_nodes(cycle, index_to_name):
    names = set()
    for edge in cycle["edges"]:
        if edge["weight"] == 0:
            continue
        for idx in edge["simplex"]:
            name = index_to_name.get(idx)
            if isinstance(name, str) and "hyb" in name:
                names.add(name)
    return names


def _drawn_edges_have_hybrid(drawn_edges, index_to_name):
    for simplex, _weight in drawn_edges:
        for idx in simplex:
            name = index_to_name.get(idx)
            if isinstance(name, str) and "hyb" in name:
                return True
    return False


def _select_thresholds(cycle_log, index_to_name, min_cycle_length):
    thresholds = []
    seen_hyb = set()
    for c in cycle_log["harmonic_cycles"]:
        if sum(1 for edge in c["edges"] if edge["weight"] != 0) <= min_cycle_length:
            continue
        hyb_nodes = _hyb_nodes(c, index_to_name)
        new_hyb = hyb_nodes - seen_hyb
        if len(new_hyb) == 0:
            continue
        thresholds.append(c["birth"])
        seen_hyb |= hyb_nodes
    return thresholds


def find_cycles(G, populated_header_fn="populated_headers.txt", which_nodes="all_nodes", sim_label="", min_cycle_length=0):
    output_path = os.path.join(HERE, os.pardir, "Outputs", "proc_phylo_outputs", sim_label)
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    cycle_output_path = os.path.join(HERE, os.pardir, "Outputs", "cycle_outputs", sim_label)
    vis_suffix = "all_nodes" if which_nodes == "all_nodes" else "leaf_nodes"
    vis_output_path = os.path.join(cycle_output_path, vis_suffix)
    if os.path.exists(vis_output_path):
        shutil.rmtree(vis_output_path)

    dp = Data_Prep(G=G, log_path=output_path, headers=False)
    pe = Populate_Edge(G=dp.G, log_path=output_path, headers=False, populated_header_fn=populated_header_fn,max_node_per_edge=1)
#
    dist_matrix = pe.populate_edges()

    if not os.path.exists(cycle_output_path):
        os.makedirs(cycle_output_path)
    hc = harmonic_cycle(dist_matrix, cycle_dim=1, sim_log=True, log_path=os.path.join(cycle_output_path,"rip.json"))
    simplices, appears_at = hc.rips_filtration()
    hc.compute_harmonics(simplices, appears_at)
    hc.run_harmonics()
    hc.save_log()

    with open(hc.log_path, "r") as f:
        cycle_log = json.load(f)

    thresholds = _select_thresholds(cycle_log, pe.index_to_name, min_cycle_length)
#
    os.makedirs(vis_output_path)
    plotter = tda_visual_from_jason(
        jason_path=hc.log_path,
        thresholds=thresholds,
        index_to_name=pe.index_to_name,
        log_path=vis_output_path,
        cycle_qualify=lambda drawn_edges: (
            len(drawn_edges) > min_cycle_length
            and _drawn_edges_have_hybrid(drawn_edges, pe.index_to_name)
        ),
    )
    plotter.cycle_plot()

    for cycle in cycle_log["harmonic_cycles"]:
        edges = cycle["edges"]
        for edge in edges:
            if abs(edge["weight"]) <= 0.01:
                continue
            u, v = edge["simplex"]
            label_u = pe.index_to_name[u]
            label_v = pe.index_to_name[v]
