import os
import json
import numpy as np
from network_lab_tda.data_prep.Data_Prep import Data_Prep
from network_lab_tda.data_prep.Populate_Edge import Populate_Edge
from network_lab_tda.tda_analysis import harmonic_cycle
from network_lab_tda.tda_visualisation.tda_visual import tda_visual_from_jason


HERE = os.path.dirname(os.path.abspath(__file__))


def main(headers=True, header_fn="header.txt", populated_header_fn="populated_headers.txt", which_nodes="all_nodes"):
    if which_nodes == "all_nodes":
        input_txt = "distance_matrix_all_nodes.txt"
    elif which_nodes == "no_hyb":
        input_txt = "distance_matrix_no_hyb.txt"
    else:
        input_txt = "distance_matrix.txt"
    input_path = os.path.join(HERE, os.pardir, "Outputs", "phylo_outputs", input_txt)
    output_path = os.path.join(HERE, os.pardir, "Outputs","proc_phylo_outputs")
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    dp = Data_Prep(filepath=input_path, log_path=output_path, headers=headers, header_fn=header_fn)
    pe = Populate_Edge(G=dp.G, log_path=output_path, headers=headers, header_fn=header_fn, populated_header_fn=populated_header_fn)
    print(f"\nEpsilon (quarter): {pe.epsilon:.3f}")
    print(f"Original node count:     {pe.original_node_count}")
#
    dist_matrix = pe.populate_edges()
#
    print(f"Populated node count:    {pe.max_index}")
    print(f"Populated distance matrix shape: {dist_matrix.shape}")

    cycle_output_path = os.path.join(HERE, os.pardir, "Outputs", "cycle_outputs")
    if not os.path.exists(cycle_output_path):
        os.makedirs(cycle_output_path)
    hc = harmonic_cycle(dist_matrix, cycle_dim=1, sim_log=True, log_path=os.path.join(cycle_output_path,"rip.json"))
    simplices, appears_at = hc.rips_filtration()
    hc.compute_harmonics(simplices, appears_at)
    hc.run_harmonics()
    hc.save_log()

    plotter = tda_visual_from_jason(
        jason_path=hc.log_path,
        index_to_name=pe.index_to_name,
        log_path=cycle_output_path,
    )
    plotter.cycle_plot()



    log_path = os.path.join(cycle_output_path, "rip.json")
    with open(log_path, "r") as f:
        cycle_log = json.load(f)

    print("\nAll edges per harmonic cycle:")
    for cycle in cycle_log["harmonic_cycles"]:
        edges = cycle["edges"]
        print(
            f"  cycle {cycle['cycle_index']:>3}  "
            f"birth={cycle['birth']:.4f}  death={cycle['death']:.4f}"
        )
        for edge in edges:
            if abs(edge["weight"]) <= 0.01:
                continue
            u, v = edge["simplex"]
            label_u = pe.index_to_name[u]
            label_v = pe.index_to_name[v]
            print(f"    edge=({label_u}, {label_v})  weight={edge['weight']:.6g}")


main(which_nodes = "no_hyb")
