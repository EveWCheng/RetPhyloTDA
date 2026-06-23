import os
import numpy as np
from network_lab_tda.data_prep.Data_Prep import Data_Prep
from network_lab_tda.data_prep.Populate_Edge import Populate_Edge
from network_lab_tda.tda_analysis import harmonic_cycle
from network_lab_tda.tda_visualisation.tda_visual import tda_visual_from_jason


HERE = os.path.dirname(os.path.abspath(__file__))


def main(headers=True, header_fn="header.txt", populated_header_fn="populated_headers.txt"):
    input_path = os.path.join(HERE, os.pardir, "Outputs", "phylo_outputs", "distance_matrix.txt")
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


main()
