import os
import warnings
from network_lab_tda.tda_analysis import harmonic_cycle, harmonic_cycle_snapshot


class CycleFinder:
    def __init__(self, G, threshold_mode, output_dir, dist_matrix, index_to_name, populated_header_fn="populated_headers.txt", sim_label="", thresholds=None, rips_threshold=float('inf'), harmonic_histories=False):
        self.G = G.to_undirected() if G.is_directed() else G
        self.harmonic_histories = harmonic_histories
        self.dist_matrix = dist_matrix
        self.index_to_name = index_to_name
        self.index_to_node = dict(enumerate(self.G.nodes()))
        self.populated_header_fn = populated_header_fn
        #the limit on the threshold for rips filtration (the last epsilon)
        self.rips_threshold = rips_threshold
        self.sim_label = sim_label
        # list of threshold-selection strategies; only "fixed" affects find_cycles() itself (selects harmonic_cycle_snapshot over harmonic_cycle below). Therest ("cyclelength", "marker") are FilterCycle-only concerns.
        self.threshold_mode = threshold_mode
        self.thresholds = thresholds if thresholds is not None else []

        self.output_path = os.path.join(output_dir, "proc_phylo_outputs", sim_label)
        self.cycle_output_path = os.path.join(output_dir, "cycle_outputs", sim_label)

    def _prepare_dirs(self):
        if not os.path.exists(self.output_path):
            os.makedirs(self.output_path)
        if not os.path.exists(self.cycle_output_path):
            os.makedirs(self.cycle_output_path)

    def find_cycles(self):
        self._prepare_dirs()

        log_path = os.path.join(self.cycle_output_path, "rip.json")

        if self.thresholds and self.threshold_mode == ["fixed"]:
            hc = harmonic_cycle_snapshot(self.dist_matrix, thresholds=self.thresholds, cycle_dim=1, sim_log=True, log_path=log_path)
            hc.run_snapshot(save=False)
        else:
            hc = harmonic_cycle(self.dist_matrix, cycle_dim=1, sim_log=True, log_path=log_path, histories=self.harmonic_histories)
            hc.run_harmonics(threshold=self.rips_threshold, save=True)
        self.cycle_log = hc.log

        return self.cycle_log
