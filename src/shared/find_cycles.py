import os
import json
import warnings
from network_lab_tda.tda_analysis import harmonic_cycle, harmonic_cycle_snapshot


class CycleFinder:
    def __init__(self, G, threshold_mode, output_dir, dist_matrix, index_to_name, populated_header_fn="populated_headers.txt", sim_label="", thresholds=None, rips_threshold=float('inf')):
        self.G = G.to_undirected() if G.is_directed() else G
        self.dist_matrix = dist_matrix
        self.index_to_name = index_to_name
        self.index_to_node = dict(enumerate(self.G.nodes()))
        self.populated_header_fn = populated_header_fn
        # cap on the Rips filtration passed to harmonic_cycle.run_harmonics() in the
        # non-"fixed" threshold_mode branch; float('inf') (default) reproduces the old
        # uncapped behavior (builds the complete simplex on every point). Callers can
        # pass e.g. the longest single edge length in G to prune the complex.
        self.rips_threshold = rips_threshold
        self.sim_label = sim_label
        # list of threshold-selection strategies; only "fixed" affects find_cycles()
        # itself (selects harmonic_cycle_snapshot over harmonic_cycle below). The
        # rest ("cyclelength", "marker") are FilterCycle-only concerns.
        self.threshold_mode = threshold_mode
        self.thresholds = thresholds if thresholds is not None else []
        # running log accumulated throughout find_cycles(); flushed as JSON to
        # self.output_path/find_cycles_log.json at the end of find_cycles()
        self.log = {}

        self.output_path = os.path.join(output_dir, "proc_phylo_outputs", sim_label)
        self.cycle_output_path = os.path.join(output_dir, "cycle_outputs", sim_label)

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
        if not os.path.exists(self.cycle_output_path):
            os.makedirs(self.cycle_output_path)

    def find_cycles(self):
        self._prepare_dirs()

        log_path = os.path.join(self.cycle_output_path, "rip.json")

        if self.thresholds and self.threshold_mode == ["fixed"]:
            hc = harmonic_cycle_snapshot(self.dist_matrix, thresholds=self.thresholds, cycle_dim=1, sim_log=True, log_path=log_path)
            hc.run_snapshot(save=False)
        else:
            hc = harmonic_cycle(self.dist_matrix, cycle_dim=1, sim_log=True, log_path=log_path)
            hc.run_harmonics(threshold=self.rips_threshold, save=True)
        self.cycle_log = hc.log

        self._write_log()
        return self.cycle_log
