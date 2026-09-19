import os
import shutil
from network_lab_tda.tda_visualisation.tda_visual import tda_visual_from_jason


class FilterCycle:
    """Selects, qualifies and visualizes the cycles a CycleFinder found."""

    WEIGHT_ZERO_TOL = 0.0

    def __init__(self, cycle_finder, cycle_qualify_mode, which_nodes="all_nodes", min_cycle_length=0, vis=True, max_plot_cycles=None):
        self.cycle_finder = cycle_finder
        # "all_nodes": keep every non-extinct node
        # "no_hyb_nodes": all_nodes, with internal hybrid-junction nodes collapsed
        self.which_nodes = which_nodes
        self.min_cycle_length = min_cycle_length
        self.vis = vis
        # list of cycle-qualification strategies, ANDed together per cycle (getattr'd as qualifying_cycle_<mode>)
        self.cycle_qualify_mode = cycle_qualify_mode
        self.qualifying_cycle_keys = []
        self.seen_markers = set()
        # cap on how many qualifying cycles get rendered in visualize; None means no cap
        self.max_plot_cycles = max_plot_cycles

        vis_suffix = "all_nodes" if which_nodes == "all_nodes" else "leaf_nodes"
        self.vis_output_path = os.path.join(self.cycle_finder.cycle_output_path, vis_suffix)

    # -- data pulled from the CycleFinder that produced these cycles --

    @property
    def G(self):
        return self.cycle_finder.G

    @property
    def index_to_name(self):
        return self.cycle_finder.index_to_name

    @index_to_name.setter
    def index_to_name(self, value):
        self.cycle_finder.index_to_name = value

    @property
    def index_to_node(self):
        return self.cycle_finder.index_to_node

    @index_to_node.setter
    def index_to_node(self, value):
        self.cycle_finder.index_to_node = value

    @property
    def cycle_log(self):
        return self.cycle_finder.cycle_log

    @cycle_log.setter
    def cycle_log(self, value):
        self.cycle_finder.cycle_log = value

    @property
    def cycle_output_path(self):
        return self.cycle_finder.cycle_output_path

    @property
    def threshold_mode(self):
        return self.cycle_finder.threshold_mode

    @property
    def thresholds(self):
        return self.cycle_finder.thresholds

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

    def _resolve_cycle_qualify(self,cycle):
        for mode in self.cycle_qualify_mode:
            if not getattr(self, f"qualifying_cycle_{mode}")(cycle):
                return False
        return True

    def visualize(self):
        self.generate_threshold_cycle_keys()
        if os.path.exists(self.vis_output_path):
            shutil.rmtree(self.vis_output_path)
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
