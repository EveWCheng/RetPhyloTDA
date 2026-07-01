"""
Birth-Death-Hybridization network simulator.
Python port of sim.bdh.age.help.2_update.R
"""

from __future__ import annotations

import numpy as np
import networkx as nx
from dataclasses import dataclass
from typing import Optional, Callable


@dataclass
class PhyloNetwork:
    """Simulated reticulate phylogeny.

    Everything lives in a single DiGraph G.

    Node attributes:
        is_leaf       -- bool
        timecreation  -- float
        extinct       -- bool: tip went extinct (not a hybrid source)
        hyb_source    -- bool: tip is a hybrid donor
        label         -- str: assigned at output time

    Edge attributes:
        edge_type     -- 'tree' or 'reticulation'
        length        -- float: genetic-distance-based branch length
        time_length   -- float: time-based branch length
    """
    G: nx.DiGraph
    inheritance: np.ndarray     # one value per hybridization event
    nleaves: int                # number of extant (non-extinct) leaves
    tip_states: Optional[list] = None
    time_in_n: Optional[np.ndarray] = None

    @property
    def leaves(self) -> list[int]:
        return [n for n in self.G if self.G.nodes[n]['is_leaf']]

    @property
    def extinct(self) -> list[int]:
        return [n for n in self.G if self.G.nodes[n].get('extinct', False)]

    @property
    def hyb_tips(self) -> list[int]:
        return [n for n in self.G if self.G.nodes[n].get('hyb_source', False)]

    @property
    def Nnode(self) -> int:
        return sum(1 for n in self.G if not self.G.nodes[n]['is_leaf'])


class SimState:
    """Mutable state for one BDH simulation run."""

    def __init__(self, mrca: bool):
        self.G = nx.DiGraph()   # all edges; edge_type='tree' or 'reticulation'

        self.time = 0.0
        self.timestep = 0.0
        self.tree_extinct = False

        self.leaves: set[int] = set()   # IDs of currently active leaves
        self.inheritance: list[float] = []
        self.time_in_n: list[float] = []
        self.trait_state: list = []
        self.ext_trait_state: list = []

        self._id = 0   # node ID counter; increments with every new node

        root  = self._new_node(is_leaf=False)
        leaf1 = self._new_node(is_leaf=True)
        self.G.add_edge(root, leaf1, edge_type='tree', length=0.0, time_length=0.0)
        self.leaves.add(leaf1)

        if mrca:
            leaf2 = self._new_node(is_leaf=True)
            self.G.add_edge(root, leaf2, edge_type='tree', length=0.0, time_length=0.0)
            self.leaves.add(leaf2)

    def _new_node(self, is_leaf: bool) -> int:
        """Allocate a new node ID and register it in G with its attributes."""
        self._id += 1
        self.G.add_node(self._id, is_leaf=is_leaf, timecreation=self.time)
        return self._id

    def _leaf_index(self, node: int) -> int:
        """Position of node in the ordered active-leaf list (for trait_state indexing)."""
        return sorted(self.leaves).index(node)

    def _seal_edge(self, parent: int, child: int):
        """Finalize both length fields on the edge (parent, child)."""
        elapsed = self.time - self.G.nodes[child]['timecreation']
        self.G[parent][child]['length'] = elapsed
        self.G[parent][child]['time_length'] = elapsed

    def _seal_incoming(self, node: int) -> int:
        """Find the parent of node, seal its incoming edge, return the parent."""
        parent = next(self.G.predecessors(node))
        self._seal_edge(parent, node)
        return parent

    def speciation(self, species: int, trait_model, is_null_trait: bool):
        self._seal_incoming(species)

        self.G.nodes[species]['is_leaf'] = False
        self.leaves.discard(species)

        child1 = self._new_node(is_leaf=True)
        child2 = self._new_node(is_leaf=True)
        self.G.add_edge(species, child1, edge_type='tree', length=0.0, time_length=0.0)
        self.G.add_edge(species, child2, edge_type='tree', length=0.0, time_length=0.0)
        self.leaves.update({child1, child2})

        if not is_null_trait:
            idx = self._leaf_index(species)
            new_traits = list(trait_model['spec_fxn']([self.trait_state[idx]]))
            self.trait_state = (self.trait_state[:idx]
                                + new_traits
                                + self.trait_state[idx + 1:])

    def extinction(self, species: int, is_null_trait: bool):
        self._seal_incoming(species)
        self.G.nodes[species]['is_leaf'] = False
        self.G.nodes[species]['extinct'] = True
        self.leaves.discard(species)

        if not is_null_trait:
            idx = self._leaf_index(species)
            self.ext_trait_state.append(self.trait_state.pop(idx))

    def _edge_length(self, u: int, v: int) -> float:
        """Current length of edge (u, v): sealed value, or elapsed time if still open."""
        stored = self.G[u][v]['length']
        if stored == 0.0 and v in self.leaves:
            return self.time - self.G.nodes[v]['timecreation']
        return stored

    def tip_distance(self, tip1: int, tip2: int) -> float:
        """Current genetic distance between two active leaves."""
        path = nx.shortest_path(self.G.to_undirected(), tip1, tip2)
        total = 0.0
        for u, v in zip(path, path[1:]):
            if self.G.has_edge(u, v):
                total += self._edge_length(u, v)
            else:
                total += self._edge_length(v, u)
        return total

    def _hyb_setup(self, sp1: int, sp2: int, inher: float, d12: float | None = None):
        """Shared opening for all hybridization events.

        Returns (primary, secondary, primary_inher, secondary_inher, d12).
        Seals both parent edges and removes them from active leaves.
        """
        d12 = d12 if d12 is not None else self.tip_distance(sp1, sp2)

        primary   = sp1 if (1 - inher) > 0.5 else sp2
        secondary = sp2 if primary == sp1 else sp1
        primary_inher   = (1 - inher) if primary == sp1 else inher
        secondary_inher = 1 - primary_inher

        for sp in (sp1, sp2):
            self._seal_incoming(sp)
            self.G.nodes[sp]['is_leaf'] = False
            self.leaves.discard(sp)

        self.inheritance.append(inher)
        return primary, secondary, primary_inher, secondary_inher, d12

    def hyb_generating(self, sp1: int, sp2: int, inher: float, d12: float | None = None):
        """LG: both parents continue; new hybrid node + leaf added."""
        primary, secondary, primary_inher, secondary_inher, d12 = self._hyb_setup(sp1, sp2, inher, d12)

        leaf_p   = self._new_node(is_leaf=True)   # primary parent continuation
        leaf_s   = self._new_node(is_leaf=True)   # secondary parent continuation
        hyb_node = self._new_node(is_leaf=False)  # internal hybrid node
        hyb_leaf = self._new_node(is_leaf=True)   # new hybrid leaf

        self.G.add_edge(sp1, leaf_p,  edge_type='tree', length=0.0, time_length=0.0)
        self.G.add_edge(sp2, leaf_s,  edge_type='tree', length=0.0, time_length=0.0)
        self.G.add_edge(primary,   hyb_node, edge_type='tree',length=secondary_inher * d12, time_length=0.0)
        self.G.add_edge(hyb_node,  hyb_leaf, edge_type='tree', length=0.0, time_length=0.0)
        self.G.add_edge(secondary, hyb_node, edge_type='reticulation',length=primary_inher * d12, time_length=0.0)

        self.leaves.update({leaf_p, leaf_s, hyb_leaf})

    def hyb_degenerative(self, sp1: int, sp2: int, inher: float, d12: float | None = None):
        """LD: secondary parent absorbed; primary continues as the hybrid leaf."""
        primary, secondary, primary_inher, secondary_inher, d12 = self._hyb_setup(sp1, sp2, inher, d12)

        hyb_leaf = self._new_node(is_leaf=True)

        self.G.add_edge(primary,   hyb_leaf, edge_type='tree', length=0.0, time_length=0.0)
        self.G.add_edge(secondary, primary,  edge_type='reticulation',length=secondary_inher * d12, time_length=0.0)

        self.G.nodes[secondary]['hyb_source'] = True
        self.leaves.add(hyb_leaf)

    def hyb_neutral(self, sp1: int, sp2: int, inher: float, d12: float | None = None):
        """LN: both parents continue; primary's lineage carries the hybrid genome."""
        primary, secondary, primary_inher, secondary_inher, d12 = \
            self._hyb_setup(sp1, sp2, inher, d12)

        hyb_leaf   = self._new_node(is_leaf=True)  # primary continuation (hybrid genome)
        donor_leaf = self._new_node(is_leaf=True)  # secondary continuation (unchanged)

        self.G.add_edge(primary,   hyb_leaf,   edge_type='tree', length=0.0, time_length=0.0)
        self.G.add_edge(secondary, donor_leaf, edge_type='tree', length=0.0, time_length=0.0)
        self.G.add_edge(secondary, primary,    edge_type='reticulation',
                        length=primary_inher * d12, time_length=0.0)

        self.leaves.update({hyb_leaf, donor_leaf})


# ── Gillespie loop ─────────────────────────────────────────────────────────────

def _nchoose2(n: int) -> int:
    return n * (n - 1) // 2


def sim_bdh_age(age: float, numbsim: int,
                lambda_: float, mu: float, nu: float,
                hybprops: list[float], hyb_inher_fxn: Callable,
                mrca: bool = False,
                hyb_rate_fxn: Optional[Callable] = None) -> list[dict]:
    """Run numbsim independent BDH simulations and return a list of result dicts."""
    return [_sim_one(age, lambda_, mu, nu, hybprops, hyb_inher_fxn,mrca, hyb_rate_fxn) for _ in range(numbsim)]


def _sim_one(age, lambda_, mu, nu, hybprops, hyb_inher_fxn, mrca, hyb_rate_fxn) -> dict:
    """Run one BDH simulation; return {phy: PhyloNetwork | 0, distance: dict}."""
    state = SimState(mrca=mrca)

    while True:
        n = len(state.leaves)
        if n == 0:
            return {'phy': 0, 'distance': None}

        spec_rate  = n * lambda_
        ext_rate   = n * mu
        hyb_rate   = _nchoose2(n) * nu
        total_rate = spec_rate + ext_rate + hyb_rate

        state.timestep = np.random.exponential(1.0 / total_rate)
        state.time    += state.timestep

        if state.time >= age:
            state.time = age
            break

        sp1, sp2 = _sample_lineages(state.leaves)
        rand = np.random.uniform()

        if rand < spec_rate / total_rate:
            state.speciation(sp1, trait_model=None, is_null_trait=True)
        elif rand < (spec_rate + ext_rate) / total_rate:
            state.extinction(sp1, is_null_trait=True)
        else:
            _try_hybridization(state, sp1, sp2, hybprops, hyb_inher_fxn, hyb_rate_fxn)

    _finalize_pendant_edges(state)
    return _build_output(state)


def _sample_lineages(leaves: set[int]) -> tuple[int, int]:
    """Sample 2 random active leaves."""
    sp1, sp2 = np.random.choice(list(leaves), size=2, replace=False)
    return int(sp1), int(sp2)


def _try_hybridization(state: SimState, sp1, sp2,
                       hybprops, hyb_inher_fxn, hyb_rate_fxn):
    inher = hyb_inher_fxn()
    d12 = None

    if hyb_rate_fxn is not None:
        d12 = state.tip_distance(sp1, sp2)
        if np.random.uniform() > hyb_rate_fxn(d12):
            return

    total = sum(hybprops)
    rand  = np.random.uniform()

    if rand < hybprops[0] / total:
        state.hyb_generating(sp1, sp2, inher, d12)
    elif rand < sum(hybprops[:2]) / total:
        state.hyb_degenerative(sp1, sp2, inher, d12)
    else:
        state.hyb_neutral(sp1, sp2, inher, d12)


def _finalize_pendant_edges(state: SimState):
    """Seal all edges leading to currently active leaves."""
    for leaf in state.leaves:
        parent = next(state.G.predecessors(leaf))
        state._seal_edge(parent, leaf)


def _build_output(state: SimState) -> dict:
    """Package the final graph into a PhyloNetwork and compute the distance matrix."""
    phy = PhyloNetwork(
        G=state.G,
        inheritance=np.array(state.inheritance),
        nleaves=len(state.leaves),
    )
    distance = dict(nx.all_pairs_dijkstra_path_length(state.G.to_undirected(),
                                                       weight='length'))
    return {'phy': phy, 'distance': distance}
