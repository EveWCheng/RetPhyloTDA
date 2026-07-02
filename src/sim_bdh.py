"""
Birth-Death-Hybridization network simulator.
Python port of sim.bdh.age.help.2_update.R
"""

from __future__ import annotations

import itertools

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
        trait         -- present only if a trait_model was supplied

    Edge attributes:
        edge_type     -- 'tree' or 'reticulation'
        length        -- float: genetic-distance-based branch length
        time_length   -- float: time-based branch length
        inher_weight  -- float, only on edges entering a reticulation node (in-degree 2):
                         probability that a gene's true history takes this edge
    """
    G: nx.DiGraph
    nleaves: int                # number of extant (non-extinct) leaves
    tip_states: Optional[dict] = None   # {species label: trait value}, only if a trait_model was used

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

    def gene_tree(self, gene_index: int) -> nx.DiGraph:
        edges = [(u, v) for u, v, attrs in self.G.edges(data=True) if gene_index in attrs.get('genes', set())]
        return self.G.edge_subgraph(edges)

    def enumerate_gene_trees(self) -> list[tuple[nx.DiGraph, float]]:
        """Every possible resolved gene tree topology, with its exact probability.

        Equivalent to the R port's Ngene=0 mode: at each reticulation node
        (in-degree 2) a gene takes exactly one of the two incoming edges, so
        the full topology space is the Cartesian product of those choices
        across all reticulation nodes, weighted by 'inher_weight'.

        Have not checked **Eve**
        """
        retic_nodes = [n for n in self.G if self.G.in_degree(n) == 2]
        choices = [list(self.G.in_edges(n, data='inher_weight')) for n in retic_nodes]

        results = []
        for combo in itertools.product(*choices):
            chosen = {(u, v) for u, v, _ in combo}
            dropped = {(u, v) for edges in choices for u, v, _ in edges} - chosen
            weight = 1.0
            for _, _, w in combo:
                weight *= w
            kept = [(u, v) for u, v in self.G.edges() if (u, v) not in dropped]
            results.append((self.G.edge_subgraph(kept), weight))
        return results



TRAIT_MODEL_KEYS = {'initial', 'time_fxn', 'spec_fxn', 'hyb_event_fxn', 'hyb_compatibility_fxn'}


class SimState:
    """Mutable state for one BDH simulation run."""

    def __init__(self, mrca: bool, Ngene: int = 0, trait_model: Optional[dict] = None):
        self.G = nx.DiGraph()   # all edges; edge_type='tree' or 'reticulation'
        self.time = 0.0
        self.timestep = 0.0
        self.tree_extinct = False
        self.Ngene = Ngene
        self.trait_model = trait_model

        self.leaves: set[int] = set()   # IDs of currently active leaves
        self._id = 0   # node ID counter; increments with every new node

        if trait_model is not None:
            missing = TRAIT_MODEL_KEYS - trait_model.keys()
            if missing:
                raise ValueError(f"trait_model is missing required keys: {sorted(missing)}")
            initial = trait_model['initial']
            n_start = 2 if mrca else 1
            if len(initial) != n_start:
                raise ValueError(f"trait_model['initial'] must supply {n_start} starting value(s), got {len(initial)}")

        root  = self._new_node(is_leaf=False)
        leaf1 = self._new_node(is_leaf=True)
        self.G.add_edge(root, leaf1, edge_type='tree', length=0.0, time_length=0.0, genes=set(range(Ngene)))
        self.leaves.add(leaf1)
        if trait_model is not None:
            self.G.nodes[leaf1]['trait'] = trait_model['initial'][0]

        if mrca:
            leaf2 = self._new_node(is_leaf=True)
            self.G.add_edge(root, leaf2, edge_type='tree', length=0.0, time_length=0.0, genes=set(range(Ngene)))
            self.leaves.add(leaf2)
            if trait_model is not None:
                self.G.nodes[leaf2]['trait'] = trait_model['initial'][1]

    def _new_node(self, is_leaf: bool) -> int:
        """Allocate a new node ID and register it in G with its attributes."""
        self._id += 1
        self.G.add_node(self._id, is_leaf=is_leaf, timecreation=self.time)
        return self._id

    def _evolve_traits(self, timestep: float):
        """Advance every active leaf's trait by one timestep via trait_model['time_fxn']."""
        time_fxn = self.trait_model['time_fxn']
        for leaf in self.leaves:
            self.G.nodes[leaf]['trait'] = time_fxn(self.G.nodes[leaf]['trait'], timestep)

    def _seal_edge(self, parent: int, child: int):
        """Finalize both length fields on the edge (parent, child)."""
        elapsed = self.time - self.G.nodes[child]['timecreation']
        self.G[parent][child]['length'] = elapsed
        self.G[parent][child]['time_length'] = elapsed

    def _seal_incoming(self, node: int) -> int:
        """Seal node's incoming edge and mark it no longer an active leaf. Returns the parent."""
        parent = next(self.G.predecessors(node))
        self._seal_edge(parent, node)
        self.G.nodes[node]['is_leaf'] = False
        return parent

    def speciation(self, species: int):
        self._seal_incoming(species)
        self.leaves.discard(species)

        child1 = self._new_node(is_leaf=True)
        child2 = self._new_node(is_leaf=True)
        all_genes = set(range(self.Ngene))
        self.G.add_edge(species, child1, edge_type='tree', length=0.0, time_length=0.0, genes=all_genes)
        self.G.add_edge(species, child2, edge_type='tree', length=0.0, time_length=0.0, genes=all_genes)
        self.leaves.update({child1, child2})

        if self.trait_model is not None:
            child1_trait, child2_trait = self.trait_model['spec_fxn'](self.G.nodes[species]['trait'])
            self.G.nodes[child1]['trait'] = child1_trait
            self.G.nodes[child2]['trait'] = child2_trait

    def extinction(self, species: int):
        self._seal_incoming(species)
        self.G.nodes[species]['extinct'] = True
        self.leaves.discard(species)

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

        Returns (primary, secondary, primary_inher, secondary_inher, d12,
        primary_genes, secondary_genes, parent_of_primary).
        Seals both parent edges and removes them from active leaves.
        parent_of_primary is primary's pre-existing ancestor, needed by
        hyb_degenerative/hyb_neutral to retroactively tag that edge's
        inher_weight once it becomes one of two incoming edges at a
        reticulation node.
        """
        d12 = d12 if d12 is not None else self.tip_distance(sp1, sp2)

        primary   = sp1 if (1 - inher) > 0.5 else sp2
        secondary = sp2 if primary == sp1 else sp1
        primary_inher   = (1 - inher) if primary == sp1 else inher
        secondary_inher = 1 - primary_inher

        parent_of_primary = None
        for sp in (sp1, sp2):
            parent = self._seal_incoming(sp)
            self.leaves.discard(sp)
            if sp == primary:
                parent_of_primary = parent

        if self.Ngene > 0:
            k = np.random.binomial(self.Ngene, primary_inher)
            primary_genes = set(np.random.choice(self.Ngene, k, replace=False))
            secondary_genes = set(range(self.Ngene)) - primary_genes
        else:
            primary_genes = set()
            secondary_genes = set()

        return (primary, secondary, primary_inher, secondary_inher, d12,
                primary_genes, secondary_genes, parent_of_primary)

    def hyb_generating(self, sp1: int, sp2: int, inher: float, d12: float | None = None,hyb_trait=None):
        """LG: both parents continue; new hybrid node + leaf added."""
        if self.trait_model is not None:
            sp1_trait, sp2_trait = self.G.nodes[sp1]['trait'], self.G.nodes[sp2]['trait']
        primary, secondary, primary_inher, secondary_inher, d12, primary_genes, secondary_genes, _ = self._hyb_setup(sp1, sp2, inher, d12)

        leaf_p   = self._new_node(is_leaf=True)   # primary parent continuation
        leaf_s   = self._new_node(is_leaf=True)   # secondary parent continuation
        hyb_node = self._new_node(is_leaf=False)  # internal hybrid node
        hyb_leaf = self._new_node(is_leaf=True)   # new hybrid leaf

        all_genes = set(range(self.Ngene))
        self.G.add_edge(sp1, leaf_p,  edge_type='tree', length=0.0, time_length=0.0, genes=all_genes)
        self.G.add_edge(sp2, leaf_s,  edge_type='tree', length=0.0, time_length=0.0, genes=all_genes)
        self.G.add_edge(primary,   hyb_node, edge_type='tree', length=secondary_inher * d12, time_length=0.0, genes=primary_genes, inher_weight=primary_inher)
        self.G.add_edge(hyb_node,  hyb_leaf, edge_type='tree', length=0.0, time_length=0.0, genes=all_genes)
        self.G.add_edge(secondary, hyb_node, edge_type='reticulation', length=primary_inher * d12, time_length=0.0, genes=secondary_genes, inher_weight=secondary_inher)

        self.leaves.update({leaf_p, leaf_s, hyb_leaf})

        if self.trait_model is not None:
            self.G.nodes[leaf_p]['trait'] = sp1_trait
            self.G.nodes[leaf_s]['trait'] = sp2_trait
            self.G.nodes[hyb_leaf]['trait'] = hyb_trait

    def hyb_degenerative(self, sp1: int, sp2: int, inher: float, d12: float | None = None,
                         hyb_trait=None):
        """LD: secondary parent absorbed; primary continues as the hybrid leaf."""
        primary, secondary, primary_inher, secondary_inher, d12, primary_genes, secondary_genes, parent_of_primary = \
            self._hyb_setup(sp1, sp2, inher, d12)

        hyb_leaf = self._new_node(is_leaf=True)

        self.G.add_edge(primary,   hyb_leaf, edge_type='tree', length=0.0, time_length=0.0, genes=primary_genes)
        self.G.add_edge(secondary, primary,  edge_type='reticulation', length=secondary_inher * d12, time_length=0.0, genes=secondary_genes, inher_weight=secondary_inher)
        self.G[parent_of_primary][primary]['inher_weight'] = primary_inher

        self.G.nodes[secondary]['hyb_source'] = True
        self.leaves.add(hyb_leaf)

        if self.trait_model is not None:
            self.G.nodes[hyb_leaf]['trait'] = hyb_trait

    def hyb_neutral(self, sp1: int, sp2: int, inher: float, d12: float | None = None,
                    hyb_trait=None):
        """LN: both parents continue; primary's lineage carries the hybrid genome."""
        primary, secondary, primary_inher, secondary_inher, d12, primary_genes, secondary_genes, parent_of_primary = \
            self._hyb_setup(sp1, sp2, inher, d12)

        if self.trait_model is not None:
            secondary_trait = self.G.nodes[secondary]['trait']

        hyb_leaf   = self._new_node(is_leaf=True)  # primary continuation (hybrid genome)
        donor_leaf = self._new_node(is_leaf=True)  # secondary continuation (unchanged)

        all_genes = set(range(self.Ngene))
        self.G.add_edge(primary,   hyb_leaf,   edge_type='tree', length=0.0, time_length=0.0, genes=primary_genes)
        self.G.add_edge(secondary, donor_leaf, edge_type='tree', length=0.0, time_length=0.0, genes=all_genes)
        self.G.add_edge(secondary, primary, edge_type='reticulation', length=primary_inher * d12, time_length=0.0, genes=secondary_genes, inher_weight=secondary_inher)
        self.G[parent_of_primary][primary]['inher_weight'] = primary_inher

        self.leaves.update({hyb_leaf, donor_leaf})

        if self.trait_model is not None:
            self.G.nodes[hyb_leaf]['trait'] = hyb_trait
            self.G.nodes[donor_leaf]['trait'] = secondary_trait


# ── Gillespie loop ─────────────────────────────────────────────────────────────

def _nchoose2(n: int) -> int:
    return n * (n - 1) // 2


def sim_bdh_age(age: float, numbsim: int,
                lambda_: float, mu: float, nu: float,
                hybprops: list[float], hyb_inher_fxn: Callable,
                mrca: bool = False,
                hyb_rate_fxn: Optional[Callable] = None,
                Ngene: int = 0,
                trait_model: Optional[dict] = None) -> list[dict]:
    """Run numbsim independent BDH simulations and return a list of result dicts.

    trait_model, if given, must supply callables under these keys:
        initial                  -- sequence of starting trait value(s), one per root leaf
        time_fxn(trait, dt)      -- evolve one lineage's trait over elapsed time dt
        spec_fxn(trait)          -- trait -> (child1_trait, child2_trait) at speciation
        hyb_event_fxn(t1, t2, inher)         -- parents' traits + inheritance -> hybrid trait
        hyb_compatibility_fxn(t1, t2, hyb_trait) -- bool: whether the hybridization can occur
    """
    return [_sim_one(age, lambda_, mu, nu, hybprops, hyb_inher_fxn, mrca, hyb_rate_fxn, Ngene, trait_model) for _ in range(numbsim)]


def _sim_one(age, lambda_, mu, nu, hybprops, hyb_inher_fxn, mrca, hyb_rate_fxn, Ngene, trait_model=None) -> dict:
    """Run one BDH simulation; return {phy: PhyloNetwork | 0, distance: dict}."""
    state = SimState(mrca=mrca, Ngene=Ngene, trait_model=trait_model)

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

        if state.trait_model is not None:
            state._evolve_traits(state.timestep)

        sp1 = _sample_one(state.leaves)
        rand = np.random.uniform()

        if rand < spec_rate / total_rate:
            state.speciation(sp1)
        elif rand < (spec_rate + ext_rate) / total_rate:
            state.extinction(sp1)
        else:
            _try_hybridization(state, sp1, hybprops, hyb_inher_fxn, hyb_rate_fxn)

    _finalize_pendant_edges(state)
    return _build_output(state)


def _sample_one(leaves: set[int]) -> int:
    return int(np.random.choice(list(leaves)))


def _try_hybridization(state: SimState, sp1,
                       hybprops, hyb_inher_fxn, hyb_rate_fxn):
    sp2 = _sample_one(state.leaves - {sp1})
    inher = hyb_inher_fxn()
    d12 = None
    hyb_trait = None

    if state.trait_model is not None:
        sp1_trait, sp2_trait = state.G.nodes[sp1]['trait'], state.G.nodes[sp2]['trait']
        hyb_trait = state.trait_model['hyb_event_fxn'](sp1_trait, sp2_trait, inher)
        if not state.trait_model['hyb_compatibility_fxn'](sp1_trait, sp2_trait, hyb_trait):
            return

    if hyb_rate_fxn is not None:
        d12 = state.tip_distance(sp1, sp2)
        if np.random.uniform() > hyb_rate_fxn(d12):
            return

    total = sum(hybprops)
    rand  = np.random.uniform()

    if rand < hybprops[0] / total:
        state.hyb_generating(sp1, sp2, inher, d12, hyb_trait)
    elif rand < sum(hybprops[:2]) / total:
        state.hyb_degenerative(sp1, sp2, inher, d12, hyb_trait)
    else:
        state.hyb_neutral(sp1, sp2, inher, d12, hyb_trait)


def _finalize_pendant_edges(state: SimState):
    """Seal all edges leading to currently active leaves."""
    for leaf in state.leaves:
        parent = next(state.G.predecessors(leaf))
        state._seal_edge(parent, leaf)


def _assign_labels(G: nx.DiGraph):
    for n, attrs in G.nodes(data=True):
        if attrs.get('extinct'):
            label = f"ext{n}"
        elif attrs.get('hyb_source'):
            label = f"hyb{n}"
        elif attrs['is_leaf']:
            label = f"sp{n}"
        else:
            label = f"-{n}"
        G.nodes[n]['label'] = label


def _build_output(state: SimState) -> dict:
    _assign_labels(state.G)

    tip_states = None
    if state.trait_model is not None:
        tip_states = {
            attrs['label']: attrs['trait']
            for n, attrs in state.G.nodes(data=True)
            if attrs['is_leaf'] or attrs.get('extinct') or attrs.get('hyb_source')
        }

    phy = PhyloNetwork(
        G=state.G,
        nleaves=len(state.leaves),
        tip_states=tip_states,
    )
    distance = dict(nx.all_pairs_dijkstra_path_length(state.G.to_undirected(),
                                                       weight='length'))
    return {'phy': phy, 'distance': distance}
