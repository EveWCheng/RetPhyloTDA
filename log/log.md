***2026-08-31 — Which reticulations a TDA-detected cycle captures, and why some strong ones are missed

## What we did

Ran the sim_phylo pipeline (`main.py`, seed 42, 71 non-extinct sims) with
`use_data_prep=False` (no phantom-node subdivision — see the separate
`max_node_per_edge` note), then analysed, per detected harmonic 1-cycle, which
reticulation edge each cycle's closing edge points back to.

Along the way we:

1. Found and fixed a shortest-path tie bug (below).
2. Retested the "cycle appears iff a tip pair is shortened past 50%" rule on the
   full 71-sim sample.
3. Characterised every reticulation that shortens a tip pair >50% but is **not**
   captured by any cycle, and explained all of them.

## The tie bug

The network has **138 `length == 0.0` tree edges** across 69 of 71 sims. With
zero-length edges, shortest paths through the network are frequently non-unique
(148 of 1048 leaf–leaf cycle edges have >1 shortest path; 14 cross *different*
reticulation edges depending on the route).

Two places picked a single arbitrary `nx.shortest_path` / `nx.single_source_dijkstra`:

- `main_result_analysis.top_cycle_edge_paths` — decides `on_closing_path`.
- `find_cycles.dress_distance_matrix_with_choice` — builds
  `reticulation_edges_on_path` in the dress log, the sole source of the
  `max_shortening` / `n_pairs` columns.

**Fix:** both now enumerate `nx.all_shortest_paths` (undirected, `weight="length"`)
and take the union of edges on *any* co-shortest path.
`main_result_analysis` gained `shortest_path_edges(net, u, v)` for this.

**Other changes made:**

- `top_k` in `main_result_analysis.main()` raised `1 → 3` (trace the 3 largest
  |harmonic-weight| cycle edges, not just the closing edge).
- Reticulate-edge table gained a `max_shortening_pair` column — the tip pair that
  achieved `max_shortening`. Dress log now stores tip **labels** in `i`/`j`
  (were numeric ids).
- `tests/test_dress_distance_matrix.py`: dropped 2 phantom-node tests that the
  `use_data_prep=False` switch made obsolete. Full suite: 54 passed.

Effect of the fixes on the "off closing path" count for reticulations shortening
>50%: **8 → 7** (tie fix recovered sim6 `hyb24-hyb32`, a genuine co-shortest
route) **→ 6** (`top_k=3` recovered sim23 `hyb32-hyb36`, rank 2–3 of its cycle).

## The 50% rule (retested, 71 sims)

Not a hard cutoff — a soft transition:

| closing-path shortening | outcome |
|---|---|
| stronger than ~ -0.52 | cycle every time (0 exceptions in 71 sims), persistent |
| -0.36 to -0.52 | usually a cycle, but faint; not guaranteed |
| weaker than ~ -0.36 | no cycle |

Sims **with** a cycle: strongest shortening -0.43 … -0.82. Sims **without**:
0.00 … -0.46. Small overlap in the -0.43…-0.46 band (e.g. sim35 has a cycle at
-0.43, sim31 has none at -0.46). The `1/2` is not fitted — it falls out of the
Vietoris–Rips filtration (proof below); node discreteness and non-round detours
blur the empirical edge to slightly weaker than `-0.5`.

## Why the threshold is exactly 1/2 — proof

### Setup

A detected cycle, as a metric space, is the geodesic circle of circumference
`L = P + r`:

- a finely-sampled **detour** arc `a = v_0, v_1, …, v_k = b`, segment lengths
  `ℓ_i`, `Σ ℓ_i = P`, all `ℓ_i ≤ δ ≪ r`;
- one **shortcut** edge `{a, b}` of length `r` (`0 < r < P`) — no sample points
  along its arc.

Metric `d(x, y)` = shorter of the two ways round, so `d(a, b) = r`. Filtration:
`VR_t` = simplices of diameter `≤ t`. `γ` = the fundamental loop
`[v_0 v_1] + … + [v_{k-1} v_k] + [v_k v_0]`.

`d_obs = r` is the dressed tip–tip distance; `d_tree = P` is the tree distance.
Shortening ratio `ρ = (d_obs − d_tree) / d_tree = r/P − 1`.

**Claim.** `H_1(VR_t)` carries a bar born at `t = r`, of positive length **iff
`r < P/2`**, i.e. **iff `ρ < −1/2`**.

### Lemma (three points cluster)

If `x, y, z ∈ S¹_L` have pairwise geodesic distance `< L/3`, they lie in a common
arc of length `< L/3`.

*Proof.* Order them; let `g_1, g_2, g_3 ≥ 0` be the gaps, `g_1 + g_2 + g_3 = L`.
From `d(x, z) = min(g_3, g_1 + g_2) < L/3`: either `g_1 + g_2 < L/3` (done), or
`g_3 < L/3`, forcing `g_1 + g_2 > 2L/3`. In that case `d(x, y) = min(g_1, L−g_1)`
and `d(y, z) = min(g_2, L−g_2)` are both `< L/3`; `g_1, g_2` cannot both be
`< L/3` (sum `> 2L/3`), so one exceeds `2L/3`, say `g_1`, and then `y, z, x` sit
in the arc of length `g_2 + g_3 = L − g_1 < L/3`. ∎

### Lower bound: a bar exists when `r < P/2`

Define the winding cochain `ω(x, y)` = signed length of the shorter arc from `x`
to `y` (well-defined whenever `d(x, y) < L/2`, and `ω(y, x) = −ω(x, y)`).

On `VR_t` with `t < L/3`: every edge has `d < L/3 < L/2`, and every triangle's
three vertices lie in a common arc `< L/3` (Lemma), on which the shorter arcs
telescope. So `ω(x, y) + ω(y, z) + ω(z, x) = 0` — `ω` is a **1-cocycle**.

Evaluate on `γ`: each detour edge contributes `+ℓ_i`, summing to `P`; the edge
`[b, a]` with `d(b, a) = r < L/3` continues forward by `r` (position
`P ↦ P + r ≡ 0`), contributing `+r`. Hence

```
⟨ω, γ⟩  =  P + r  =  L  ≠  0,
```

so `[γ] ≠ 0` in `H_1(VR_t; ℝ)` for every `r ≤ t < L/3`. `γ` first becomes a
cycle at `t = r` (its longest edge), so the bar is `[r, ≥ L/3)`, nonempty iff
`r < L/3`.

### Upper bound: no bar when `r ≥ P/2`

Let `m` be the detour midpoint, `d(a, m) = d(m, b) = P/2` (the shortcut gives no
shorter route to `m`). For `t ≥ max(r, P/2)`:

- `{a, m, b} ∈ VR_t` — all three pairwise distances (`P/2`, `P/2`, `r`) are `≤ t`;
- for each `v_j` between `a` and `m`, `d(a, v_j) = pos(v_j) ≤ P/2 ≤ t`, so
  `{a, v_j, v_{j+1}} ∈ VR_t` — these triangles make the sub-path `a → … → m`
  homologous rel endpoints to the single edge `[a, m]`; likewise `m → … → b`.

Therefore `γ ~ [a, m] + [m, b] + [b, a] = ∂{a, m, b} ~ 0`. When `r ≥ P/2`, `γ`
exists only for `t ≥ r ≥ P/2`, where `[γ] = 0` — no bar.

### Threshold

`r < L/3  ⟺  3r < P + r  ⟺  r < P/2`. Combining the bounds: a positive `H_1`
bar exists **iff `r < P/2`**, equivalently **`ρ < −1/2`**. The two bounds pin the
death to `[L/3, P/2]`; its exact value is `L/3`, so the barcode is
`birth = r`, `death = (P + r)/3`, `length = (P − 2r)/3`.

### What is and isn't already in the literature

- **The engine** — a Vietoris–Rips complex of the geodesic circle `S¹_L` is
  homotopy equivalent to `S¹` for `t ∈ (0, L/3]` (and to `S³, S⁵, …`, all with
  `H_1 = 0`, above that): **Adamaszek & Adams, "The Vietoris–Rips complexes of a
  circle," Pacific J. Math. 290 (2017), 1–40.** Gives the exact death `L/3`.
- **The general principle** — a geodesic loop of length `L` is seen by Rips
  exactly for `t < L/3`, with `H_1` persistence set by a shortest homology
  basis: **Virk, "1-dimensional intrinsic persistence of geodesic spaces"
  (2019)** and *"Approximations of 1-dimensional intrinsic persistence …"*.
- **The winding-number / degree-1 cocycle** (lower-bound technique here):
  **de Silva & Ghrist, "Coverage in sensor networks via persistent homology"
  (2007)**.
- **Metric graphs**: Gasparovic et al. (2018) characterise the intrinsic
  **Čech** `H_1` diagram of metric graphs (a length-`L` loop → `[0, L/4)` for
  Čech; `[0, L/3)` for Rips).

Not in the literature: the corollary in the form used here — loop = detour `P`
plus one bridging edge `r` from the dressed distance matrix, born at `r`, hence
detected iff `r < P/2`. It is a one-line specialization of Adamaszek–Adams /
Virk, not a new theorem, and the dressed-tip-distance framing is particular to
this project.

## Main result — the 6 strong reticulations that are missed

Of **92** reticulations that shorten some tip pair by >50%:

- **86** lie on a top-3 cycle closing path (captured). **0 false positives** in
  111 cycles — every captured reticulation is real.
- **6 missed:**

| sim | reticulation | n_pairs | max_short | most-shortened pair | inher_w |
|---|---|--:|--:|---|--:|
| 9  | hyb24 – hyb28 | 3 | -0.550 | sp25 – sp29  | 0.368 |
| 18 | -6 – hyb16    | 8 | -0.593 | sp31 – sp35  | 0.463 |
| 18 | -22 – hyb28   | 7 | -0.591 | hyb28 – sp35 | 0.495 |
| 27 | -18 – hyb26   | 4 | -0.624 | sp27 – sp28  | 0.481 |
| 33 | -24 – hyb30   | 5 | -0.550 | sp32 – sp35  | 0.492 |
| 42 | -11 – hyb22   | 6 | -0.596 | sp24 – sp32  | 0.394 |

### Why all 6 are missed

For every one, the `max_shortening_pair` shares **exactly one tip** with a
**captured** reticulation's `max_shortening_pair`:

| sim | missed pair | captured pair (= the cycle's closing edge) | shared tip |
|---|---|---|---|
| 9  | sp25 – **sp29**  | sp29 – hyb32                | sp29 |
| 18 | sp31 – **sp35**  | hyb34 – sp35                | sp35 |
| 18 | hyb28 – **sp35** | hyb34 – sp35                | sp35 |
| 27 | **sp27** – sp28  | sp27 – hyb30                | sp27 |
| 33 | **sp32** – **sp35** | sp32 – sp36  and  sp35 – hyb38 | sp32, sp35 |
| 42 | **sp24** – sp32  | sp24 – sp28                | sp24 |

The structure is consistent:

1. A captured reticulation's `max_shortening_pair` **is** the cycle's closing
   edge (the weight-1.0 edge). The closing edge is the most-shortened tip pair,
   and its shortcut runs through that reticulation.
2. The missed reticulation's `max_shortening_pair` shares one endpoint with that
   closing edge. So one tip (sp29, sp35, sp27, sp24) is a **pivot** carrying two
   strong shortcuts — one wins the closing edge, the other is a chord anchored on
   the same vertex.

A harmonic cycle is a global object: it passes through both endpoints of *every*
strong shortcut (checked — all 92 have both endpoints among the ~all-vertex
harmonic support). But only **one** tip pair per cycle becomes the closing edge —
the one whose late filtration appearance births the class. A second shortcut
anchored on an existing loop vertex does not open a new independent H1 class, so
it gets no separate bar and no separate closing edge. It is absorbed into the
same class, nudging the harmonic weights.

sim33's `-24-hyb30` is the extreme case: wedged between two cycles, sharing
`sp32` with cycle 1's closing edge and `sp35` with cycle 0's.

## Takeaway

`detected independent cycles  ≈  homologically independent strong shortcuts
                             ≤  reticulations shortening a tip pair >50%`

The "on closing path" metric is structurally **one reticulation per cycle**, so
multi-shortcut regions (reticulation chains, or fans anchored on a shared tip)
will always under-count by exactly the redundant shortcuts. Precision is 100%;
recall against all 376 reticulations is ~30% (most reticulations don't shorten
any pair enough), and against the 92 strong shortcuts it is 86/92 — the 6-edge
gap is fully explained.

To attribute the redundant shortcuts you would need to leave the closing-edge
heuristic: leave-one-out (drop reticulation R from `original_G` before dressing,
re-run the TDA, check whether a bar moves) is the honest per-reticulation test.

## What TDA actually buys here — and a simpler algorithm

The proof above collapses to a plain distance test. A cycle is born at the
shortened tip distance `d_obs` and dies at `(d_tree + d_obs)/3`, so it exists
iff `d_obs < d_tree / 2` — the tip pair is shortened by >50%. That is the whole
signal; the homology is bookkeeping on top of it.

Checked directly on the dress logs (all 71 sims): flag every leaf pair with
`(d_obs - d_tree) / d_tree <= -0.5`.

| threshold | flagged pairs | reticulations recovered | precision |
|---|--:|--:|--:|
| -0.30 | 474 | 199 / 376 (53%) | 100% |
| -0.40 | 311 | 145 / 376 (39%) | 100% |
| **-0.50** | 161 | **92 / 376 (24%)** | 100% |
| -0.60 | 56 | 40 / 376 (11%) | 100% |

Precision is 100% at every threshold and this is not a coincidence — it is
exact: a reticulation-aware distance can only beat the tree distance by using a
reticulation edge, so `d_obs < d_tree` **provably** implies a reticulation on the
i–j path. The 50% is not a detection threshold, it is a noise / strength cutoff.
The -0.5 test recovers all 92 "strong" reticulations, the 6 TDA merges included;
the ~half never recovered are minor events (low inheritance weight, short
reticulation branch) that barely move any distance.

So a defensible, simpler pipeline: estimate the tree -> per leaf pair compute
`(d_obs - d_tree) / d_tree` -> keep pairs past a strength cutoff (~-0.5 for
strong evidence) -> cluster kept pairs by shared path segment -> report the
reticulation edges. Counting and coherence-based noise rejection, which the
homology gives for free, are a few lines around the kept pairs.

**What TDA still adds:** essentially the derivation of the ~50% cutoff (why 50
and not 30 or 70), and one genuine capability — when there is **no single
reference tree** to difference against (a bag of discordant gene trees, the
`tree_main.py` framing), TDA finds loops in the merged distance structure
without first committing to a topology. In the species-tree + pairwise-distance
setting it is largely a repackaging of the shortening test.

**Bottom line:** the contribution is a validated, interpretable criterion — a
leaf pair whose realized distance is under half its tree distance is strong
evidence of a reticulation on that path — plus a mechanistic account of why the
threshold sits at 50%. Caveat for real data: `d_tree` there is estimated from
the same sequences, and reticulations bias that estimate toward smaller
divergences, so the observed shortening is attenuated below the true value
(affects any tree-differencing method, not just this one).
